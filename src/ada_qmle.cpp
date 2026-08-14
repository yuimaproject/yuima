#include <RcppArmadillo.h>
#include <cmath>
#include <vector>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

namespace {

bool inverse_spd(const arma::mat& a, arma::mat& inverse, double& logdet) {
    double sign = 0.0;
    if (!arma::log_det(logdet, sign, a) || sign <= 0.0 || !std::isfinite(logdet)) {
        return false;
    }
    return arma::inv_sympd(inverse, a) && inverse.is_finite();
}

arma::mat matrix_term(const NumericMatrix& term, int observation, int dimension) {
    arma::mat out(dimension, dimension);
    for (int column = 0; column < dimension; ++column) {
        for (int row = 0; row < dimension; ++row) {
            out(row, column) = term(observation, row + dimension * column);
        }
    }
    return out;
}

arma::vec vector_term(const NumericMatrix& term, int observation, int dimension) {
    arma::vec out(dimension);
    for (int row = 0; row < dimension; ++row) {
        out(row) = term(observation, row);
    }
    return out;
}

} // namespace

// Evaluate arbitrary scalar model terms at all observation points.  Keeping this
// alongside the contrast kernels avoids repeatedly allocating one R environment
// per expression in the R-level optimizer.
// [[Rcpp::export]]
NumericMatrix adaEvalTermsCpp(ExpressionVector terms,
                              CharacterVector state,
                              const arma::mat& data,
                              CharacterVector timeVariable,
                              NumericVector time,
                              Environment env) {
    const int observations = data.n_rows;
    const int states = state.size();
    NumericMatrix result(observations, terms.size());

    if (data.n_cols != static_cast<unsigned int>(states)) {
        stop("The number of state variables does not match the data dimension.");
    }
    const bool hasTime = timeVariable.size() == 1;
    if (hasTime && time.size() != observations) {
        stop("The time vector does not match the number of observations.");
    }

    for (int observation = 0; observation < observations; ++observation) {
        for (int j = 0; j < states; ++j) {
            env.assign(as<std::string>(state[j]), data(observation, j));
        }
        if (hasTime) {
            env.assign(as<std::string>(timeVariable[0]), time[observation]);
        }
        for (int j = 0; j < terms.size(); ++j) {
            SEXP value = Rf_eval(terms[j], env);
            result(observation, j) = as<double>(value);
        }
    }
    return result;
}

// C++ kernel for the initial V/W contrasts, the plug-in V contrasts and the
// full U contrast.  meanTerms contains the coefficients of h^j, j >= 2, in
// the conditional mean. momentTerms contains the coefficients of h^j,
// j >= 2, in E[(X_{t+h}-X_t)(X_{t+h}-X_t)'].
// [[Rcpp::export]]
double adaContrastCpp(const arma::mat& increments,
                      const arma::mat& drift,
                      const arma::mat& diffusion,
                      List meanTerms,
                      List momentTerms,
                      double h,
                      std::string contrast,
                      int order) {
    const int observations = increments.n_rows;
    const int dimension = increments.n_cols;
    if (h <= 0.0 || !std::isfinite(h)) {
        return R_PosInf;
    }
    if (drift.n_rows != increments.n_rows || drift.n_cols != increments.n_cols ||
        diffusion.n_rows != increments.n_rows ||
        diffusion.n_cols % dimension != 0) {
        stop("Incompatible dimensions in adaQmle contrast inputs.");
    }

    double result = 0.0;
    const int k = order / 2;

    for (int observation = 0; observation < observations; ++observation) {
        const arma::vec dx = increments.row(observation).t();
        const arma::vec b = drift.row(observation).t();
        const int wienerDimension = diffusion.n_cols / dimension;
        arma::mat sigma(dimension, wienerDimension);
        for (int column = 0; column < wienerDimension; ++column) {
            for (int row = 0; row < dimension; ++row) {
                sigma(row, column) = diffusion(observation, row + dimension * column);
            }
        }
        const arma::mat a = sigma * sigma.t();

        if (contrast == "initial_least_squares_diffusion") {
            const arma::mat residual = dx * dx.t() - h * a;
            result += arma::accu(arma::square(residual)) / (h * h);
            continue;
        }
        if (contrast == "initial_least_squares_drift") {
            const arma::vec residual = dx - h * b;
            result += arma::dot(residual, residual) / h;
            continue;
        }

        arma::mat inverseA;
        double logdetA = 0.0;
        if (!inverse_spd(a, inverseA, logdetA)) {
            return R_PosInf;
        }

        if (contrast == "initial_gaussian_diffusion") {
            result += logdetA + arma::as_scalar(dx.t() * inverseA * dx) / h;
            continue;
        }
        if (contrast == "initial_gaussian_drift") {
            const arma::vec residual = dx - h * b;
            result += arma::as_scalar(residual.t() * inverseA * residual) / h;
            continue;
        }
        if (contrast == "plugin_diffusion") {
            arma::mat adjusted = dx * dx.t();
            const int highestMoment = (order + 1) / 2;
            for (int j = 2; j <= highestMoment; ++j) {
                NumericMatrix term = momentTerms[j - 2];
                adjusted -= std::pow(h, j) * matrix_term(term, observation, dimension);
            }
            result += logdetA + arma::trace(inverseA * adjusted) / h;
            continue;
        }
        if (contrast == "plugin_drift") {
            arma::vec residual = dx - h * b;
            const int highestMean = order / 2;
            for (int j = 2; j <= highestMean; ++j) {
                NumericMatrix term = meanTerms[j - 2];
                residual -= std::pow(h, j) * vector_term(term, observation, dimension);
            }
            result += arma::as_scalar(residual.t() * inverseA * residual) / h;
            continue;
        }
        if (contrast != "coupled") {
            stop("Unknown adaQmle contrast mode.");
        }

        // gamma[j] is the coefficient of h^(j+1) in the conditional
        // covariance. gamma[0] = a.  The mean coefficient of h is b.
        std::vector<arma::vec> mean(k + 1);
        mean[1] = b;
        for (int j = 2; j <= k; ++j) {
            NumericMatrix term = meanTerms[j - 2];
            mean[j] = vector_term(term, observation, dimension);
        }

        std::vector<arma::mat> gamma(k + 1);
        gamma[0] = a;
        for (int j = 2; j <= k + 1; ++j) {
            NumericMatrix term = momentTerms[j - 2];
            arma::mat covarianceCoefficient = matrix_term(term, observation, dimension);
            for (int left = 1; left < j; ++left) {
                const int right = j - left;
                if (left <= k && right <= k) {
                    covarianceCoefficient -= mean[left] * mean[right].t();
                }
            }
            gamma[j - 1] = covarianceCoefficient;
        }

        // Coefficients of the formal inverse of
        // a + h gamma[1] + ... + h^k gamma[k].
        std::vector<arma::mat> inverseCoefficient(k + 1);
        inverseCoefficient[0] = inverseA;
        for (int degree = 1; degree <= k; ++degree) {
            arma::mat convolution(dimension, dimension, arma::fill::zeros);
            for (int j = 1; j <= degree; ++j) {
                convolution += gamma[j] * inverseCoefficient[degree - j];
            }
            inverseCoefficient[degree] = -inverseA * convolution;
        }

        arma::mat inverseExpansion(dimension, dimension, arma::fill::zeros);
        for (int degree = 0; degree <= k; ++degree) {
            inverseExpansion += std::pow(h, degree) * inverseCoefficient[degree];
        }

        // d/dh log det(C(h)) = tr(C(h)^(-1) C'(h)); integrating
        // its formal series gives the required log-determinant coefficients.
        double logdetExpansion = logdetA;
        for (int degree = 1; degree <= k; ++degree) {
            double derivativeCoefficient = 0.0;
            for (int inverseDegree = 0; inverseDegree < degree; ++inverseDegree) {
                const int covarianceDegree = degree - inverseDegree;
                derivativeCoefficient += covarianceDegree *
                    arma::trace(inverseCoefficient[inverseDegree] * gamma[covarianceDegree]);
            }
            logdetExpansion += std::pow(h, degree) * derivativeCoefficient / degree;
        }

        arma::vec residual = dx;
        for (int j = 1; j <= k; ++j) {
            residual -= std::pow(h, j) * mean[j];
        }
        const double contribution = logdetExpansion +
            arma::as_scalar(residual.t() * inverseExpansion * residual) / h;
        if (!std::isfinite(contribution)) {
            return R_PosInf;
        }
        result += contribution;
    }

    return std::isfinite(result) ? result : R_PosInf;
}
