poest = function(yuima, estimator = "qmle", start, prior=NULL, lower, upper, penalty, method = "L-BFGS-B", ...) {
  ### Explanation of Parameters
  # yuima:                  a yuima object.
  # estimator:              estimator used before and after LASSO step. Currently, 'qmle' and 'adaBayes' are available.
  # start:                  initial value to be passed to optimizer in first QMLE.
  # prior:                  a list of prior distributions for the parameters specified by 'code'. If 'estimator' is 'adaBayes', 'prior' is required.
  #                         Currently, dunif(z, min, max), dnorm(z, mean, sd), dbeta(z, shape1, shape2), dgamma(z, shape, rate) are available.
  # lower:                  a named list for specifying lower bounds of parameters.
  # upper:                  a named list for specifying upper bounds of parameters.
  # penalty:                list of 2 sub-lists named drift and diffusion. each sub-list must contain q , gamma and r.
  # ...:                    parameters for "qmle" function and "optim" function.
  
  ### Step0: validation of arguments
    if(missing(yuima)) {
      yuima.stop("yuima object is missing.")
    }

    if(missing(penalty)) {
      penalty <- list(
        drift = list(q = 1, gamma = 1, r = 1.5),
        diffusion = list(q = 1, gamma = 1, r = 1.5)
      )
    } else {
      if(!is.element("drift", names(penalty))) {
        penalty$drift <- list(q = 1, gamma = 1, r = 1.5)
      } else {
        if(!is.element("q", names(penalty$drift))) {
          penalty$drift$q <- 1
        }
        if(!is.element("gamma", names(penalty$drift))) {
          penalty$drift$gamma <- 1
        }
        if(!is.element("r", names(penalty$drift))) {
          penalty$drift$r <- (3 - penalty$drift$q + penalty$drift$gamma) / 2
        }

        # check the conditions of q and gamma.
        if(penalty$drift$q <= 0 || 1 < penalty$drift$q) {
          yuima.stop("Parameter \"penalty$drift$q\" for LASSO step must satisfy: 0 < q <= 1.")
        }
        if(penalty$drift$gamma <= -(1 - penalty$drift$q)) {
          yuima.stop("Parameter \"penalty$drift$gamma\" for LASSO step must satisfy: penalty$drift$gamma > -(1 - penalty$drift$q).")
        }
      }
      
      if(!is.element("diffusion", names(penalty))) {
        penalty$diffusion <- list(q = 1, gamma = 1, r = 1.5)
      } else {
        if(!is.element("q", names(penalty$diffusion))) {
          penalty$diffusion$q <- 1
        }
        if(!is.element("gamma", names(penalty$diffusion))) {
          penalty$diffusion$gamma <- 1
        }
        if(!is.element("r", names(penalty$diffusion))) {
          penalty$diffusion$r <- (3 - penalty$diffusion$q + penalty$diffusion$gamma) / 2
        }
      
        # check the conditions of q and gamma.
        if(penalty$diffusion$q <= 0 || 1 < penalty$diffusion$q) {
          yuima.stop("Parameter \"penalty$diffusion$q\" for LASSO step must satisfy: 0 < q <= 1.")
        }
        if(penalty$diffusion$gamma <= -(1 - penalty$diffusion$q)) {
          yuima.stop("Parameter \"penalty$diffusion$gamma\" for LASSO step must satisfy: penalty$diffusion$gamma > -(1 - penalty$diffusion$q).")
        }

      }
    }
    
  ### Step1: estimate parameters using QMLE
    if(estimator == "qmle") {
      init.est <- qmle(yuima, start = start, lower = lower, upper = upper, method = method, ...)
    } else if(estimator == "adaBayes") {
      if(is.null(prior)) {
        yuima.stop("When \"estimator\" is \"adaBayes\", \"prior\" is required.")
      }
      init.est <- adaBayes(yuima, prior = prior, start = start, lower = lower, upper = upper, method = method, ...)
    } else {
      yuima.stop("At present, \"estimator\" must be \"adaBayes\" or \"qmle\"")
    }
    init.est.coef <- init.est@coef
    # separate the diffusion and drift parameters.
    diffusion.init.est <- init.est.coef[is.element(names(init.est.coef), yuima@model@parameter@diffusion)]
    drift.init.est <- init.est.coef[is.element(names(init.est.coef), yuima@model@parameter@drift)]

    class(lower) <- "numeric"
    diffusion.lower <- lower[is.element(names(lower), yuima@model@parameter@diffusion)]
    drift.lower <- lower[is.element(names(lower), yuima@model@parameter@drift)]

    class(upper) <- "numeric"
    diffusion.upper <- upper[is.element(names(upper), yuima@model@parameter@diffusion)]
    drift.upper <- upper[is.element(names(upper), yuima@model@parameter@drift)]
    
    # check that there is no duplicate in diffusion parameters and drift parameters.
    if(length(diffusion.init.est) + length(drift.init.est) != length(yuima@model@parameter@all)) {
      yuima.stop("There is duplication in dissusion paramters and drift parameters.")
    }


  ### Step2: select non-zero components using LASSO
    
    # function to check the condition of parameters: q, gamma, alpha, and r. returns alpha.
    get.alpha <- function(params, n) {
      q <- params$q
      gamma <- params$gamma
      r <- params$r

      # validation is done before STEP1
      # # check the conditions of q and gamma.
      # if(q <= 0 || 1 < q) {
      #   yuima.stop("Parameter \"q\" for LASSO step must satisfy: 0 < q <= 1.")
      # }
      # if(gamma <= -(1 - q)) {
      #   yuima.stop("Parameter \"gamma\" for LASSO step must satisfy: gamma > -(1 - q).")
      # }
      
      alpha <- n^(-r/2)
      return(alpha)
    }
    
    # function to calculate directly when q = 1.
    f1 <- function(x) {

      theta <- x["theta"]
      kappa <- x["kappa"]
      lower <- x["lower"]
      upper <- x["upper"]

      if(theta >  kappa / 2) {
        res <- (theta - kappa / 2)
      } else if(theta < -kappa / 2) {
        res <- (theta + kappa / 2)
      } else {
        res <- 0
      }

      if(res < lower) {
        res <- lower
      }
      if(upper < res) {
        res <- upper
      }

      return(res) 
    }
    
    # number of observations
    n <- yuima@sampling@n[1]

    # length of observation time
    t <- n * yuima@sampling@delta


    # function to optimize when q < 1.
    f2 <- function(x, alpha, gamma, q, n, ...) {
      # objective function to be used when q < 1.
      par <- x["par"]
      lower <- x["lower"]
      upper <- x["upper"]

      # Q <- function(x) {
      #   kappa <- alpha * abs(par) ^ (- gamma)
      #   (x - par)^2 + kappa * abs(x) ^ q
      # }

      # res.optim <- optim(
      #   # par = par,
      #   par = 0,
      #   fn = Q, 
      #   lower = lower, upper = upper, 
      #   method = method, ...
      # )
      # return(res.optim$par)

      kappa <- alpha * abs(par) ^ (- gamma)
      tmp <- (2/(2-q)) * (2*(1-q)/(2-q))^(1-q)*abs(par)^(2-q)
      if(all(kappa > tmp, lower <= 0, upper >= 0)) return(0)
      else return(par)
    }
    
    
    ### calculation of LASSO
    
    # diffusion term.
    if(length(yuima@model@parameter@diffusion) > 0) {
      # print("Diffusion Lasso Start")
      diffusion.alpha <- get.alpha(params = penalty$diffusion, n =n)
      if(penalty$diffusion$q == 1) {
        # calculate directly when q = 1.
        df <- data.frame(
          theta = diffusion.init.est,
          kappa = diffusion.alpha * abs(diffusion.init.est) ^ (- penalty$diffusion$gamma),
          lower = diffusion.lower,
          upper = diffusion.upper
        )
        theta.diffusion.lasso <- apply(
          df, MARGIN = 1, FUN = f1
        )
      } else {
        # use "optim" function for each cordinate when q < 1.
        df <- data.frame(
          par = diffusion.init.est,
          lower = diffusion.lower,
          upper = diffusion.upper
        )
        
        theta.diffusion.lasso <- apply(
          df, MARGIN = 1, FUN = f2, 
          alpha = diffusion.alpha, gamma = penalty$diffusion$gamma, q = penalty$diffusion$q, n = n
        )
      }
    } else {
      theta.diffusion.lasso <- NULL
    }

    
    # drift term.
    if(length(yuima@model@parameter@drift) > 0) {
      # print("Drift Lasso Start")
      drift.alpha = get.alpha(params = penalty$drift, n = t)
      if(penalty$drift$q == 1) {
        # calculate directly when q = 1.
        df <- data.frame(
          theta = drift.init.est,
          kappa = drift.alpha * abs(drift.init.est) ^ (- penalty$drift$gamma),
          lower = drift.lower,
          upper = drift.upper
        )
        theta.drift.lasso <- apply(
          df, MARGIN = 1, FUN = f1
        )
      } else {
        # use "optim" function for each cordinate when q < 1.
        df <- data.frame(
          par = drift.init.est,
          lower = drift.lower,
          upper = drift.upper
        )
        
        theta.drift.lasso <- apply(
          df, MARGIN = 1, FUN = f2, 
          alpha = drift.alpha, gamma = penalty$drift$gamma, q = penalty$drift$q, n = t
        )
      }
    } else {
      theta.drift.lasso <- NULL
    }
    
    # result of LASSO
    theta.lasso <- c(
      theta.diffusion.lasso, 
      theta.drift.lasso
    )

  ### Step3: estimate parameters using QMLE again

    fixed <- as.list(theta.lasso[theta.lasso == 0])
    start = as.list(init.est.coef)
    if(estimator == "qmle") {
      final.est <- qmle(yuima, start = start, lower = lower, upper = upper, method = method, fixed = fixed, ...)
    } else if(estimator == "adaBayes") {
      final.est <- adaBayes(yuima, prior = prior, start = start, lower = lower, upper = upper, method = method, fixed = fixed, ...)
    }
    
    est <- new("yuima.poest", 
      initial_estimator = init.est,
      final_estimator = final.est,
      model = yuima@model,
      estimator = estimator,
      call = match.call(), 
      coef = final.est@fullcoef,
      vcov = final.est@vcov,
      selected_coef = names(theta.lasso[theta.lasso != 0])
      )

    return(est)

}

setClass("yuima.poest",
  representation(
    initial_estimator="ANY",
    final_estimator="ANY",
    model="yuima.model",
    estimator = "character",
    call="call",
    coef="numeric",
    vcov="matrix",
    selected_coef="character"
  )
)
