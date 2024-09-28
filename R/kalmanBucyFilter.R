kalmanBucyFilter <- function(yuima, params, mean_init, vcov_init = NULL, delta.vcov.solve = 0.001, are = FALSE, explicit = FALSE, env = globalenv()) {
    # are : flag if use algebraic Riccati equation or not
    # validate input
    if(!inherits(yuima@model, "yuima.linear_state_space_model")) {
        yuima.stop("model must be yuima.linear_state_space_model")
    }
    
    ## distinguish observed/unobserved variables
    is.observed.equation <- yuima@model@is.observed
    observed.variables <- yuima@model@state.variable[is.observed.equation]
    unobserved.variables <- yuima@model@state.variable[!is.observed.equation]

    if(are) {
      if (!is.null(vcov_init)) {
        yuima.warn("vcov init is given, but ignored because are = TRUE")
      }
      vcov_init = matrix()
    } else {
      if(is.null(vcov_init)) {
        yuima.stop("'vcov_init' is required if are = FALSE")
      }
      ## if length of observed.variables is 1 and numeric of length 1 is given to vcov_init, convert vcov_init to matrix.
      if(inherits(vcov_init, "numeric")) {
        if(length(vcov_init) != 1) {
            yuima.stop("Invalid value for 'vcov_init', must be square matrix of size length(unobserved.variables) * length(unobserved.variables) or numeric of length 1.")
        }
        if(length(unobserved.variables) > 1) {
          yuima.stop("Invalid value for 'vcov_init', must be square matrix of size length(unobserved.variables) * length(unobserved.variables) when length(observed.variables) > 1.")
        }
        vcov_init <- matrix(vcov_init)
      } else if(inherits(vcov_init, "matrix")){
        if(!all(dim(vcov_init) == length(unobserved.variables))) {
            yuima.stop("Invalid value for 'vcov_init', must be square matrix of size length(unobserved.variables) * length(unobserved.variables) or numeric of length 1.")
        }
      } else {
          yuima.stop("Invalid value for 'vcov_init', must be square matrix of size length(unobserved.variable) * length(unobserved.variable) or numeric of length 1.")
      }
    }
    
    # create coefficient matrix of diffusion term of observed/unobserved variables
    col_num = length(yuima@model@diffusion[[1]])
    is.observed.column = rep(TRUE, col_num)
    is.unobserved.column = rep(TRUE, col_num)
    for(i in 1:col_num) {
      is.observed.column[i]   = all( yuima@model@is.observed | sapply(yuima@model@diffusion, function(row) as.character(row[i])=="(0)"))
      is.unobserved.column[i] = all(!yuima@model@is.observed | sapply(yuima@model@diffusion, function(row) as.character(row[i])=="(0)"))
    }
    if(!all(is.observed.column | is.unobserved.column)) {
      yuima.stop('Invalid diffusion matrix. Cannnot divide columns to observed/unobserved.')
    }
    observed.diffusion.expr   <- lapply(yuima@model@diffusion[is.observed.equation],  function(x) x[is.observed.column])
    unobserved.diffusion.expr <- lapply(yuima@model@diffusion[!is.observed.equation], function(x) x[is.unobserved.column])

    # get coefficient matrix of drift term of observed/unobserved variables
    observed.drift.slope.expr <- yuima@model@drift_slope[is.observed.equation]
    observed.drift.intercept.expr <- yuima@model@drift_intercept[is.observed.equation]
    unobserved.drift.slope.expr <- yuima@model@drift_slope[!is.observed.equation]
    unobserved.drift.intercept.expr <- yuima@model@drift_intercept[!is.observed.equation]

    
    tmp.env <- new.env(parent = env)
    for(i in 1:length(params)){
        assign(names(params)[i],params[[i]], envir=tmp.env)
    }
   
    if(are) {
        # Use algebraic Ricatti equation. a, b, c, sigma is independent of t. so returns evaluated values in 2-dim matrix.
        delta <- yuima@sampling@delta
        eval_exp <- function(expr) {
            env=tmp.env
            nrow = length(expr)
            ncol = length(expr[[1]])
            res = matrix(nrow=nrow, ncol=ncol)
            for( r in 1:nrow ) {
                for( c in 1:ncol ) {
                    res[r,c] <- eval(expr[[r]][c], envir=tmp.env)
                }
            }
            return(res)
        }
    } else {
         # evaluate each matrix expression
        delta_smaller_than = delta.vcov.solve
        K <- ceiling( yuima@sampling@delta / delta_smaller_than )
        delta <- yuima@sampling@delta / K
        n <- yuima@sampling@n[1]* K 
        time.points = as.matrix((0:n)*delta)
        # Use Ricatti equation. a, b, c, sigma is dependent of t. so returns evaluated values in 3-dim matrix.
        eval_exp <- function(expr) {
            variable=yuima@model@time.variable
            data=time.points
            env=tmp.env
            vec <- diffusionTermCpp(expr, variable, data, env)
            res <- array(vec, dim = c(length(expr),length(expr[[1]]), length(data)))
            return(res)
        }
    }
    unobserved.drift.slope <- eval_exp(unobserved.drift.slope.expr)
    unobserved.drift.intercept <- eval_exp(unobserved.drift.intercept.expr)
    unobserved.diffusion <- eval_exp(unobserved.diffusion.expr)
    observed.drift.slope <- eval_exp(observed.drift.slope.expr)
    observed.diffusion <- eval_exp(observed.diffusion.expr)
    observed.drift.intercept <- eval_exp(observed.drift.intercept.expr)

    # create other variables
    delta.observed.variable <- array(dim = c(length(observed.variables), length(yuima@data@zoo.data[[1]]) - 1), dimnames=list(observed.variables))
    for(variable in observed.variables) {
        delta.observed.variable[variable,] <- diff(matrix(yuima@data@zoo.data[[which(yuima@model@state.variable== variable)]]))
    }

    # calculate vcov and mean

    if(are) {
        ## solve algebraic Ricatti equation
        # solve vcov
        if(sum(yuima@model@is.observed) == 1 && sum(!yuima@model@is.observed) == 1) {
          res <- (observed.diffusion %*% t(observed.diffusion)) * abs(unobserved.drift.slope) / (observed.drift.slope^2) * (sqrt(1 + ((unobserved.diffusion %*% t(unobserved.diffusion)) * observed.drift.slope ^ 2)/((observed.diffusion %*% t(observed.diffusion)) * unobserved.drift.slope ^ 2)) - 1)
        } else {
          res <- calc_filter_vcov_are(
                  unobserved.drift.slope,
                  unobserved.diffusion,
                  observed.drift.slope,
                  observed.diffusion)
        }
        # solve mean
        if(explicit) {
            vcov <- res
            mean <- calc_filter_mean_explicit(
                unobserved.drift.slope,
                unobserved.drift.intercept,
                observed.drift.slope,
                observed.drift.intercept,
                observed.diffusion,
                vcov,
                mean_init,
                delta,
                delta.observed.variable
            )
            vcov <- array(vcov, dim=c(dim(vcov),yuima@sampling@n[1] + 1)) # vcov to return
        } else {
            arrayize <- function(mat) {
                array(mat, dim=c(dim(mat),yuima@sampling@n[1] + 1))
            }
            vcov <- arrayize(res)
            inv_sq_ob_diff <- arrayize(solve(observed.diffusion %*% t(observed.diffusion)))
            unobserved.drift.intercept <- arrayize(unobserved.drift.intercept)
            unobserved.drift.slope <- arrayize(unobserved.drift.slope)
            unobserved.diffusion <- arrayize(unobserved.diffusion)
            observed.drift.intercept <- arrayize(observed.drift.intercept)
            observed.drift.slope <- arrayize(observed.drift.slope)
            observed.diffusion <- arrayize(observed.diffusion)
            mean <- calc_filter_mean(
                unobserved.drift.slope, 
                unobserved.drift.intercept,
                unobserved.diffusion,
                observed.drift.slope,
                observed.drift.intercept,
                observed.diffusion,
                vcov, 
                inv_sq_ob_diff,
                mean_init,
                yuima@sampling@delta, 
                delta.observed.variable
            )
        }
    } else {
	    res <- calc_filter_vcov(
            unobserved.drift.slope, 
            unobserved.diffusion,
            observed.drift.slope,
            observed.diffusion,
            vcov_init, delta)
        vcov <- res$vcov
        inv_sq_ob_diff <- res$inv_sq_ob_diff
        
        # subsampling for solving m
        get.subsampling <- function(array.obj) {
            res <- array.obj[,,0:yuima@sampling@n[1] * K + 1]
            res <- array(res, dim=c(dim(array.obj)[1],dim(array.obj)[2],yuima@sampling@n[1] + 1))
            return(res)
        }
        unobserved.drift.slope <- get.subsampling(unobserved.drift.slope)
        unobserved.drift.intercept <- get.subsampling(unobserved.drift.intercept)
        unobserved.diffusion <- get.subsampling(unobserved.diffusion)
        
        observed.drift.slope <- get.subsampling(observed.drift.slope)
        observed.drift.intercept <- get.subsampling(observed.drift.intercept)
        observed.diffusion <- get.subsampling(observed.diffusion)
        
        vcov <- get.subsampling(vcov)
        inv_sq_ob_diff <- get.subsampling(inv_sq_ob_diff)
        
        # solve m
        mean <- calc_filter_mean(
            unobserved.drift.slope, 
            unobserved.drift.intercept,
            unobserved.diffusion,
            observed.drift.slope,
            observed.drift.intercept,
            observed.diffusion,
            vcov, 
            inv_sq_ob_diff,
            mean_init,
            yuima@sampling@delta, 
            delta.observed.variable)
    }
    rownames(mean) <- unobserved.variables
    ts.mean <- ts(t(mean), start = start(yuima@data@zoo.data[[1]]), frequency = frequency(yuima@data@zoo.data[[1]]))
    res <- new(
        "yuima.kalmanBucyFilter",
        model=yuima@model,
        mean=ts.mean,
        vcov=vcov,
        mean.init=mean_init,
        vcov.init=vcov_init,
        delta=yuima@sampling@delta,
        data=yuima@data
    )
    return(res)
}
