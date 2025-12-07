#' Class for the Asymptotic Expansion of Diffusion Processes
#' 
#' The \code{yuima.ae} class is used to describe the output of the functions \code{\link{ae}} and \code{\link{aeMarginal}}.
#' Two generic methods are provided for the class: \code{\link{initialize,yuima.ae-method}} and \code{\link{plot,yuima.ae-method}}.
#' 
#' @slot order integer. The order of the expansion.
#' @slot var character. The state variables.
#' @slot u.var character. The variables of the characteristic function.
#' @slot eps.var character. The perturbation variable.
#' @slot characteristic expression. The characteristic function.
#' @slot density expression. The probability density function.
#' @slot Z0 numeric. The solution to the deterministic process obtained by setting the perturbation to zero.
#' @slot Mu numeric. The drift vector for the representation of Z1.
#' @slot Sigma matrix. The diffusion matrix for the representation of Z1.
#' @slot c.gamma list. The coefficients of the Hermite polynomials.
#' @slot h.gamma list. Hermite polynomials.
#' 
setClass("yuima.ae", slots = c(
  order = "integer",
  var = "character",
  u.var = "character",
  eps.var = "character",
  characteristic = "expression",
  density = "expression",
  Z0 = "numeric",
  Mu = "numeric",
  Sigma = "matrix",
  c.gamma = "list",
  h.gamma = "list"
))

#' @title Constructor for \code{yuima.ae} Class
#' 
#' @description Construct an object of class \code{\link{yuima.ae-class}}.
#' 
#' @aliases initialize,yuima.ae-method
#' @aliases initialize,yuima.ae,ANY-method
#' 
#' @param .Object an object of class \code{\link{yuima.ae-class}}.
#' @param order integer. The order of the expansion.
#' @param var character. The state variables.
#' @param u.var character. The variables of the characteristic function.
#' @param eps.var character. The perturbation variable.
#' @param characteristic expression. The characteristic function.
#' @param Z0 numeric. The solution to the deterministic process obtained by setting the perturbation to zero.
#' @param Mu numeric. The drift vector for the representation of Z1.
#' @param Sigma matrix. The diffusion matrix for the representation of Z1.
#' @param c.gamma list. The coefficients of the Hermite polynomials.
#' @param verbose logical. Print on progress? Default \code{FALSE}.
#' 
setMethod("initialize", "yuima.ae", function(.Object, order, var, u.var, eps.var, characteristic, Z0, Mu, Sigma, c.gamma, verbose){
  
  ######################################################## 
  # Set density                                          #
  ########################################################
  if(verbose) {
    cat(paste0('Computing distribution...'))
    time <- Sys.time()
  }
  
  tmp <- calculus::hermite(var = var, sigma = solve(Sigma), order = 3*order)
  h.gamma <- lapply(tmp, function(x) x$f)
  .Object@density <- sapply(0:order, function(m) {
    
    if(m>0){
      pdf <- sapply(1:m, function(m){
        g <- c.gamma[[m]][sapply(c.gamma[[m]], function(x) x!=0)]
        if(length(g)==0) return(0)
        h <- h.gamma[names(g)]
        p <- paste(calculus::wrap(g), calculus::wrap(h), sep = "*", collapse = " + ")
        paste(paste0(eps.var, "^", m), calculus::wrap(p), sep = " * ")
      })
      pdf <- paste(pdf, collapse = " + ")
      pdf <- paste0(1, " + ", pdf)
    }
    else {
      pdf <- 1
    }
    
    kernel <- sprintf('%s * exp(%s)', ((2*pi)^(-length(var)/2)*(det(Sigma)^(-0.5))), (-0.5 * solve(Sigma)) %inner% (var %outer% var))
    
    pdf <- paste0(kernel, " * (", pdf, ")")
    
    for(i in 1:length(var)) {
      z <- sprintf("(((%s - %s)/%s) - %s)", var[i], Z0[var[i]], eps.var, Mu[i])
      pdf <- gsub(x = pdf, pattern = paste0("\\b",var[i],"\\b"), replacement = z) 
      pdf <- sprintf("(%s) / abs(%s)", pdf, eps.var)
    }
    
    return(parse(text = pdf))
    
  })
  
  if(verbose) {
    cat(sprintf(' (%s sec)\n', difftime(Sys.time(), time, units = "secs")[[1]]))
    time <- Sys.time()
  }
  
  ######################################################## 
  # Set object                                           #
  ########################################################
  .Object@order <- order
  .Object@var <- var
  .Object@u.var <- u.var
  .Object@eps.var <- eps.var
  .Object@characteristic <- characteristic
  .Object@Z0 <- Z0
  .Object@Mu <- Mu
  .Object@Sigma <- Sigma
  .Object@c.gamma <- c.gamma
  .Object@h.gamma <- h.gamma
  
  # return 
  return(.Object)
  
})

#' @title Plot Method for \code{yuima.ae} Class
#' 
#' @description Plot an object of class \code{\link{yuima.ae-class}}.
#' 
#' @aliases plot,yuima.ae-method
#' @aliases plot,yuima.ae,ANY-method
#' 
#' @param x an object of class \code{\link{yuima.ae-class}}.
#' @param grids list. A named list of vectors specifying the grid to evaluate the density. The names must match the state variables.
#' @param eps numeric. The intensity of the perturbation.
#' @param order integer. The expansion order. If \code{NULL} (default), it uses the maximum order used in \code{ae}.
#' @param ... additional arguments passed to the plot function.
#' 
setMethod("plot", signature(x = "yuima.ae"), function(x, grids = list(), eps = 1, order = NULL, ...){
  
  n <- length(x@var)
  
  for(z in x@var){
    margin <- aeMarginal(ae = x, var = z)
    grid <- grids[z]
    if(is.null(grid[[1]])){
      mu <- aeMean(ae = margin, eps = eps, order = order)[[1]]
      sd <- aeSd(ae = margin, eps = eps, order = order)[[1]]
      grid <- list(seq(mu-5*sd, mu+5*sd, length.out = 1000))
      names(grid) <- z
    }
    dens <- do.call('aeDensity', c(grid, list(ae = margin, eps = eps, order = order)))
    plot(x = grid[[1]], y = dens, type = 'l', xlab = z, ylab = 'Density', ...)
  }
  
})

#' Asymptotic Expansion
#' 
#' Asymptotic expansion of uni-dimensional and multi-dimensional diffusion processes.
#' 
#' @param model an object of \code{\link{yuima-class}} or \code{\link{yuima.model-class}}.
#' @param xinit initial value vector of state variables.
#' @param order integer. The asymptotic expansion order. Higher orders lead to better approximations but longer computational times.
#' @param true.parameter named list of parameters.
#' @param sampling a \code{\link{yuima.sampling-class}} object.
#' @param eps.var character. The perturbation variable.
#' @param solver the solver for ordinary differential equations. One of \code{"rk4"} (more accurate) or \code{"euler"} (faster).
#' @param verbose logical. Print on progress? Default \code{FALSE}.
#' 
#' @return An object of \code{\link{yuima.ae-class}}
#' 
#' @details 
#' If \code{sampling} is not provided, then \code{model} must be an object of \code{\link{yuima-class}} with non-empty \code{sampling}.
#' 
#' if \code{eps.var} does not appear in the model specification, then it is internally added in front of the diffusion matrix to apply the asymptotic expansion scheme.
#' 
#' @author 
#' Emanuele Guidotti <emanuele.guidotti@unine.ch>
#' 
#' @examples
#' \dontrun{
#' # model
#' gbm <- setModel(drift = 'mu*x', diffusion = 'sigma*x', solve.variable = 'x')
#' 
#' # settings
#' xinit <- 100
#' par <- list(mu = 0.01, sigma = 0.2)
#' sampling <- setSampling(Initial = 0, Terminal = 1, n = 1000)
#' 
#' # asymptotic expansion
#' approx <- ae(model = gbm, sampling = sampling, order = 4, true.parameter = par, xinit = xinit)
#' 
#' # exact density
#' x <- seq(50, 200, by = 0.1)
#' exact <- dlnorm(x = x, meanlog = log(xinit)+(par$mu-0.5*par$sigma^2)*1, sdlog = par$sigma*sqrt(1))
#' 
#' # compare
#' plot(x, exact, type = 'l', ylab = "Density")
#' lines(x, aeDensity(x = x, ae = approx, order = 1), col = 2)
#' lines(x, aeDensity(x = x, ae = approx, order = 2), col = 3)
#' lines(x, aeDensity(x = x, ae = approx, order = 3), col = 4)
#' lines(x, aeDensity(x = x, ae = approx, order = 4), col = 5)
#' }
#' @importFrom calculus %dot%
#' @importFrom calculus %inner%
#' @importFrom calculus %mx%
#' @importFrom calculus %outer%
#' @importFrom calculus %prod%
#' @importFrom calculus %sum%
#' @importFrom calculus %gradient%
#' 
#' @export
#' 
ae <- function(model, xinit, order = 1L, true.parameter = list(), sampling = NULL, eps.var = 'eps', solver = "rk4", verbose = FALSE){
  
  obj.class <- class(model)
  
  if(!is.null(sampling)){
    
    if(obj.class=='yuima'){
      stop('model must be of class yuima.model when sampling is provided')
    }
    else if(obj.class=='yuima.model'){
      model <- setYuima(model = model, sampling = sampling)
    }
    else stop('model must be of class yuima or yuima.model')
    
  }
  else if(obj.class!='yuima'){
    stop('model must be of class yuima when sampling is not provided')
  }
  
  if(!model@sampling@regular)
    stop('ae needs regular sampling')
  
  if(length(sampling@grid)!=1)
    stop('ae needs a unidimensional sampling grid')
  
  if(length(model@model@jump.coeff)>0)
    stop('ae does not support jump processes')
  
  if(model@model@hurst!=0.5)
    stop('ae does not support fractional processes')
  
  # temporary disable calculus.auto.wrap 
  calculus.auto.wrap <- options(calculus.auto.wrap = FALSE)
  on.exit(options(calculus.auto.wrap), add = TRUE)
  
  
  ######################################################## 
  # dz                                                   #
  ########################################################
  dz <- function(i){
    
    paste0('d', i, '.', AE$z)
    
  }
  
  dz.nu <- function(i, nu){
    
    dz <- dz(i)
    
    z <- sapply(1:AE$d, function(i) rep(dz[i], nu[i]), simplify = F)
    z <- paste(unlist(z), collapse = ' * ') 
    
    if(z=='') z <- '1'
    return(z)
    
  }
  
  ######################################################## 
  # Z.I                                                  #
  ########################################################
  Z.I <- function(I, bar = FALSE) {
    
    tmp <- NULL
    for(i in I){
      
      z.i <- dz(i = i)
      if(bar) {
        z.i <- c( z.i, ifelse(i==1, 1, 0) )
      }
      z.i <- array(z.i)
      
      if(is.null(tmp)) tmp <- z.i
      else tmp <- tmp %outer% z.i
      
    }
    
    return(tmp)
    
  }

  # Hermite Valued Expectation
  HVE <- function(nu, K){
    
    # convert to label
    nu <- label(nu)
    
    # for each nu' 
    H <- sapply(names(AE$c.nu[[nu]]), function(nu.prime) {
      c(AE$c.nu[[nu]][[nu.prime]]) * as.numeric(AE$Ez.T[AE$Ez.K[[label(K)]][[nu.prime]]])
    })
    
    if(is.null(dim(H))) 
      H <- sum(H)
    else 
      H <- rowSums(H)
    
    return(H)
  }

  # Tensor Valued Expectation
  TVE <- function(K){
    
    K <- K + 1
    
    nu <- AE$nu[[sum(K)]]
    
    E <- apply(nu, 2, function(nu){
      c   <- (1i)^sum(nu) / prod(factorial(nu)) * HVE(nu = nu, K = K)
      c   <- c/prod(factorial(K))
      paste0(AE$u,"^",nu, collapse = "*") %prod% array(calculus::wrap(c))
    })
    
    if(is.null(dim(E))) E <- cpp_collapse(E, ' + ')
    else E <- apply(E, 1, function(x) cpp_collapse(x, ' + ')) 
    
    E <- array(E, dim = rep(AE$d, length(K)))
    
    return(E)
  }
  
  # Characteristic function
  psi <- function(m){
    
    martingale <- sprintf('exp(%s)', (calculus::wrap(1i*AE$Mu) %inner% AE$u) %sum% ((-0.5 * AE$Sigma) %inner% (AE$u %outer% AE$u)))
    
    if(m>0){
      psi <- cpp_collapse(paste0(AE$eps.var, "^", (1:m)) %prod% calculus::wrap(AE$P.m[1:m]), " + ")
      psi <- 1 %sum% psi
    }
    else {
      psi <- 1
    }
    
    psi <- paste0(martingale," * (", psi, ")")
    
    for(i in 1:AE$d)
      psi <- gsub(x = psi, pattern = paste0("\\b",AE$u[i],"\\b"), replacement = calculus::wrap(paste(AE$eps.var,"*",AE$u[i])))
    
    psi <- paste0("exp((",1i,") * (",AE$u %inner% AE$Ez.T[AE$z],")) * (", psi, ")")
    
    return(parse(text = psi))
    
  }
  
  # label for partitions 
  label <- function(I){
    
    paste(I, collapse = ',')
    
  }
  
  ######################################################## 
  # AE environment                                       #
  ########################################################
  AE          <- new.env()
  
  # expansion order
  AE$m        <- as.integer(order)
  
  # expansion variable
  AE$eps.var  <- eps.var
  
  # model parameters
  AE$par      <- true.parameter
  AE$par[[AE$eps.var]] <- 0
  
  # model dimension
  AE$d        <- model@model@equation.number
  
  # noise dimension
  AE$r        <- model@model@noise.number + 1
  
  # solve variables 
  AE$z        <- model@model@solve.variable
  
  # extended solve variables 
  AE$z.bar    <- c(AE$z, AE$eps.var)
  
  # characteristic function variables
  AE$u        <- paste0('u', 1:AE$d)
  
  # simulation initial value
  AE$xinit    <- xinit
  
  # timing if verbose
  if(verbose) time <- Sys.time()
  
  ######################################################## 
  # V: V[,1] -> V_{0}; V[,2] -> V_{1}                    #
  ########################################################
  AE$V <- NULL
  
  # drift 
  drift     <- calculus::e2c(model@model@drift)
  
  # diffusion 
  diffusion <- unlist(lapply(model@model@diffusion, calculus::e2c))
  diffusion <- as.array(matrix(diffusion, nrow = AE$d, ncol = AE$r-1, byrow = TRUE))
  
  # expansion coefficient
  is.eps.diffusion <- any(grepl(x = diffusion, pattern = paste0('\\b',AE$eps.var,'\\b')))
  is.eps.drift     <- any(grepl(x = drift, pattern = paste0('\\b',AE$eps.var,'\\b')))
  
  if(!is.eps.diffusion){
    
    diffusion <- AE$eps.var %prod% calculus::wrap(diffusion)
    
  } else {
    
    test <- parse(text = diffusion)

    env  <- AE$par
    env[[AE$eps.var]] <- 0
    for(z in AE$z)
      env[[z]] <- runif(n = 100, min = -999, max = 999)

    is.ok <- suppressWarnings(sapply(test, function(expr){
      all(eval(expr, env)==0, na.rm = TRUE)
    }))
    
    if(!all(is.ok)){
      stop('diffusion must vanish when evaluated at epsilon = 0')
    }
    
  }
  
  # building V...
  AE$V <- array(c(drift, diffusion), c(AE$d, AE$r))
  
  ######################################################## 
  # dV                                                   #
  ########################################################
  AE$dV       <- list()
  if(verbose) cat('Computing dV...')
  
  # parse v only once
  expr <- array(parse(text = AE$V), dim = dim(AE$V))
  
  # for j up to twice the expansion order...
  for(j in 1:(2*AE$m)) {
    
    # differentiate expression 
    expr <- calculus::derivative(expr, var = AE$z.bar, deparse = FALSE)
    
    # convert to char
    tmp <- calculus::e2c(expr)
    
    # break if dV vanishes
    if(all(tmp=="0")) break
    
    # evaluate at eps = 0
    if(!is.eps.diffusion & !is.eps.drift){
      
      # auto: eps * diffusion -> boost performance
      zero <- grepl(x = tmp, pattern = paste0('\\b',AE$eps.var,'\\b'))
      tmp[zero] <- "0"
      
    } 
    else {
      
      # user defined: gsub
      tmp[] <- gsub(x = tmp, pattern = paste0('\\b',AE$eps.var,'\\b'), replacement = "0")  
      
    }
    
    # drop white spaces (!important)
    tmp[] <- gsub(x = tmp, pattern = ' ', replacement = '', fixed = T)
    
    # store dV
    AE$dV[[j]] <- tmp
    
  }
  
  if(verbose) {
    cat(sprintf(' %s derivatives (%s sec)\n', length(unlist(AE$dV)), difftime(Sys.time(), time, units = "secs")[[1]]))
    time <- Sys.time()
  }
  
  ######################################################## 
  # dZ                                                   #
  ########################################################
  AE$dZ       <- list()
  if(verbose) cat('Computing dZ...')
  
  # for j up to twice the expnasion order...
  for(k in 1:(2*AE$m)) {
    
    # get partitions of k
    I.set <- calculus::partitions(n = k)
    
    # building dZ...
    AE$dZ[[k]] <- "0"
    
    # for each I in I.set
    lapply(I.set, function(I){
      
      # length I
      j   <- length(I)
      
      # if dV[[j]] does not vanish
      if(j <= length(AE$dV)){
        
        # compute U
        U <- (factorial(sum(I))/prod(factorial(c( I, table(I) )))) %prod% calculus::wrap(AE$dV[[j]])
        
        # isolate coefficients
        U[] <- paste0('{',U,'}')
        
        # add to dZ
        AE$dZ[[k]] <- AE$dZ[[k]] %sum% ( U %dot% Z.I(I = I, bar = TRUE) )
        
      }
      
      # void
      return(NULL)
      
    })
    
  }
  
  if(verbose) {
    cat(sprintf(' (%s sec)\n', difftime(Sys.time(), time, units = "secs")[[1]]))
    time <- Sys.time()
  }
  
  ######################################################## 
  # K.set                                                #
  ########################################################
  AE$K.set <- NULL
  
  # for n up to 4 times the expansion order...
  for(n in 1:(4*AE$m)){
    
    # store K.set
    AE$K.set <- c(AE$K.set, calculus::partitions(n = n, max = 2*AE$m))
    
  }
  
  ######################################################## 
  # Z.K                                                  #
  ########################################################
  AE$Z.K        <- lapply(AE$K.set, function(K) Z.I(I = K))
  names(AE$Z.K) <- lapply(AE$K.set, label)
  
  ########################################################
  
  if(verbose) cat('Computing Ito ODE system...')
  
  ito <- cpp_ito(K_set = AE$K.set, dZ = AE$dZ, Z_K = AE$Z.K, d = AE$d, r = AE$r)
  AE$ito.lhs     <- ito$lhs
  AE$ito.rhs     <- ito$rhs
  AE$ito.rhs.var <- ito$rhs.var
  
  if(verbose) {
    cat(sprintf(' %s equations (%s sec)\n', length(AE$ito.rhs)+AE$d, difftime(Sys.time(), time, units = "secs")[[1]]))
    time <- Sys.time()
  }
  
  ########################################################
  
  if(verbose) cat('Reducing Ito ODE system...')
  
  AE$nu     <- list()
  AE$Ez.K   <- list()
  
  # nu indices
  for(k in 1:(2*AE$m)) 
    AE$nu[[k]] <- calculus::partitions(n = k, length = AE$d, fill = TRUE, perm = TRUE, equal = FALSE)
  
  # find needed terms
  for(i in 1:AE$m) {
    for(l in 1:i){
      
      K.set <- calculus::partitions(n = i, length = l, perm = TRUE)
      
      lapply(K.set, function(K) {
        
        K <- K+1
        
        AE$Ez.K[[label(K)]] <- list()
        
        apply(AE$nu[[sum(K)]], 2, function(nu) {
          # store
          AE$Ez.K[[label(K)]][[label(nu)]] <- cpp_E( dz.nu(i = 1, nu = nu) %prod% Z.I(I = K) )
          # void
          return(NULL)
        })
        
        # void
        return(NULL)
      })
    }
  }
  
  ########################################################
  
  AE$Ez.T   <- list()
  AE$Ez     <- unique(unlist(AE$Ez.K))
  
  # drop duplicated equations
  idx <- which(!duplicated(AE$ito.lhs))
  AE$ito.lhs     <- AE$ito.lhs[idx]
  AE$ito.rhs     <- AE$ito.rhs[idx]
  AE$ito.rhs.var <- AE$ito.rhs.var[idx]
  
  while(TRUE){
    # needed equation id
    idx <- which(AE$ito.lhs %in% AE$Ez)
    # additional needed var
    add <- unique(unlist(AE$ito.rhs.var[idx]))
    add <- add[!(add %in% AE$Ez)]
    # break or store
    if(length(add)==0) break
    AE$Ez <- c(AE$Ez, add)
  }
  
  # drop equations not needed
  idx <- which(AE$ito.lhs %in% AE$Ez)
  AE$ito.lhs     <- AE$ito.lhs[idx]
  AE$ito.rhs     <- AE$ito.rhs[idx]
  AE$ito.rhs.var <- AE$ito.rhs.var[idx]
  
  # plug zero if empty equation
  while(TRUE){
    
    idx  <- which(AE$ito.rhs=='')
    if(length(idx)==0) break
    
    zero <- AE$ito.lhs[idx] 
    AE$Ez.T[zero] <- 0
    AE$ito.lhs <- AE$ito.lhs[-idx]
    AE$ito.rhs <- AE$ito.rhs[-idx]
    
    pattern <- gsub(zero, pattern = '.', replacement = '\\.', fixed = T)
    pattern <- gsub(pattern, pattern = '_', replacement = '\\_', fixed = T)
    pattern <- paste0(pattern, collapse = '|')
    pattern <- paste0(' \\* (',pattern,')\\b')
    
    AE$ito.rhs <- unlist(lapply(strsplit(AE$ito.rhs, split = ' + ', fixed = T), function(x){
      is.zero <- grepl(x = x, pattern = pattern)
      x <- x[!is.zero]
      if(length(x)==0) return("")
      else return(paste(x, collapse = ' + '))
    }))
    
  }
  
  if(verbose) {
    cat(sprintf(' %s equations (%s sec)\n', length(AE$ito.rhs)+AE$d, difftime(Sys.time(), time, units = "secs")[[1]]))
    time <- Sys.time()
  }
  
  ########################################################
  
  if(verbose) cat('Solving Ito ODE system...')
  
  # Ito ODE
  lhs <- c(AE$z, AE$ito.lhs)
  rhs <- c(AE$V[,1], AE$ito.rhs)
  
  # Initial values
  xinit <- c(rep(AE$xinit, length.out = AE$d), rep(0, length(lhs)-AE$d))
  names(xinit) <- lhs
  
  # Solve
  Ez.T <- calculus::ode(f = rhs, var = xinit, times = sampling@grid[[1]], params = AE$par, timevar = model@model@time.variable, drop = TRUE, method = solver)
  AE$Ez.T <- c(as.list(Ez.T), AE$Ez.T)
  
  if(verbose) {
    cat(sprintf(' (%s sec)\n', difftime(Sys.time(), time, units = "secs")[[1]]))
    time <- Sys.time()
  }
  
  ########################################################
  
  if(verbose) cat('Computing Sigma matrix...')
  
  y.lhs <- array(paste('y', gsub(AE$z %outer% AE$z, pattern = " * ", replacement = "_", fixed = T), sep = '.'), dim = rep(AE$d, 2))
  y.rhs <- calculus::wrap(AE$V[,1] %gradient% AE$z) %mx% y.lhs
  y.0 <- diag(AE$d)
  
  dim(y.lhs) <- NULL
  dim(y.rhs) <- NULL
  dim(y.0) <- NULL
  
  y.lhs <- c(AE$z, y.lhs)
  y.rhs <- c(AE$V[,1], y.rhs)
  
  y.0 <- c(rep(AE$xinit, length.out = AE$d), y.0)
  names(y.0) <- y.lhs
  
  times <- sampling@grid[[1]]
  n.times <- length(times)
  
  x <- calculus::ode(f = y.rhs, var = y.0, times = times, params = AE$par, timevar = model@model@time.variable, method = solver)
  y <- x[, -c(1:AE$d), drop = FALSE]
  z <- as.data.frame(x[, c(1:AE$d), drop = FALSE])
  if(!is.null(model@model@time.variable))
    z[[model@model@time.variable]] <- times

  y.s <- lapply(1:n.times, function(i) matrix(y[i,], nrow = AE$d))
  y_inv.s <- lapply(y.s, solve)
  
  dV   <- AE$V %gradient% AE$eps.var
  dV.s <- calculus::evaluate(dV, as.data.frame(c(z, AE$par)))
  dV.s <- lapply(1:n.times, function(i) array(dV.s[i,], dim = dim(dV)))
  if(length(dV.s)==1)
    dV.s = rep(dV.s, n.times)
  
  # Mu
  if(all(dV[,1,]=="0")){
    
    AE$Mu <- array(0, dim = AE$d)
    
  }
  else {
    
    Mu <- sapply(1:n.times, function(i){
      y_inv.s[[i]] %*% dV.s[[i]][,1,]
    })
    
    if(is.null(dim(Mu))) {
      n <- length(Mu)
      Mu <- ( sum(Mu[c(-1,-n)]) + sum(Mu[c(1,n)])/2 ) * sampling@delta
    }
    else {
      n <- ncol(Mu)
      Mu <- ( rowSums(Mu[,c(-1,-n)]) + rowSums(Mu[,c(1,n)])/2 ) * sampling@delta
    }
    
    AE$Mu <- array(y.s[[n.times]] %*% Mu)
    
  }
  
  # Sigma
  S <- sapply(1:n.times, function(i){
    a <- y.s[[n.times]] %*% y_inv.s[[i]] %*% dV.s[[i]][,-1,]
    a %*% t(a)    
  })
  
  if(is.null(dim(S))) {
    n <- length(S)
    S <- ( sum(S[c(-1,-n)]) + sum(S[c(1,n)])/2 ) * sampling@delta
  }
  else {
    n <- ncol(S)
    S <- ( rowSums(S[,c(-1,-n)]) + rowSums(S[,c(1,n)])/2 ) * sampling@delta
  }
  
  AE$Sigma <- array(S, dim = c(AE$d,AE$d))
  
  if(verbose) {
    cat(sprintf(' (%s sec)\n', difftime(Sys.time(), time, units = "secs")[[1]]))
    time <- Sys.time()
  }
  
  ########################################################
  
  # Hermite coefficients  
  if(verbose) cat(paste0('Computing Hermite...'))
  
  tmp <- calculus::hermite(var = AE$z, sigma = AE$Sigma, order = 2*AE$m, 
                 transform = solve(AE$Sigma) %dot% calculus::wrap(AE$z %sum% -AE$Mu))
  
  AE$c.nu <- lapply(tmp, function(x) {
    coef <- as.list(x$terms$coef)
    names(coef) <- rownames(x$terms)
    return(coef)
  })
  AE$h.nu <- lapply(tmp, function(x) x$f)
  
  if(verbose) {
    cat(sprintf(' (%s sec)\n', difftime(Sys.time(), time, units = "secs")[[1]]))
    time <- Sys.time()
  }
  
  ########################################################
  
  AE$ul <- list(array(AE$u))
  if(AE$m > 1) for(l in 2:AE$m){
    AE$ul[[l]] <- AE$ul[[l-1]] %outer% AE$u
  }
  AE$ul <- lapply(AE$ul, function(ul){
    array(unlist(lapply(strsplit(ul, split = " * ", fixed = T), function(u) {
      u <- table(u)
      paste0(names(u),"^",u, collapse = "*")
    })), dim = dim(ul))
  })
  
  ########################################################
  
  if(verbose) cat('Computing characteristic function...')
  
  AE$psi <- list()
  for(m in 1:AE$m) {
    
    AE$psi[[m]] <- list()
    
    for(l in 1:m){
      
      K.set <- calculus::partitions(n = m, length = l, perm = TRUE)
      
      psi.m.l <- unlist(lapply(K.set, function(K){
        calculus::wrap((1i)^l) %prod% calculus::wrap((calculus::wrap(TVE(K = K)) %inner% AE$ul[[l]]))
      }))
      
      expr <- (1/factorial(l)) %prod% calculus::wrap(cpp_collapse(psi.m.l, ' + '))
      AE$psi[[m]][[l]] <- calculus::taylor(expr, var = AE$u, order = m+2*l)$f
      
    }
    
  }
  
  AE$P.m = sapply(AE$psi, function(p.m.l) cpp_collapse(unlist(p.m.l), " + "))
  AE$c.gamma <- lapply(1:AE$m, function(m) {
    p <- calculus::taylor(AE$P.m[m], var = AE$u, order = 3*m)
    coef <- Re(p$terms$coef/(1i)^p$terms$degree)
    coef <- as.list(coef)
    names(coef) <- rownames(p$terms)
    return(coef)
  })
  
  AE$PSI <- sapply(0:AE$m, psi)
  
  if(verbose) {
    cat(sprintf(' (%s sec)\n', difftime(Sys.time(), time, units = "secs")[[1]]))
    time <- Sys.time()
  }
  
  ########################################################
  
  return(new(
    "yuima.ae",
    order = AE$m,
    var = AE$z,
    u.var = AE$u,
    eps.var = AE$eps.var,
    characteristic = AE$PSI,
    Z0 = unlist(AE$Ez.T[AE$z]),
    Mu = as.numeric(AE$Mu),
    Sigma = AE$Sigma,
    c.gamma = AE$c.gamma,
    verbose = verbose
  ))
  
}

#' Asymptotic Expansion - Density
#' 
#' @param ... named argument, data.frame, list, or environment specifying the grid to evaluate the density. See examples.
#' @param ae an object of class \code{\link{yuima.ae-class}}.
#' @param eps numeric. The intensity of the perturbation.
#' @param order integer. The expansion order. If \code{NULL} (default), it uses the maximum order used in \code{ae}.
#' 
#' @return Probability density function evaluated on the given grid.
#' 
#' @examples
#' \dontrun{
#' # model
#' gbm <- setModel(drift = 'mu*x', diffusion = 'sigma*x', solve.variable = 'x')
#' 
#' # settings
#' xinit <- 100
#' par <- list(mu = 0.01, sigma = 0.2)
#' sampling <- setSampling(Initial = 0, Terminal = 1, n = 1000)
#' 
#' # asymptotic expansion
#' approx <- ae(model = gbm, sampling = sampling, order = 4, true.parameter = par, xinit = xinit)
#' 
#' # The following are all equivalent methods to specify the grid via ....
#' # Notice that the character 'x' corresponds to the solve.variable of the yuima model.
#' 
#' # 1) named argument  
#' x <- seq(50, 200, by = 0.1)
#' density <- aeDensity(x = x, ae = approx, order = 4)
#' # 2) data frame
#' df <- data.frame(x = seq(50, 200, by = 0.1))
#' density <- aeDensity(df, ae = approx, order = 4)
#' # 3) environment
#' env <- new.env()
#' env$x <- seq(50, 200, by = 0.1)
#' density <- aeDensity(env, ae = approx, order = 4)
#' # 4) list
#' lst <- list(x = seq(50, 200, by = 0.1))
#' density <- aeDensity(lst, ae = approx, order = 4)
#' 
#' # exact density
#' exact <- dlnorm(x = x, meanlog = log(xinit)+(par$mu-0.5*par$sigma^2)*1, sdlog = par$sigma*sqrt(1))
#' 
#' # compare
#' plot(x = exact, y = density, xlab = "Exact", ylab = "Approximated")
#' }
#' @export
#' 
aeDensity <- function(..., ae, eps = 1, order = NULL){
  
  if(is.null(order)) 
    order <- ae@order
  
  pdf <- ae@density[order+1]
  
  env <- list(...)
  if(length(env)==1) 
    if(is.list(env[[1]]) | is.environment(env[[1]]))
      env <- env[[1]]
  
  env[[ae@eps.var]] <- eps
  
  return(eval(expr = pdf, envir = env))
  
}

#' Asymptotic Expansion - Marginals
#' 
#' @param ae an object of class \code{\link{yuima.ae-class}}.
#' @param var variables of the marginal distribution to compute.
#' 
#' @return An object of \code{\link{yuima.ae-class}}
#' 
#' @examples
#' \dontrun{
#' # multidimensional model
#' gbm <- setModel(drift = c('mu*x1','mu*x2'), 
#'                 diffusion = matrix(c('sigma1*x1',0,0,'sigma2*x2'), nrow = 2), 
#'                 solve.variable = c('x1','x2'))
#' 
#' # settings
#' xinit <- c(100, 100)
#' par <- list(mu = 0.01, sigma1 = 0.2, sigma2 = 0.1)
#' sampling <- setSampling(Initial = 0, Terminal = 1, n = 1000)
#' 
#' # asymptotic expansion
#' approx <- ae(model = gbm, sampling = sampling, order = 3, true.parameter = par, xinit = xinit)
#' 
#' # extract marginals
#' margin1 <- aeMarginal(ae = approx, var = "x1")
#' margin2 <- aeMarginal(ae = approx, var = "x2")
#' 
#' # compare with exact solution for marginal 1
#' x1 <- seq(50, 200, by = 0.1)
#' exact <- dlnorm(x = x1, meanlog = log(xinit[1])+(par$mu-0.5*par$sigma1^2), sdlog = par$sigma1)
#' plot(x1, exact, type = 'p', ylab = "Density")
#' lines(x1, aeDensity(x1 = x1, ae = margin1, order = 3), col = 2)
#' 
#' # compare with exact solution for marginal 2
#' x2 <- seq(50, 200, by = 0.1)
#' exact <- dlnorm(x = x2, meanlog = log(xinit[2])+(par$mu-0.5*par$sigma2^2), sdlog = par$sigma2)
#' plot(x2, exact, type = 'p', ylab = "Density")
#' lines(x2, aeDensity(x2 = x2, ae = margin2, order = 3), col = 2)
#' 
#' }
#' @export
#' 
aeMarginal <- function(ae, var){
  
  # init
  keep <- ae@var %in% var
  if(sum(!keep)==0) return(ae)
  
  # vanish marginal coefficients
  c.gamma <- ae@c.gamma
  for(i in 1:length(c.gamma)) for(j in 1:length(c.gamma[[i]])){
    o <- as.numeric(unlist(strsplit(x = names(c.gamma[[i]][j]), split = ',')))
    if(any(o[!keep]>0)){
      c.gamma[[i]][[j]] <- 0
    }
    else {
      names(c.gamma[[i]])[j] <- paste(o[keep], collapse = ',')
    }
  }

  # characteristic
  characteristic <- sapply(calculus::e2c(ae@characteristic), function(psi) {
    for(u in ae@u.var[!keep]) 
      psi <- gsub(x = psi, pattern = paste0("\\b",u,"\\b"), replacement = 0)
    return(calculus::c2e(psi))
  })  
  
  # return
  return(new(
    "yuima.ae",
    order = ae@order,
    var = ae@var[keep],
    u.var = ae@u.var[keep],
    eps.var = ae@eps.var,
    characteristic = characteristic,
    Z0 = ae@Z0[keep],
    Mu = ae@Mu[keep],
    Sigma = ae@Sigma[keep, keep, drop = FALSE],
    c.gamma = c.gamma,
    verbose = FALSE
  ))
  
}

#' Asymptotic Expansion - Functionals
#' 
#' Compute the expected value of functionals. 
#' 
#' @param f character. The functional.
#' @param bounds named list of integration bounds in the form \code{list(x = c(xmin, xmax), y = c(ymin, ymax), ...)} 
#' @param ae an object of class \code{\link{yuima.ae-class}}.
#' @param eps numeric. The intensity of the perturbation.
#' @param order integer. The expansion order. If \code{NULL} (default), it uses the maximum order used in \code{ae}.
#' @param ... additional arguments passed to \code{\link[cubature]{cubintegrate}}.
#' 
#' @return return value of \code{\link[cubature]{cubintegrate}}. The expectation of the functional provided. 
#' 
#' @examples 
#' \dontrun{
#' # model
#' gbm <- setModel(drift = 'mu*x', diffusion = 'sigma*x', solve.variable = 'x')
#' 
#' # settings
#' xinit <- 100
#' par <- list(mu = 0.01, sigma = 0.2)
#' sampling <- setSampling(Initial = 0, Terminal = 1, n = 1000)
#' 
#' # asymptotic expansion
#' approx <- ae(model = gbm, sampling = sampling, order = 4, true.parameter = par, xinit = xinit)
#' 
#' # compute the mean via integration
#' aeExpectation(f = 'x', bounds = list(x = c(0,1000)), ae = approx)
#' 
#' # compare with the mean computed by differentiation of the characteristic function
#' aeMean(approx)
#' }
#' @export
#' 
aeExpectation <- function(f, bounds, ae, eps = 1, order = NULL, ...){
  
  f      <- calculus::c2e(f)
  var    <- names(bounds)
  ae     <- aeMarginal(ae = ae, var = var)
  lower  <- sapply(bounds, function(x) x[1])
  upper  <- sapply(bounds, function(x) x[2])
  
  args <- list(...)
  args$f <- function(x) {
    names(x) <- var
    x <- as.list(x)
    eval(f, envir = x) * aeDensity(x, ae = ae, eps = eps, order = order)
  }
  
  args$lower <- lower
  args$upper <- upper
  
  if(is.null(args$method))
    args$method <- 'hcubature'
  
  return(do.call("cubintegrate", args))
  
}

#' Asymptotic Expansion - Characteristic Function
#' 
#' @param ... named argument, data.frame, list, or environment specifying the grid to evaluate the characteristic function. See examples.
#' @param ae an object of class \code{\link{yuima.ae-class}}.
#' @param eps numeric. The intensity of the perturbation.
#' @param order integer. The expansion order. If \code{NULL} (default), it uses the maximum order used in \code{ae}.
#' 
#' @return Characteristic function evaluated on the given grid.
#' 
#' @examples 
#' \dontrun{
#' # model
#' gbm <- setModel(drift = 'mu*x', diffusion = 'sigma*x', solve.variable = 'x')
#' 
#' # settings
#' xinit <- 100
#' par <- list(mu = 0.01, sigma = 0.2)
#' sampling <- setSampling(Initial = 0, Terminal = 1, n = 1000)
#' 
#' # asymptotic expansion
#' approx <- ae(model = gbm, sampling = sampling, order = 4, true.parameter = par, xinit = xinit)
#' 
#' # The following are all equivalent methods to specify the grid via ....
#' # Notice that the character 'u1' corresponds to the 'u.var' of the ae object.
#' approx@u.var
#' 
#' # 1) named argument  
#' u1 <- seq(0, 1, by = 0.1)
#' psi <- aeCharacteristic(u1 = u1, ae = approx, order = 4)
#' # 2) data frame
#' df <- data.frame(u1 = seq(0, 1, by = 0.1))
#' psi <- aeCharacteristic(df, ae = approx, order = 4)
#' # 3) environment
#' env <- new.env()
#' env$u1 <- seq(0, 1, by = 0.1)
#' psi <- aeCharacteristic(env, ae = approx, order = 4)
#' # 4) list
#' lst <- list(u1 = seq(0, 1, by = 0.1))
#' psi <- aeCharacteristic(lst, ae = approx, order = 4)
#' }
#' @export
#' 
aeCharacteristic <- function(..., ae, eps = 1, order = NULL){
  
  if(is.null(order)) 
    order <- ae@order
  
  psi <- ae@characteristic[order+1]
  
  env <- list(...)
  if(length(env)==1) 
    if(is.list(env[[1]]) | is.environment(env[[1]]))
      env <- env[[1]]
  
  env[[ae@eps.var]] <- eps
  
  return(eval(expr = psi, envir = env))
  
}

#' Asymptotic Expansion - Moments
#' 
#' @param ae an object of class \code{\link{yuima.ae-class}}.
#' @param m integer. The moment order. In case of multidimensional processes, it is possible to compute cross-moments by providing a vector of the same length as the state variables.
#' @param eps numeric. The intensity of the perturbation.
#' @param order integer. The expansion order. If \code{NULL} (default), it uses the maximum order used in \code{ae}.
#' 
#' @return numeric.
#' 
#' @examples 
#' \dontrun{
#' # model
#' gbm <- setModel(drift = 'mu*x', diffusion = 'sigma*x', solve.variable = 'x')
#' 
#' # settings
#' xinit <- 100
#' par <- list(mu = 0.01, sigma = 0.2)
#' sampling <- setSampling(Initial = 0, Terminal = 1, n = 1000)
#' 
#' # asymptotic expansion
#' approx <- ae(model = gbm, sampling = sampling, order = 4, true.parameter = par, xinit = xinit)
#' 
#' # second moment, expansion order max
#' aeMoment(ae = approx, m = 2)
#' 
#' # second moment, expansion order 3
#' aeMoment(ae = approx, m = 2, order = 3)
#' 
#' # second moment, expansion order 2
#' aeMoment(ae = approx, m = 2, order = 2)
#' 
#' # second moment, expansion order 1
#' aeMoment(ae = approx, m = 2, order = 1)
#' }
#' @export
#' 
aeMoment <- function(ae, m = 1, eps = 1, order = NULL){
  
  if(is.null(order)) 
    order <- ae@order
  
  psi <- ae@characteristic[order+1]
  
  env <- c()
  env[ae@u.var] <- 0
  env[ae@eps.var] <- eps
  
  return(Re((-1i)^sum(m) * calculus::evaluate(calculus::derivative(psi, var = ae@u.var, order = m, deparse = FALSE), env)))
  
}

#' Asymptotic Expansion - Mean
#' 
#' @param ae an object of class \code{\link{yuima.ae-class}}.
#' @param eps numeric. The intensity of the perturbation.
#' @param order integer. The expansion order. If \code{NULL} (default), it uses the maximum order used in \code{ae}.
#' 
#' @return numeric.
#' 
#' @examples 
#' \dontrun{
#' # model
#' gbm <- setModel(drift = 'mu*x', diffusion = 'sigma*x', solve.variable = 'x')
#' 
#' # settings
#' xinit <- 100
#' par <- list(mu = 0.01, sigma = 0.2)
#' sampling <- setSampling(Initial = 0, Terminal = 1, n = 1000)
#' 
#' # asymptotic expansion
#' approx <- ae(model = gbm, sampling = sampling, order = 4, true.parameter = par, xinit = xinit)
#' 
#' # expansion order max
#' aeMean(ae = approx)
#' 
#' # expansion order 1
#' aeMean(ae = approx, order = 1)
#' }
#' @export
#' 
aeMean <- function(ae, eps = 1, order = NULL){
  
  return(aeMoment(ae = ae, m = 1, eps = eps, order = order))
  
}

#' Asymptotic Expansion - Standard Deviation
#' 
#' @param ae an object of class \code{\link{yuima.ae-class}}.
#' @param eps numeric. The intensity of the perturbation.
#' @param order integer. The expansion order. If \code{NULL} (default), it uses the maximum order used in \code{ae}.
#' 
#' @return numeric.
#' 
#' @examples 
#' \dontrun{
#' # model
#' gbm <- setModel(drift = 'mu*x', diffusion = 'sigma*x', solve.variable = 'x')
#' 
#' # settings
#' xinit <- 100
#' par <- list(mu = 0.01, sigma = 0.2)
#' sampling <- setSampling(Initial = 0, Terminal = 1, n = 1000)
#' 
#' # asymptotic expansion
#' approx <- ae(model = gbm, sampling = sampling, order = 4, true.parameter = par, xinit = xinit)
#' 
#' # expansion order max
#' aeSd(ae = approx)
#' 
#' # expansion order 1
#' aeSd(ae = approx, order = 1)
#' }
#' @export
#' 
aeSd <- function(ae, eps = 1, order = NULL){
  
  m1 <- aeMoment(ae = ae, m = 1, eps = eps, order = order)
  m2 <- aeMoment(ae = ae, m = 2, eps = eps, order = order)
  
  return(sqrt(m2-m1^2))
  
}

#' Asymptotic Expansion - Skewness
#' 
#' @param ae an object of class \code{\link{yuima.ae-class}}.
#' @param eps numeric. The intensity of the perturbation.
#' @param order integer. The expansion order. If \code{NULL} (default), it uses the maximum order used in \code{ae}.
#' 
#' @return numeric.
#' 
#' @examples 
#' \dontrun{
#' # model
#' gbm <- setModel(drift = 'mu*x', diffusion = 'sigma*x', solve.variable = 'x')
#' 
#' # settings
#' xinit <- 100
#' par <- list(mu = 0.01, sigma = 0.2)
#' sampling <- setSampling(Initial = 0, Terminal = 1, n = 1000)
#' 
#' # asymptotic expansion
#' approx <- ae(model = gbm, sampling = sampling, order = 4, true.parameter = par, xinit = xinit)
#' 
#' # expansion order max
#' aeSkewness(ae = approx)
#' 
#' # expansion order 1
#' aeSkewness(ae = approx, order = 1)
#' }
#' @export
#' 
aeSkewness <- function(ae, eps = 1, order = NULL){
  
  m1 <- aeMoment(ae = ae, m = 1, eps = eps, order = order)
  m2 <- aeMoment(ae = ae, m = 2, eps = eps, order = order)
  m3 <- aeMoment(ae = ae, m = 3, eps = eps, order = order)
  
  return((m3-3*m1*m2+2*m1^3)/sqrt(m2-m1^2)^3)
  
}

#' Asymptotic Expansion - Kurtosis
#' 
#' @param ae an object of class \code{\link{yuima.ae-class}}.
#' @param eps numeric. The intensity of the perturbation.
#' @param order integer. The expansion order. If \code{NULL} (default), it uses the maximum order used in \code{ae}.
#' 
#' @return numeric.
#' 
#' @examples
#' \dontrun{
#' # model
#' gbm <- setModel(drift = 'mu*x', diffusion = 'sigma*x', solve.variable = 'x')
#' 
#' # settings
#' xinit <- 100
#' par <- list(mu = 0.01, sigma = 0.2)
#' sampling <- setSampling(Initial = 0, Terminal = 1, n = 1000)
#' 
#' # asymptotic expansion
#' approx <- ae(model = gbm, sampling = sampling, order = 4, true.parameter = par, xinit = xinit)
#' 
#' # expansion order max
#' aeKurtosis(ae = approx)
#' 
#' # expansion order 1
#' aeKurtosis(ae = approx, order = 1)
#' }
#' @export
#' 
aeKurtosis <- function(ae, eps = 1, order = NULL){
  
  m1 <- aeMoment(ae = ae, m = 1, eps = eps, order = order)
  m2 <- aeMoment(ae = ae, m = 2, eps = eps, order = order)
  m3 <- aeMoment(ae = ae, m = 3, eps = eps, order = order)
  m4 <- aeMoment(ae = ae, m = 4, eps = eps, order = order)
  
  return((m4-4*m1*m3+6*m1^2*m2-3*m1^4)/sqrt(m2-m1^2)^4)
  
}

