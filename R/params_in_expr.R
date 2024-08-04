# check if given params are included in given expression. return vector of logical
params_in_expr <- function(params, expr) {
  return(is.element(params, all.names(expr)))
}

# check if given params are included in any element of given list of expressions. return vector of logical
params_in_exprs <- function(params, exprs) {
  is.element.res <- rep(0, length(params))
  for (i in 1:length(exprs)) {
    is.element.res <- is.element.res + as.numeric(params_in_expr(params, exprs[[i]]))
  }
  res <- is.element.res > 0
  return(res)
}
