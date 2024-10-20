# check if given expression is linear with given variables
is.linear <- function(drift, variables) {
  variables.num <- length(variables)

  # try symbolic differentiation to check linearity
  hessians <- try(derivative(f = derivative(f = drift, var = variables), var = variables), silent = T)
  if (!inherits(hessians, "try-error")) {
    return(all(hessians == "0"))
  }

  # if symbolic differentiation filed,raise warning and return TRUE
  yuima.warn("Symbolic differentiation of drift term failed. Could not validate that drift term is linear with given variables.")
  return(TRUE)
}
