# partially evaluate given expression. returns expression.
partial.eval <- function(expr, env) {
  partial_eval_call <- function(call, env) {
    do.call("substitute", args = list(call, env = env))
  }
  as.expression(lapply(expr, partial_eval_call, env = env))
}
