library(distrEx)
library(parallel)

#' Marginal survival function for exponential
#' (possibly conditionally exponential)
#' random variable.
#' 
#' @export
#' @param t A number or vector, the survival times
#' @param a A number, the shape parameter of Gamma frailty distribution
#' @param b A number, the rate parameter of Gamma frailty distribution
#' @param parameters A list, with an element named lambda, the rate
#' parameter of the exponential distribution
#' @param correlated If TRUE, random variable is conditionally exponential
#' on the frailty parameter, that has Gamma(shape=a, rate=b) distribution.
#' @return Survival function evaluated at t
#' @examples
#' S.Y.Exponential(t=0, a=1, b=1, parameters=list(lambda=1), correlated=FALSE)
#' S.Y.Exponential(t=c(1, 2, 3), a=1, b=1, parameters=list(lambda=1), correlated=FALSE)
S.Y.Exponential <- function(
  t,
  a,
  b,
  parameters,
  correlated){
  if(correlated){
    return((b / (b + parameters$lambda * t))^a)
  } else {
    return(exp(-t * parameters$lambda))
  }
}
