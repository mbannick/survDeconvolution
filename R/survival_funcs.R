library(distrEx)
library(parallel)
library(stats)

#' Marginal survival function for exponential
#' (possibly conditionally exponential)
#' (possibly conditionally exponential)
#' random variable.
#' 
#' @export
#' @importFrom stats dweibull
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

#' Marginal survival function for the sum
#' of two exponential, potentially correlated random variables.
#' 
#' \eqn{X} and \eqn{Y}.
#' If \code{correlated}, then the distributions are
#' correlated through a shared frailty term, \code{k},
#' where \eqn{k \sim Gamma(a, b)}.
#' 
#' @export
#' @param t A number or vector, the survival times
#' @param a A number, the shape parameter of Gamma frailty distribution
#' @param b A number, the rate parameter of Gamma frailty distribution
#' @param parameters_x A list, with an element named lambda, the rate
#' parameter of the exponential distribution
#' @param parameters_y A list, with an element named lambda, the rate
#' parameter of the exponential distribution
#' @param correlated If TRUE, random variables that make up the sum
#' are each conditionally exponential on the frailty parameter,
#' that has Gamma(shape=a, rate=b) distribution.
#' @return Survival function evaluated at t
#' @examples
#' S.Z.Exponential(t=2, a=1, b=1, parameters_x=list(lambda=1),
#' parameters_y=list(lambda=2), correlated=TRUE)
S.Z.Exponential <- function(
  t,
  a,
  b,
  parameters_x,
  parameters_y,
  correlated){
  if(correlated){
    return(
      b^(a) / (parameters_x$lambda - parameters_y$lambda) *
        (parameters_x$lambda/(b + t * parameters_y$lambda)^(a) -
           parameters_y$lambda/(b + t * parameters_x$lambda)^(a))
    )
  } else {
    return(
      (parameters_x$lambda * exp(-(parameters_y$lambda * t)) -
         parameters_y$lambda * exp(-(parameters_x$lambda * t))) /
        (parameters_x$lambda - parameters_y$lambda)
    )
  }
}

#' Marginal density function for a Weibull
#' with scale parameter a function of the scalar k.
#' 
#' @export 
#' @param x A vector
#' @param k The scalar that influences the scale
#' parameter for the Weibull. Consider a Weibull
#' with rate parameter \eqn{p} and \eqn{\lambda}.
#' If the distribution is re-parameterized such that
#' the rate parameter is \eqn{\lambda \cdot k^{p}},
#' then \eqn{k} proportionally scales the hazard function
#' for the distribution.
#' @param parameters A list of the two parameters
#' of the Weibull distribution \eqn{p} and \eqn{\lambda}
#' before the re-parameterization described for \eqn{k}.
#' @return The density function evaluated at x
#' @examples
#' f.Weibull(x=c(0, 0.5, 1, 5), k=0.5, parameters=list(p=1, lambda=0.5))
#' f.Weibull(x=c(0, 0.5, 1, 5), k=1, parameters=list(p=1, lambda=0.5))
#' f.Weibull(x=c(0, 0.5, 1, 5), k=4, parameters=list(p=1, lambda=0.5))
f.Weibull <- function(
  x,
  k,
  parameters){
  dweibull(
    x,
    shape=parameters$p,
    scale=1/(parameters$lambda*k^(1/parameters$p)))
}
