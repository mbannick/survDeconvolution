library(distrEx)
library(parallel)
library(stats)

#' Marginal survival function for exponential
#' random variable \eqn{Y}, with frailty parameter \eqn{K}.
#' 
#' If \code{!correlated}: \eqn{Y ~ Exponential(\lambda)}
#' If \code{correlated}: \eqn{Y ~ Exponential(\lambda \cdot k) \quad K ~ Gamma(a, b)}
#' 
#' @importFrom stats dweibull
#' @param t A number or vector, the survival times
#' @param a A number, the shape parameter of Gamma frailty distribution
#' @param b A number, the rate parameter of Gamma frailty distribution
#' @param parameters A list, with an element named lambda, the rate
#' parameter of the exponential distribution
#' @param correlated If TRUE, random variable is conditionally exponential
#' on the frailty parameter, that has \code{Gamma(shape=a, rate=b)} distribution.
#' @return Survival function evaluated at t
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
#' of two exponential random variables, \eqn{Z=X+Y}, with shared frailty parameter \eqn{K}.
#' 
#' If \code{!correlated}: \eqn{X \sim Exponential(\lambda_{x}) \quad Y \sim Exponential(\lambda_{y})} \cr
#' If \code{correlated}: \eqn{X \sim Exponential(\lambda_{x} \cdot k) \quad Y \sim Exponential(\lambda_{y} \cdot k) \quad K \sim Gamma(a, b)}
#' 
#' @param t A number or vector, the survival times
#' @param a A number, the shape parameter of Gamma frailty distribution
#' @param b A number, the rate parameter of Gamma frailty distribution
#' @param parameters_x A list, with an element named \code{lambda} (\eqn{\lambda_x}), the rate
#' parameter of the exponential distribution
#' @param parameters_y A list, with an element \code{lambda} (\eqn{\lambda_y}), the rate
#' parameter of the exponential distribution
#' @param correlated If TRUE, random variables that make up the sum
#' are each conditionally exponential on the frailty parameter,
#' that has Gamma(shape=a, rate=b) distribution.
#' @return Survival function evaluated at t
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

#' Conditional density function for Weibull random variable \eqn{X} with
#' scalar multiplying the hazard \eqn{k}. 
#' 
#' \eqn{X ~ Weibull(shape=p, rate=\lambda \cdot k^{p})}
#' 
#' @param x A vector
#' @param k The scalar that influences the scale
#' parameter for the Weibull. Consider a Weibull
#' with rate parameter \eqn{p} and \eqn{\lambda}.
#' If the distribution is re-parameterized such that
#' the rate parameter is \eqn{\lambda \cdot k^{p}},
#' then \eqn{k} proportionally scales the hazard function
#' for the distribution.
#' @param parameters A list of the two parameters
#' of the Weibull distribution p and lambda
#' before the re-parameterization described for k.
#' @return The density function evaluated at x
f.Weibull <- function(
  x,
  k,
  parameters){
  dweibull(
    x,
    shape=parameters$p,
    scale=1/(parameters$lambda*k^(1/parameters$p)))
}

#' Marginal density function for conditional Weibull random variable \eqn{Y|K}
#' marginalized over \eqn{K}.
#' 
#' If \code{!correlated}: \eqn{Y \sim Weibull(shape=p, rate=\lambda)} \cr
#' If \code{correlated}: \eqn{Y \sim Weibull(shape=p, rate=\lambda \cdot k^{p}) \quad K \sim Gamma(a, b)}
#' 
#' @importFrom stats dweibull
#' @param y A vector at which to evaluate the density function
#' @param parameters A list of the two parameters
#' of the Weibull distribution p and lambda
#' before the re-parameterization described for \code{k} in \code{f.Weibull}.
#' @param gamma_shape The shape parameter for gamma distribution
#' where in \code{f.Weibull}, \eqn{k ~ Gamma(gamma_shape, gamma_rate)}.
#' @param gamma_rate The shape parameter for gamma distribution
#' where in \code{f.Weibull}, \eqn{k ~ Gamma(gamma_shape, gamma_rate)}.
#' @param correlated If \code{FALSE}, then this is just a regular Weibull
#' density that you could get with \code{dweibull}. If \code{TRUE}, then include
#' the frailty parameterization marginalizing over \code{k} in \code{f.Weibull}.
f.Y.Weibull <- function(
  y,
  parameters,
  gamma_shape,
  gamma_rate,
  correlated){
  lambda <- parameters$lambda
  p <- parameters$p
  a <- gamma_shape
  b <- gamma_rate
  if(correlated){
    result <- sapply(y, function(y_i)
      lambda^p * p * y_i^(p-1) * b^a * a /
        (b + (lambda * y_i)^p)^(a+1))
  } else {
    result <- dweibull(y, shape=parameters$p,
                       scale=1/parameters$lambda)
  }
  return(result)
}

#' Joint density function for two conditionally Weibull
#' random variables with shared frailty \eqn{k}.
#' 
#' \eqn{X \sim Weibull(shape=p_{x}, rate=\lambda_{x} \cdot k^{p_{x}})} \cr
#' \eqn{Y \sim Weibull(shape=p_{y}, rate=\lambda_{y} \cdot k^{p_{y}})} \cr
#' \eqn{Z := X + Y \implies X = Z - Y} \cr
#' \eqn{K \sim Gamma(shape=gamma_shape, rate=gamma_rate)}
#' 
#' @importFrom stats dgamma
#' @param k A number or vector, the frailty term, \eqn{K = k}
#' @param y A number, value of \eqn{Y = y}
#' @param z A number, value of \eqn{Z = z}
#' @param gamma_shape A number, the shape parameter of Gamma frailty distribution
#' @param gamma_rate A number, the rate parameter of Gamma frailty distribution
#' @param parameters_x A list, with an element named \code{lambda} (\eqn{\lambda_x}), and an element named \code{p} (\eqn{p_x}).
#' @param parameters_y A list, with an element \code{lambda} (\eqn{\lambda_y}), and an element named \code{p} (\eqn{p_y}).
#' @return joint density function evaluated at \code{K=k, Y=y, and X=z-y}.
f.YZK.Weibull <- function(
  k,
  y,
  z,
  parameters_x,
  parameters_y,
  gamma_shape,
  gamma_rate){
  f.Weibull(y, k, parameters_y) *
    f.Weibull(z-y, k, parameters_x) *
    dgamma(k, shape=gamma_shape, rate=gamma_rate)
}

#' Joint density function for two conditionally Weibull
#' random variables, marginalized over the frailty term.
#' 
#' If \code{correlated}, the function computes the marginal joint density
#' of Y and Z by marginalizing over the support of K. For computational efficiency,
#' it integrates from \code{epsilon} to \code{int.upper}, where \code{epsilon} can be 0
#' and \code{int.upper} can be \code{Inf}. Small deviations from 0 and \code{Inf} will
#' approximate the true integral.
#' 
#' \eqn{X \sim Weibull(shape=p_{x}, rate=\lambda_{x} \cdot k^{p_{x}})} \cr
#' \eqn{Y \sim Weibull(shape=p_{y}, rate=\lambda_{y} \cdot k^{p_{y}})} \cr
#' \eqn{Z := X + Y \implies X = Z - Y} \cr
#' \eqn{K \sim Gamma(shape=gamma_shape, rate=gamma_rate)}
#' 
#' @importFrom stats dweibull
#' @importFrom distrEx distrExIntegrate
#' @param y A number at which to evaluate the density \eqn{Y = y}
#' @param z A number at which to evaluate the density \eqn{X = z-y}
#' @param gamma_shape A number, the shape parameter of Gamma frailty distribution
#' @param gamma_rate A number, the rate parameter of Gamma frailty distribution
#' @param parameters_x A list, with an element named \code{lambda} (\eqn{\lambda_x}), and an element named \code{p} (\eqn{p_x}).
#' @param parameters_y A list, with an element \code{lambda} (\eqn{\lambda_y}), and an element named \code{p} (\eqn{p_y}).
#' @param epsilon A float, the lower bound for support of K
#' @param int.upper A float, the upper bound for support of K
#' @param control_list A list, with an element named \code{rel.tol} to be passed to `distrExIntegrate` as the \code{rel.tol} argument (see \code{?distrExIntegrate})
#' @return joint density function evaluated at \code{Y=y, and X=z-y}.
#' @param correlated If \code{FALSE}, then \eqn{X} and \eqn{Y} are independent. If \code{TRUE}, then
#' there is a shared frailty term \eqn{K} that induces correlation between \eqn{X} and \eqn{Y}.
f.YZ.Weibull <- function(
  y,
  z,
  parameters_x,
  parameters_y,
  gamma_shape,
  gamma_rate,
  epsilon,
  control_list,
  int.upper, correlated=TRUE){
  if(correlated){
    result <- mapply(function(y_i, z_i) distrExIntegrate(
      f.YZK.Weibull, lower=epsilon, upper=int.upper,
      y=y_i, z=z_i,
      parameters_x=parameters_x,
      parameters_y=parameters_y,
      gamma_shape=gamma_shape,
      gamma_rate=gamma_rate,
      subdivisions=control_list$subdivisions,
      rel.tol=control_list$rel.tol), y, z)
  } else {
    result <- dweibull(z-y, shape=parameters_x$p, scale=1/parameters_x$lambda) *
      dweibull(y, shape=parameters_y$p, scale=1/parameters_y$lambda)
  }
  return(result)
}

#' Density function for the sum of two conditionally Weibull random variables
#' with a shared frailty term.
#' 
#' If \code{correlated}, the function computes the density of Z by
#' first marginalizing over the support of K, then over Y from 0 to z. For computational efficiency,
#' it integrates K from \code{epsilon} to \code{int.upper} and Y from \code{epsilon} to z, where \code{epsilon} can be 0
#' and \code{int.upper} can be \code{Inf}. Small deviations from 0 and \code{Inf} will
#' approximate the true integral.
#' 
#' \eqn{X \sim Weibull(shape=p_{x}, rate=\lambda_{x} \cdot k^{p_{x}})} \cr
#' \eqn{Y \sim Weibull(shape=p_{y}, rate=\lambda_{y} \cdot k^{p_{y}})} \cr
#' \eqn{Z := X + Y} \cr
#' \eqn{K \sim Gamma(shape=gamma_shape, rate=gamma_rate)}
#' 
#' @importFrom stats dweibull
#' @importFrom distrEx distrExIntegrate
#' @param z A number or vector at which to evaluate the density function
#' @param gamma_shape A number, the shape parameter of Gamma frailty distribution
#' @param gamma_rate A number, the rate parameter of Gamma frailty distribution
#' @param parameters_x A list, with an element named \code{lambda} (\eqn{\lambda_x}), and an element named \code{p} (\eqn{p_x}).
#' @param parameters_y A list, with an element \code{lambda} (\eqn{\lambda_y}), and an element named \code{p} (\eqn{p_y}).
#' @param epsilon A float, the lower bound for support of K and of Y
#' @param int.upper A float, the upper bound for support of K
#' @param control_list A list, with an element named \code{rel.tol} to be passed
#' to \code{distrExIntegrate} as the \code{rel.tol} argument,
#' and an element named \code{subdivisions} which is the maximum number of subintervals \code{subdivisions}
#' argument in \code{distrExIntegrate} (see \code{?distrExIntegrate}).
#' @param correlated If \code{FALSE}, then \eqn{X} and \eqn{Y} are independent. If \code{TRUE}, then
#' there is a shared frailty term \eqn{K} that induces correlation between \eqn{X} and \eqn{Y}.
#' @return density function evaluated at \code{Z=z}.
f.Z.Weibull <- function(
  z,
  parameters_x,
  parameters_y,
  gamma_shape,
  gamma_rate,
  epsilon,
  int.upper,
  control_list, correlated){
  cat("+")
  sapply(z, function(z_i) distrExIntegrate(
    f.YZ.Weibull, lower=epsilon, upper=z_i,
    z=z_i, parameters_x=parameters_x,
    parameters_y=parameters_y,
    gamma_shape=gamma_shape, gamma_rate=gamma_rate,
    epsilon=epsilon,
    int.upper=int.upper,
    correlated=correlated,
    control_list=control_list,
    subdivisions=control_list$subdivisions,
    rel.tol=control_list$rel.tol
  ))
}

#' Survival function for conditionally Weibull random variable \eqn{Y|K}
#' marginalized over \eqn{K}.
#' 
#' If \code{!correlated}: \eqn{Y \sim Weibull(shape=p, rate=\lambda)} \cr
#' If \code{correlated}: \eqn{Y \sim Weibull(shape=p, rate=\lambda \cdot k^{p}) \quad K \sim Gamma(a, b)}
#' 
#' @importFrom stats dweibull
#' @param t A number or vector at which to evaluate the density function
#' @param parameters A list of the two parameters
#' of the Weibull distribution p and lambda
#' before the re-parameterization described for \code{k} in \code{f.Weibull}.
#' @param a The shape parameter for gamma distribution
#' where in \code{f.Weibull}, \eqn{k ~ Gamma(gamma_shape, gamma_rate)}.
#' @param b The shape parameter for gamma distribution
#' where in \code{f.Weibull}, \eqn{k ~ Gamma(gamma_shape, gamma_rate)}.
#' @param epsilon A float, the lower bound for support of Y
#' @param correlated If \code{FALSE}, then this is just a regular Weibull
#' density that you could get with integrating \code{dweibull} or exact survival function.
#' If \code{TRUE}, then include the frailty parameterization
#' marginalizing over \code{k} in \code{f.Weibull}.
S.Y.Weibull <- function(
  t,
  parameters,
  a, b,
  epsilon,
  correlated){
  if(correlated){
    int_result <- sapply(t, function(t_i) 1 - distrExIntegrate(
      f.Y.Weibull, lower=epsilon, upper=t_i,
      parameters=parameters,
      gamma_shape=a, gamma_rate=b,
      correlated=correlated))
  } else {
    int_result <- exp(-(t*parameters$lambda)^parameters$p)
  }
  return(int_result)
}

#' Survival function for the sum of two conditionally Weibull random variables
#' with a shared frailty term.
#' 
#' If \code{correlated}, the function computes the density of Z by
#' first marginalizing over the support of K, then over Y from 0 to z. For computational efficiency,
#' it integrates K from \code{epsilon} to \code{int.upper} and Y from \code{epsilon} to z, where \code{epsilon} can be 0
#' and \code{int.upper} can be \code{Inf}. Small deviations from 0 and \code{Inf} will
#' approximate the true integral.
#' 
#' \eqn{X \sim Weibull(shape=p_{x}, rate=\lambda_{x} \cdot k^{p_{x}})} \cr
#' \eqn{Y \sim Weibull(shape=p_{y}, rate=\lambda_{y} \cdot k^{p_{y}})} \cr
#' \eqn{Z := X + Y} \cr
#' \eqn{K \sim Gamma(shape=a, rate=b)}
#' 
#' @importFrom distrEx distrExIntegrate
#' @importFrom parallel mclapply
#' @param t A number or vector at which to evaluate the survival function
#' @param a A number, the shape parameter of Gamma frailty distribution
#' @param b A number, the rate parameter of Gamma frailty distribution
#' @param parameters_x A list, with an element named \code{lambda} (\eqn{\lambda_x}), and an element named \code{p} (\eqn{p_x}).
#' @param parameters_y A list, with an element \code{lambda} (\eqn{\lambda_y}), and an element named \code{p} (\eqn{p_y}).
#' @param epsilon A float, the lower bound for support of K, Y, and Z
#' @param int.upper A float, the upper bound for support of K
#' @param control_list A list, with an element named \code{rel.tol} to be passed
#' to \code{distrExIntegrate} as the \code{rel.tol} argument,
#' and an element named \code{subdivisions} which is the maximum number of subintervals \code{subdivisions}
#' argument in \code{distrExIntegrate} (see \code{?distrExIntegrate}).
#' @param correlated If \code{FALSE}, then \eqn{X} and \eqn{Y} are independent. If \code{TRUE}, then
#' there is a shared frailty term \eqn{K} that induces correlation between \eqn{X} and \eqn{Y}.
#' @return density function evaluated at \code{Z=z}.
S.Z.Weibull <- function(
  t,
  parameters_x,
  parameters_y,
  a, b,
  epsilon,
  int.upper,
  control_list,
  correlated=TRUE){
  cat("\nIntegrating S.Z with parameter values",
      "\nlambda: ", parameters_x$lambda,
      ", mu: ", parameters_y$lambda,
      ", p: ", parameters_x$p,
      ", q: ", parameters_y$p,
      ", b: ", b,
      "\n")
  s <- unlist(mclapply(t, function(t_i) 1 - distrExIntegrate(
    f.Z.Weibull, lower=epsilon, upper=t_i,
    parameters_x=parameters_x,
    parameters_y=parameters_y,
    gamma_shape=a, gamma_rate=b,
    epsilon=epsilon,
    int.upper=int.upper,
    correlated=correlated,
    control_list=control_list,
    subdivisions=control_list$subdivisions,
    rel.tol=control_list$rel.tol
  )))
  return(s)
}
