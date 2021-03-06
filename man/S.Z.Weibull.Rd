% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/survival_funcs.R
\name{S.Z.Weibull}
\alias{S.Z.Weibull}
\title{Survival function for the sum of two conditionally Weibull random variables
with a shared frailty term.}
\usage{
S.Z.Weibull(t, parameters_x, parameters_y, a, b, epsilon, int.upper,
  control_list, correlated = TRUE)
}
\arguments{
\item{t}{A number or vector at which to evaluate the survival function}

\item{parameters_x}{A list, with an element named \code{lambda} (\eqn{\lambda_x}), and an element named \code{p} (\eqn{p_x}).}

\item{parameters_y}{A list, with an element \code{lambda} (\eqn{\lambda_y}), and an element named \code{p} (\eqn{p_y}).}

\item{a}{A number, the shape parameter of Gamma frailty distribution}

\item{b}{A number, the rate parameter of Gamma frailty distribution}

\item{epsilon}{A float, the lower bound for support of K, Y, and Z}

\item{int.upper}{A float, the upper bound for support of K}

\item{control_list}{A list, with an element named \code{rel.tol} to be passed
to \code{distrExIntegrate} as the \code{rel.tol} argument,
and an element named \code{subdivisions} which is the maximum number of subintervals \code{subdivisions}
argument in \code{distrExIntegrate} (see \code{?distrExIntegrate}).}

\item{correlated}{If \code{FALSE}, then \eqn{X} and \eqn{Y} are independent. If \code{TRUE}, then
there is a shared frailty term \eqn{K} that induces correlation between \eqn{X} and \eqn{Y}.}
}
\value{
density function evaluated at \code{Z=z}.
}
\description{
If \code{correlated}, the function computes the density of Z by
first marginalizing over the support of K, then over Y from 0 to z. For computational efficiency,
it integrates K from \code{epsilon} to \code{int.upper} and Y from \code{epsilon} to z, where \code{epsilon} can be 0
and \code{int.upper} can be \code{Inf}. Small deviations from 0 and \code{Inf} will
approximate the true integral.
}
\details{
\eqn{X \sim Weibull(shape=p_{x}, rate=\lambda_{x} \cdot k^{p_{x}})} \cr
\eqn{Y \sim Weibull(shape=p_{y}, rate=\lambda_{y} \cdot k^{p_{y}})} \cr
\eqn{Z := X + Y} \cr
\eqn{K \sim Gamma(shape=a, rate=b)}
}
