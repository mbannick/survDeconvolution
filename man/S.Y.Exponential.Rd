% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/survival_funcs.R
\name{S.Y.Exponential}
\alias{S.Y.Exponential}
\title{Marginal survival function for exponential
random variable \eqn{Y}, with frailty parameter \eqn{K}.}
\usage{
S.Y.Exponential(t, a, b, parameters, correlated)
}
\arguments{
\item{t}{A number or vector, the survival times}

\item{a}{A number, the shape parameter of Gamma frailty distribution}

\item{b}{A number, the rate parameter of Gamma frailty distribution}

\item{parameters}{A list, with an element named lambda, the rate
parameter of the exponential distribution}

\item{correlated}{If TRUE, random variable is conditionally exponential
on the frailty parameter, that has \code{Gamma(shape=a, rate=b)} distribution.}
}
\value{
Survival function evaluated at t
}
\description{
If \code{!correlated}: \eqn{Y ~ Exponential(\lambda)}
If \code{correlated}: \eqn{Y ~ Exponential(\lambda \cdot k) \quad K ~ Gamma(a, b)}
}
