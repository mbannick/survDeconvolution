% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/optimization.R
\name{fit.deconvolution}
\alias{fit.deconvolution}
\title{Fit a deconvolution model}
\usage{
fit.deconvolution(data, correlated = TRUE,
  parametric_distribution = "Weibull", par.init = c(1, 1, 1.1, 1, 1),
  par.lower = c(0.01, 0.01, 0.01, 0.11, 0.01), par.upper = c(5, 5, 5,
  5, 20), epsilon = 1e-05, itnmax = 100, int.upper = Inf,
  optim.rel.tol = 0.001, control_list = list(subdivisions = 100,
  rel.tol = 0.01), method = c("nlminb"))
}
\arguments{
\item{data}{A list of input data with named elements (see \sQuote{Details}).}

\item{correlated}{If \code{TRUE}, then \code{data$y} and \code{x = data$z - data$y} are correlated through a shared
frailty term. If \code{FALSE}, then they are independent and the function
ignores \code{beta}.}

\item{parametric_distribution}{Which parametric distribution to use: 'Weibull' or 'Exponential'}

\item{par.init}{A vector of initial parameter values}

\item{par.lower}{A vector of lower bounds for parameters}

\item{par.upper}{A vector of upper bounds for parameters}

\item{epsilon}{The smallest value from which to integrate}

\item{itnmax}{A number specifying the maximum number of iterations to \code{optimx}}

\item{int.upper}{The largest value from which to integrate}

\item{optim.rel.tol}{Relative tolerance for \code{optimx} optimization}

\item{control_list}{A named list of options to be passed to \code{distrExIntegrate},
called \code{subdivisions} and \code{rel.tol}.}

\item{method}{A string specifying the optimization method to use in \code{optimx}}
}
\description{
Fits a deconvolution model with either Exponential or Weibull
distributions, and with or without correlation induced by adding
a frailty term. If fitting a Weibull model, need 4-element initial/lower/upper parameter vectors.
If fitting an exponential model, need 2-element initial/lower/upper parameter vectors. If \code{correlated},
one additional parameter is needed in both Weibull and exponential models.
}
\details{
The \code{data} argument is a named list with four elements that represent
the observed input data to be fit for the deconvolution model. Its elements are:
\describe{
  \item{\code{y}}{Observed survival at times \code{t_y}}
  \item{\code{t_y}}{Survival times from intermediate time point}
  \item{\code{z}}{Observed survival at times \code{t_z}}
  \item{\code{t_z}}{Survival times from initial time point}
}
\code{y} and \code{t_y} must be of the same length, as do \code{z} and \code{z_y}.
}
