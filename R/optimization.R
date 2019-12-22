#' Fit a deconvolution model
#' 
#' Fits a deconvolution model with either Exponential or Weibull
#' distributions, and with or without correlation induced by adding
#' a frailty term. If fitting a Weibull model, need 4-element initial/lower/upper parameter vectors.
#' If fitting an exponential model, need 2-element initial/lower/upper parameter vectors. If \code{correlated},
#' one additional parameter is needed in both Weibull and exponential models.
#' 
#' @importFrom optimx optimx
#' 
#' @export
#' @param data A list of input data with named elements (see \sQuote{Details}).
#' @param par.init A vector of initial parameter values
#' @param par.lower A vector of lower bounds for parameters
#' @param par.upper A vector of upper bounds for parameters
#' @param epsilon The smallest value from which to integrate
#' @param int.upper The largest value from which to integrate
#' @param optim.rel.tol Relative tolerance for \code{optimx} optimization
#' @param control_list A named list of options to be passed to \code{distrExIntegrate},
#' called \code{subdivisions} and \code{rel.tol}.
#' @param correlated If \code{TRUE}, then \code{data$y} and \code{x = data$z - data$y} are correlated through a shared
#' frailty term. If \code{FALSE}, then they are independent and the function
#' ignores \code{beta}.
#' 
#' @details The \code{data} argument is a named list with four elements that represent
#' the observed input data to be fit for the deconvolution model. Its elements are:
#' \describe{
#'   \item{\code{y}}{Observed survival at times \code{t_y}}
#'   \item{\code{t_y}}{Survival times from intermediate time point}
#'   \item{\code{z}}{Observed survival at times \code{t_z}}
#'   \item{\code{t_z}}{Survival times from initial time point}
#' }
#' \code{y} and \code{t_y} must be of the same length, as do \code{z} and \code{z_y}.
#' 
fit.deconvolution <- function(
  data,
  correlated=TRUE,
  parametric_distribution='Weibull',
  par.init=c(1, 1, 1.1, 1, 1),
  par.lower=c(0.01, 0.01, 0.01, 0.11, 0.01),
  par.upper=c(5, 5, 5, 5, 20),
  epsilon=1e-5,
  itnmax=100,
  int.upper=Inf,
  optim.rel.tol=1e-3,
  control_list=list(
    subdivisions=100,
    rel.tol=1e-2),
  method=c("nlminb")){

  if(length(data$y) != length(data$t_y)) stop("y and t_y are not of the same length")
  if(length(data$z) != length(data$t_z)) stop("z and t_z are not of the same length")
  
  if(parametric_distribution == 'Weibull'){
    if(correlated){
      if(length(par.init) != 5) stop("Initial parameters must be of length 5 for correlated Weibull")
      if(length(par.lower) != 5) stop("Parameter lower bound must be of length 5 for correlated Weibull")
      if(length(par.upper) != 5) stop("Parameter upper bound must be of length 5 for correlated Weibull")
    }
    if(!correlated){
      if(length(par.init) != 4) stop("Initial parameters must be of length 4 for uncorrelated Weibull")
      if(length(par.lower) != 4) stop("Parameter lower bound must be of length 4 for uncorrelated Weibull")
      if(length(par.upper) != 4) stop("Parameter upper bound must be of length 4 for uncorrelated Weibull")
    }
  } else if(parametric_distribution == 'Exponential'){
    if(correlated){
      if(length(par.init) != 3) stop("Initial parameters must be of length 3 for correlated exponential")
      if(length(par.lower) != 3) stop("Parameter lower bound must be of length 3 for correlated exponential")
      if(length(par.upper) != 3) stop("Parameter upper bound must be of length 3 for correlated exponential")
    }
    if(!correlated){
      if(length(par.init) != 2) stop("Initial parameters must be of length 2 for uncorrelated exponential")
      if(length(par.lower) != 2) stop("Parameter lower bound must be of length 2 for uncorrelated exponential")
      if(length(par.upper) != 2) stop("Parameter upper bound must be of length 2 for uncorrelated exponential")
    }
  } else {
    stop("Unrecognized parametric distribution, can take do 'Exponential' or 'Weibull'.")
  }
  
  if(class(method) != "character") stop("method argument must be character")
  if(class(epsilon) != "numeric") stop("epsilon argument must be numeric")
  if(class(int.upper) != "numeric") stop("int.upper argument must be numeric")
  if(class(optim.rel.tol) != "numeric") stop("optim.rel.tol argument must be numeric")
  if(class(itnmax) != "numeric") stop("itnmax argument must be numeric")
  
  if(class(par.init) != "numeric") stop("initial parameters must be numeric")
  if(class(par.lower) != "numeric") stop("parameter lower bound must be numeric")
  if(class(par.upper) != "numeric") stop("parameter upper bound must be numeric")
  
  if(class(control_list) != "list") stop("control_list argument must be a list")
  if(class(control_list$subdivisions) != "numeric") stop("control list must have numeric subdivisions element")
  if(class(control_list$rel.tol) != "numeric") stop("control list must have numeric rel.tol element")
  
  # different objective functions for the parametric distributions
  if(parametric_distribution == 'Weibull'){
    obj.func <- objective.function.Weibull
  } else if(parametric_distribution == 'Exponential'){
    obj.func <- objective.function.Exponential
  }
  
  opt <- optimx(
    par=par.init,
    fn=obj.func,
    control=list(
      rel.tol=optim.rel.tol),
    y=data$y, z=data$z,
    t_y=data$t_y, t_z=data$t_z,
    lower=par.lower,
    upper=par.upper,
    epsilon=epsilon,
    int.upper=int.upper,
    correlated=correlated,
    control_list=control_list,
    itnmax=itnmax,
    method=method)
  
  return(opt)
}
