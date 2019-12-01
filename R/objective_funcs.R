#' Objective function for the Weibull deconvolution
#' 
#' Returns the objective function evaluated with the current parameters and observed data.
#' Used to perform non-linear least squares. This objective function
#' is specific to Weibull deconvolution, so it requires 5 parameters
#' (two for each Weibull and one for the frailty).
#' 
#' @param y A vector of observations representing survival at \code{t_y}
#' @param z A vector of observations representing survival at \code{t_z}
#' @param t_y A vector of times at which survival is observed for \code{y}
#' @param t_z A vector of times at which survival is observed for \code{z}
#' @param params A list of parameters in the order \code{lambda, p, mu, q, beta}
#' @param alpha A number, usually fixed to 1, that is the shape parameter of the frailty distribution
#' @param control_list A named list of options to be passed to \code{distrExIntegrate},
#' called \code{subdivisions} and \code{rel.tol}.
#' @param epsilon The smallest value from which to integrate, passed to integration functions, for approximations
#' @param int.upper The largest value from which to integrate (in place of \code{Inf}) for approximations
#' @param correlated If \code{TRUE}, then \code{y} and \code{x = z - y} are correlated through a shared
#' frailty term with rate parameter \code{beta}. If \code{FALSE}, then they are independent and the function
#' ignores \code{beta}.
#' @return A float, the objective function evaluated at the current parameters
objective.function.Weibull <- function(
  params,
  y, t_y,
  z, t_z,
  control_list,
  alpha=1,
  epsilon=epsilon,
  correlated=TRUE,
  int.upper=Inf){
  parameters <- list(
    lambda=params[1],
    p=params[2],
    mu=params[3],
    q=params[4]
  )
  if(correlated){
    parameters[["beta"]] <- params[5]
  } else {
    parameters[["beta"]] <- 1 # throw away value - won't be used
  }
  y.int.value <- S.Y.Weibull(
    t=t_y,
    parameters=list(
      lambda=parameters$mu,
      p=parameters$q),
    a=alpha,
    b=parameters$beta,
    epsilon=epsilon,
    correlated=correlated)
  
  z.int.value <- S.Z.Weibull(
    t=t_z,
    parameters_x=list(
      lambda=parameters$lambda,
      p=parameters$p),
    parameters_y=list(
      lambda=parameters$mu,
      p=parameters$q),
    a=alpha,
    b=parameters$beta,
    epsilon=epsilon,
    correlated=correlated,
    int.upper=int.upper,
    control_list=control_list)
  
  SS.y <- sum((y-y.int.value)^2)
  SS.z <- sum((z-z.int.value)^2)
  SS <- SS.y + SS.z
  cat("\nSum of Squares Y: ", SS.y)
  cat("\nSum of Squares Z: ", SS.z)
  cat("\nSum of Squares: ", SS)
  return(SS)
}

#' Objective function for the Exponential deconvolution
#' 
#' Returns the objective function evaluated with the current parameters and observed data.
#' Used to perform non-linear least squares. This objective function
#' is specific to Exponential deconvolution, so it requires 3 parameters
#' (one for each Weibull and one for the frailty).
#' 
#' @param y A vector of observations representing survival at \code{t_y}
#' @param z A vector of observations representing survival at \code{t_z}
#' @param t_y A vector of times at which survival is observed for \code{y}
#' @param t_z A vector of times at which survival is observed for \code{z}
#' @param params A list of parameters in the order \code{lambda, mu, beta}
#' @param alpha A number, usually fixed to 1, that is the shape parameter of the frailty distribution
#' @param control_list A named list of options to be passed to \code{distrExIntegrate},
#' called \code{subdivisions} and \code{rel.tol}.
#' @param epsilon The smallest value from which to integrate, passed to integration functions, for approximations
#' @param int.upper The largest value from which to integrate (in place of \code{Inf}) for approximations
#' @param correlated If \code{TRUE}, then \code{y} and \code{x = z - y} are correlated through a shared
#' frailty term with rate parameter \code{beta}. If \code{FALSE}, then they are independent and the function
#' ignores \code{beta}.
#' @return A float, the objective function evaluated at the current parameters
objective.function.Exponential <- function(
  params,
  y, t_y,
  z, t_z,
  control_list,
  alpha=1,
  epsilon=epsilon,
  correlated=TRUE,
  int.upper=Inf){
  parameters <- list(
    lambda=params[1],
    mu=params[2]
  )
  if(correlated){
    parameters[["beta"]] <- params[3]
  } else {
    parameters[["beta"]] <- 1 # throw away value - won't be used
  }
  y.int.value <- S.Y.Exponential(
    t=t_y,
    parameters=list(lambda=parameters$mu),
    a=alpha,
    b=parameters$beta,
    correlated=correlated
  )
  
  z.int.value <- S.Z.Exponential(
    t=t_z,
    parameters_x=list(lambda=parameters$lambda),
    parameters_y=list(lambda=parameters$mu),
    a=alpha,
    b=parameters$beta,
    correlated=correlated
  )
  
  SS.y <- sum((y-y.int.value)^2)
  SS.z <- sum((z-z.int.value)^2)
  SS <- SS.y + SS.z
  cat("\nSum of Squares Y: ", SS.y)
  cat("\nSum of Squares Z: ", SS.z)
  cat("\nSum of Squares: ", SS)
  return(SS)
}
