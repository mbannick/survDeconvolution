library(data.table)
library(magrittr)
library(survival)

one.simulate <- function(
  parameters,
  form='Weibull',
  correlated=TRUE){
  if(correlated){
    k <- rgamma(
      n=1, shape=parameters$alpha,
      rate=parameters$beta)
  }
  if(form=='Weibull'){
    if(correlated){
      scale1 <- 1/(parameters$lambda*k^(1/parameters$p))
      scale2 <- 1/(parameters$mu*k^(1/parameters$q))
    } else {
      scale1 <- 1/parameters$lambda
      scale2 <- 1/parameters$mu
    }
    
    x <- rweibull(
      n=1,
      scale=scale1,
      shape=parameters$p)
    y <- rweibull(
      n=1,
      scale=scale2,
      shape=parameters$q)
    
  } else if(form == 'Exponential'){
    if(correlated){
      rate1 <- (parameters$lambda*k)
      rate2 <- (parameters$mu*k)
    } else {
      rate1 <- (parameters$lambda)
      rate2 <- (parameters$mu)
    }
    x <- rexp(
      n=1,
      rate=rate1)
    y <- rexp(
      n=1,
      rate=rate2)
  }
  
  z <- x + y
  event <- 1
  return(c(x, y, z, event))
}

data.simulate <- function(
  n,
  parameters,
  form,
  correlated,
  times=seq(0, 20, by=0.01)){
  
  sims <- replicate(
    n,
    one.simulate(
      parameters=parameters,
      form=form,
      correlated=correlated)) %>% t %>% as.data.table
  setnames(sims, c("x", "y", "z", "event"))
  
  km_x <- extract.KM(time=sims$x, event=sims$event, times=times)
  km_y <- extract.KM(time=sims$y, event=sims$event, times=times)
  km_z <- extract.KM(time=sims$z, event=sims$event, times=times)
  
  df_x <- data.table(
    time=sims$x,
    event=sims$event
  )
  df_y <- data.table(
    time=sims$y,
    event=sims$event
  )
  df_z <- data.table(
    time=sims$z,
    event=sims$event
  )
  
  result <- list(x=list(df=df_x, km=km_x),
                 y=list(df=df_y, km=km_y),
                 z=list(df=df_z, km=km_z))
  return(result)
}

extract.KM <- function(time, event, times=seq(0, 20, by=0.1), parametric_distribution=NA){
  surv <- Surv(time=time, event=event)
  if(is.na(parametric_distribution)){
    fit <- survfit(surv ~ 1, conf.type='log-log')
    km <- summary(fit, times=times)
    result <- data.table(
      time=km$time,
      n.risk=km$n.risk,
      n.event=km$n.event,
      n.censor=km$n.censor,
      surv=km$surv,
      std.err=km$std.err
    )
  } else {
    if(parametric_distribution == 'weibull'){
      fit <- survreg(surv ~ 1, dist='weibull')
      shape <- 1 / fit$scale
      scale <- (exp(coef(fit)))
      result <- list(
        surv=exp(-(times/scale)^shape),
        shape=shape,
        scale=scale
      )
    } else if(parametric_distribution == 'exponential'){
      fit <- survreg(surv ~ 1, dist='exponential')
      scale <- exp(coef(fit))
      result <- list(
        surv=exp(-(times/scale)),
        scale=scale
      )
    }
  }
  return(list(result=result, fit=fit))
}