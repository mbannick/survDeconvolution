context("Optimization")

test_that("Exponential deconvolution objective function", {
  data <- list(
    y=c(1, 0.25, 0),
    z=c(1, 0.5, 0.25),
    t_y=c(0, 1, 2),
    t_z=c(0, 1, 2)
  )
  fit <- fit.deconvolution(
    data,
    par.init=c(1, 1.1),
    par.lower=c(0.01, 0.01),
    par.upper=c(5, 5),
    epsilon=1e-5,
    int.upper=Inf,
    control_list=list(
      subdivisions=100,
      rel.tol=1e-2),
    correlated=FALSE,
    parametric_distribution='Exponential',
    method=c("nlminb"),
    optim.rel.tol=1e-3,
    itnmax=1
  )
  expect_equal(c(fit$p1, fit$p2, fit$value), c(1.1582, 1.34284, 0.02784572), tolerance=1e-3)
})

test_that("Weibull deconvolution objective function", {
  data <- list(
    y=c(1, 0.25, 0),
    z=c(1, 0.5, 0.25),
    t_y=c(0, 1, 2),
    t_z=c(0, 1, 2)
  )
  fit <- fit.deconvolution(
    data,
    par.init=c(1, 1, 1.1, 1),
    par.lower=c(0.01, 0.01, 0.01, 0.11),
    par.upper=c(5, 5, 5, 5),
    epsilon=1e-5,
    int.upper=Inf,
    control_list=list(
      subdivisions=100,
      rel.tol=1e-2),
    correlated=FALSE,
    parametric_distribution='Weibull',
    method=c("nlminb"),
    optim.rel.tol=1e-3,
    itnmax=1
  )
  expect_equal(c(fit$p1, fit$p2, fit$p3, fit$p4, fit$value),
               c(1.158205, 0.9741616, 1.341286, 1.032461, 0.02683126), tolerance=1e-3)
})

test_that("Exponential deconvolution objective function error", {
  data <- list(
    y=c(1, 0.25, 0),
    z=c(1, 0.5, 0.25),
    t_y=c(0, 1, 2),
    t_z=c(0, 1, 2)
  )
  expect_error(
    fit <- fit.deconvolution(
      data,
      par.init=c(1, 1.1),
      par.lower=c(0.01, 0.01),
      par.upper=c(5, 5),
      epsilon=1e-5,
      int.upper=Inf,
      control_list=list(
        subdivisions=100,
        rel.tol=1e-2),
      correlated=TRUE,
      parametric_distribution='Exponential',
      method=c("nlminb"),
      optim.rel.tol=1e-3,
      itnmax=1
    ),
    "Initial parameters must be of length 3 for correlated exponential"
  )
})

test_that("Weibull deconvolution objective function error", {
  data <- list(
    y=c(1, 0.25, 0),
    z=c(1, 0.5, 0.25),
    t_y=c(0, 1, 2),
    t_z=c(0, 1, 2)
  )
  expect_error(
    fit <- fit.deconvolution(
      data,
      par.init=c(1, 1),
      par.lower=c(0.11, 0.01),
      par.upper=c(5, 20),
      epsilon=1e-5,
      int.upper=Inf,
      control_list=list(
        subdivisions=100,
        rel.tol=1e-2),
      correlated=TRUE,
      parametric_distribution='Weibull',
      method=c("nlminb"),
      optim.rel.tol=1e-3,
      itnmax=1
    ),
    "Initial parameters must be of length 5 for correlated Weibull"
  )
})