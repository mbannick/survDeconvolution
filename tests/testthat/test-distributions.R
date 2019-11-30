context("Probability distributions")

test_that("exponential survival", {
  expect_equal(S.Y.Exponential(t=0, a=1, b=1, parameters=list(lambda=1), correlated=F), 1)
  expect_equal(S.Y.Exponential(t=0, a=1, b=1, parameters=list(lambda=1), correlated=T), 1)
  expect_equal(S.Y.Exponential(t=4, a=1, b=1, parameters=list(lambda=1), correlated=F), exp(-4))
})

test_that("sum of exponential survival", {
  expect_equal(S.Z.Exponential(t=0, a=1, b=1,
                               parameters_x=list(lambda=1), parameters_y=list(lambda=2), correlated=F), 1)
  expect_equal(S.Z.Exponential(t=0, a=1, b=1,
                               parameters_x=list(lambda=2), parameters_y=list(lambda=1), correlated=T), 1)
})

test_that("density for Weibull", {
  tol <- 1e-5
  x <- seq(0, 100, by=tol)
  expect_equal(sum(tol * f.Weibull(x, k=0.5, parameters=list(p=1, lambda=0.5))),
               1, tolerance=tol)
})

test_that("density for conditional Weibull", {
  f1 <- f.Y.Weibull(y=c(0, 1), parameters=list(lambda=1, p=1),
                    gamma_shape=1, gamma_rate=1, correlated=TRUE)
  f2 <- f.Y.Weibull(y=c(0, 1), parameters=list(lambda=1, p=1),
                    gamma_shape=1, gamma_rate=1, correlated=FALSE)
  expect_equal(f1, c(1.0, 0.25))
  expect_equal(f2, c(1.0, 0.3678794))
})

test_that("joint Z Y K density Weibull", {
  f <- f.YZK.Weibull(k=c(0, 1, 2, 3), z=1, y=1, parameters_x=list(lambda=1, p=1),
                     parameters_y=list(lambda=0.5, p=1),
                     gamma_shape=1, gamma_rate=1)
  expect_equal(f, c(0.00000000, 0.11156508, 0.09957414, 0.04999048))
})

test_that("joint Z Y density Weibull", {
  # Integrate to infinity
  f <- f.YZ.Weibull(z=1, y=1, parameters_x=list(lambda=1, p=1), parameters_y=list(lambda=0.5, p=1),
                    gamma_shape=1, gamma_rate=1, epsilon=1e-3, control_list=list(rel.tol=1e-3), int.upper=Inf)
  expect_equal(f, 0.2962963)
  # Lower upper bound
  f <- f.YZ.Weibull(z=1, y=1, parameters_x=list(lambda=1, p=1), parameters_y=list(lambda=0.5, p=1),
                    gamma_shape=1, gamma_rate=1, epsilon=1e-3, control_list=list(rel.tol=1e-3), int.upper=5)
  expect_equal(f, 0.2902943)
  f <- f.YZ.Weibull(z=1, y=1, parameters_x=list(lambda=1, p=1), parameters_y=list(lambda=0.5, p=1),
                    gamma_shape=1, gamma_rate=1, epsilon=1e-3, control_list=list(rel.tol=1e-3), int.upper=5,
                    correlated=FALSE)
  expect_equal(f, 0.3032653)
})

test_that("Z density Weibull", {
  f <- f.Z.Weibull(z=1, parameters_x=list(lambda=1, p=1), parameters_y=list(lambda=0.5, p=1), gamma_shape=1, gamma_rate=1,
                   epsilon=1e-3, control_list=list(rel.tol=1e-3, subdivisions=100), int.upper=Inf, correlated=TRUE)
  expect_equal(f, 0.1943194)
})
# ADD IN TOLERANCE
test_that("Z survival Weibull", {
  f <- S.Z.Weibull(t=c(0, 1, Inf), parameters_x=list(lambda=1, p=1),
                   parameters_y=list(lambda=0.5, p=1),
                   a=1, b=1, epsilon=0, control_list=list(rel.tol=1e-3, subdivisions=100),
                   int.upper=Inf, correlated=TRUE)
  expect_equal(f, c(1.000000e+00, 8.333333e-01, 3.090206e-06))
  f <- S.Z.Weibull(t=c(0, 1, Inf), parameters_x=list(lambda=1, p=1),
                   parameters_y=list(lambda=0.5, p=1),
                   a=1, b=1, epsilon=1e-4, control_list=list(rel.tol=1e-3, subdivisions=100),
                   int.upper=Inf, correlated=FALSE)
  expect_equal(f, c(1.000000e+00, 8.452135e-01, 4.999911e-05))
})

test_that("Y survival Weibull", {
  f1 <- S.Y.Weibull(t=c(0), parameters=list(lambda=1, p=1), a=1, b=1, correlated=TRUE, epsilon=0)
  f2 <- S.Y.Weibull(t=c(0), parameters=list(lambda=1, p=1), a=1, b=1, correlated=FALSE, epsilon=0)
  f3 <- S.Y.Weibull(t=c(Inf), parameters=list(lambda=1, p=1), a=1, b=1, correlated=TRUE, epsilon=0)
  f4 <- S.Y.Weibull(t=c(Inf), parameters=list(lambda=1, p=1), a=1, b=1, correlated=TRUE, epsilon=1e-3)
  
  expect_equal(f1, 1)
  expect_equal(f2, 1)
  expect_equal(f3, 0)
  expect_equal(f4, 0, tolerance=1e-3)
})