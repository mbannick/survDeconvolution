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