context("Objective functions")

test_that("Weibull deconvolution objective function", {
  obj <- objective.function.Weibull(params=c(1, 2, 1, 2, 1), y=c(1), t_y=c(0), z=c(1), t_z=c(0),
                                    control_list=list(rel.tol=1e-10, subdivisions=100), epsilon=1e-10,
                                    correlated=FALSE, int.upper=100)
  expect_equal(obj, 0)
  obj <- objective.function.Weibull(params=c(1, 2, 1, 2, 1), y=c(0), t_y=c(0), z=c(0), t_z=c(0),
                                    control_list=list(rel.tol=1e-10, subdivisions=100), epsilon=1e-10,
                                    correlated=FALSE, int.upper=100)
  expect_equal(obj, 2)
})

test_that("Exponential deconvolution objective function", {
  obj <- objective.function.Exponential(params=c(1, 2, 1), y=c(1), t_y=c(0), z=c(1), t_z=c(0),
                                    control_list=list(rel.tol=1e-10, subdivisions=100), epsilon=1e-10,
                                    correlated=FALSE, int.upper=100)
  expect_equal(obj, 0)
  obj <- objective.function.Exponential(params=c(1, 2, 1), y=c(0), t_y=c(0), z=c(0), t_z=c(0),
                                    control_list=list(rel.tol=1e-10, subdivisions=100), epsilon=1e-10,
                                    correlated=FALSE, int.upper=100)
  expect_equal(obj, 2)
})

