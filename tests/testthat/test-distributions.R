context("Probability distributions")

test_that("exponential survival function", {
  expect_equal(S.Y.Exponential(t=0, a=1, b=1, parameters=list(lambda=1), correlated=F), 1)
  expect_equal(S.Y.Exponential(t=0, a=1, b=1, parameters=list(lambda=1), correlated=T), 1)
  expect_equal(S.Y.Exponential(t=4, a=1, b=1, parameters=list(lambda=1), correlated=F), exp(-4))
})