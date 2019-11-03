context("Probability distributions")

test_that("equal qnorm", {
  expect_equal(qnorm(0.5), 0)
})