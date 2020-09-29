context("Bootstrap")

test_that("Data simulate", {
  set.seed(1000)
  df <- data.simulate(
    n=10,
    parameters=list(lambda=1, mu=1, p=1, beta=1, alpha=1, q=1),
    correlated=TRUE,
    form='Weibull',
    times=seq(0, 10, by=1)
  )
  expect_equal(nrow(df$x$df), 10)
  expect_equal(nrow(df$x$km$result), 11)
  expect_equal(nrow(df$y$df), 10)
  expect_equal(nrow(df$y$km$result), 11)
  expect_equal(nrow(df$z$df), 10)
  expect_equal(nrow(df$z$km$result), 11)
})

