context("COmputation of log sum of exponentials")

x <- 1:10
test_that("", {
  expect_that(getLogSumOfExponentials(x), equals(log(sum(exp(x)))))
})