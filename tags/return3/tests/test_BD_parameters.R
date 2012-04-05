context("Bayesian Dirichlet prior parameters")

cardinalities <- c(3, 4, 2)
node <- 3
parents <- c(1, 2)

equivalentSampleSize <- 4
expectedAlphas <- matrix(4/(2*3*4), nrow=(3*4), ncol=2)
functBDeuParams <- createBDeuPriorParamsProvider(cardinalities, equivalentSampleSize)
test_that("Correct BDeu prior parameter matrix", {
  expect_that(functBDeuParams(node, parents), equals(expectedAlphas))
})


alpha <- 1
expectedAlphas <- matrix(1, nrow=(3*4), ncol=2)
functBDK2Params <- createBDK2PriorParamsProvider(cardinalities, alpha)
test_that("Correct K2 type prior parameter matrix", {
  expect_that(functBDK2Params(node, parents), equals(expectedAlphas))
})