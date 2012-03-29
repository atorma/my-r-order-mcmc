context("Bayesian Dirichlet Equivalent parameters")

cardinalities <- c(3, 4, 2)
equivSampleSize <- 4
node <- 3
parents <- c(1, 2)
expectedAlphas <- matrix(4/(2*3*4), nrow=(3*4), ncol=2)
test_that("getBDEuParams returns expected parameter matrix", {
  expect_that(getBDEuParams(node, parents, cardinalities, equivSampleSize), equals(expectedAlphas))
})
