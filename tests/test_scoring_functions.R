context("Scoring functions")

vOrder <- c(2, 1, 3)
test_that("Log local structure prior given order produces expected result when no actual parents", {
  expect_that(getLogLocalStructurePrior(2, integer(0), vOrder), equals(log(1)))
}) 

test_that("Log local structure prior given order produces expected result when no actual parents", {
  expect_that(getLogLocalStructurePrior(1, integer(0), vOrder), equals(log(1)))
}) 

expected <- log(getNChooseKPrior(2, 1))
test_that("Log local structure prior given order produces expected result when two possible parents and one actual", {
  expect_that(getLogLocalStructurePrior(3, 1, vOrder), equals(expected))
}) 

test_that("Log local structure prior is -Inf for impossible cases", {
  expect_that(getLogLocalStructurePrior(2, c(1,3), vOrder), equals(-Inf))
  expect_that(getLogLocalStructurePrior(100, c(1,3), vOrder), equals(-Inf))
})



node <- 1
parents <- vector()
mockGetAlphas = function(node, parents) {
  if (node != 1 || length(parents) != 0) stop("Wrong node or parents")
  
  matrix(1/3, 1, 3)
}
mockGetSuffStats = function(node, parents) {
  if (node != 1 || length(parents) != 0) stop("Wrong node or parents")
  
  retVal <- matrix(0, nrow=1, ncol=3)
  retVal[1,] <- c(0, 3, 2)
  return(retVal)
}
logScore <- (lgamma(3*1/3) - lgamma(3*1/3 + 0+3+2)) + (lgamma(1/3+0) + lgamma(1/3+3) + lgamma(1/3+2) - 3*lgamma(1/3))

test_that("Log likelihood of a node without parents computed correctly", {
  expect_that(getLogLocalDataLikelihood(node, parents, mockGetAlphas, mockGetSuffStats), equals(logScore))
})


node <- 1
parents <- 2
mockgetAlphas = function(node, parents) {
  if (node != 1 || parents != 2) stop("Wrong node or parents")
  
  matrix(1/9, nrow=3, ncol=3)
}
mockGetSuffStats = function(node, parents) {
  if (node != 1 || parents != 2) stop("Wrong node or parents")
  
  retVal <- matrix(0, nrow=3, ncol=3)
  retVal[1,] <- c(0, 1, 1)
  retVal[2,] <- c(0, 2, 1)
  retVal[3,] <- c(0, 0, 0)
  return(retVal)
}
logScore <- ( (lgamma(3/9) - lgamma(3/9 + 0+1+1)) + (lgamma(1/9 + 0) + lgamma(1/9 + 1) + lgamma(1/9 + 1) - 3*lgamma(1/9)) )  + ( (lgamma(3/9) - lgamma(3/9 + 0+2+1)) + (lgamma(1/9+0) + lgamma(1/9+2) + lgamma(1/9+1) - 3*lgamma(1/9)) )  +  ( (lgamma(3/9) - lgamma(3/9 + 0+0+0)) + (3*lgamma(3/9+0) - 3*lgamma(3/9)) )

test_that("Log likelihood of a node with 1 parent computed correctly", {
  expect_that(getLogLocalDataLikelihood(node, parents, mockgetAlphas, mockGetSuffStats), equals(logScore))
})



# Mock function
functLogLocalOrderScore <- function(node, vOrder) {
  if (node == 1) return(-1)
  else if (node == 2) return(-2)
  else if (node == 3) return(-3)
  else if (node == 4) return(-4)
}
test_that("All order scores computed as expected", {
  expect_that(getLogLocalOrderScores(vOrder=c(4,3,2,1), functLogLocalOrderScore), equals(c(-1,-2,-3,-4)))
})




# For future reference: these are the structure and parameters from which mObs are drawn in random
cardinalities <- c(3, 3, 3)
maxParents <- 2

mAdj <- matrix(FALSE, 3, 3)
mAdj[1,2] <- TRUE
mAdj[1,3] <- TRUE

arrTheta <- array(NA, c(3, 3, 3))
arrTheta[1,,1] <- c(0.05, 0.75, 0.2)
arrTheta[1,,2] <- c(0.25, 0.25, 0.5)
arrTheta[2,,2] <- c(0.40, 0.40, 0.2)
arrTheta[3,,2] <- c(0.80, 0.15, 0.05)
arrTheta[1,,3] <- c(0.55, 0.05, 0.4)
arrTheta[2,,3] <- c(0.60, 0.10, 0.3)
arrTheta[3,,3] <- c(0.20, 0.20, 0.60)

mObs <- matrix(NA, 5, 3)
mObs[1,] <- c(3, 1, 3)
mObs[2,] <- c(2, 3, 3)
mObs[3,] <- c(2, 1, 3)
mObs[4,] <- c(2, 2, 1)
mObs[5,] <- c(2, 1, 3)

# This order score assumes using prior P(G | <) = (n-1 choose |Pa(Xi)|)
logScoreOrder123 <- -14.83126557412244928003

functBDPriorParams <- createBDeuPriorParamsProvider(cardinalities, equivalentSampleSize=1)
functSuffStats <- createSufficientStatsProvider(cardinalities, mObs)
functLogLocalStructureScore <- createLogLocalStructureScoringFunction(cardinalities, functBDPriorParams, functSuffStats)
functLogLocalOrderScore <- createCustomLogLocalOrderScoringFunction(maxParents=2, functLogLocalStructureScore)
test_that("Order score computed as expected when scoring function created using factory function", {
  expect_that(sum(getLogLocalOrderScores(vOrder=c(1,2,3), functLogLocalOrderScore)), equals(logScoreOrder123))
})

