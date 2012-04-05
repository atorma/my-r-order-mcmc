context("Computing sufficient statistics")

cardinalities <- c(2, 2, 3)
mObs <- matrix(NA, 4, 3)
mObs[1,] <- c(1, 2, 2)
mObs[2,] <- c(2, 1, 1)
mObs[3,] <- c(1, 2, 2) # this observed now twice
mObs[4,] <- c(1, 1, 3)


getSuffStat <- createSufficientStatsProvider(cardinalities, mObs)

# Case when node has parents
node <- 3
parents <- c(1,2)

expectedSuffStat <- matrix(0, prod(cardinalities[parents]), cardinalities[node])
expectedSuffStat[3,2] <- 2
expectedSuffStat[2,1] <- 1
expectedSuffStat[1,3] <- 1
test_that("Sufficient statistics by each parent configuration reported as matrix when node has parents", {
  expect_that(getSuffStat(node, parents), equals(expectedSuffStat))
})


# Case when node has no parents
node <- 3
parents <- integer(0)

expectedSuffStat <- matrix(0, prod(cardinalities[parents]), cardinalities[node])
expectedSuffStat[1,1] <- 1
expectedSuffStat[1,2] <- 2
expectedSuffStat[1,3] <- 1

test_that("Sufficient statistics by each parent configuration reported as matrix when node has no parents", {
  expect_that(getSuffStat(node, parents), equals(expectedSuffStat))
})

test_that("Counting node's sufficient statistics throws error if a parent index falls outside of range", {
  expect_that(getSuffStat(node, parents=100), throws_error())
})


