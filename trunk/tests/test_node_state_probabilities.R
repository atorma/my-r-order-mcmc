context("Computation of node state posterior probabilities")

cardinalities <- c(2, 2, 3)

mObs <- matrix(NA, 4, 3)
mObs[1,] <- c(1, 2, 2)
mObs[2,] <- c(2, 1, 1)
mObs[3,] <- c(1, 2, 2) 
mObs[4,] <- c(1, 1, 2)
# note that state 3 of node 3 never observed

# function under test f(node, nodeState, parents, parentStates)
# should return "E[theta | alpha, D] = (N_ijk + alpha_ijk)/(N_ij + alpha_ij)",
# alpha_ijk = equivSampleSize/(r_i*q_i)
getStateProb <- createStateProbabilityFunction(cardinalities, mObs, equivalentSampleSize=1)

test_that("State probability when node has no parents", {
  expect_that( getStateProb(node=1, nodeState=1), equals( (3 + 1/2)/(4 + 2*1/2) ) ) 
  expect_that( getStateProb(node=3, nodeState=3), equals( (0 + 1/3)/(4 + 3*1/3) ) )
})

test_that("State probability when node has 1 parent", {
  expect_that( getStateProb(node=3, nodeState=1, parents=1, parentStates=1), 
               equals( (0 + 1/(3*2))/(3 + 3*1/(3*2)) ) ) 
  expect_that( getStateProb(node=3, nodeState=2, parents=1, parentStates=1), 
               equals( (3 + 1/(3*2))/(3 + 3*1/(3*2)) ) )
})

test_that("State probability when node has 2 parents", {
  expect_that( getStateProb(node=3, nodeState=1, parents=c(1,2), parentStates=c(1,2)), 
               equals( (0 + 1/(3*2*2))/(2 + 3*1/(3*2*2)) ) ) 
  expect_that( getStateProb(node=3, nodeState=2, parents=c(1,2), parentStates=c(1,2)), 
               equals( (2 + 1/(3*2*2))/(2 + 3*1/(3*2*2)) ) )
  expect_that( getStateProb(node=3, nodeState=2, parents=c(1,2), parentStates=c(1,1)), 
               equals( (1 + 1/(3*2*2))/(1 + 3*1/(3*2*2)) ) )
})

test_that("State probability not sensitive to parent ordering", {
  expect_that( getStateProb(node=3, nodeState=1, parents=c(2,1), parentStates=c(2,1)), 
               equals( (0 + 1/(3*2*2))/(2 + 3*1/(3*2*2)) ) ) 
})


test_that("Error if node number outside range defined by length of cardinality vector", {
  expect_that( getStateProb(node=100, nodeState=1), throws_error() ) 
})
test_that("Error if node state outside node's cardinality", {
  expect_that( getStateProb(node=1, nodeState=100), throws_error() ) 
})
test_that("Error if parent vector given but not parent states", {
  expect_that( getStateProb(node=3, nodeState=1, parents=c(1,2)), throws_error() ) 
})
test_that("Error if parent state vector given but not parent vector", {
  expect_that( getStateProb(node=3, nodeState=1, parentStates=c(1,2)), throws_error() ) 
})
test_that("Error if a parent state outside its cardinality", {
  expect_that( getStateProb(node=3, nodeState=1, parents=c(1,2), parentStates=c(100,2)), throws_error() ) 
})
