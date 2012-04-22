context("Computation of node state posterior probabilities")

cardinalities <- c(2, 2, 3)

mObs <- matrix(NA, 4, 3)
mObs[1,] <- c(1, 2, 2)
mObs[2,] <- c(2, 1, 1)
mObs[3,] <- c(1, 2, 2) 
mObs[4,] <- c(1, 1, 2)
# note that state 3 of node 3 never observed

# function under test f(node, nodeState, parents, parentStates)
# should return "E[theta_ijk | alpha_ij, D] = (N_ijk + alpha_ijk)/(N_ij + alpha_ij)",
# alpha_ijk = equivSampleSize/(r_i*q_i)
functBDeuPriorParams <- createBDeuPriorParamsProvider(cardinalities, equivalentSampleSize=1)
getStateProb <- createStateProbabilityFunction(cardinalities, mObs, functBDeuPriorParams)

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



# Like above, but now the output is matrix with parameters for all states and parent configurations
# Parent configuration index is obtained by sorting parents in node order and varying parent states the quicker
# the higher the parent is on the sorted list (like function getIndexFromConfig() does). 

thetas.node3.noparents    <- matrix(c((1 + 1/3)/(4 + 3/3), (3 + 1/3)/(4 + 3/3), (0 + 1/3)/(4 + 3/3)), 
                                    nrow=1, ncol=3, byrow=TRUE)

thetas.node3.parent1      <- matrix(c((0 + 1/6)/(3 + 3/6), (3 + 1/6)/(3 + 3/6), (0 + 1/6)/(3 + 3/6), 
                                      (1 + 1/6)/(1 + 3/6), (0 + 1/6)/(1 + 3/6), (0 + 1/6)/(1 + 3/6)), 
                                    nrow=2, ncol=3, byrow=TRUE)

thetas.node3.parents1And2 <- matrix(c((0 + 1/12)/(1 + 3/12), (1 + 1/12)/(1 + 3/12), (0 + 1/12)/(1 + 3/12),  #(1, 1)
                                      (1 + 1/12)/(1 + 3/12), (0 + 1/12)/(1 + 3/12), (0 + 1/12)/(1 + 3/12),  #(2, 1)
                                      (0 + 1/12)/(2 + 3/12), (2 + 1/12)/(2 + 3/12), (0 + 1/12)/(2 + 3/12),  #(1, 2)
                                      (0 + 1/12)/(0 + 3/12), (0 + 1/12)/(0 + 3/12), (0 + 1/12)/(0 + 3/12)), #(2, 2)
                                    nrow=4, ncol=3, byrow=TRUE)

getStateProbMatrix <- createStateProbabilityMatrixFunction(cardinalities, mObs, functBDeuPriorParams)

test_that("Parameter table for each state and parent configuration computed", {
  expect_that(getStateProbMatrix(node=3),                  equals(thetas.node3.noparents))
  expect_that(getStateProbMatrix(node=3, parents=1),       equals(thetas.node3.parent1))
  expect_that(getStateProbMatrix(node=3, parents=c(1, 2)), equals(thetas.node3.parents1And2))
})