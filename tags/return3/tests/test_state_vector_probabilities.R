context("Computation of state vector posterior probabilities")

cardinalities <- c(3, 2, 2)
inputOrder <- c(3, 2, 1)
inputState <- c(3, 1, 1)

mockFamiliesAndLogStructureScores <- function(node, vOrder) {
  #cat("Node", node, "order", vOrder, fill=T)
  
  if (!all(vOrder == inputOrder)) stop("Invalid order")
  
  if (node == 1) {
    return(list(
      parentSets = list(integer(0), 3, 2, c(2, 3)),
      scores = c(-1, -2, -3, -4)
    ))
  } else if (node == 2) {
    return(list(
      parentSets = list(integer(0), 3),
      scores = c(-1, -2)
    ))
  } else if (node == 3) {
    return(list(
      parentSets = list(integer(0)),
      scores = -1
    ))
  } 
  
  stop("Invalid node or order")
}

mockStateProbability <- function(node, nodeState, vParents, vParentStates, parentsSorted) {
  #cat("Node", node, "state", nodeState, "parents", vParents, "config", vParentStates, fill=T)
  
  if (node == 1 && nodeState == 3) {
    if (length(vParents) == 0) return(1/2)
    else if (all(vParents == 3) && all(vParentStates == 1)) return(1/3)
    else if (all(vParents == 2) && all(vParentStates == 1)) return(1/4)
    else if (all(vParents == c(2,3)) && all(vParentStates == c(1,1))) return(1/5)
             
  } else if (node == 2 && nodeState == 1) {
    if (length(vParents) == 0) return(1/2)
    else if (all(vParents == 3) && all(vParentStates == 1)) return(1/3)
    
  } else if (node == 3 && nodeState == 1) {
    if (length(vParents) == 0) return(1/2)
  }
  
  stop("Invalid node, nodeState, parent set, or parent states")
}

# Terms of denominator P(D | <) (not in ln scale!)
orderScores <- c(exp(-1) + exp(-2) + exp(-3) + exp(-4),
                 exp(-1) + exp(-2),
                 exp(-1))

# Terms of numerator (not in ln scale!)
stateProbsTimesScores <- c(1/2*exp(-1) + 1/3*exp(-2) + 1/4*exp(-3) + 1/5*exp(-4),
                           1/2*exp(-1) + 1/3*exp(-2),
                           1/2*exp(-1))   
# Expected P(x | D, <)
expectedProb <- prod(stateProbsTimesScores)/prod(orderScores)

# Computed P(x | D, <)
computedProb <- getStateVectorProbability(inputState, inputOrder, mockStateProbability, mockFamiliesAndLogStructureScores)
test_that("Posterior probability of state vector computed as expected", {
  expect_that(computedProb, equals(expectedProb))
})

# Computed P(x | D, <)
computedProb <- getStateVectorProbability(inputState, inputOrder, mockStateProbability, mockFamiliesAndLogStructureScores)
test_that("Posterior probability of state vector computed as expected when order scores given", {
  expect_that(computedProb, equals(expectedProb))
})