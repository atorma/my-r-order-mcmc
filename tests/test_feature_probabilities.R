context("Computation of edge probabilities from given node orders")


inputOrder <- c(3,2,1)

mockFamiliesAndLogStructureScores <- function(node, vOrder) {
  if (!all(vOrder == inputOrder)) stop("wrong order")
  
  if (node == 1) {
    return(list(
      parentSets = list(c(2,3), 2, 3, integer(0)),
      scores = c(-1, -2, -3, -4)
    ))    
  } else if (node == 2) {
    return(list(
      parentSets = list(integer(0), 3),
      scores = c(-1, -2)
    ))
  }
  
  stop("Unexpected node")
}

logNodeScore2 <- log( exp(-1) + exp(-2) )
logEdgeScore3to2 <- log( exp(-2) )
probEdge3to2 <- exp(logEdgeScore3to2 - logNodeScore2)

logNodeScore1 <- log( exp(-4) + exp(-3) + exp(-2) + exp(-1) )
logEdgeScore3to1 <- log( exp(-3) + exp(-1) )
logEdgeScore2to1 <- log( exp(-2) + exp(-1) )
probEdge2to1 <- exp(logEdgeScore2to1 - logNodeScore1)
probEdge3to1 <- exp(logEdgeScore3to1 - logNodeScore1) 

# All other edge probabilities given order 3-2-1 are zero
mAllEdgeProbs <- matrix(0, 3, 3)
mAllEdgeProbs[2,1] <- probEdge2to1
mAllEdgeProbs[3,1] <- probEdge3to1
mAllEdgeProbs[3,2] <- probEdge3to2

edge2to1 <- c(2,1)

test_that("Probability of a single edge computed correctly", {
  expect_that(getEdgeProbability(edge2to1, inputOrder, mockFamiliesAndLogStructureScores), equals(probEdge2to1))
})

test_that("If edge is not consistent with order then probability is zero", {
  expect_that(getEdgeProbability(c(1,3), inputOrder, mockFamiliesAndLogStructureScores), equals(0))
})

test_that("Probabilities of all edges computed correctly from a single order", {
  expect_that(getEdgeProbabilities(inputOrder, mockFamiliesAndLogStructureScores), equals(mAllEdgeProbs))
})





mInputOrders <- matrix(c(inputOrder,
                         2, 1, 3), 
                       nrow=2, byrow=TRUE)
mockFamiliesAndLogStructureScores <- function(node, vOrder) {
  #cat("Node", node, "order", vOrder, fill=T)
  
  if (node == 1 && all(vOrder == mInputOrders[1,])) {
    return(list(
      parentSets = list(c(2,3), 2, 3, integer(0)),
      scores = c(-1, -2, -3, -4)
      ))
    
  } else if (node == 1 && all(vOrder == mInputOrders[2,])) {
    return(list(
      parentSets = list(2, integer(0)),
      scores = c(-1, -5)
      ))
  }
  
  stop("Invalid node or order")
}
probsEdge2to1 <- numeric(2)
probsEdge2to1[1] <- probEdge2to1
probsEdge2to1[2] <- exp( log( exp(-1) ) - log( exp(-5) + exp(-1) ) ) 
meanProbEdge2to1 <- mean(probsEdge2to1)

test_that("Probability of an edge is average probability over all input orders", {
  expect_that(getEdgeProbability(edge2to1, mInputOrders, mockFamiliesAndLogStructureScores), equals(meanProbEdge2to1))
})