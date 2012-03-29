context("Computation of edge probabilities from given node orders")


inputOrder <- c(3,2,1)
maxParents <- 2

mockLogLocalStructureScore <- function(node, vParents, vOrder) {
  #cat("Node", node, "parents", vParents, "order", vOrder, fill=T)
  
  if (node == 1 && all(vOrder == inputOrder)) {
    if (length(vParents) == 0) return(-4)
    else if (all(vParents == 3)) return(-3)
    else if (all(vParents == 2)) return(-2)
    else if (all(vParents == c(2,3))) return(-1)
    
  } else if (node == 2 && all(vOrder == inputOrder)) {
    if (length(vParents) == 0) return(-1)
    else if (all(vParents == 3)) return(-2)    
    
  }  
  
  stop("Invalid node, order, or parent set")
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
  expect_that(getEdgeProbability(edge2to1, inputOrder, maxParents, mockLogLocalStructureScore), equals(probEdge2to1))
})

test_that("If edge is not consistent with order then probability is zero", {
  expect_that(getEdgeProbability(c(1,3), inputOrder, maxParents, mockLogLocalStructureScore), equals(0))
})

test_that("Probabilities of all edges computed correctly from a single order", {
  expect_that(getEdgeProbabilities(inputOrder, maxParents, mockLogLocalStructureScore), equals(mAllEdgeProbs))
})



mockLogLocalStructureScore <- function(node, vParents, vOrder) {
  #cat("Node", node, "parents", vParents, "order", vOrder, fill=T)
  
  if (node == 1 && all(vOrder == inputOrder)) {
    # Now that the order score of the node is provided, the function
    # under test can only request parent sets including node 2, the
    # source of the edge 2 -> 1
    if (all(vParents == 2)) return(-2)
    else if (all(vParents == c(2,3))) return(-1)
  }
  
  stop("Invalid node, order, or parent set")
}

test_that("Probability of an edge computed correctly when target's order score available", {
  expect_that(getEdgeProbability(edge2to1, inputOrder, maxParents, mockLogLocalStructureScore, targetLogLocalOrderScore=logNodeScore1), equals(probEdge2to1))
})



mInputOrders <- matrix(c(inputOrder,
                         2, 1, 3), 
                       nrow=2, byrow=TRUE)
mockLogLocalStructureScore <- function(node, vParents, vOrder) {
  #cat("Node", node, "parents", vParents, "order", vOrder, fill=T)
  
  if (node == 1 && all(vOrder == mInputOrders[1,])) {
    if (length(vParents) == 0) return(-4)
    else if (all(vParents == 3)) return(-3)
    else if (all(vParents == 2)) return(-2)
    else if (all(vParents == c(2,3))) return(-1)
    
  } else if (node == 1 && all(vOrder == mInputOrders[2,])) {
    if (length(vParents) == 0) return(-5)
    else if (all(vParents == 2)) return(-1)
  }
  
  stop("Invalid node, order, or parent set")
}
probsEdge2to1 <- numeric(2)
probsEdge2to1[1] <- probEdge2to1
probsEdge2to1[2] <- exp( log( exp(-1) ) - log( exp(-5) + exp(-1) ) ) 
meanProbEdge2to1 <- mean(probsEdge2to1)

test_that("Probability of an edge is average probability over all input orders", {
  expect_that(getEdgeProbability(edge2to1, mInputOrders, maxParents, mockLogLocalStructureScore), equals(meanProbEdge2to1))
})