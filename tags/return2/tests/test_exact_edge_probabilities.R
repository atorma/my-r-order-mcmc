context("Computation of exact edge probabilities")

numNodes <- 3
maxParents <- 1
# All orderings of nodes 1, 2, 3
mOrders <- matrix(c(1, 2, 3,
                    2, 1, 3,
                    3, 2, 1,
                    1, 3, 2,
                    2, 3, 1,
                    3, 1, 2),
                    nrow=6, byrow=T)

mockLogLocalStructureScore <- function(node, vParents, vOrder) {
  #cat("Node", node, "parents", vParents, "order", vOrder, fill=T)
  
  if (length(vParents) == 0) return(-1)
  
  if (all(vOrder == mOrders[1,])) {
    if (node == 2) {
      if (vParents == 1) return(-2)    
      
    } else if (node == 3) {
      if (vParents == 1) return(-3)
      else if (vParents == 2) return(-4)
      
    }
  }
  
  if (all(vOrder == mOrders[2,])) {
    
    if (node == 1) {
      if (vParents == 2) return(-2)    
      
    } else if (node == 3) {
      if (vParents == 2) return(-3)
      else if (vParents == 1) return(-4)
      
    }
  }
  
  if (all(vOrder == mOrders[3,])) {
    
    if (node == 2) {
      if (vParents == 3) return(-2)    
      
    } else if (node == 1) {
      if (vParents == 3) return(-3)
      else if (vParents == 2) return(-4)
      
    }
  }
  
  if (all(vOrder == mOrders[4,])) {
    
    if (node == 3) {
      if (vParents == 1) return(-2)    
      
    } else if (node == 2) {
      if (vParents == 1) return(-3)
      else if (vParents == 3) return(-4)
      
    }
  }
  
  if (all(vOrder == mOrders[5,])) {
    
    if (node == 3) {
      if (vParents == 2) return(-2)    
      
    } else if (node == 1) {
      if (vParents == 2) return(-3)
      else if (vParents == 3) return(-4)
      
    }
  }
  
  if (all(vOrder == mOrders[6,])) {
    
    if (node == 1) {
      if (vParents == 3) return(-2)    
      
    } else if (node == 2) {
      if (vParents == 3) return(-3)
      else if (vParents == 1) return(-4)
      
    }
  }
  
  stop("Invalid node, order, or parent set")
}

# all orders have the same score
logOrderScore <- -1 + log(exp(-1) + exp(-2)) + log(exp(-1) + exp(-3) + exp(-4))
# all orders have the same probability after normalization
orderProb <- exp(logOrderScore - log(6*exp(logOrderScore)) )

# probability of edge from position 1 to position 2 always same (nodes differ between orders)
probEdge1to2 <- exp(-2 - log(exp(-1) + exp(-2)) )
# probability of edge from position 1 to position 3 always same (nodes differ between orders)
probEdge1to3 <- exp(-3 - log(exp(-1) + exp(-3) + exp(-4)) )
# probability of edge from position 2 to position 3 always same (nodes differ between orders)
probEdge2to3 <- exp(-4 - log(exp(-1) + exp(-3) + exp(-4)) )

# All edges appear in three different orders, xy*, x*y, *xy
# Then P(edge | D) = sum( P(edge | D)P(< | D), <) is 
edgeProb <- orderProb*(probEdge1to2 + probEdge1to3 + probEdge2to3)
mExactEdgeProbs <- matrix(edgeProb, nrow=numNodes, ncol=numNodes)
# Self-refering edges obviously impossible
mExactEdgeProbs[1,1] <- 0
mExactEdgeProbs[2,2] <- 0
mExactEdgeProbs[3,3] <- 0

mComputedEdgeProbs <- getExactEdgeProbabilities(numNodes, maxParents, mockLogLocalStructureScore)
test_that("Exact edge probabilitites computed correctly when all edges have the same probability", {
  expect_that(mComputedEdgeProbs, equals(mExactEdgeProbs))
})