context("Computation of edge probabilities")

inputEdge <- c(2,1) # 2 -> 1
inputOrder <- c(3,2,1)
maxParents <- 2

mockLogLocalStructureScore <- function(node, vParents, vOrder) {
  #cat("Node", node, "parents", vParents, "order", vOrder, fill=T)
  
  if (node == 1 && all(vOrder == inputOrder)) {
    if (length(vParents) == 0) return(-4)
    else if (all(vParents == 3)) return(-3)
    else if (all(vParents == 2)) return(-2)
    else if (all(vParents == c(2,3))) return(-1)
  }
  
  stop("Invalid node, order, or parent set")
}
expectedPr <- exp( log( exp(-2) + exp(-1) ) - log( exp(-4) + exp(-3) + exp(-2) + exp(-1) ) )

test_that("Probability of an edge computed correctly", {
  expect_that(getEdgeProbability(inputEdge, inputOrder, maxParents, mockLogLocalStructureScore), equals(expectedPr))
})
test_that("If edge is not possible in order then probability is zer", {
  expect_that(getEdgeProbability(c(1,3), inputOrder, maxParents, mockLogLocalStructureScore), equals(0))
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
logLogalOrderScore <- log( exp(-4) + exp(-3) + exp(-2) + exp(-1) )

test_that("Probability of an edge computed correctly when target's order score available", {
  expect_that(getEdgeProbability(inputEdge, inputOrder, maxParents, mockLogLocalStructureScore, targetLogLocalOrderScore=logLogalOrderScore), equals(expectedPr))
})





