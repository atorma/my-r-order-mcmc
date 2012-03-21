# Computes P(edge | D, <), the posterior probability of an edge given data and a node ordering.
#
# edge: a vector of two nodes where edge[1] is the parent and edge[2] is the child
getEdgeProbability <- function(edge, vOrder, maxParents, functLogLocalStructureScore) {
  if (length(edge) != 2) stop("Edge does not contain exactly two nodes")
  
  # Compute denominator
  functLogLocalOrderScore <- createCustomLogLocalOrderScoringFunction(maxParents, functLogLocalStructureScore)
  logNodeOrderScore <- functLogLocalOrderScore(edge[2], vOrder)
  
  # Log sum of scores over sets with given parent
  parentSets <- getParentSetsIncludingParent(edge[2], vOrder, 0:maxParents, edge[1])
  familyScores <- vector()
  for (p in 1:length(parentSets)) {
    familyScores[p] <- functLogLocalStructureScore(edge[2], parentSets[[p]], vOrder)
  }
  logEdgeOrderScore <- getLogSumOfExponentials(familyScores)
  
  return(exp(logEdgeOrderScore - logNodeOrderScore))
}

# Computes P(edge | D, <) for all edges.
getEdgeProbabilities <- function(vOrder, maxParents, functLogLocalStructureScore) {
  numNodes <- length(vOrder)
  mEdgeProb <- matrix(0, numNodes, numNodes)
  for (i in 1:numNodes) {
    for (j in i:numNodes) {
      if (i != j) {
        parent <- vOrder[i]
        child <- vOrder[j]
        mEdgeProb[parent, child] <- getEdgeProbability(c(parent,child), vOrder, maxParents, functLogLocalStructureScore)
      }
    }
  }
  return(mEdgeProb)
}