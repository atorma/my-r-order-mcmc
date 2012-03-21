# Computes P(edge | D, <), the posterior probability of an edge given data and a node ordering.
#
# edge: a vector of two nodes where edge[1] is the parent and edge[2] is the child
# vOrder: an ordering of nodes
# maxParents: max number of parents a node can have
# functLogLocalStructureScore: function returning log P(X, Pa(X) | D, <) for each node X
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

# TODO reuse order log scores when computing edge probability
getExactEdgeProbabilities <- function(numNodes, functLogLocalStructureScore) {
  
  functLogLocalOrderScore <- createCustomLogLocalOrderScoringFunction(maxParents, functLogLocalStructureScore)
  functLogOrderScore <- function(vOrder) {
    return( sum(getLogLocalOrderScores(vOrder, functLogLocalOrderScore)))
  }
  
  allOrders <- permutations(numNodes, numNodes, 1:numNodes)
  # Compute log P(D | <)
  allOrderLogScores <- apply(as.matrix(allOrders), 1, functLogOrderScore)
  # Compute P(< | D) assuming P(<) = 1
  logNormalizer <- getLogSumOfExponentials(allOrderLogScores)
  allOrderProbs <- exp(allOrderLogScores - logNormalizer)
  
  # Since we're not sampling from posterior, we must do (the normal)
  # P(e | D) = sum(P(e | D, <)P(< | D), <) 
  mExactEdgeProb <- matrix(0, numNodes, numNodes)
  for (s in 1:nrow(allOrders)) { 
    mExactEdgeProb <- mExactEdgeProb + allOrderProbs[s] * getEdgeProbabilities(allOrders[s,], maxParents, functLogLocalStructureScore) 
  }
  
  return(mExactEdgeProb)
}
