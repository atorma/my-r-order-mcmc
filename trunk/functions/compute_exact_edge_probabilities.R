# TODO reuse order log scores when computing edge probability
computeExactEdgeProbabilities <- function(numNodes, functLogLocalStructureScore) {
  
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
