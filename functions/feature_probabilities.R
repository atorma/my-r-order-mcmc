# Computes the average P(edge | D, <) over given node orderings. If the orderings
# are samples from P(< | D), then the result will approximate P(edge | D).
#
# edge: a vector of two nodes where edge[1] is the parent and edge[2] is the child
# mOrder: a matrix where each row is an ordering of nodes
# maxParents: max number of parents a node can have
# functLogLocalStructureScore: function returning log P(X, Pa(X) | D, <) for each node X
# targetLogLocalOrderScore (optional): pre-calculated log order scores of the target of the edge for each row in mOrder 
getEdgeProbability <- function(edge, mOrders, maxParents, functLogLocalStructureScore, targetLogLocalOrderScores=NULL) {
  if (!is.matrix(mOrders)) {
    mOrders <- matrix(mOrders, nrow=1)
  }
  
  if (length(mOrders) == 0) return(NA)
  if (length(edge) != 2) stop("Edge does not contain exactly two nodes")
  if (length(targetLogLocalOrderScores) > 0 && (nrow(mOrders) != length(targetLogLocalOrderScores))) stop("Number of orders inconsistent with number of input node scores")
  
  if (is.null(targetLogLocalOrderScores)) {
    functLogLocalOrderScore <- createCustomLogLocalOrderScoringFunction(maxParents, functLogLocalStructureScore)
  }
  
  # To compute the probability from one order
  getPr <- function(orderIndex) {
    vOrder <- mOrders[orderIndex,]
    targetLogLocalOrderScore <- targetLogLocalOrderScores[orderIndex]
    
    # Impossible edge given order
    if (which(vOrder == edge[1]) > which(vOrder == edge[2])) return(0)
    
    # Compute denominator, sum of scores over all families of target edge, if needed
    if (is.null(targetLogLocalOrderScore)) {
      targetLogLocalOrderScore <- functLogLocalOrderScore(edge[2], vOrder)
    }
    
    # Log sum of scores over families with given parent
    parentSets <- getParentSetsIncludingParent(edge[2], vOrder, 0:maxParents, edge[1])
    familyScores <- numeric(length(parentSets))
    for (p in 1:length(parentSets)) {
      familyScores[p] <- functLogLocalStructureScore(edge[2], parentSets[[p]], vOrder)
    }
    logEdgeOrderScore <- getLogSumOfExponentials(familyScores)
    
    return(exp(logEdgeOrderScore - targetLogLocalOrderScore))
  }
  
  probabilities <- sapply(1:nrow(mOrders), getPr)
  return(mean(probabilities))
}


# Computes mean P(edge | D, <) for all edges from a given order.
#
# mOrders: a matrix where each row is a node order
# maxParents: max number of parents a node can have
# functLogLocalStructureScore: function returning log P(X, Pa(X) | D, <) for each node X
# mLogLocalOrderScores (optional): a matrix where each column contains pre-calculated log order scores of a node and there are as many rows as inn mOrder
getEdgeProbabilities <- function(mOrders, maxParents, functLogLocalStructureScore, mLogLocalOrderScores=NULL) {
  if (!is.matrix(mOrders)) {
    mOrders <- matrix(mOrders, nrow=1)
  }
  if (!is.null(mLogLocalOrderScores) && !is.matrix(mLogLocalOrderScores)) {
    mLogLocalOrderScores <- matrix(mLogLocalOrderScores, nrow=1)
  }
  
  numNodes <- ncol(mOrders)
  mEdgeProb <- matrix(0, numNodes, numNodes)
  edges <- permutations(numNodes, 2)
  for (i in 1:nrow(edges)) {
    edge <- edges[i,]
    parent <- edge[1]
    child <- edge[2]
    mEdgeProb[parent, child] <- getEdgeProbability(edge, mOrders, maxParents, functLogLocalStructureScore, mLogLocalOrderScores[,child])
  }
  return(mEdgeProb)
}

# Computes the exact posterior probability P(edge | D) for all edges. 
# Warning: very slow for more than about 6 nodes!
getExactEdgeProbabilities <- function(numNodes, maxParents, functLogLocalStructureScore) {
  
  functLogLocalOrderScore <- createCustomLogLocalOrderScoringFunction(maxParents, functLogLocalStructureScore)
  functLogOrderScores <- function(vOrder) {
    return( getLogLocalOrderScores(vOrder, functLogLocalOrderScore) )
  }
  
  allOrders <- permutations(numNodes, numNodes, 1:numNodes)
  allLogLocalOrderScores <- t(apply(as.matrix(allOrders), 1, functLogOrderScores)) # score for each node (col) in each order (row)
  allLogOrderScores <- rowSums(allLogLocalOrderScores)
  
  # Compute P(< | D) assuming P(<) = 1
  logNormalizer <- getLogSumOfExponentials(allLogOrderScores)
  allOrderProbs <- exp(allLogOrderScores - logNormalizer)
  
  # Marginalize P(e | D) = sum P(e | D, <)P(< | D) over < 
  mExactEdgeProb <- matrix(0, numNodes, numNodes)
  for (s in 1:nrow(allOrders)) { 
    mExactEdgeProb <- mExactEdgeProb + allOrderProbs[s] * getEdgeProbabilities(allOrders[s,], maxParents, functLogLocalStructureScore, allLogLocalOrderScores[s,]) 
  }
  
  return(mExactEdgeProb)
}
