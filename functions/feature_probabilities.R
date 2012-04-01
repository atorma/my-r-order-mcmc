# Computes the average P(edge | D, <) over given node orderings. If the orderings
# are samples from P(< | D), then the result will approximate P(edge | D).
#
# edge: a vector of two nodes where edge[1] is the parent and edge[2] is the child
# mOrders: a matrix where each row is an ordering of nodes
# functFamiliesAndLogStructureScores 
getEdgeProbability <- function(edge, mOrders, functFamiliesAndLogStructureScores) {
  if (!is.matrix(mOrders)) {
    mOrders <- matrix(mOrders, nrow=1)
  }
  
  if (length(mOrders) == 0) return(NA)
  if (length(edge) != 2) stop("Edge does not contain exactly two nodes")
  

  # To compute the probability from one order
  getProbability <- function(orderIndex) {
    vOrder <- mOrders[orderIndex,]

    # Impossible edge given order
    if (which(vOrder == edge[1]) > which(vOrder == edge[2])) return(0)
    
    # For the denominator we need sum of structure scores over all consistent families of 
    # of the target node
    parentSetsAndLogStructureScores <- functFamiliesAndLogStructureScores(edge[2], vOrder)
    parentSets <- parentSetsAndLogStructureScores$parentSets # all parent sets are sorted ascending
    logStructureScores <- parentSetsAndLogStructureScores$scores
    targetLogLocalOrderScore <- getLogSumOfExponentials(logStructureScores)
    
    # For the numerator we need sum of structure scores of all consistent families of
    # the target node that include the source node
    logNumeratorScores <- numeric(length(logStructureScores)) # shrink it later
    index <- 0
    
    for (j in 1:length(parentSets)) {
      pSet <- parentSets[[j]]
      if (edge[1] %in% pSet) {
        index <- index + 1
        logNumeratorScores[index] <- logStructureScores[j]
      }
    }
    
    # We didn't find any parent sets containing the source node consistent with the order,
    # so probability of the edge given the order is zero
    if (index == 0) {
      return(0)
    }
    
    logNumeratorScores <- logNumeratorScores[1:index]
    logEdgeOrderScore <- getLogSumOfExponentials(logNumeratorScores)
    
    return(exp(logEdgeOrderScore - targetLogLocalOrderScore))
  }
  
  probabilities <- sapply(1:nrow(mOrders), getProbability)
  return(mean(probabilities))
}

# Computes mean P(edge | D, <) for all edges from given node orders.
#
# mOrders: a matrix where each row is a node order
# functFamiliesAndLogStructureScores: 
getEdgeProbabilities <- function(mOrders, functFamiliesAndLogStructureScores) {
  if (!is.matrix(mOrders)) {
    mOrders <- matrix(mOrders, nrow=1)
  }
  
  numNodes <- ncol(mOrders)
  mEdgeProb <- matrix(0, numNodes, numNodes)
  edges <- permutations(numNodes, 2)
  for (i in 1:nrow(edges)) {
    edge <- edges[i,]
    parent <- edge[1]
    child <- edge[2]
    mEdgeProb[parent, child] <- getEdgeProbability(edge, mOrders, functFamiliesAndLogStructureScores)
    #cat("Edge", edge, "probability", mEdgeProb[parent, child], fill=T)
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
  
  # Marginalize: P(e | D) = sum P(e | D, <)P(< | D) over < 
  
  functFamiliesAndLogStructureScores <- function(node, vOrder) {
    parentSets <- getParentSets(node, vOrder, 0:maxParents)
    logStructureScores <- numeric(length(parentSets))
    for (j in 1:length(parentSets)) {
      logStructureScores[j] <- functLogLocalStructureScore(node, parentSets[[j]], vOrder)
    }
    ordering <- sort.list(logStructureScores, decreasing=TRUE, method="quick", na.last=NA)
    logStructureScores <- logStructureScores[ordering]
    parentSets <- parentSets[ordering]
    
    return(list(scores=logStructureScores, parentSets=parentSets))
  }
  
  mExactEdgeProb <- matrix(0, numNodes, numNodes)
  for (s in 1:nrow(allOrders)) { 
    mExactEdgeProb <- mExactEdgeProb + allOrderProbs[s] * getEdgeProbabilities(allOrders[s,], functFamiliesAndLogStructureScores) 
  }
  
  return(mExactEdgeProb)
}
