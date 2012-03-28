# Creates a function f(node, parents) 
# such that that f returns a matrix theta containing 
# the expected values of multinomial distribution 
# parameters of the node for each parent configuration:
# theta[j, k] is the parameter for node's state k when
# parent configuration is j.
#
# Parent configurations are indexed by sorting the parents
# in asceding order and then varying states the quicker 
# the higher the parent is on the list. E.g.
# (1, 1) -> 1
# (2, 1) -> 2
# (1, 2) -> 3
# (2, 2) -> 4
#
# Uses Bayesian Dirichlet Equivalent uniform (BDEu) prior 
# with given equivalent sample size.
#
# The returned function caches parameters tables.
#
# See functions getIndexFromConfig() and getConfigFromIndex()
createStateProbabilityMatrixFunction <- function(cardinalities, mObs, equivalentSampleSize=1) {
  getSuffStats <- createSufficientStatsHelper(cardinalities, mObs)
  
  cache <- hash()
  getKey <- function(node, parents) {
    paste(c(parents, node), collapse=" ")
  }
  
  return(function(node, parents=integer(0), parentsSorted=FALSE) {
    #cat("Requested node", node, "parents", parents, fill=T)
    
    # For indexation to be insentitive to parent ordering, we must sort parents ascending,
    # unless we know they're sorted, of course
    if (!parentsSorted) {
      parents <- sort(parents)
    }
    
    key <- getKey(node, parents)
    thetasExpected <- cache[[key]]
    
    if (is.null(thetasExpected)) {
      # assumption: getSuffStats() returns parent states 
      # i.e. rows of matrix are in order as described above
      suffStats <- getSuffStats(node, parents, parentsSorted)
      alphas <- getBDEuParams(node, parents, cardinalities, equivalentSampleSize)
      posteriorCounts <- suffStats + alphas
      
      thetasExpected <- posteriorCounts/rowSums(posteriorCounts)
      
      cache[[key]] <- thetasExpected
    }
    
    return(thetasExpected)
  })
}


# Creates a function f(node, nodeState, parents, parentStates) 
# such that that f returns the posterior probility
# of the node's state given parents and their states.
#
# Uses Bayesian Dirichlet Equivalent uniform (BDEu) prior 
# with given equivalent sample size.
createStateProbabilityFunction <- function(cardinalities, mObs, equivalentSampleSize=1) {
  
  getThetaMatrix <- createStateProbabilityMatrixFunction(cardinalities, mObs, equivalentSampleSize)
  
  return(function(node, nodeState, parents=integer(0), parentStates=integer(0), parentsSorted=FALSE) {
    #cat("Requested node", node, "state", nodeState, "parents", parents, "configuration", parentStates, fill=T)
    
    if (length(parents) != length(parentStates)) stop("Inconsistent length of parent and parent state vectors")
    if (length(nodeState) != 1 || length(node) != 1) stop("Only one node and node state allowed")

    # for indexation to be insentitive to parent ordering, we must sort (parents, parentStates) 
    # by parent node number
    if (!parentsSorted) {
      sortedPerm <- order(parents)
      parents <- parents[sortedPerm]
      parentStates <- parentStates[sortedPerm]
    }
    parentConfigIndex <- getIndexFromConfig(parentStates, cardinalities[parents])
    
    thetasExpected <- getThetaMatrix(node, parents, parentsSorted=TRUE)
    thetaExpected <- thetasExpected[parentConfigIndex, nodeState]
    
    return(thetaExpected)
  })
}

createStateProbabilityFunctionOld <- function(cardinalities, mObs, equivalentSampleSize=1) {
  getSuffStats <- createSufficientStatsHelper(cardinalities, mObs)
  
  cache <- hash()
  getKey <- function(node, nodeState, parents, parentStates) {
    paste(c("n:", node, "s:", nodeState, "p:", parents, "c:", parentStates), collapse=" ")
  }
  
  return(function(node, nodeState, parents=integer(0), parentStates=integer(0), parentsSorted=FALSE) {
    #cat("Requested node", node, "state", nodeState, "parents", parents, "configuration", parentStates, fill=T)
    
    if (length(parents) != length(parentStates)) stop("Inconsistent length of parent and parent state vectors")
    if (length(nodeState) != 1 || length(node) != 1) stop("Only one node and node state allowed")
    
    # for indexation to be insentitive to parent ordering, we must sort (parents, parentStates) 
    # by parent node number
    if (!parentsSorted) {
      orderIndex <- order(parents)
      parents <- parents[orderIndex]
      parentStates <- parentStates[orderIndex]
    }
    
    key <- getKey(node, nodeState, parents, parentStates)
    thetaExpected <- cache[[key]]
    
    if (is.null(thetaExpected)) {
      # assumption: parent states i.e. rows of matrix are in order obtained 
      # by varying states the quicker the smaller the number of the node
      suffStats <- getSuffStats(node, parents)
      alphas <- getBDEuParams(node, parents, cardinalities, equivalentSampleSize)
      posteriorCounts <- suffStats + alphas
      
      parentConfigIndex <- getIndexFromConfig(parentStates, cardinalities[parents])
      thetaExpected <- posteriorCounts[parentConfigIndex, nodeState]/sum(posteriorCounts[parentConfigIndex,])
      
      cache[[key]] <- thetaExpected
    }
    
    return(thetaExpected)
  })
}



# Computes the probabilities of state vectors given node orderings.
# If multiple orders are inputed, then the result is the mean probability of 
# the state over the orders.
#
# mStates: vStates[s, i] is the state of node i in state vector s
# mOrders: mOrders[o, p] is the node at position p of the o:th order
# maxParents: maximum number of parents any node can have
# functNodeStateProbability: function f(node, nodeState, parents, parentState) that returns the probability of a node state given a parent configuration
# functLogLocalStructureScore: function f(node, parents, order) that returns the natural logarithm of score(X, Pa(X) | D, <)
# mlogLocalOrderScores (optional): mlogLocalOrderScores[o, node] is the log score contribution of node to the o:th order score
getStateVectorProbability <- function(mStates, mOrders, maxParents, functNodeStateProbability, functLogLocalStructureScore, mlogLocalOrderScores=NULL) {
  if (!is.matrix(mStates)) {
    mStates <- matrix(mStates, nrow=1)
  }
  if (!is.matrix(mOrders)) {
    mOrders <- matrix(mOrders, nrow=1)
  }
  if (length(mlogLocalOrderScores) > 0 && !is.matrix(mlogLocalOrderScores)) {
    mlogLocalOrderScores <- matrix(mlogLocalOrderScores, nrow=1)
  }
  
  if (length(mOrders) == 0) return(NA)
  
  if (length(mlogLocalOrderScores) > 0 && (nrow(mOrders) != nrow(mlogLocalOrderScores))) stop("Number of orders inconsistent with number of order score vectors")
  if (ncol(mStates) != ncol(mOrders)) stop("state vectors must have as many elements as order vectors") 
  
  
  if (is.null(mlogLocalOrderScores)) {
    functLogLocalOrderScore <- createCustomLogLocalOrderScoringFunction(maxParents, functLogLocalStructureScore)
  }
  
  # In log scale, P(X | D, <) is computed as sum of node specific terms. 
  # Each term requires summing over all possible parent sets of the node,
  # which is expensive. 
  # If we rearrange state vector X in the same order as nodes, we can
  # reuse earlier computations.
  # If the order is e.g. (1, 2, 3, *, *), the state is e.g.
  # (1, 3, 1, *, *), and we're computing the term of node 3,
  # we can use previously computed term of node 3 in its state 1 with 
  # _possible_ parents {1, 2} in states {1, 3} regardless of what nodes
  # and state come after node 3 in the order.
  
  cache <- hash()

  getCacheKey <- function(node, predecessorNodes, vStates) {
    nodeState <- vStates[node]
    predecessorNodes <- sort(predecessorNodes)
    predecessorStates <- vStates[predecessorNodes]
    
    return(paste(c("n:", node, "ns:", nodeState, "p:", predecessorNodes, "ps:", predecessorStates), collapse=" "))
  }
  
  # log P(Xi | D, <) for given node in a given state vector.
  # This is an expensive operation.
  getLogNodeStateProbability <- function(node, vStates, orderIndex) {
    vOrder <- mOrders[orderIndex,]
    
    # numerator in log scale
    # an expensive operation
    parentSets <- getParentSets(node, vOrder, 0:maxParents) # parent sets are sorted ascending
    familyScores <- numeric(length(parentSets))
    for (j in 1:length(parentSets)) {
      parents <- parentSets[[j]]
      familyScores[j] <- log(functNodeStateProbability(node, vStates[node], parents, vStates[parents], parentsSorted=TRUE)) + functLogLocalStructureScore(node, parents, vOrder)
    }
    logWeighedProbabilities <- getLogSumOfExponentials(familyScores)
    
    # denominators in log scale
    logSumOfWeights <- mlogLocalOrderScores[orderIndex, node]
    if (is.null(logSumOfWeights)) {
      # an expensive operation
      logSumOfWeights <- functLogLocalOrderScore(node, vOrder)
    }
    
    return(logWeighedProbabilities - logSumOfWeights)
  }
  
  
  # P(X | D, <) for given order index an state vector
  getProbabilityGivenOrder <- function(vStates, orderIndex) {
    vOrder <- mOrders[orderIndex,]
 
    logNodeStateProbabilities <- numeric(length(vStates))
    for (i in 1:length(vOrder)) {
      node <- vOrder[i]
      
      if (i == 1) {
        predecessorNodes <- integer(0)
      } else {
        predecessorNodes <- vOrder[1:(i-1)]
      }
      cacheKey <- getCacheKey(node, predecessorNodes, vStates)
      
      logNodeStateProbability <- cache[[cacheKey]]
      if (is.null(logNodeStateProbability)) {
        logNodeStateProbability <- getLogNodeStateProbability(node, vStates, orderIndex)
        cache[[cacheKey]] <- logNodeStateProbability
      }
      
      logNodeStateProbabilities[node] <- logNodeStateProbability
    }
    
    return(exp(sum(logNodeStateProbabilities)))
  }
  
  # Mean P(X | D) over all the input orders 
  getMeanProbabilityOverOrders <- function(vStates) {
    probabilities <- sapply(1:nrow(mOrders), function(orderIndex) getProbabilityGivenOrder(vStates, orderIndex))
    return(mean(probabilities))
  }
  
  # Mean P(X | D) for all input states
  probabilities <- apply(mStates, 1, getMeanProbabilityOverOrders)
  return(probabilities)
}