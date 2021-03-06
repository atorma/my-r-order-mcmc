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
# cardinalities: 
#   cardinalities[i] is the cardinality of the domain of node i
# functBDPriorParams: 
#   functBDPriorParams(node, parents) returns the BD prior 
#   parameter matrix for a family (node, parents)
# functSuffStats
#   functSuffStats(node, parents[, parentsSorted]) returns the 
#   sufficent stats matrix for a family (node, parents).
#   Each row corresponds to a parent configuration and rows must 
#   have the same order as in matrix from functBDPriorParams. 
#   Optional boolean parameter parentsSorted is for optimisation 
#   and indicates whether parents is already sorted ascending 
#   by node number.
# useCache (optional, default FALSE):
#   Whether to cache the results for reuse or not. Caching may 
#   cause high memory consumption with little benefit if sufficient
#   if sufficient stats are already cached.
#
# See also functions getIndexFromConfig() and getConfigFromIndex() 
# for how parent configurations are expected to be indexed.
createStateProbabilityMatrixFunction <- function(functBDPriorParams, functSuffStats, useCache=TRUE) {

  if (useCache) {
    cache <- hash()
    getKey <- function(node, parents) {
      paste(c(parents, node), collapse=" ")
    }
  }
  
  return(function(node, parents=integer(0), parentsSorted=FALSE) {
    #cat("Requested node", node, "parents", parents, fill=T)
    
    # For indexation to be insentitive to parent ordering, we must sort parents ascending,
    # unless we know they're sorted, of course
    if (!parentsSorted) {
      parents <- sort(parents)
    }
    
    if (useCache) {
      key <- getKey(node, parents)
      thetasExpected <- cache[[key]]
    } else {
      thetasExpected <- NULL
    }
    
    if (is.null(thetasExpected)) {
      # assumption: functSuffStats() returns parent states 
      # i.e. rows of matrix are in order as described above
      suffStats <- functSuffStats(node, parents, parentsSorted)
      alphas <- functBDPriorParams(node, parents)
      posteriorCounts <- suffStats + alphas
      
      thetasExpected <- posteriorCounts/rowSums(posteriorCounts)
      
      if (useCache) cache[[key]] <- thetasExpected
    }
    
    return(thetasExpected)
  })
}


# Creates a function f(node, nodeState, parents, parentStates) 
# such that that f returns the posterior probility
# of the node's state given parents and their states.
#
# cardinalities: 
#   cardinalities[i] is the cardinality of the domain of node i
# functBDPriorParams: 
#   functBDPriorParams(node, parents) returns the BD prior 
#   parameter matrix for a family (node, parents)
# functSuffStats
#   functSuffStats(node, parents[, parentsSorted]) returns the 
#   sufficent stats matrix for a family (node, parents).
#   Each row corresponds to a parent configuration and rows must 
#   have the same order as in matrix from functBDPriorParams. 
#   Optional boolean parameter parentsSorted is for optimisation 
#   and indicates whether parents is already sorted ascending 
#   by node number.
# useCache (optional, default FALSE):
#   Whether to cache the results for reuse or not. Caching may 
#   cause high memory consumption with little benefit if sufficient
#   if sufficient stats are already cached.
#
# See also functions getIndexFromConfig() and getConfigFromIndex() for 
# how parent configurations are expected to be indexed.
createStateProbabilityFunction <- function(cardinalities, functBDPriorParams, functSuffStats, useCache=FALSE) {
  
  getThetaMatrix <- createStateProbabilityMatrixFunction(functBDPriorParams, functSuffStats, useCache)
  
  return(function(node, nodeState, parents=integer(0), parentStates=integer(0), parentsSorted=FALSE) {
    #cat("Requested node", node, "state", nodeState, "parents", parents, "configuration", parentStates, fill=T)
    
    if (length(parents) != length(parentStates)) stop("Inconsistent length of parent and parent state vectors")
    if (length(nodeState) != 1 || length(node) != 1) stop("Only one node and node state allowed")

    # for indexation to be insentitive to parent ordering, we must sort (parents, parentStates) by parent node number
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

# Computes the probabilities of state vectors given node orderings.
# If multiple orders are inputed, then the result is the mean probability of 
# the state over the orders.
#
# mStates: 
#   vStates[s, i] is the state of node i in state vector s
# mOrders: 
#   mOrders[o, p] is the node at position p of the o:th order
# functNodeStateProbability: 
#   function f(node, nodeState, parents, parentState) that returns the 
#   probability of a node state given a parent configuration
# functFamiliesAndLogStructureScores:
# useCache (default FALSE):
#   Whether to cache results for each (node, nodeState predecessors, predecessorStates) 
#   This reduce computation time the more the similar input states and orders are. It can
#   also be a huge memory hog.
getStateVectorProbability <- function(mStates, mOrders, functNodeStateProbability, functFamiliesAndLogStructureScores, useCache=FALSE) {
  if (!is.matrix(mStates)) {
    mStates <- matrix(mStates, nrow=1)
  }
  if (!is.matrix(mOrders)) {
    mOrders <- matrix(mOrders, nrow=1)
  }
  
  if (length(mOrders) == 0) return(NA)
  if (ncol(mStates) != ncol(mOrders)) stop("state vectors must have as many elements as order vectors") 
    
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
  if (useCache) {
    cache <- hash()
  
    getCacheKey <- function(node, predecessorNodes, vStates) {
      nodeState <- vStates[node]
      predecessorNodes <- sort(predecessorNodes)
      predecessorStates <- vStates[predecessorNodes]
      
      return(paste(c("n:", node, "ns:", nodeState, "p:", predecessorNodes, "ps:", predecessorStates), collapse=" "))
    }
  }
  
  # log P(Xi | D, <) for given node in a given state vector.
  # This is an expensive operation.
  getLogNodeStateProbability <- function(node, vStates, orderIndex) {
    vOrder <- mOrders[orderIndex,]
    
    parentSetsAndLogStructureScores <- functFamiliesAndLogStructureScores(node, vOrder)
    
    # numerator in log scale
    # an expensive operation
    parentSets <- parentSetsAndLogStructureScores$parentSets # all parent sets are sorted ascending
    logStructureScores <- parentSetsAndLogStructureScores$scores
    
    familyScores <- numeric(length(parentSets))
    for (j in 1:length(parentSets)) {
      parents <- parentSets[[j]]
      familyScores[j] <- log(functNodeStateProbability(node, vStates[node], parents, vStates[parents], parentsSorted=TRUE)) + logStructureScores[j]
    }
    logWeighedProbabilities <- getLogSumOfExponentials(familyScores)
    
    # denominators in log scale
    logSumOfWeights <- getLogSumOfExponentials(logStructureScores)
    
    #cat("Computed P(node=", node, "state=", vStates[node], "X=", vStates, "order=", vOrder, fill=T)
    
    return(logWeighedProbabilities - logSumOfWeights)
  }
  
  #For printing only
  stateIndex <- 1
  
  # P(X | D, <) for given order index an state vector
  getProbabilityGivenOrder <- function(vStates, orderIndex) {
    vOrder <- mOrders[orderIndex,]
 
    logNodeStateProbabilities <- numeric(length(vStates))
    for (i in 1:length(vOrder)) {
      node <- vOrder[i]
      
      if (useCache) {
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
        
      } else {
        logNodeStateProbability <- getLogNodeStateProbability(node, vStates, orderIndex)
      }
      
      logNodeStateProbabilities[node] <- logNodeStateProbability
    }
    
    #cat("Computed probability of state ", stateIndex, "given order ", orderIndex, fill=T)
    
    return(exp(sum(logNodeStateProbabilities)))
  }
  
  # Mean P(X | D) over all the input orders 
  getMeanProbabilityOverOrders <- function(vStates) {
    probabilities <- sapply(1:nrow(mOrders), function(orderIndex) getProbabilityGivenOrder(vStates, orderIndex))
    if (stateIndex %% 100 == 0) {
      cat("Computed mean probability of state", stateIndex, "over all orders", fill=T)
    }
    stateIndex <<- stateIndex + 1
    return(mean(probabilities))
  }
  
  # Mean P(X | D) for all input states
  probabilities <- apply(mStates, 1, getMeanProbabilityOverOrders)
  return(probabilities)
}