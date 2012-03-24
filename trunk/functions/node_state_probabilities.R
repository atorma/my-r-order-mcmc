
# Creates a function f(node, nodeState, parents, parentStates) 
# such that that f returns the posterior probility
# of the node's state given parents and their states.
#
# Uses Bayesian Dirichlet Equivalent uniform (BDEu) prior 
# with given equivalent sample size.
createStateProbabilityFunction <- function(cardinalities, mObs, equivalentSampleSize=1) {
  getSuffStats <- createSufficientStatsHelper(cardinalities, mObs)
  
  cache <- hash()
  getKey <- function(node, nodeState, parents, parentStates) {
    paste(c("n:", node, "s:", nodeState, "p:", parents, "c:", parentStates), collapse=" ")
  }
  
  return(function(node, nodeState, parents=integer(0), parentStates=integer(0)) {
    #cat("Requested node", node, "state", nodeState, "parents", parents, "configuration", parentStates, fill=T)
    
    if (length(parents) != length(parentStates)) stop("Inconsistent length of parent and parent state vectors")
    if (length(nodeState) != 1 || length(node) != 1) stop("Only one node and node state allowed")
    
    # for indexation to be insentitive to parent ordering, we must sort (parents, parentStates) 
    # by parent node number
    orderIndex <- order(parents)
    parents <- parents[orderIndex]
    parentStates <- parentStates[orderIndex]
    
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

# Computes the probability of a state of the network given a node ordering.
# If multple orders are inputed, then the result is the mean probability of 
# the state over the orders.
#
# vStates: vStates[i] is the state of node i
# mOrders: mOrders[o, p] is the node at position p of the o:th order
# maxParents: maximum number of parents any node can have
# functStateProbability: function f(node, nodeState, parents, parentState) that returns the probability of a node state given a parent configuration
# functLogLocalStructureScore: function f(node, parents, order) that returns the natural logarithm of score(X, Pa(X) | D, <)
# mlogLocalOrderScores (optional): mlogLocalOrderScores[o, node] is the log score contribution of node to the o:th order score
getStateVectorProbability <- function(vStates, mOrders, maxParents, functStateProbability, functLogLocalStructureScore, mlogLocalOrderScores=NULL) {
  if (!is.matrix(mOrders)) {
    mOrders <- matrix(mOrders, nrow=1)
  }
  if (length(mlogLocalOrderScores) > 0 && !is.matrix(mlogLocalOrderScores)) {
    mlogLocalOrderScores <- matrix(mlogLocalOrderScores, nrow=1)
  }
  
  if (length(mOrders) == 0) return(NA)
  if (length(mlogLocalOrderScores) > 0 && (nrow(mOrders) != nrow(mlogLocalOrderScores))) stop("Number of orders inconsistent with number of order score vectors")
  
  
  if (is.null(mlogLocalOrderScores)) {
    functLogLocalOrderScore <- createCustomLogLocalOrderScoringFunction(maxParents, functLogLocalStructureScore)
  }
  
  # Pr(X | D, <) for given order index
  getProbability <- function(orderIndex) {
    vOrder <- mOrders[orderIndex,]
    
    # denominator in log scale
    logLocalOrderScores <- mlogLocalOrderScores[orderIndex,]
    if (is.null(logLocalOrderScores)) {
      logLocalOrderScores <- getLogLocalOrderScores(vOrder, functLogLocalOrderScore)
    }  
  
    # log Pr(node | D, <) for given node 
    getLogNodeStateScore <- function(node) {
      # numerator
      parentSets <- getParentSets(node, vOrder, 0:maxParents)
      familyScores <- numeric(length(parentSets))
      for (j in 1:length(parentSets)) {
        parents <- parentSets[[j]]
        familyScores[j] <- log(functStateProbability(node, vStates[node], parents, vStates[parents])) + functLogLocalStructureScore(node, parents, vOrder)
      }
      logStateScore <- getLogSumOfExponentials(familyScores)
      
      return(logStateScore - logLocalOrderScores[node])
    }
    
    logStateScores <- numeric(length(vStates))
    for (node in 1:length(vStates)) {
      logStateScores[node] <- getLogNodeStateScore(node)
    }
    return(exp(sum(logStateScores)))
  }
  
  probabilities <- sapply(1:nrow(mOrders), getProbability)
  return(mean(probabilities))
}