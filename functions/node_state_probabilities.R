
# Creates a function f(node, nodeState, parents, parentStates) 
# such that that f returns the posterior probility
# of the node's state given parents and their states.
#
# Uses Bayesian Dirichlet Equivalent uniform (BDEu) prior 
# with given equivalent sample size.
createStateProbabilityFunction <- function(cardinalities, mObs, equivalentSampleSize=1) {
  getSuffStats <- createSufficientStatsHelper(cardinalities, mObs)
  
  return(function(node, nodeState, parents=integer(0), parentStates=integer(0)) {
    #cat("Requested node", node, "state", nodeState, "parents", parents, "configuration", parentStates, fill=T)
    
    if (length(parents) != length(parentStates)) stop("Inconsistent length of parent and parent state vectors")
    
    # assumption: parent states i.e. rows of matrix are in order obtained 
    # by varying states the quicker the smaller the number of the node
    suffStats <- getSuffStats(node, parents)
    alphas <- getBDEuParams(node, parents, cardinalities, equivalentSampleSize)
    posteriorCounts <- suffStats + alphas
    
    # for indexation to be insentitive to parent ordering, we must sort (parents, parentStates) 
    # by parent node number
    orderIndex <- order(parents)
    parents <- parents[orderIndex]
    parentStates <- parentStates[orderIndex]
    parentConfigIndex <- getIndexFromConfig(parentStates, cardinalities[parents])
    
    thetaExpected <- posteriorCounts[parentConfigIndex, nodeState]/sum(posteriorCounts[parentConfigIndex,])
    return(thetaExpected)
  })
}

# Computes the probability of a state of the network given a node ordering.
#
# vStates: vStates[i] is the state of node i
# vOrder: vOrder[p] is the node at position p of the order
# maxParents: maximum number of parents any node can have
# functStateProbability: function f(node, nodeState, parents, parentState) that returns the probability of a node state given a parent configuration
# functLogLocalStructureScore: function f(node, parents, order) that returns the natural logarithm of score(X, Pa(X) | D, <)
getStateVectorProbability <- function(vStates, vOrder, maxParents, functStateProbability, functLogLocalStructureScore) {
  # denominator in log scale
  functLogLocalOrderScore <- createCustomLogLocalOrderScoringFunction(maxParents, functLogLocalStructureScore)
  logLocalOrderScores <- getLogLocalOrderScores(vOrder, functLogLocalOrderScore)
  
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