
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