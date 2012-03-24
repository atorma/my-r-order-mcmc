# Returns an unnormalized local graph structure prior probability.
getNChooseKPrior <- function(cPossibleParents, cActualParents) {
  cParentSets <- choose(cPossibleParents, cActualParents)
  if (cParentSets > 0) {
    return(1/cParentSets)
  } else {
    return(0)
  }
}

# Computes the term of a node in log(P(G | <))
getLogLocalStructurePrior <- function(node, vParents, vOrder) {
  cPossibleParents <- max(which(vOrder == node) - 1, 0) # max handles the case when vParents is inconsistent with node and vOrder
  cActualParents <- length(vParents)
  score <- log(getNChooseKPrior(cPossibleParents, cActualParents))
  return(score)
}


# Computes the term of a node in log(P(D | G))
getLogLocalDataLikelihood <- function(node, vParents, paramsAndData) {
  alphas <- paramsAndData$getAlphas(node, vParents)
  suffStat <- paramsAndData$countSuffStats(node, vParents)
  
  score <- sum( lgamma(rowSums(alphas)) - lgamma(rowSums(alphas) + rowSums(suffStat)) + rowSums(lgamma(alphas + suffStat) - lgamma(alphas)) )
  
  return(score)
}


# Precomputes family (parent set) related scores for all families up to maximum family size for each node
# using the given scoring function functScore(node, parents).
#
# Returns an accessor function f to the cache such that f(node, parents) returns score
# of the family.
#
# This implementation uses a hash map where each key is a string representation of a family: 
# first parents in ascending order and then their child. The cache is only very slightly worse performing
# than using matrices with column indexes computed from parent configuration (function getIndexFromConfig)
# and saves memory.
#
# Note that this doesn't scale up very well, but if the initial computation
# investment and memory limits are ok, then it's fast to get a family score.
createFamilyScoreCache <- function(functScore, numberOfNodes, maxParents) {
  
  nodes <- 1:numberOfNodes
  familyToScoreMap <- hash()
  
  getFamilyKey <- function(node, parents) {
    # must not be sensitive to parent order
    parents <- sort(parents) 
    # must have a separator between node numbers, otherwise 12 3 is same as 1 2 3 (duh, yeah. but it fooled me.)
    familyKey <- paste(c(parents, node), collapse=' ')
    return(familyKey)
  }
  
  setCachedScore <- function(node, parents) {
    familyKey <- getFamilyKey(node, parents)
    
    if (!is.null(familyToScoreMap[[familyKey]])) {
      stop(paste("Score for node", node, "parents", paste(parents, collapse=' '), "is being overwritten", collapse=' '))
    }
    
    familyToScoreMap[[familyKey]] <- functScore(node, parents)
  }
  
  # Must handle families of one node separately because combinations can't handle the case
  for (node in nodes) {
    setCachedScore(node, integer(0))
  }
  
  if (maxParents > 0) {
    for (s in 1:maxParents) {
      for (node in nodes) {
        otherVars <- nodes[nodes != node]
        parentCombinations <- combinations(length(otherVars), s, otherVars) 
        for (p in 1:nrow(parentCombinations)) {
          parents <- parentCombinations[p,]
          setCachedScore(node, parents)
        }
      }
    }
  }
  
  cacheAccessor <- function(node, parents) {   
    familyKey <- getFamilyKey(node, parents)
    score <- familyToScoreMap[[familyKey]]
    if (is.null(score)) {
      stop(paste("Score for node", node, "parents", paste(parents, collapse=' '), "not in cache", collapse=' '))
    }
    return(score)
  }
  return(cacheAccessor)
}


# Computes log(score(Xi, Pa(Xi) | D, <)) using given local prior and data likelihood scoring functions.
#
#
# This is the same scoring function as when searching over network structures, 
# but the structure prior depends on node order.
getLogLocalStructureScore <- function(node, vParents, vOrder, functLogLocalStructurePrior, functLogLocalDataLikelihood) {
  logScore <- functLogLocalStructurePrior(node, vParents, vOrder) + functLogLocalDataLikelihood(node, vParents)
  return(logScore)
}


# Computes the term of a node in log(P(D | <)) using the given 
# function f(node, parents, order) to compute log(score(Xi, Pa(Xi) | D, <))
#
# 
getLogLocalOrderScore <- function(node, vOrder, maxParents, functLogLocalStructureScore) {
  families <- getParentSets(node, vOrder, 0:maxParents)
  familyScores <- vector() # TODO do not grow vector
  for (f in 1:length(families)) {
    familyScores[f] <- functLogLocalStructureScore(node, families[[f]], vOrder)
  }
  
  return(getLogSumOfExponentials(familyScores))
}


# Computes all terms in log(P(D | <)) as a vector using the given scoring function
getLogLocalOrderScores <- function(vOrder, functLogLocalOrderScore) {
  logScores <- vector() # TODO do not grow vector
  for (node in 1:length(vOrder)) {
    logScores[node] <- functLogLocalOrderScore(node, vOrder)
  }
  return(logScores)
}


# Factory function that returns a function f(node, parents)
# for computing log local order score.
#
# This is for comparison only. It's inefficient.
createLogLocalOrderScoringFunction <- function(cardinalities, mObs, maxParents) {
 
  logLocalStructurePrior <- function(node, vParents, vOrder) {
    getLogLocalStructurePrior(node, vParents, vOrder)
  }
  
  paramsAndData <- wrapParamsAndData(cardinalities, mObs, maxParents)
  logLocalDataLikelihood <- function(node, vParents) {
    getLogLocalDataLikelihood(node, vParents, paramsAndData)
  }
  
  logLocalStructureScore <- function(node, vParents, vOrder) {
    getLogLocalStructureScore(node, vParents, vOrder, logLocalStructurePrior, logLocalDataLikelihood)
  }
  
  logLogLocalOrderScore <- function(node, vOrder) {
    getLogLocalOrderScore(node, vOrder, maxParents, logLocalStructureScore)
  }
  
  return(logLogLocalOrderScore)
}


# Factory that returns a function f(node, parents, order)
# for computing log local data likelhood so that likelihoods
# are cached for each family of each node.
createCachedLogLocalStructureScoringFunction <- function(cardinalities, mObs, maxParents) {
  
  logLocalStructurePrior <- function(node, vParents, vOrder) {
    getLogLocalStructurePrior(node, vParents, vOrder)
  }
  
  paramsAndData <- wrapParamsAndData(cardinalities, mObs, maxParents)
  logLocalDataLikelihood <- function(node, vParents) {
    getLogLocalDataLikelihood(node, vParents, paramsAndData)
  }
  cache <- createFamilyScoreCache(logLocalDataLikelihood, paramsAndData$numberOfVariables, maxParents)
  
  
  logLocalStructureScore <- function(node, vParents, vOrder) {
    getLogLocalStructureScore(node, vParents, vOrder, logLocalStructurePrior, cache)
  }
  
  return(logLocalStructureScore)
}

# Factory that returns a function f(node, order)
# for computing log local order scores. 
#
# maxParents: the maximum number of parents any node can have
# functLogLocalStructureScore: function f(node, parents, order) that returns the natural logarithm of score(X, Pa(X) | D, <)
createCustomLogLocalOrderScoringFunction <- function(maxParents, functLogLocalStructureScore) {
  
  logLogLocalOrderScore <- function(node, vOrder) {
    getLogLocalOrderScore(node, vOrder, maxParents, functLogLocalStructureScore)
  }
  
  return(logLogLocalOrderScore)
}
