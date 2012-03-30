# Returns an unnormalized local graph structure prior probability
# computed as (numberOfPossibleParents choose numberOfActualParents)^(-1)
getNChooseKPrior <- function(cPossibleParents, cActualParents) {
  cParentSets <- choose(cPossibleParents, cActualParents)
  if (cParentSets > 0) {
    return(1/cParentSets)
  } else {
    return(0)
  }
}

# Computes the term of a node in log P(G | <)
getLogLocalStructurePrior <- function(node, vParents, vOrder) {
  cPossibleParents <- max(which(vOrder == node) - 1, 0) # max handles the case when vParents is inconsistent with node and vOrder
  cActualParents <- length(vParents)
  score <- -1*lchoose(cPossibleParents, cActualParents)
  if (score < Inf) return(score)
  else return(-Inf)
}


# Computes the term of a node in log P(D | G)
getLogLocalDataLikelihood <- function(node, vParents, paramsAndData) {
  alphas <- paramsAndData$getAlphas(node, vParents)
  suffStat <- paramsAndData$countSuffStats(node, vParents)
  
  score <- sum( lgamma(rowSums(alphas)) - lgamma(rowSums(alphas) + rowSums(suffStat)) + rowSums(lgamma(alphas + suffStat) - lgamma(alphas)) )
  
  return(score)
}


# Computes the term of a node in log P(G | D, <), i.e. log score(Xi, Pa(Xi) | D, <)
# using given the local prior and data likelihood scoring functions.
#
# This is the same scoring function as when searching over network structures, 
# but the structure prior depends on node order.
#
# functLogLocalStructurePrior: f(node, parents, order)
# functLogLocalDataLikelihood: g(node, parents)
getLogLocalStructureScore <- function(node, vParents, vOrder, functLogLocalStructurePrior, functLogLocalDataLikelihood) {
  logScore <- functLogLocalStructurePrior(node, vParents, vOrder) + functLogLocalDataLikelihood(node, vParents)
  return(logScore)
}


# Precomputes family (node + parent set) specific scores for all families up to a maximum 
# family size for each node using the given scoring function functScore(node, parents).
#
# Returns an accessor function f(node, parents) such that f returns the score of the family.
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

  # assumption: same parents always sorted in same order
  getFamilyKey <- function(node, parents, parentsSorted) {
    if (!parentsSorted) {
      parents <- sort(parents)
    }
    
    # must have a separator between node numbers, otherwise 12 3 is same as 1 2 3 (duh, yeah. but it fooled me.)
    familyKey <- paste(c(parents, node), collapse=' ')
    return(familyKey)
  }
  
  setCachedScore <- function(node, parents) {
    familyKey <- getFamilyKey(node, parents, parentsSorted=TRUE)
    
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
        parentCombinations <- combinations(length(otherVars), s, otherVars) # each combination sorted ascending even if otherVars isn't
        for (p in 1:nrow(parentCombinations)) {
          parents <- parentCombinations[p,]
          setCachedScore(node, parents)
        }
      }
    }
  }
  
  cacheAccessor <- function(node, parents, parentsSorted=FALSE) {   
    familyKey <- getFamilyKey(node, parents, parentsSorted)
    score <- familyToScoreMap[[familyKey]]
    if (is.null(score)) {
      stop(paste("Score for node", node, "parents", paste(parents, collapse=' '), "not in cache", collapse=' '))
    }
    return(score)
  }
  return(cacheAccessor)
}


# Precomputes family (node + parent set) specific scores for all families up to a maximum 
# family size for each node using the given scoring function functScore(node, parents).
#
# Returns a list object with the following elements
# getScore: 
#   function f(node, parents) that returns the cached score for given family
# mSortedScores: 
#   mSortedScores[f, i] is the score of the f:th best family for node i
# sortedFamilies: 
#   list l1 such that l1[[i]] contains another list l2 where l2[[f]] is the f:th best parent set (vector) for node i
#
# Note that this doesn't scale up very well, but if the initial computation
# investment and memory limits are ok, then it's fast to get a family score.
computeFamilyScores <- function(functScore, numberOfNodes, maxParents) {
  
  nodes <- 1:numberOfNodes
  familyToScoreMap <- hash()
  mSortedScores <- matrix(NA, nrow=getNumParentSets(numberOfNodes, nodes, 0:maxParents), ncol=numberOfNodes)
  sortedFamilies <- vector("list", numberOfNodes)
  
  # assumption: same parents always sorted in same order
  getFamilyKey <- function(node, parents, parentsSorted) {
    if (!parentsSorted) {
      parents <- sort(parents)
    }
    familyKey <- paste(c(parents, node), collapse=' ')
    return(familyKey)
  }
  
  computeAndCacheScore <- function(node, parents) {
    familyKey <- getFamilyKey(node, parents, parentsSorted=TRUE)
    
    if (!is.null(familyToScoreMap[[familyKey]])) {
      stop(paste("Score for node", node, "parents", paste(parents, collapse=' '), "is being overwritten", collapse=' '))
    }
    
    score <- functScore(node, parents)
    familyToScoreMap[[familyKey]] <- score
    return(score)
  }
  
  
  for (node in nodes) {
    # Compute scores
    orderWithOurNodeLast <- c(nodes[nodes != node], node)
    parentSets <- getParentSets(node, orderWithOurNodeLast, 0:maxParents)
    
    for (p in 1:length(parentSets)) {
      parents <- parentSets[[p]]
      score <- computeAndCacheScore(node, parents)
      mSortedScores[p, node] <- score
    }
    
    # Sort scores and families
    ordering <- sort.list(mSortedScores[,node], decreasing=TRUE, method="quick", na.last=NA)
    mSortedScores[,node] <- mSortedScores[ordering, node]
    sortedFamilies[[node]] <- parentSets[ordering]
  }

  cacheAccessor <- function(node, parents, parentsSorted=FALSE) {   
    familyKey <- getFamilyKey(node, parents, parentsSorted)
    score <- familyToScoreMap[[familyKey]]
    if (is.null(score)) {
      stop(paste("Score for node", node, "parents", paste(parents, collapse=' '), "not in cache", collapse=' '))
    }
    return(score)
  }
  
  return(list(getScore=cacheAccessor, mSortedScores=mSortedScores, sortedFamilies=sortedFamilies))
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
  logLocalDataLikelihoodCached <- createFamilyScoreCache(logLocalDataLikelihood, paramsAndData$numberOfVariables, maxParents)
  
  
  logLocalStructureScore <- function(node, vParents, vOrder) {
    getLogLocalStructureScore(node, vParents, vOrder, logLocalStructurePrior, logLocalDataLikelihoodCached)
  }
  
  return(logLocalStructureScore)
}

# Factory that returns a list object with the following elements
# functLogLocalStructureScore: 
#   function f(node, parents, order) that computes the log score term of 
#   node having given parents in a given order
# mSortedScores: 
#   mSortedScores[f, i] is the log structure score of the f:th best family for node i
# sortedFamilies: 
#   list l1 such that l1[[i]] contains another list l2 where l2[[f]] is the f:th best parent set for node i
#
# 1. Allows reuse of a sufficient statistics function (e.g. a cache)
# 2. P(G | <) = P(G) = product of (numNodes-1 choose numParents)^(-1) regardless of order
# 3. Caches terms of P(G | D, <) instead of P(D | G, <) because P(G) independent of order
createCachedLogLocalStructureScoringFunction2 <- function(cardinalities, maxParents, functSufficientStats) {
  
  numNodes <- length(cardinalities)
  
  logLocalStructurePrior <- function(node, vParents) {
    -1*lchoose(numNodes-1, length(vParents))
  }
  
  paramsAndData <- list(
    getAlphas = function(node, parents) {
      getBDEuParams(node, parents, cardinalities, equivalentSampleSize=1)
    },
    countSuffStats = functSufficientStats
  )
  logLocalDataLikelihood <- function(node, vParents) {
    getLogLocalDataLikelihood(node, vParents, paramsAndData)
  }
  
  logLocalStructureScore <- function(node, vParents) {
    logLocalStructurePrior(node, vParents) + logLocalDataLikelihood(node, vParents)
  }
  logLocalStructureScoreCached <- createFamilyScoreCache(logLocalStructureScore, numNodes, maxParents)
  
  
  # The interface expected by clients has vOrder, even though we don't use it here
  return(function(node, vParents, vOrder) logLocalStructureScoreCached(node, vParents) )
}



# Computes the term of a node in log P(D | <) using the given 
# function f(node, parents, order) returning log score(Xi, Pa(Xi) | D, <)
#
# Optimally, f = functLogLocalStructureScore is a cache of 
# family scores.
getLogLocalOrderScore <- function(node, vOrder, maxParents, functLogLocalStructureScore) {
  families <- getParentSets(node, vOrder, 0:maxParents)
  familyScores <- numeric(length(families)) 
  for (f in 1:length(families)) {
    familyScores[f] <- functLogLocalStructureScore(node, families[[f]], vOrder)
  }
  
  return(getLogSumOfExponentials(familyScores))
}


# Factory that returns a function f(node, order)
# for computing node specific terms of log P(D | <).
#
# maxParents: the maximum number of parents any node can have
# functLogLocalStructureScore: function f(node, parents, order) that returns the natural logarithm of score(X, Pa(X) | D, <)
createCustomLogLocalOrderScoringFunction <- function(maxParents, functLogLocalStructureScore) {
  
  logLogLocalOrderScore <- function(node, vOrder) {
    getLogLocalOrderScore(node, vOrder, maxParents, functLogLocalStructureScore)
  }
  
  return(logLogLocalOrderScore)
}


# A convenience function to compute all terms in log P(D | <) as a vector 
# 
# vOrder: an order
# functLogLocalOrderScore: a function f(node, order) returning node specific term of log P(D | <)
getLogLocalOrderScores <- function(vOrder, functLogLocalOrderScore) {
  logScores <- numeric(length(vOrder)) 
  for (node in 1:length(vOrder)) {
    logScores[node] <- functLogLocalOrderScore(node, vOrder)
  }
  return(logScores)
}