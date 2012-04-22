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


# Computes the term of a node in log P(D | G) using BDe scoring.
#
# node:
#   node index
# vParents:
#   indexes of parent nodes
# functGetPriorParamMatrix
#   a function f(node, parents) that returns a matrix alpha of Dirichlet parameters such that
#   alpha[j, k] is the parameter related to j:th configuration of parents and k:th node state
# functGetSuffStatsMatrix
#   a function g(node, parents) that returns a matrix N of sufficient statistics such that
#   N[j, k] is the number of times node has been in state k and when parents were in j:th configuration
getLogLocalDataLikelihood <- function(node, vParents, functGetPriorParamMatrix, functGetSuffStatsMatrix) {
  alphas <- functGetPriorParamMatrix(node, vParents)
  suffStat <- functGetSuffStatsMatrix(node, vParents)
  
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


# Factory that returns a function log score(Xi, Pa(X) | D, <)
#  
# cardinalities:
#   cardinalities[i] is the cardinality of domain of node i
# functSufficientStats:
#   f(node, parents) that returns a matrix N so that N[j, k] is the count of observations
#   where node has been observed in state k when parents were in j:th configuration
# functBDPriorParams:
#   f(node, parents) that returns matrix alpha so that alpha[j, k] is the Dirichlet parameter 
#   (prior count + 1) when node is observed in state k and its parents are in j:th configuration
#
# See getIndexFromConfig() for configuration to index mapping
createLogLocalStructureScoringFunction <- function(cardinalities, functBDPriorParams, functSufficientStats) {
  
  numNodes <- length(cardinalities)
  
  logLocalStructurePrior <- function(node, vParents) {
    -1*lchoose(numNodes-1, length(vParents))
  }
  
  logLocalDataLikelihood <- function(node, vParents) {
    getLogLocalDataLikelihood(node, vParents, functBDPriorParams, functSufficientStats)
  }
  
  logLocalStructureScore <- function(node, vParents) {
    logLocalStructurePrior(node, vParents) + logLocalDataLikelihood(node, vParents)
  }
  
  
  # The interface expected by clients has parameter vOrder, even though we don't use it here
  return(function(node, vParents, vOrder=NA) logLocalStructureScore(node, vParents) )
}


# Precomputes family (node + parent set) specific scores for all families up to a maximum 
# family size for each node using the given scoring function functScore(node, parents).
#
# Returns a list object with the following elements
# getFamilyScore: 
#   f(node, parents) that returns the cached score for the given family. 
# getFamiliesAndScores:
#   f(node, order[, diff]) that returns a list where list$scores contains cached family scores 
#   of the node consistent with the order. The scores are sorted descending and the last 
#   score is smallest score with (bestScore - score) <= diff. The default diff is Inf.
#   list$parentSets contains the parent sets producing the scores. 
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
  
  getFamilyKey <- function(node, parents, parentsSorted) {
    if (!parentsSorted) {
      parents <- sort(parents)
    }
    familyKey <- paste(c(parents, node), collapse=' ')
    return(familyKey)
  }
  
  computeAndCacheScore <- function(node, parents, parentsSorted) {
    familyKey <- getFamilyKey(node, parents, parentsSorted)
    
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
      score <- computeAndCacheScore(node, parents, parentsSorted=TRUE)
      mSortedScores[p, node] <- score
    }
    
    # Sort scores and families
    ordering <- sort.list(mSortedScores[,node], decreasing=TRUE, method="quick", na.last=NA)
    mSortedScores[,node] <- mSortedScores[ordering, node]
    sortedFamilies[[node]] <- parentSets[ordering]
  }
  
  getFamilyScore <- function(node, parents, parentsSorted=FALSE) {   
    familyKey <- getFamilyKey(node, parents, parentsSorted)
    score <- familyToScoreMap[[familyKey]]
    if (is.null(score)) {
      stop(paste("Score of node", node, "with parents", paste(parents, collapse=' '), "not in cache", collapse=' '))
    }
    return(score)
  }
  
  getFamiliesAndScores <- function(node, vOrder, pruningDiff=Inf) {
    nodePos <- which(vOrder == node) 
    if (!length(nodePos)) stop("Node not found in order or found more than once")
    
    possibleParents <- vOrder[0:(nodePos-1)] # may be empty vector
    
    allParentSetsSorted <- sortedFamilies[[node]]
    allScoresSorted <- mSortedScores[,node]
    
    # Initialize output vectors. Better to reduce size later than to grow.
    resultParentSets <- vector("list", length(allParentSetsSorted))
    resultScores <- numeric(length(allScoresSorted))
    
    bestFamilyScore <- NA
    resultFamilyIndex <- 0 # start at 0 because the result may be empty
    for (f in 1:length(allParentSetsSorted)) {
      pSet <- allParentSetsSorted[[f]]
      if (all(pSet %in% possibleParents)) {
        score <- allScoresSorted[f]
        if (is.na(bestFamilyScore)) {
          bestFamilyScore <- score
        }
        if ((bestFamilyScore - score) <= pruningDiff) {
          resultFamilyIndex <- resultFamilyIndex + 1 
          resultParentSets[[resultFamilyIndex]] <- pSet
          resultScores[[resultFamilyIndex]] <- score
        } else {
          #cat("Pruning threshold exceeded at all family index", f, ", result family index", resultFamilyIndex, fill=T)
          break; # only low scoring families left
        }
      }
    }
    
    return(list(scores=resultScores[0:resultFamilyIndex], parentSets=resultParentSets[0:resultFamilyIndex]))
  }
  
  return(list(
    getFamilyScore = getFamilyScore, 
    getFamiliesAndScores = getFamiliesAndScores,
    mSortedScores = mSortedScores, 
    sortedFamilies = sortedFamilies
    ))
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