# Don't use, numerically unstable! In MCMC can cause reporting of unrealistic scores!
# 
# Computes the new order and log local order scores for each node 
# when the new order is obtained by switching two nodes in the input 
# order.
#
# This version outperforms flipNodesInOrder2 when the number of nodes is large.
# However, it suffers from numerical instability in a step log(exp(x1) - exp(x2) + exp(x3))
# where x:s are can be large negative.
flipNodesInOrder <- function(vOrder, logLocalOrderScores, placesToSwitch, maxParents, functLogLocalStructureScore) {
  if (length(placesToSwitch) != 2) {
    stop("There must be exactly two places in order to switch")
  }
  if (min(placesToSwitch) < 1 || max(placesToSwitch) > length(vOrder)) {
    stop("Places to switch outside the range of the order")
  }
  if (placesToSwitch[1] == placesToSwitch[2]) {
    return(list(newOrder=vOrder, newLogOrderScores=logLocalOrderScores))
  }
  
  # Convenience
  functLogLocalOrderScore <- createCustomLogLocalOrderScoringFunction(maxParents, functLogLocalStructureScore)
  
  # Ensure i < j, where i=placesToSwitch[1] and j=placesToSwitch[2]
  placesToSwitch <- sort(placesToSwitch)
  
  # When switching nodes at positions i and j, i < j, order scores of nodes at positions
  # * ... x < i or x > j stay the same
  # * ... i and j must be computed from scratch
  # * .. i < x < j is "score in old order" - "terms in old order where parent set contains i" + "terms in new order where parent set contains j"
  
  vOldOrder <- vOrder
  
  vNewOrder <- vOrder
  tempNode <- vNewOrder[placesToSwitch[1]]
  vNewOrder[placesToSwitch[1]] <- vNewOrder[placesToSwitch[2]]
  vNewOrder[placesToSwitch[2]] <- tempNode
  
  for (pos in placesToSwitch) {
    logLocalOrderScores[vNewOrder[pos]] <- functLogLocalOrderScore(vNewOrder[pos], vNewOrder)
  } 
  
  if (placesToSwitch[2] - placesToSwitch[1] > 1) { # if there's any nodes between the switched ones
    for (pos in (placesToSwitch[1] + 1):(placesToSwitch[2] - 1)) {
      node <- vNewOrder[pos]
      
      # Note: getParentSetsIncludingParent computes parent sets from scratch twice. This could be
      # optimized to use the same base parent sets (those without switched nodes) for both the old
      # and the new case---at a cost of less code reuse.
      
      vOldLogScores <- vector()
      parentSets <- getParentSetsIncludingParent(node, vOldOrder, 0:maxParents, vOldOrder[placesToSwitch[1]])
      for (s in 1:length(parentSets)) {
        vOldLogScores[s] <- functLogLocalStructureScore(node, parentSets[[s]], vOldOrder)
      }
      logScoreToSubtract <- getLogSumOfExponentials(vOldLogScores)
      
      vNewLogScores <- vector()
      parentSets <- getParentSetsIncludingParent(node, vNewOrder, 0:maxParents, vNewOrder[placesToSwitch[1]])
      for (s in 1:length(parentSets)) {
        vNewLogScores[s] <- functLogLocalStructureScore(node, parentSets[[s]], vNewOrder)
      }
      logScoreToAdd <- getLogSumOfExponentials(vNewLogScores)
      
      #cat("Node", node, "Old score", logLocalOrderScores[node], "Subtract", logScoreToSubtract, "Add", logScoreToAdd, fill=T)
      
      # For numerical stability. 
      #m <- max(c(logLocalOrderScores[node], logScoreToSubtract, logScoreToAdd))
      #print(c(logLocalOrderScores[node], logScoreToSubtract, logScoreToAdd))
      #sumExp <- exp(logLocalOrderScores[node]-m) - exp(logScoreToSubtract-m) + exp(logScoreToAdd-m)
      #cat("sum exp", sumExp, fill=T)
      # Despite our attempt, the subtraction sometimes produces a negative sum close to 0. Believing it's because of
      # numerical instability, truncate to zero in such a case.
      #sumExp <- max(0, sumExp) 
      #logLocalOrderScores[node] <- m + log(sumExp)
      
      if (logLocalOrderScores[node] == logScoreToSubtract) {
        logLocalOrderScores[node] <- logScoreToAdd
      } else if (logScoreToAdd == logScoreToSubtract) {
        # nothing changes
      } else {
        m <- max(c(logLocalOrderScores[node], logScoreToSubtract, logScoreToAdd))
        sumExp <- exp(logLocalOrderScores[node]-m) - exp(logScoreToSubtract-m) + exp(logScoreToAdd-m)
        logLocalOrderScores[node] <- m + log(sumExp)
      }
      
      #sumExp <- exp(logLocalOrderScores[node]) - exp(logScoreToSubtract) + exp(logScoreToAdd)
      #cat("sum exp", sumExp, fill=T)
      # The subtraction sometimes produces a negative sum close to 0. Believing it's because of
      # numerical instability, truncate the sum to zero in such a case.
      #sumExp <- max(0, sumExp)
      #logLocalOrderScores[node] <- log(sumExp)
      
      #cat("Node", node, "score", logLocalOrderScores[node], fill=T)
    }
  }
  
  return(list(newOrder=vNewOrder, newLogOrderScores=logLocalOrderScores))
}


# Computes the new order and log local order scores for each node 
# when the new order is obtained by switching two nodes in the input 
# order.
#
# This version outperforms flipNodesInOrder when the number of nodes is small
# and local data likelihoods are fast to compute (e.g. scores of all families cached).
# Also, this version does not suffer from numeric instability when subracting scores.
flipNodesInOrder2 <- function(vOrder, logLocalOrderScores, placesToSwitch, maxParents, functLogLocalStructureScore) {
  if (length(placesToSwitch) != 2) {
    stop("There must be exactly two places in order to switch")
  }
  if (min(placesToSwitch) < 1 || max(placesToSwitch) > length(vOrder)) {
    stop("Places to switch outside the range of the order")
  }
  if (placesToSwitch[1] == placesToSwitch[2]) {
    return(list(newOrder=vOrder, newLogOrderScores=logLocalOrderScores))
  }
  
  functLogLocalOrderScore <- createCustomLogLocalOrderScoringFunction(maxParents, functLogLocalStructureScore)
  
  # Ensure i < j, where i=placesToSwitch[1] and j=placesToSwitch[2]
  placesToSwitch <- sort(placesToSwitch)

  vNewOrder <- vOrder
  tempNode <- vNewOrder[placesToSwitch[1]]
  vNewOrder[placesToSwitch[1]] <- vNewOrder[placesToSwitch[2]]
  vNewOrder[placesToSwitch[2]] <- tempNode
  
  # When switching nodes at positions i and j, i < j, only scores of nodes at positions i <= p <= j are affected.
  affectedNodes <- vOrder[placesToSwitch[1]:placesToSwitch[2]]
  #cat("Affected", affectedNodes, fill=T)
  for (node in affectedNodes) {
    logLocalOrderScores[node] <- functLogLocalOrderScore(node, vNewOrder)
  } 
 
  return(list(newOrder=vNewOrder, newLogOrderScores=logLocalOrderScores))
}



runOrderMCMC <- function(numNodes, maxParents, functLogLocalStructureScore, maxSteps) {
  mcmcSamples <- matrix(NA, nrow=0, ncol=numNodes)
  mcmcLogScores <- matrix(NA, nrow=0, ncol=numNodes)
  recordStep <- function(vOrder, logOrderScores) {
    mcmcSamples <<- rbind(mcmcSamples, vOrder)
    mcmcLogScores <<- rbind(mcmcLogScores, logOrderScores)
  }
  
  # Initialize
  vOrder <- permute(1:numNodes) # start with a random permutation
  functLogLocalOrderScore <- createCustomLogLocalOrderScoringFunction(maxParents, functLogLocalStructureScore)
  logOrderScores <- getLogLocalOrderScores(vOrder, functLogLocalOrderScore)
  recordStep(vOrder, logOrderScores)

  step <- 2
  places <- 1:numNodes
  while (step <= maxSteps) {
    # switch places of two nodes at random and recalculate score of this suggested new order
    #cat("Current order", vOrder, fill=T)
    #cat("Current scores", logOrderScores, fill=T)
    placesToSwitch <- sample(places, 2)
    #cat("Places to switch", placesToSwitch, fill=T)
    switchResult <- flipNodesInOrder2(vOrder, logOrderScores, placesToSwitch, maxParents, functLogLocalStructureScore)
    #print("Result of flip algorithm")
    #print(switchResult)
    
    # Test without order switching "improvement"
#     placesToSwitch <- sort(placesToSwitch)
#     proposedOrder <- vOrder
#     tempNode <- proposedOrder[placesToSwitch[1]]
#     proposedOrder[placesToSwitch[1]] <- proposedOrder[placesToSwitch[2]]
#     proposedOrder[placesToSwitch[2]] <- tempNode
#     proposedLogOrderScores <- getLogLocalOrderScores(proposedOrder, functLogLocalOrderScore)
#     switchResult <- list(newOrder=proposedOrder, newLogOrderScores=proposedLogOrderScores)
    
    #print("Result without optimization")
    #print(switchResult)
    
    accept <- FALSE
    currentLogOrderScore <- sum(logOrderScores)
    proposedLogOrderScore <- sum(switchResult$newLogOrderScores)
    if (proposedLogOrderScore > currentLogOrderScore) {
      accept <- TRUE
    } else if (proposedLogOrderScore > -Inf) {
      probAccept <- exp(proposedLogOrderScore - currentLogOrderScore)
      
      if (runif(n=1, min=0, max=1) <= probAccept) {
        accept <- TRUE
      }
    } else if (proposedLogOrderScore == -Inf) {
      cat("-Inf score at step", step, fill=T)
    }
    
    if (accept) {
      logOrderScores <- switchResult$newLogOrderScores
      vOrder <- switchResult$newOrder
    }
    
    recordStep(vOrder, logOrderScores)
    step <- step + 1
  }
  
  return(list(logScores=mcmcLogScores, samples=mcmcSamples))
}