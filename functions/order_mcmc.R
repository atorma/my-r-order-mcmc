# Computes the new order and log local order scores for each node 
# when the new order is obtained by switching two nodes in the input 
# order.
flipNodesInOrder <- function(vOrder, logLocalOrderScores, placesToSwitch, maxParents, functLogLocalOrderScore) {
  if (length(placesToSwitch) != 2) {
    stop("There must be exactly two places in order to switch")
  }
  if (min(placesToSwitch) < 1 || max(placesToSwitch) > length(vOrder)) {
    stop("Places to switch outside the range of the order")
  }
  if (placesToSwitch[1] == placesToSwitch[2]) {
    return(list(newOrder=vOrder, newLogOrderScores=logLocalOrderScores))
  }
  
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



runOrderMCMC <- function(numNodes, maxParents, functLogLocalOrderScore, maxSteps) {
  mcmcSamples <- matrix(NA, nrow=0, ncol=numNodes)
  mcmcLogScores <- matrix(NA, nrow=0, ncol=numNodes)
  recordStep <- function(vOrder, logOrderScores) {
    mcmcSamples <<- rbind(mcmcSamples, vOrder)
    mcmcLogScores <<- rbind(mcmcLogScores, logOrderScores)
  }
  
  numAccepted <- 0
  numRejected <- 0
  
  # Initialize
  vOrder <- permute(1:numNodes) # start with a random permutation
  logOrderScores <- getLogLocalOrderScores(vOrder, functLogLocalOrderScore)
  recordStep(vOrder, logOrderScores)

  step <- 2
  places <- 1:numNodes
  while (step <= maxSteps) {
    placesToSwitch <- sample(places, 2)
    switchResult <- flipNodesInOrder(vOrder, logOrderScores, placesToSwitch, maxParents, functLogLocalOrderScore)
    
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
    }
    
    if (accept) {
      logOrderScores <- switchResult$newLogOrderScores
      vOrder <- switchResult$newOrder
      numAccepted <- numAccepted + 1
    } else {
      numRejected <- numRejected + 1
    }
    
    recordStep(vOrder, logOrderScores)
    step <- step + 1
  }
  
  return(list(logScores=mcmcLogScores, samples=mcmcSamples, numAccepted=numAccepted, numRejected=numRejected))
}