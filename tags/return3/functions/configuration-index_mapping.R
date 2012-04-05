# Computes a unique index for a state configuration.
#
# The index is computed by variying the states 
# the quicker the smaller a node's index in the 
# order determined by vector stateConfig is.
#
# vector stateConfig: state vector indexed by node indexes
# cardinalities: vector of numbers of possible states for each node
getIndexFromConfig <- function(stateConfig, cardinalities) {
  cNodes <- length(stateConfig)
  
  if (cNodes == 0) {
    return(1)
  }
  
  if (cNodes == 1) {
    return(stateConfig)
  }
  
  if (length(cardinalities) == 1) {
    cardinalities <- rep(cardinalities, cNodes)
  }
  
  if (cNodes != length(cardinalities)) {
    stop("Number of states in configuration does not match number of possible states for each variable")
  }
  
  index <- stateConfig[1]
  stateCombs <- cardinalities[1]
  for (s in 2:cNodes) {
    index <- index + stateCombs*(stateConfig[s] - 1)
    stateCombs <- stateCombs*cardinalities[s]
  }
  return(index)
}


# Computes the configuration corresponding to a unique index assuming
# the index has been computed by variying the states 
# the quicker the smaller a node's index in the 
# order determined by a configuration vector is.
# 
# index: the index of a configuration
# possibleStates: vector of numbers of possible states for each node in the configuration represented by index
getConfigFromIndex <- function(index, cardinalities) {
  if (index > prod(cardinalities)) {
    stop("Index greater than number of configurations")
  }
  
  cVars <- length(cardinalities)
  
  stateCombs <- prod(cardinalities)
  stateConfig <- vector("integer")
  for (i in cVars:1) {
    stateCombs <- stateCombs/cardinalities[i]
    stateConfig[i] <- ((index - 1) %/% stateCombs) %% cardinalities[i] + 1
  }
  return(stateConfig)
}