
# Returns a list where list$data are counts of observed variable configurations 
# indexed using getIndexFromConfig() and list$getCount(vars, config) is a function that counts
# the number of times a given state of given variables was observed
countObservations <- function(mObs, cardinalities) {
  counts <- rep(0, prod(cardinalities))
  for (s in 1:nrow(mObs)) {
    index <- getIndexFromConfig(mObs[s,], cardinalities)
    counts[index] <- counts[index] + 1
  }
  
  return(list(counts=counts))
}


# Counts the number of times each state in a given observation matrix
# has been observed. 
#
# mObs: a matrix where each row is an observed state and columns are variables.
# returns: a matrix where each row is an observed state, columns are variables except the last column
# contains the count of the state
countStates <- function(mObs) {
  cVar <- ncol(mObs)
  cObs <- nrow(mObs)
  
  counted <- matrix(NA, nrow=1, ncol=cVar+1)  
  counted[1,] <- c(mObs[1,], 1)
  for (o in 2:cObs) {
    testMatrix <- matrix(mObs[o,], nrow=nrow(counted), ncol=cVar, byrow=T)
    rowIndex <- which(apply(counted[,1:cVar] == testMatrix, 1, all))
    if (length(rowIndex) == 0) { # new state
      counted <- rbind(counted, c(mObs[o,], 1))
    } else {
      counted[rowIndex, cVar+1] <- counted[rowIndex, cVar+1] + 1
    }
  }
  return(counted)
}

# Returns the sufficient statistics of a local multinomial distribution where
# each row contains counts given a parent configuration and each column corresponds
# to a state of their child node.
#
# mObsCounts: a matrix where each row is a state, the last column contains the number of times each
# state has been observed, and other columns contain states of variables
countSufficientStats <- function(node, vParents, cardinalities, mObsCounts) {
  cVar <- ncol(mObsCounts)-1
  if (length(cardinalities) == 1) {
    cardinalities <- rep(cardinalities, cVar) 
  }
  
  vParents <- sort(vParents) # nodes must be in ascending order
  parentCards <- cardinalities[vParents]
  cParentConfigs <- prod(parentCards)
  mCounts <- matrix(0, cParentConfigs, cardinalities[node])
  for (s in 1:nrow(mObsCounts)) {
    parentConfig <- mObsCounts[s, vParents]
    iParentConfig <- getIndexFromConfig(parentConfig, parentCards)
    iNodeState <- mObsCounts[s, node]
    mCounts[iParentConfig, iNodeState] <- mCounts[iParentConfig, iNodeState] + mObsCounts[s, cVar+1]
  }
  return(mCounts)
}


# Returns a function f initialized with given observations and 
# cardinalities such that f(node, parents) returns the 
# sufficient statistics matrix for the given node with the given
# parents.
#
# Makes client code needing sufficient statistics cleaner and 
# allows improving the efficiency of computing the sufficent stats
# later.
createSufficientStatsHelper <- function(cardinalities, mObs) {
  mObsCounts <- countStates(mObs)
  return(function(node, parents) {
    countSufficientStats(node, parents, cardinalities, mObsCounts)
  })
}