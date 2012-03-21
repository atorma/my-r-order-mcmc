library("MCMCpack")
library("gtools")


sampleDiscrete <- function(prob, size) {
  sample(x=1:length(prob), size=size, replace=TRUE, prob=prob)
}

# Generates a random DAG. Node indexes are sorted in total order
# such that if i is a parent of j, then i < j.
genRandomDag <- function(numNodes, maxParents) {
  iNodes <- 1:numNodes;
  mAdj <- matrix(F, numNodes, numNodes);
  for (i in 1:numNodes) {
    cParents <- sample(0:min(i-1, maxParents), 1);
    iParents <- sample(1:(i-1), min(i-1, cParents));
    mAdj[iParents, i] <- TRUE;
  }
  return(mAdj);
}

getAlphaBDEu <- function(cPriorSamples, cPossibleStates, cParents) {
  alpha <- cPriorSamples/(cPossibleStates^(1+cParents));
  alpha[alpha == Inf] <- 1;
  return(alpha)
}

genMultinomParams <- function(mAdj, cPriorSamples, cStates) {
  cParents <- colSums(mAdj);
  
  cNodes = nrow(mAdj);
  arrTheta <- array(NA, c(cStates^max(cParents), cStates, cNodes) ); 
  cParentStates <- cStates^cParents;
  for (i in 1:cNodes) {
	alpha <- getAlphaBDEu(cPriorSamples, cStates, cParents[i])
    arrTheta[1:cParentStates[i],,i] <- rdirichlet(cParentStates[i], rep(alpha, cStates));
  }
	
  return(arrTheta);
}	

getChildren <- function(iNode, mAdj) { 
  which(mAdj[iNode,]);
}

getParents <- function(iNode, mAdj) { 
  which(mAdj[,iNode]);
}

getLeaves <- function(mAdj) { 
  which(apply(mAdj, 1, max) == 0); 
}

getRoots <- function(mAdj) { 
  which(apply(mAdj, 2, max) == 0);
}

getAllStates <- function(iNodes, cStates) {
  if (length(iNodes) == 0) {
    return(matrix(1))
  }

  allStates <- 1:cStates
  if (length(iNodes) == 1) {
    return(matrix(allStates))
  }
	
  i <- 2
  while (i <= length(iNodes)) {
    allStates <- matrix(allStates, ncol=i-1)
    newStates <- {}
    for (v in 1:cStates) {
      newStates <- rbind(newStates, cbind(allStates, v))
    }
    allStates <- matrix(newStates, ncol=i)
    i <- i+1
  }
  return(allStates)
}

# In matrix of all parent configurations, states vary the quicker
# the smaller the parent's index in the total order is.
getParentConfigIndex <- function(parentStates, cStates) {
  cParentStates <- length(parentStates)

  if (cParentStates == 0) {
    return(1)
  }

  index <- 1
  for (s in 1:cParentStates) {
    index <- index + cStates^(s - 1) * (parentStates[s] - 1)
  }
  return(index)
}

# Generates random samples from the given model,
# assuming nodes are indexed in a topological order
genModelSamples <- function(mAdj, arrTheta, size) {
  cStates <- dim(arrTheta)[2]
  cNodes <- nrow(mAdj)
  
  mSamples <- matrix(NA, size, cNodes)
  for (s in 1:size) {
    vSample <- rep(NA, cNodes)
	for (i in 1:cNodes) {
	  iConfig <- getParentConfigIndex(vSample[getParents(i, mAdj)], cStates)
	  vSample[i] <- which(rmultinom(1, 1, arrTheta[iConfig,,i]) == 1)	
	}	
	mSamples[s,] <- vSample
  }
  
  return(mSamples)
}

# mObs is a matrix where each column contains observed states of a node.
# Retutrns a matrix where rows are nodes and each column contains 
# the number of times a state was observed in mObs.
countObservedStates <- function(mObs, cPossibleStates) {
  cNodes <- ncol(mObs)
  cSamples <- nrow(mObs)
  mCounts <- matrix(0, cNodes, cPossibleStates)
  for (s in 1:cSamples) {
    mIndex <- cbind(1:cNodes, mObs[s,])
    mCounts[mIndex] <- mCounts[mIndex] + 1
  }
  return(mCounts)
}

getNChooseKPrior <- function(cPossibleParents, cActualParents) {
  1/choose(cPossibleParents, cActualParents)
}

# vOrder is an a vector where elements are nodes
getLogLocalStructurePrior <- function(node, vParents, vOrder) {
  cPossibleParents <- which(vOrder == node) - 1
  cActualParents <- length(vParents)
  score <- log(getNChooseKPrior(cPossibleParents, cActualParents))
  return(score)
}

# Sufficient statistics of a local multinomial distribution.
# Each row contains counts given a parent configuration.
countSufficientStats <- function(node, vParents, cPossibleStates, mObs) {
  vParents <- sort(vParents) # nodes must be in ascending order
  cParentConfigs <- cPossibleStates^(length(vParents))
  mCounts <- matrix(0, cParentConfigs, cPossibleStates)
  for (s in 1:nrow(mObs)) {
    parentConfig <- mObs[s, vParents]
    iParentConfig <- getParentConfigIndex(parentConfig, cPossibleStates)
    iNodeState <- mObs[s, node]
    mCounts[iParentConfig, iNodeState] <- mCounts[iParentConfig, iNodeState] + 1
  }
  return(mCounts)
}


# node and vParents are the names of nodes. 
# vOrder is a vector of node names in a total order.
# mObs is a matrix of observed vectors on rows.
getLogLocalStructureBelief <- function(node, cPossibleStates, vParents, mObs) {
  # BDEu prior for parameters with equivalent sample size 1
  # This is the scalar hyperparameter value for each (parent config, node state)
  cActualParents <- length(vParents)
  alpha <- getAlphaBDEu(1, cPossibleStates, cActualParents)
  
  suffStat <- countSufficientStats(node, vParents, cPossibleStates, mObs)
  sumsOverNodeStates <- rowSums(suffStat)
  score <- sum( lgamma(cPossibleStates*alpha) - lgamma(cPossibleStates*alpha + sumsOverNodeStates) + rowSums(lgamma(alpha + suffStat) - lgamma(alpha)) )

  return(score)
}

# Computes log score of a node given parents
getLogLocalScore <- function(node, cPossibleStates, vParents, mObs, vOrder) {
  logScore <- getLogLocalStructurePrior(node, vParents, vOrder) + getLogLocalStructureBelief(node, cPossibleStates, vParents, mObs)
  return(logScore)
}

# Computes log score of a node over all parent sets consistent with given order
getLogLocalOrderScore <- function(node, vOrder, cPossibleStates, maxParents, mObs) {
  # When summing over probabilities (marginalizing), we can't use log scores. 
  # With large number of observations log scores are large negative and cause underflow
  # when computing exp(x). Transformation 
  # log(sum(exp(x)) = m + log(sum(exp(x - m))),
  # where m = min(x), improves numerical stability.
  
  # Initialise with family size 0.
  orderOfNode <- which(vOrder == node)
  vFamilyLogScores <- getLogLocalScore(node, cPossibleStates, integer(0), mObs, vOrder)
  if (orderOfNode > 1) {
	  for (k in 1:min(orderOfNode-1, maxParents)) {
		mFamiliesOfSizeK <- combinations(orderOfNode-1, k, vOrder[1:(orderOfNode-1)])
		for (f in 1:nrow(mFamiliesOfSizeK)) {
		   vFamilyLogScores[length(vFamilyLogScores)+1] <- getLogLocalScore(node, cPossibleStates, mFamiliesOfSizeK[f,], mObs, vOrder)
		}	  
	  }
  }
  m <- min(vFamilyLogScores)

  logScore <- m + log(sum(exp(vFamilyLogScores - m)))
  return(logScore)
}
  
# Computes getLogLocalOrderScore for each node in an order
getLogLocalOrderScores <- function(vOrder, cPossibleStates, maxParents, mObs) {
  logScores <- vector("numeric")
  for (i in 1:length(vOrder)) {
    logScores[i] <- getLogLocalOrderScore(vOrder[i], vOrder, cPossibleStates, maxParents, mObs)
  }
  return(logScores)
}

getLogOrderScore <- function(vOrder, cPossibleStates, maxParents, mObs) {
  logScore <- getLogLocalScore(vOrder[1], cPossibleStates, integer(0), mObs, vOrder)
  for (i in 2:length(vOrder)) {
    # When summing over probabilities, we can't use log scores. This is bad for numerical stability.
	# Transformation log(sum(exp(x)) = m + log(sum(exp(x - m))),
	# where m = min(x), improves numerical stability
  
    vFamilyLogScores <- getLogLocalScore(vOrder[i], cPossibleStates, integer(0), mObs, vOrder)
    for (k in 1:min(i-1, maxParents)) {
      mFamiliesOfSizeK <- combinations(i-1, k, vOrder[1:(i-1)])
      for (f in 1:nrow(mFamiliesOfSizeK)) {
         vFamilyLogScores[length(vFamilyLogScores)+1] <- getLogLocalScore(vOrder[i], cPossibleStates, mFamiliesOfSizeK[f,], mObs, vOrder)
      }	  
    }
	m <- min(vFamilyLogScores)
	
	logScore <- logScore + m + log(sum(exp(vFamilyLogScores - m)))
  }
  return(logScore)
}


cPossibleStates <- 3
mAdj <- matrix(FALSE, 3, 3)
mAdj[1,2] <- TRUE
mAdj[1,3] <- TRUE

arrTheta <- array(NA, c(3, 3, 3))
arrTheta[1,,1] <- c(0.05, 0.75, 0.2)
arrTheta[1,,2] <- c(0.25, 0.25, 0.5)
arrTheta[2,,2] <- c(0.40, 0.40, 0.2)
arrTheta[3,,2] <- c(0.80, 0.15, 0.05)
arrTheta[1,,3] <- c(0.55, 0.05, 0.4)
arrTheta[2,,3] <- c(0.60, 0.10, 0.3)
arrTheta[3,,3] <- c(0.20, 0.20, 0.60)

mObs <- matrix(NA, 5, 3)
mObs[1,] <- c(3, 1, 3)
mObs[2,] <- c(2, 3, 3)
mObs[3,] <- c(2, 1, 3)
mObs[4,] <- c(2, 2, 1)
mObs[5,] <- c(2, 1, 3)

getLogOrderScore(c(1,2,3), cPossibleStates, 2, mObs)