
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

getAllStates <- function(cardinalities) {
  cVars <- length(cardinalities)

  if (cVars == 0) {
    return(matrix(1))
  }

  allStates <- 1:cardinalities[1]
  if (cVars == 1) {
    return(matrix(allStates))
  }
	
  i <- 2
  while (i <= cVars) {
    allStates <- matrix(allStates, ncol=i-1)
    newStates <- {}
    for (v in 1:cardinalities[i]) {
      newStates <- rbind(newStates, cbind(allStates, v))
    }
    allStates <- matrix(newStates, ncol=i)
    i <- i+1
  }
  return(allStates)
}


# Generates a random DAG. Node indexes are sorted in total order
# such that if i is a parent of j, then i < j.
generateRandomDag <- function(numNodes, maxParents) {
  iNodes <- 1:numNodes;
  mAdj <- matrix(F, numNodes, numNodes);
  for (i in 1:numNodes) {
    cParents <- sample(0:min(i-1, maxParents), 1);
    iParents <- sample(1:(i-1), min(i-1, cParents));
    mAdj[iParents, i] <- TRUE;
  }
  return(mAdj);
}


# Draws random multinomial parameters for nodes in network
# using the equivalent sample size method. This makes the 
# distribution of multinomial parameters consistent with
# the BDE prior.
generateMultinomialParams <- function(mAdj, cardinalities, cPriorSamples=1) {
  cParents <- colSums(mAdj)
  cNodes <- nrow(mAdj)
  
  arrTheta <- array(NA, c(max(cardinalities)^max(cParents), max(cardinalities), cNodes) )
  
  for (node in 1:cNodes) {
    parents <- getParents(node, mAdj)
    cParentStates <- prod(cardinalities[parents])
    alpha <- getAlphaBDeu(cPriorSamples, cardinalities[node], cParentStates)
    arrTheta[1:cParentStates,,node] <- rdirichlet(cParentStates, rep(alpha, cardinalities[node]))
  }
  
  return(arrTheta);
}

# Draws random multinomial parameters for nodes in network using a 
# dirichlet weight parameter.
# 
# dirichletWeightParam: 1 to draw uniformly distributed parameters (default)
generateMultinomialParamsUsingWeightParam <- function(mAdj, cardinalities, dirichletWeightParam=1) {
    cParents <- colSums(mAdj)
    cNodes <- nrow(mAdj)
    
    arrTheta <- array(NA, c(max(cardinalities)^max(cParents), max(cardinalities), cNodes) )
    
    for (node in 1:cNodes) {
      parents <- getParents(node, mAdj)
      cParentStates <- prod(cardinalities[parents])
      alphas <- rep(dirichletWeightParam, cardinalities[node])
      arrTheta[1:cParentStates,,node] <- rdirichlet(cParentStates, alphas)
    }
    
    return(arrTheta);
}



# Generates random samples from the given model,
# assuming nodes are indexed in a topological order
#
# TODO allow varying node state cardinality. 
# It's now fixed to number of columns in arrTheta (max cardinality).
generateSamplesFromModel <- function(mAdj, arrTheta, size) {
  cStates <- dim(arrTheta)[2]
  cNodes <- nrow(mAdj)
  
  mSamples <- matrix(NA, size, cNodes)
  for (s in 1:size) {
    vSample <- rep(NA, cNodes)
  	for (i in 1:cNodes) {
      parentConfig <- vSample[getParents(i, mAdj)]
  	  iConfig <- getIndexFromConfig(parentConfig, cStates)
      # Draw one omultinomial sample with total number of objects = 1.
      # Only one of the states is 1 and the rest are 0. 
  	  vSample[i] <- which(rmultinom(1, 1, arrTheta[iConfig,,i]) == 1)	
  	}	
  	mSamples[s,] <- vSample
  }
  
  return(mSamples)
}


# Computes the probabilities of generated observation vectors,
# assuming nodes are indexed in a topological order.
#
# TODO allow varying node state cardinality. 
# It's now fixed to number of columns in arrTheta (max cardinality).
computeObsProbs <- function(mAdj, arrTheta, mObs) {
  cardinality <- dim(arrTheta)[2]
  
  getObsProb <- function(vStates) {
    probs <- numeric(length(vStates))
    for (node in 1:nrow(mAdj)) {
      nodeState <- vStates[node]
      parents <- getParents(node, mAdj)
      parentStates <- vStates[parents]
      
      counts <- integer(cardinality)
      counts[nodeState] <- 1
      parentConfigIndex <- getIndexFromConfig(parentStates, cardinality)
      thetas <- arrTheta[parentConfigIndex, ,node]
      probs[node] <- dmultinom(counts, size=1, prob=thetas)
    }
    return(prod(probs))
  }
  
  return( apply(mObs, 1, getObsProb) )
  
}