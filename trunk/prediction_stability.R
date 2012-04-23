# Compute edge probabilities 1
system.time(mEdgeProb1 <- getEdgeProbabilities(samples1, functFamiliesAndLogStructureScores))
rownames(mEdgeProb1) <- varNames
colnames(mEdgeProb1) <- varNames

# List of edge probabilities 1
sourceNames <- character(numNodes^2 - numNodes)
targetNames <- character(numNodes^2 - numNodes)
edgeProbs1 <- numeric(numNodes^2 - numNodes)
row <- 1
for (i in 1:numNodes) {
  for (j in 1:numNodes) {
    if (i != j) {
      sourceNames[row] <- varNames[i]
      targetNames[row] <- varNames[j]
      edgeProbs1[row] <- mEdgeProb1[i,j]
      row <- row+1
    }
  }
}


# Compute edge probabilities 2
system.time(mEdgeProb2 <- getEdgeProbabilities(samples2, functFamiliesAndLogStructureScores))
rownames(mEdgeProb2) <- varNames
colnames(mEdgeProb2) <- varNames

# List of edge probabilities 2
sourceNames <- character(numNodes^2 - numNodes)
targetNames <- character(numNodes^2 - numNodes)
edgeProbs2 <- numeric(numNodes^2 - numNodes)
row <- 1
for (i in 1:numNodes) {
  for (j in 1:numNodes) {
    if (i != j) {
      sourceNames[row] <- varNames[i]
      targetNames[row] <- varNames[j]
      edgeProbs2[row] <- mEdgeProb2[i,j]
      row <- row+1
    }
  }
}

plot(jitter(edgeProbs1), jitter(edgeProbs2), xlab="arc probs run 1", ylab="arc probs run 2")
xy <- matrix(c(0,1,0,1), 2, 2)
lines(xy)

# compute the predicted test vector probabilities using a subset of the samples from run 1
sampleSubset1 <- samples[sample(1:nrow(samples1), 5),] 
functNodeStateProb <- createStateProbabilityFunction(cardinalities, functBDPriorParams, functSuffStats)
system.time({
  vEstimatedTestProbs1 <- getStateVectorProbability(mTestObs, sampleSubset1, functNodeStateProb, functFamiliesAndLogStructureScores)
})

# compute the predicted test vector probabilities using a subset of the samples from run 2
sampleSubset2 <- samples[sample(1:nrow(samples2), 5),] 
functNodeStateProb <- createStateProbabilityFunction(cardinalities, functBDPriorParams, functSuffStats)
system.time({
  vEstimatedTestProbs2 <- getStateVectorProbability(mTestObs, sampleSubset2, functNodeStateProb, functFamiliesAndLogStructureScores)
})