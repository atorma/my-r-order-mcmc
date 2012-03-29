# Data and parameters for MCMC
mObs <- as.matrix(training_data)
varNames <- names(training_data)
numNodes <- length(varNames)
cardinalities <- rep(3, numNodes)
maxParents <- 3

# Function for log(score(Xi, Pa(Xi) | D, <))
functLogLocalStructureScore <- createCachedLogLocalStructureScoringFunction(cardinalities, mObs, maxParents)
functLogLocalOrderScore <- createCustomLogLocalOrderScoringFunction(maxParents, functLogLocalStructureScore)

# order-MCMC
numSamples <- 25000
system.time(result <- runOrderMCMC(numNodes, maxParents, functLogLocalStructureScore, numSamples))
plot(rowSums(result$logScores), type="l")

# Compute edge probabilities
sampleIdx <- seq(from=5000, to=numSamples, by=200)
samples <- result$samples[sampleIdx,]
sampleLogScores <- result$logScores[sampleIdx,]
mEdgeProb <- getEdgeProbabilities(samples, maxParents, functLogLocalStructureScore, sampleLogScores)

# List of edge probabilities
sourceNames <- character(numNodes^2 - numNodes)
targetNames <- character(numNodes^2 - numNodes)
edgeProbs <- numeric(numNodes^2 - numNodes)
row <- 1
for (i in 1:numNodes) {
  for (j in 1:numNodes) {
    if (i != j) {
      sourceNames[row] <- varNames[i]
      targetNames[row] <- varNames[j]
      edgeProbs[row] <- mEdgeProb[i,j]
      row <- row+1
    }
  }
}
edgeRanking <- data.frame(source=sourceNames, target=targetNames, probability=edgeProbs)
edgeRanking <- edgeRanking[order(edgeProbs, sourceNames, targetNames, decreasing=TRUE), ]


# compute the predicted test vector probabilities using all the samples
functNodeStateProb <- createStateProbabilityFunction(cardinalities, mObs)
system.time({
    vEstTestObsProbs <- getStateVectorProbability(mTestObs, samples, maxParents, functNodeStateProb, functLogLocalStructureScore, samples.logLocalOrderScores)
})
# normalize vEstTestObsProbs
vEstTestObsProbsNorm <- vEstTestObsProbs/sum(vEstTestObsProbs)

