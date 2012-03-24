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


# compute the predicted test vector probabilities using the high-scoring order
samples.logOrderScores <- rowSums(samples.logLocalOrderScores)
orderRanking <- order(samples.logOrderScores, decreasing=TRUE)
bestOrder <- samples[orderRanking[1], ]
bestOrderLogLocalOrderScores <- samples.logLocalOrderScores[orderRanking[1],]

functNodeStateProb <- createStateProbabilityFunction(cardinalities, mTestObs)
system.time({
  vEstTestObsProbs <- numeric(nrow(mTestObs))
  for (o in 1:nrow(mTestObs)) {
    vEstTestObsProbs[o] <- getStateVectorProbability(mTestObs[o,], bestOrder, maxParents, functNodeStateProb, functLogLocalStructureScore, bestOrderLogLocalOrderScores)
  }
})
# normalize vEstTestObsProbs
vEstTestObsProbsNorm <- vEstTestObsProbs/sum(vEstTestObsProbs)

