# Function for log(score(Xi, Pa(Xi) | D, <))
functLogLocalStructureScore <- createCachedLogLocalStructureScoringFunction(cardinalities, mObs, maxParents)
functLogLocalOrderScore <- createCustomLogLocalOrderScoringFunction(maxParents, functLogLocalStructureScore)

# Score of best order
logScoreBestOrder <- sum(getLogLocalOrderScores(1:numNodes, functLogLocalOrderScore))

# order-MCMC
numSamples <- 5000
system.time(result <- runOrderMCMC(numNodes, maxParents, functLogLocalStructureScore, numSamples))
plot(rowSums(result$logScores), type="l")

# Compute edge probabilities
sampleIdx <- seq(from=1000, to=numSamples, by=100)
samples <- result$samples[sampleIdx,]
sampleLogScores <- result$logScores[sampleIdx,]
mEdgeProb <- getEdgeProbabilities(samples, maxParents, functLogLocalStructureScore, sampleLogScores) 

# Plot ROC curve
roc <- getRocCurve(mEdgeProb, mAdj)
xy <- matrix(c(0,1,0,1), 2, 2)
plot(roc, type="l", xlim=c(0,1), ylim=c(0,1), col="blue")
lines(xy, col="red")



functNodeStateProb <- createStateProbabilityFunction(cardinalities, mObs)
system.time({
vEstimatedObsProbs <- numeric(nrow(mUniqueObs))
for (o in 1:20) {
  vEstimatedObsProbs[o] <- getStateVectorProbability(mUniqueObs[o,], samples, maxParents, functNodeStateProb, functLogLocalStructureScore, sampleLogScores)
}
})
getKLDivergence(vObsProbs, vEstimatedObsProbs)