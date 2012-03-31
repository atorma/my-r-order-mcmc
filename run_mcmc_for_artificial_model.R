# Function for sufficient stats
functSuffStats <- createSufficientStatsProvider(cardinalities, mObs)

# Local structure score function log(score(Xi, Pa(Xi) | D, <)) 
functLogLocalStructureScore <- createLogLocalStructureScoringFunction(cardinalities, maxParents, functSuffStats)

# Cache all the local structure scores 
scoreList <- computeFamilyScores(functLogLocalStructureScore, numNodes, maxParents)

# Replace the local structure scoring function with its cached version
functLogLocalStructureScore <- function(node, parents, vOrder) {
  scoreList$getFamilyScore(node, parents) 
}

# Local order score i.e term of node in log P(D | <)
pruningDiff <- exp(7) # max difference between log(bestFamilyScore) - log(familyScore)
functLogLocalOrderScore <- function(node, vOrder) {
  getLogSumOfExponentials( scoreList$getFamiliesAndScores(node, vOrder, pruningDiff)$scores )
}

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



functNodeStateProb <- createStateProbabilityFunction(cardinalities, mObs, functSuffStats=functSuffStats)
system.time({
  vEstimatedObsProbs <- getStateVectorProbability(mUniqueObs, samples, maxParents, functNodeStateProb, functLogLocalStructureScore, sampleLogScores)
})
getKLDivergence(vObsProbs, vEstimatedObsProbs)