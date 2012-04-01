# Data and parameters for MCMC
mObs <- as.matrix(training_data)
varNames <- names(training_data)
numNodes <- length(varNames)
cardinalities <- rep(3, numNodes)
maxParents <- 3

mTestObs <- as.matrix(test_data)

# Function for sufficient stats
functSuffStats <- createSufficientStatsProvider(cardinalities, mObs)

# Local structure score function log(score(Xi, Pa(Xi) | D, <)) 
functLogLocalStructureScore <- createLogLocalStructureScoringFunction(cardinalities, maxParents, functSuffStats)

# Cache all the local structure scores 
system.time(scoreList <- computeFamilyScores(functLogLocalStructureScore, numNodes, maxParents))

# Replace the local structure scoring function with its cached version
functLogLocalStructureScore <- function(node, parents, vOrder) {
  scoreList$getFamilyScore(node, parents) 
}

# Local order score i.e term of node in log P(D | <)
pruningDiff <- 7 # best family consistent with an order exp(pruningDiff) times more probable than worst included in computations
functLogLocalOrderScore <- function(node, vOrder) {
  getLogSumOfExponentials( scoreList$getFamiliesAndScores(node, vOrder, pruningDiff)$scores )
}


# order-MCMC
numSamples <- 25000
# Chain 1
system.time(result1 <- runOrderMCMC(numNodes, maxParents, functLogLocalOrderScore, numSamples))
plot(rowSums(result1$logScores), type="l")
# Chain 2
system.time(result2 <- runOrderMCMC(numNodes, maxParents, functLogLocalOrderScore, numSamples))
plot(rowSums(result2$logScores), type="l")


# Combine samples of two chains
sampleIdx <- seq(from=5000, to=numSamples, by=200)
samples1 <- result1$samples[sampleIdx,]
sampleLogScores1 <- result1$logScores[sampleIdx,]
samples2 <- result2$samples[sampleIdx,]
sampleLogScores2 <- result2$logScores[sampleIdx,]
samples <- rbind(samples1, samples2)
sampleLogScores <- rbind(sampleLogScores1, sampleLogScores2)

# Compute edge probabilities
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
functNodeStateProb <- createStateProbabilityFunction(cardinalities, mObs, functSuffStats=functSuffStats)
functFamiliesAndLogStructureScores <- function(node, vOrder) {
  scoreList$getFamiliesAndScores(node, vOrder, pruningDiff)
}
system.time({
  vEstimatedObsProbs <- getStateVectorProbability(mTestObs, samples, maxParents, functNodeStateProb, functFamiliesAndLogStructureScores)
})

# normalize vEstTestObsProbs
vEstTestObsProbsNorm <- vEstTestObsProbs/sum(vEstTestObsProbs)