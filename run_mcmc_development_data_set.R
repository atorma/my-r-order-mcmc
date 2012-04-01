# Data and parameters for MCMC
mObs <- as.matrix(devel_data)
varNames <- names(devel_data)
numNodes <- length(varNames)
cardinalities <- rep(3, numNodes)
maxParents <- 3

mTestObs <- as.matrix(devel_data)

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

# This is needed to compute edge and vector probabilities 
functFamiliesAndLogStructureScores <- function(node, vOrder) {
  scoreList$getFamiliesAndScores(node, vOrder, pruningDiff)
}


# order-MCMC
numSamples <- 25000
# Chain 1
system.time(result1 <- runOrderMCMC(numNodes, maxParents, functLogLocalOrderScore, numSamples))
plot(rowSums(result1$logScores), type="l", col="red")
# Chain 2
system.time(result2 <- runOrderMCMC(numNodes, maxParents, functLogLocalOrderScore, numSamples))
lines(rowSums(result2$logScores), type="l", col="blue")


# Combine samples of two chains
sampleIdx <- seq(from=5000, to=numSamples, by=200)
samples1 <- result1$samples[sampleIdx,]
sampleLogScores1 <- result1$logScores[sampleIdx,]
samples2 <- result2$samples[sampleIdx,]
sampleLogScores2 <- result2$logScores[sampleIdx,]
samples <- rbind(samples1, samples2)
sampleLogScores <- rbind(sampleLogScores1, sampleLogScores2)

# Compute edge probabilities
system.time(mEdgeProb <- getEdgeProbabilities(samples, functFamiliesAndLogStructureScores))
rownames(mEdgeProb) <- varNames
colnames(mEdgeProb) <- varNames

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
rownames(edgeRanking) <- NULL

# compute the predicted test vector probabilities using all the samples
sampleSubset <- samples[sample(1:nrow(samples), 50),] # assuming all orders equally probable!
functNodeStateProb <- createStateProbabilityFunction(cardinalities, mObs, functSuffStats=functSuffStats)
system.time({
  vEstimatedObsProbs <- getStateVectorProbability(mObs, sampleSubset, functNodeStateProb, functFamiliesAndLogStructureScores)
})

# normalize vEstTestObsProbs
vEstTestObsProbsNorm <- vEstTestObsProbs/sum(vEstTestObsProbs)


# known probabilities
vKnownProbs <- as.vector(devel_probs)
vKnownProbs <- vKnownProbs/sum(vKnownProbs)
getKLDivergence(vKnownProbs, vEstimatedObsProbs)