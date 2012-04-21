# Input data 
mObs <- as.matrix(training_data) # training samples
varNames <- names(training_data)
numNodes <- length(varNames)
cardinalities <- rep(3, numNodes)

mTestObs <- as.matrix(test_data) # We must compute the probabilities of these vectors

mDevelObs <- as.matrix(devel_data) # We know the probabilities of these vectors. They are a biased sample, not to be used for learning!
vDevelProbs <- devel_probs[,1] # Known vector probsl, already normalized


# Scoring and parameters

maxParents <- 4

# We use a K2 prior
functBDPriorParams <- createBDK2PriorParamsProvider(cardinalities, alpha=0.25)

# Function for sufficient stats
functSuffStats <- createSufficientStatsProvider(cardinalities, mObs)

# Local structure score function log(score(Xi, Pa(Xi) | D, <)) 
functLogLocalStructureScore <- createLogLocalStructureScoringFunction(cardinalities, functBDPriorParams, functSuffStats)

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


# Run the order-MCMC using four independent chains. They seem to converge after 1000 steps.

numSamples <- 5000
# Chain 1
system.time(result1 <- runOrderMCMC(numNodes, maxParents, functLogLocalOrderScore, numSamples))
plot(rowSums(result1$logScores), type="l", col="red")
# Chain 2
system.time(result2 <- runOrderMCMC(numNodes, maxParents, functLogLocalOrderScore, numSamples))
lines(rowSums(result2$logScores), type="l", col="blue")
# Chain 3
system.time(result3 <- runOrderMCMC(numNodes, maxParents, functLogLocalOrderScore, numSamples))
lines(rowSums(result3$logScores), type="l", col="green")
# Chain 4
system.time(result4 <- runOrderMCMC(numNodes, maxParents, functLogLocalOrderScore, numSamples))
lines(rowSums(result4$logScores), type="l", col="yellow")

# Combine samples of the chains. Burn-in 1000 samples, thinning every 100 samples
sampleIdx <- seq(from=1000, to=numSamples, by=100)
samples1 <- result1$samples[sampleIdx,]
samples2 <- result2$samples[sampleIdx,]
samples3 <- result3$samples[sampleIdx,]
samples4 <- result4$samples[sampleIdx,]
samples <- rbind(samples1, samples2, samples3, samples4)

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

# compute the predicted test vector probabilities using a subset of the samples
sampleSubset <- samples[sample(1:nrow(samples), 20),] 
functNodeStateProb <- createStateProbabilityFunction(cardinalities, functBDPriorParams, functSuffStats)
system.time({
  vEstimatedTestProbs <- getStateVectorProbability(mTestObs, sampleSubset, functNodeStateProb, functFamiliesAndLogStructureScores)
})

# normalize the state vector probabilities
vEstimatedTestProbsNorm <- vEstimatedTestProbs/sum(vEstimatedTestProbs)