# Data and parameters for MCMC
mObs <- as.matrix(training_data)
varNames <- names(training_data)
numNodes <- length(varNames)
cardinalities <- rep(3, numNodes)
maxParents <- 4

# We must compute the probabilities of these vectors
mTestObs <- as.matrix(test_data)

# Development vectors to compute probabilities for. 
mDevelObs <- as.matrix(devel_data)
# Known development vector probabilities. 
# The probabilities are based on the true model from which the training samples are from, but 
# are a biased sample!
vDevelProbs <- devel_probs[,1] # already normalized


# We use a K2 prior
functBDPriorParams <- createBDK2PriorParamsProvider(cardinalities, alpha=0.175)

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


# order-MCMC
numSamples <- 5000
# Chain 1
system.time(result1 <- runOrderMCMC(numNodes, maxParents, functLogLocalOrderScore, numSamples))
plot(rowSums(result1$logScores), type="l", col="red")
# Chain 2
system.time(result2 <- runOrderMCMC(numNodes, maxParents, functLogLocalOrderScore, numSamples))
lines(rowSums(result2$logScores), type="l", col="blue")


# Combine samples of two chains
sampleIdx <- seq(from=1000, to=numSamples, by=50)
samples1 <- result1$samples[sampleIdx,]
samples2 <- result2$samples[sampleIdx,]
samples <- rbind(samples1, samples2)


# compute the predicted test vector probabilities using a subset of the samples
sampleSubset <- samples[sample(1:nrow(samples), 5),] 
functNodeStateProb <- createStateProbabilityFunction(cardinalities, functBDPriorParams, functSuffStats)
system.time({
  vEstimatedDevelProbs <- getStateVectorProbability(mDevelObs, sampleSubset, functNodeStateProb, functFamiliesAndLogStructureScores)
})

# normalize estimated vector probabilities
vEstimatedDevelObsProbsNorm <- vEstimatedDevelProbs/sum(vEstimatedDevelProbs)

# So what's our KL-divergence?
getKLDivergence(vDevelProbs, vEstimatedDevelObsProbsNorm)