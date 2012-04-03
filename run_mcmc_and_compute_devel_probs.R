# Data and parameters for MCMC
mObs <- as.matrix(training_data)
varNames <- names(training_data)
numNodes <- length(varNames)
cardinalities <- rep(3, numNodes)

maxParents <- 3
equivalentSampleSize <- 10

# Development vectors to compute probabilities for. 
mDevelObs <- as.matrix(devel_data)
# Known development vector probabilities. 
# The probabilities are based on the true model from which the training samples are from, but 
# are a biased sample!
vDevelProbs <- devel_probs[,1] # already normalized


# Function for sufficient stats
functSuffStats <- createSufficientStatsProvider(cardinalities, mObs)

# Local structure score function log(score(Xi, Pa(Xi) | D, <)) 
functLogLocalStructureScore <- createLogLocalStructureScoringFunction(cardinalities, functSuffStats, equivalentSampleSize)

# Cache all the local structure scores 
system.time(scoreList <- computeFamilyScores(functLogLocalStructureScore, numNodes, maxParents))

# Replace the local structure scoring function with its cached version
functLogLocalStructureScore <- function(node, parents, vOrder) {
  scoreList$getFamilyScore(node, parents) 
}

# For each node, best family consistent with an order is exp(pruningDiff) 
# times more probable than worst included in computations
pruningDiff <- 7 

# Local order score i.e term of node in log P(D | <)
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

sampleIdx <- seq(from=5000, to=numSamples, by=200)
samples <- result1$samples[sampleIdx,]


# compute the predicted test vector probabilities using MCMC samples from the training model
sampleSubset <- samples[sample(1:nrow(samples), 10),] # assuming all sampled orders scored equally!
functNodeStateProb <- createStateProbabilityFunction(cardinalities, mObs, functSuffStats=functSuffStats)
system.time({
  vEstimatedObsProbs <- getStateVectorProbability(mDevelObs, sampleSubset, functNodeStateProb, functFamiliesAndLogStructureScores)
})

# normalize estimated vector probabilities
vEstimatedObsProbsNorm <- vEstimatedObsProbs/sum(vEstimatedObsProbs)

# So what's our KL-divergence
getKLDivergence(vDevelProbs, vEstimatedObsProbsNorm)