# **Create a random network**

# Parameters for creating a random network
numNodes <- 10
cardinalities <- rep(3, numNodes)
maxParents <- 4
mAdj <- generateRandomDag(numNodes, maxParents)
arrThetas <- generateMultinomialParams(mAdj, cardinalities)
mTrainingObs <- generateSamplesFromModel(mAdj, arrThetas, 500) # training samples

mTestObs <- generateSamplesFromModel(mAdj, arrThetas, 500) # more network states for testing how well probability densities are estimated
mTestObs <- as.matrix(unique.data.frame(mTestObs))
vDevelProbs <- computeObsProbs(mAdj, arrThetas, mTestObs)



# **Learning**

# Parameters for learning
maxParents <- 4
pruningDiff <- 7 # best family consistent with an order is exp(pruningDiff) times more probable than worst included in computations

# We use a BDeu prior which corresponds to how the random network was parameterized
functBDPriorParams <- createBDeuPriorParamsProvider(cardinalities, equivalentSampleSize=1)
#functBDPriorParams <- createBDK2PriorParamsProvider(cardinalities, alpha=1)

# Function for sufficient stats
functSuffStats <- createSufficientStatsProvider(cardinalities, mTrainingObs)

# Local structure score function log(score(Xi, Pa(Xi) | D, <)) 
functLogLocalStructureScore <- createLogLocalStructureScoringFunction(cardinalities, functBDPriorParams, functSuffStats)

# Cache all the local structure scores 
system.time(scoreList <- computeFamilyScores(functLogLocalStructureScore, numNodes, maxParents))

# Replace the local structure scoring function with its cached version
functLogLocalStructureScore <- function(node, parents, vOrder) {
  scoreList$getFamilyScore(node, parents) 
}

# Local order score i.e term of node in log P(D | <)
functLogLocalOrderScore <- function(node, vOrder) {
  getLogSumOfExponentials( scoreList$getFamiliesAndScores(node, vOrder, pruningDiff)$scores )
}

# This is needed to compute edge and vector probabilities 
functFamiliesAndLogStructureScores <- function(node, vOrder) {
  scoreList$getFamiliesAndScores(node, vOrder, pruningDiff)
}


# Compute the score of a highest score order since we know one (there may be several equally good)
logScoreBestOrder <- sum(getLogLocalOrderScores(1:numNodes, functLogLocalOrderScore))



# order-MCMC
numSamples <- 5000
# Chain 1
system.time(result1 <- runOrderMCMC(numNodes, maxParents, functLogLocalOrderScore, numSamples))
plot(rowSums(result1$logScores), type="l", col="red")
# Chain 2
system.time(result2 <- runOrderMCMC(numNodes, maxParents, functLogLocalOrderScore, numSamples))
lines(rowSums(result2$logScores), type="l", col="blue")

# Combine samples from all the chains
sampleIdx <- seq(from=1000, to=numSamples, by=20)
samples1 <- result1$samples[sampleIdx,]
samples2 <- result2$samples[sampleIdx,]
samples <- rbind(samples1, samples2)


# Compute edge probabilities
system.time(mEdgeProb <- getEdgeProbabilities(samples, functFamiliesAndLogStructureScores))

# Plot ROC curve and compute AUROC
roc <- getRocCurve(mEdgeProb, mAdj)
xy <- matrix(c(0,1,0,1), 2, 2)
plot(roc, type="l", xlim=c(0,1), ylim=c(0,1), col="blue")
lines(xy, col="red")
getAuroc(roc)

# Predict probability of each unique development state vector using a subset of the mcmc samples
sampleSubset <- samples[sample(1:nrow(samples), nrow(samples)),] # now using all samples
functNodeStateProb <- createStateProbabilityFunction(cardinalities, functBDPriorParams, functSuffStats)
system.time({
  vEstimatedDevelProbs <- getStateVectorProbability(mTestObs, sampleSubset, functNodeStateProb, functFamiliesAndLogStructureScores)
})

# normalize the estimated vector probabilities
vEstimatedDevelProbsNorm <- vEstimatedDevelProbs/sum(vEstimatedDevelProbs)

# So what's our KL-divergence?
getKLDivergence(vDevelProbs, vEstimatedDevelProbsNorm)