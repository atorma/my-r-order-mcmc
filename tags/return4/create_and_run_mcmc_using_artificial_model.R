# **Create a random network**

# Parameters for creating a random network
numNodes <- 10
cardinalities <- rep(3, numNodes)
maxParents <- 4
mAdj <- generateRandomDag(numNodes, maxParents)
arrThetas <- generateMultinomialParams(mAdj, cardinalities)
mTrainingObs <- generateSamplesFromModel(mAdj, arrThetas, 500) # draw training samples

mTestObs <- generateSamplesFromModel(mAdj, arrThetas, 2500) # further samples for testing KL-divergence
mTestObs <- as.matrix(unique.data.frame(mTestObs))
mTestObs <- mTestObs[1:min(200, nrow(mTestObs)),]
vTestProbs <- computeStateProbs(mAdj, arrThetas, mTestObs)



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

# Plot ROC curve and compute AUROC
roc <- getRocCurve(mEdgeProb, mAdj)
xy <- matrix(c(0,1,0,1), 2, 2)
plot(roc, type="l", xlim=c(0,1), ylim=c(0,1), col="blue")
lines(xy, col="red")
getAuroc(roc)

# Predict probability of each test state vector using a subset of the mcmc samples
sampleSubset <- samples[sample(1:nrow(samples), 10),] 
functNodeStateProb <- createStateProbabilityFunction(cardinalities, functBDPriorParams, functSuffStats)
system.time({
  vEstimatedTestProbs <- getStateVectorProbability(mTestObs, sampleSubset, functNodeStateProb, functFamiliesAndLogStructureScores)
})

# So what's our KL-divergence?
getKLDivergence(vTestProbs, vEstimatedTestProbs)