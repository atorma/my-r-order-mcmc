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