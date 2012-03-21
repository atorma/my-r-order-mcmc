library("gtools")
library("hash")
library("R.utils")
library("rbenchmark")

# This makes MCM about 13% faster than without caching
getParentSets <- createCachedParentSetsProvider(numNodes, maxParents)

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
mEdgeProb <- matrix(0, numNodes, numNodes)
for (s in 1:nrow(samples)) { 
  mEdgeProb <- mEdgeProb + getEdgeProbabilities(samples[s,], maxParents, functLogLocalStructureScore) 
}
mEdgeProb <- mEdgeProb/nrow(samples)


# Compute ROC curve
roc <- getRocCurve(mEdgeProb, mAdj)

xy <- matrix(c(0,1,0,1), 2, 2)
plot(roc, type="l", xlim=c(0,1), ylim=c(0,1), col="blue")
lines(xy, col="red")

# Exact edge probabilities. BE CAREFUL!
mExactEdgeProb <- computeExactEdgeProbabilities(numNodes, functLogLocalStructureScore)

exactVersusMcmc <- matrix(NA, nrow=length(mEdgeProb), ncol=2)
exactVersusMcmc[,1] <- mExactEdgeProb
exactVersusMcmc[,2] <- mEdgeProb
plot(exactVersusMcmc)
lines(xy, col="red")