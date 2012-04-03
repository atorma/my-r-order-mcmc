# Development vectors to compute probabilities for. 
mDevelObs <- as.matrix(devel_data)
# Known probabilities. The probabilities are based on the true model from which the training samples are from, but 
# are a biased sample!
vDevelProbs <- devel_probs[,1] # already normalized

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