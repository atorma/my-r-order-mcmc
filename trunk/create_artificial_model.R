numNodes <- 10
cardinalities <- rep(3, numNodes)
maxParents <- 3
mAdj <- generateRandomDag(numNodes, maxParents)
arrThetas <- generateMultinomialParams(mAdj, cardinalities)
mObs <- generateSamplesFromModel(mAdj, arrThetas, 500)
vObsProbs <- computeObsProbs(mAdj, arrThetas, mObs)