numNodes <- 10
cardinalities <- rep(3, numNodes)
maxParents <- 3
mAdj <- generateRandomDag(numNodes, maxParents)
arrThetas <- generateMultinomialParamsUsingWeightParameter(mAdj, cardinalities)
mObs <- generateSamplesFromModel(mAdj, arrThetas, 100)