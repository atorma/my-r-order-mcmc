library("gtools")
library("hash")
library("R.utils")
library("rbenchmark")

sourceDir <- "~/Hommia/Santsaus/Project in Probabilistic Models (582637)/my-r-order-mcmc/functions"
sourceDirectory(sourceDir)

numNodes <- 6
cardinalities <- rep(3, numNodes)
maxParents <- 3
mAdj <- generateRandomDag(numNodes, maxParents)
arrThetas <- generateMultinomialParams(mAdj, cardinalities)
mObs <- generateSamplesFromModel(mAdj, arrThetas, 50)