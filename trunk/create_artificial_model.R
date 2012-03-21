library("gtools")
library("hash")
library("R.utils")
library("rbenchmark")

wd <- "~/Hommia/Santsaus/Project in Probabilistic Models (582637)/R"
sourceDir <- "~/Hommia/Santsaus/Project in Probabilistic Models (582637)/R/sources"

setwd(wd)
sourceDirectory(sourceDir)

numNodes <- 8
cardinalities <- rep(3, numNodes)
maxParents <- 3
mAdj <- generateRandomDag(numNodes, maxParents)
arrThetas <- generateMultinomialParamsUsingWeightParam(mAdj, cardinalities)
mObs <- generateSamplesFromModel(mAdj, arrThetas, 300)