library("gtools")
library("hash")
library("R.utils")
library("rbenchmark")
library("testthat")

# Set the directory to the location of functions
sourceDir <- "~/Hommia/Santsaus/Project in Probabilistic Models (582637)/my-r-order-mcmc/functions"
sourceDirectory(sourceDir)

# Set the directory to the location of data sets. It will also be the working directory.
workingDir <- "~/Hommia/Santsaus/Project in Probabilistic Models (582637)/my-r-order-mcmc"
setwd(workingDir)

# Read the data in
training_data <- read.table("training_data.txt", header=TRUE)
test_data <- read.table("test_data.txt", header=TRUE)
devel_data <- read.table("devel_data.txt", header=TRUE)
devel_probs <- read.table("devel_probs.txt", header=FALSE)