library("gtools")
library("hash")
library("testthat")

sourceDir <- "~/Hommia/Santsaus/Project in Probabilistic Models (582637)/R/sources"
testDir <- "~/Hommia/Santsaus/Project in Probabilistic Models (582637)/R/tests"
setwd(sourceDir)
auto_test(code_path=sourceDir, test_path=testDir)