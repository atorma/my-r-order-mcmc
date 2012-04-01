context("Mappings between states and indexes")

cardinalities <- c(2,4,3)
allStates <- getAllStates(cardinalities)

test_that("Mapping from state configuration to index is correct", {
  expect_that(getIndexFromConfig(stateConfig=integer(0), cardinalities=3), equals(1))
  expect_that(getIndexFromConfig(stateConfig=1, cardinalities=3), equals(1))
  expect_that(getIndexFromConfig(stateConfig=c(1,1,1), cardinalities=3), equals(1))
  expect_that(getIndexFromConfig(allStates[1,], cardinalities), equals(1))
  expect_that(getIndexFromConfig(allStates[18,], cardinalities), equals(18)) 
  expect_that(getIndexFromConfig(allStates[18,], c(1, 2)), throws_error())
})


# You can check the expected indexed by eyeing print(allStates)
#test_that("Indexes corresponding to state of a subset of variables computed correctly", {
  #expect_that(getIndexesFromConfig(c(2,2), cardinalities, 1:100), throws_error())
  #expect_that(getIndexesFromConfig(c(2,2), cardinalities, c(1,3)), equals(c(10,12,14,16)) )
#})


test_that("Mapping from index to state configuration is correct", {
  expect_that(getConfigFromIndex(index=1, cardinalities), equals(c(1,1,1)))
  expect_that(getConfigFromIndex(index=1000000, cardinalities), throws_error())
  expect_that(getConfigFromIndex(index=18, cardinalities), equals(allStates[18,]))
  expect_that(getConfigFromIndex(index=nrow(allStates), cardinalities), equals(allStates[nrow(allStates),]))
})

