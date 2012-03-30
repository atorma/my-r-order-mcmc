context("Construction of parent sets consistent with an order")

vOrder <- c(1, 2, 3, 4)

# Requirement is that each parent set is sorted ascending
expected <- list(
  integer(0),
  1,
  2,
  3,
  c(1, 2),
  c(1, 3),
  c(2, 3),
  c(1, 2, 3)
)
test_that("Parents sets of given sizes", {
  expect_that(getParentSets(node=4, vOrder, size=2), equals(expected[5:7]))
  expect_that(getParentSets(node=4, vOrder, size=0:3), equals(expected))
  expect_that(getParentSets(node=4, vOrder, size=0:100), equals(expected)) # max parents for node 4 in order is 3
})

test_that("Parent sets sorted ascending regardless order of possible parent nodes", {
  expect_that(getParentSets(node=4, c(2, 1, 3, 4), size=0:3), equals(expected))
})

test_that("Error cases", {
  expect_that(getParentSets(node=100, vOrder, size=3), throws_error())
  expect_that(getParentSets(node=4, vOrder, size=-1), throws_error()) 
})




expected <- list(
  2,
  c(1, 2),
  c(2, 3),
  c(1, 2, 3)
)
test_that("Parent sets of given sizes containing given parent", {
  expect_that(getParentSetsIncludingParent(node=4, vOrder, size=3, parent=2), equals(expected[4]))
  expect_that(getParentSetsIncludingParent(node=4, vOrder, size=0:3, parent=2), equals(expected)) 
  expect_that(getParentSetsIncludingParent(node=4, vOrder, size=0:3, parent=100), equals(list()))
})

test_that("Parent sets containing given parent sorted ascending regardless order of possible parent nodes", {
  expect_that(getParentSetsIncludingParent(node=4, c(2, 1, 3, 4), size=0:3, parent=2), equals(expected)) 
})


test_that("Numbers of possible parent sets computed correctly", {
  expect_that(getNumParentSets(26, 1:26, 0:3), equals(1 + choose(25, 1) + choose(25, 2) + choose(25, 3)))
})