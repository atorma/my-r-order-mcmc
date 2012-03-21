context("Construction of parent sets consistent with an order")

vOrder <- c(1, 2, 3, 4)


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

test_that("Error cases", {
  expect_that(getParentSets(node=100, vOrder, size=3), throws_error())
  expect_that(getParentSets(node=4, vOrder, size=-1), throws_error()) 
})

# Using cache
getParentSetsCached <- createCachedParentSetsProvider(length(vOrder), maxParents=3);
test_that("Parents sets of given sizes", {
  expect_that(getParentSetsCached(node=4, vOrder, size=2), equals(expected[5:7]))
  expect_that(getParentSetsCached(node=4, vOrder, size=0:3), equals(expected))
  expect_that(getParentSetsCached(node=4, vOrder, size=0:100), equals(expected)) # max parents for node 4 in order is 3
})


expected <- list(
  2,
  c(1, 2),
  c(2, 3),
  c(1, 2, 3)
)
test_that("Parents sets of given sizes containing given parent", {
  expect_that(getParentSetsIncludingParent(node=4, vOrder, size=3, parent=2), equals(expected[4]))
  expect_that(getParentSetsIncludingParent(node=4, vOrder, size=0:3, parent=2), equals(expected)) 
  expect_that(getParentSetsIncludingParent(node=4, vOrder, size=0:3, parent=100), equals(list()))
})
