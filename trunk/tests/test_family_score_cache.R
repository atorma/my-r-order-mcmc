context("Caching of scores for each node family")

nodes <- 1:12
logFamilyScore <- function(node, vParents) {
  if (node == 1) {    
    if (length(vParents) == 0) return(0)
    else if (all(vParents == 2)) return(-1)
    else if (all(vParents == 3)) return(-2)
    else if (all(vParents == c(2,3))) return(-3)
    else return(-4)
  } else if (node == 2) {
    if (length(vParents) == 0) return(-5)
    else if (all(vParents == 1)) return(-6)
    else if (all(vParents == 3)) return(-7)
    else if (all(vParents == c(1,3))) return(-8)
    else return(-9)
  } else if (node == 3) {
    if (length(vParents) == 0) return(-10)
    else if (all(vParents == 1)) return(-11)
    else if (all(vParents == 2)) return(-12)
    else if (all(vParents == c(1,2))) return(-13) # 123
    else if (all(vParents == 12)) return(-14) # 123 again when family intepreted as number sequence! 
    else return(-15)
  } else if (node %in% nodes) {
    return(-16)
  }
  stop(paste("Illegal node or family", node, vParents))
}

resultList <- computeFamilyScores(logFamilyScore, numberOfNodes=length(nodes), maxParents=2)


# (node, parentSet) specific scores

cache <- resultList$getFamilyScore

test_that("Expected value retrieved from cache", {
  expect_that(cache(1, integer(0)), equals(0))
  expect_that(cache(1, 2),          equals(-1))
  expect_that(cache(1, 3),          equals(-2))
  expect_that(cache(1, c(2,3)),     equals(-3))
  expect_that(cache(1, c(10,11)),   equals(-4))
  expect_that(cache(1, c(7,8)),     equals(-4))
  
  expect_that(cache(2, integer(0)), equals(-5))
  expect_that(cache(2, 1),          equals(-6))
  expect_that(cache(2, 3),          equals(-7))
  expect_that(cache(2, c(1,3)),     equals(-8))
  expect_that(cache(2, c(10,11)),   equals(-9))
  expect_that(cache(2, c(7,8)),     equals(-9))
  
  expect_that(cache(10, c(1,2)),    equals(-16))
})

test_that("Cache is not sensitive to order in parent list", {
  expect_that(cache(1, c(2,3)), equals(cache(1, c(3,2))))
  expect_that(cache(2, c(3,1)), equals(cache(2, c(1,3))))
  expect_that(cache(3, c(2,1)), equals(cache(3, c(1,2))))
})

test_that("Cache isn't fooled by node list that gives the same number sequence", {
  expect_that(cache(3, c(1,2)), equals(-13))
  expect_that(cache(3, c(12)), equals(-14))
})

test_that("Error if family is not in cache", {
  expect_that(cache(100, integer(0)), throws_error())
  expect_that(cache(3, c(1,2,4,5,6)), throws_error())
})


# Sorted scores and families

mSortedScores <- resultList$mSortedScores
sortedFamilies <- resultList$sortedFamilies

test_that("Family scores of node 1 sorted descending", {
  expect_that(mSortedScores[1,1], equals(0))
  expect_that(mSortedScores[2,1], equals(-1))
  expect_that(mSortedScores[3,1], equals(-2))           
  expect_that(mSortedScores[4,1], equals(-3))           
})

test_that("Families of node 1 sorted descending by score", {
  expect_that(sortedFamilies[[1]][[1]], equals(integer(0)))
  expect_that(sortedFamilies[[1]][[2]], equals(2))
  expect_that(sortedFamilies[[1]][[3]], equals(3))           
  expect_that(sortedFamilies[[1]][[4]], equals(c(2,3)))
})


# Getting high scoring families consistent with an order
getFamiliesAndScores <- resultList$getFamiliesAndScores

test_that("Only empty parent set can be the high scoring when node is the first in order", {
  expect_that(getFamiliesAndScores(node=1, 1:12, pruningDiff=1)$scores, equals(0) )
  expect_that(getFamiliesAndScores(node=1, 1:12, pruningDiff=1)$parentSets, equals( list(integer(0)) ) )
})

test_that("Only scores and families withing the pruning threshold returned", {
  expect_that(getFamiliesAndScores(node=1, 12:1, pruningDiff=1)$scores, equals( c(0, -1) ) )
  expect_that(getFamiliesAndScores(node=1, 12:1, pruningDiff=1)$parentSets, equals( list(integer(0), 2) ) )
})