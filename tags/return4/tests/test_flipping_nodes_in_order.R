context("Test switching places of nodes in an order")

maxParents <- 2
oldOrder <- c(1, 2, 3, 4, 5, 6)
oldOrderScores <- c(0, 0, 0, 0, 0, 0)

placesToSwitch <- c(2, 5)
newOrder <- c(1, 5, 3, 4, 2, 6)


# Mock local order score for an algorithm that computes terms of all affected nodes from scratch.
# This way the subraction exp(oldLogScore) - exp(logTermsRelatedToOldParent) + - exp(logTermsRelatedToNewParent)
# cannot produce a negative number and log(score) = NaN. 

mockLogLocalOrderScore <- function(node, vOrder) {
  #cat("Local score requested for node", node, "parents", vParents, fill=T)
  
  if (!all(vOrder == newOrder)) stop("wrong order")
  
  if (node == 2) return(-2)
  else if (node == 3) return(-3)
  else if (node == 4) return(-4)
  else if (node == 5) return(-5)
  
  stop("Illegal node or parent set")
}

newOrderScores <- c(0, 
                    -2, 
                    -3, 
                    -4,
                    -5,
                    0)

result <- flipNodesInOrder(oldOrder, oldOrderScores, placesToSwitch, maxParents, mockLogLocalOrderScore)

test_that("Switching node orders computes only the necessary changes.", {
  expect_that(result$newOrder, equals(newOrder))
  expect_that(result$newLogOrderScores, equals(newOrderScores))
})

test_that("Switching node orders stops on invalid input", {
  expect_that(flipNodesInOrder(oldOrder, oldOrderScores, placesToSwitch=integer(0), maxParents=3, mockLogLocalOrderScore), throws_error())
  expect_that(flipNodesInOrder(oldOrder, oldOrderScores, placesToSwitch=1, maxParents=3, mockLogLocalOrderScore), throws_error())
  expect_that(flipNodesInOrder(oldOrder, oldOrderScores, placesToSwitch=1:3, maxParents=3, mockLogLocalOrderScore), throws_error())
  expect_that(flipNodesInOrder(oldOrder, oldOrderScores, placesToSwitch=c(100,1000), maxParents=3, mockLogLocalOrderScore), throws_error())
})

