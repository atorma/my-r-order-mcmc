context("Test switching places of nodes in an order")

maxParents <- 1
oldOrder <- c(5, 2, 3, 4, 1)
oldOrderScores <- c(0, 0, 0, 0, 0)
placesToSwitch <- c(2,4)
newOrder <- c(5, 4, 3, 2, 1)

# Mock scoring functions, fixed for flipping of nodes at positions 2 and 4 in the order,
# assuming maxParents.
mockLogLocalStructureScore <- function(node, vParents, vOrder) {  
  #cat("Local score requested for node", node, "parents", vParents, fill=T)
  
  # Node 2 needs to be recalculated from scratch
  if (node == 2) {
    if (length(vParents) == 0) return(-0.5)
    else if (vParents == 5) return(-5)
    else if (vParents == 4) return(-4)
    else if (vParents == 3) return(-3)
  }
  
  # For node 3, the scores associated with old parent node 2 need to be subracted
  # and scores associated with new parent node 4 need to be added
  if (node == 3) {
    if (vParents == 2) return(-2) # to be subtracted
    else if (vParents == 4) return(-4) # to be added
  }
  
  # Node 4 needs to be recalculated from scratch
  if (node == 4) {
    if (length(vParents) == 0) return(-0.5)
    else if (vParents == 5) return(-5)
  }

  stop("Illegal node or parents requested")
}

newOrderScores <- c(0, 
                    log(exp(-0.5) + exp(-5) + exp(-4) + exp(-3)), 
                    log(exp(0) - exp(-2) + exp(-4)), 
                    log(exp(-0.5) + exp(-5)), 
                    0)

result <- flipNodesInOrder(oldOrder, oldOrderScores, placesToSwitch, maxParents, mockLogLocalStructureScore)
test_that("Switching node orders computes only the necessary changes. Case max parents = 1.", {
  expect_that(result$newOrder, equals(newOrder))
  expect_that(result$newLogOrderScores, equals(newOrderScores))
})

result <- flipNodesInOrder(oldOrder, oldOrderScores, placesToSwitch=c(2,2), maxParents=3, mockLogLocalStructureScore)
test_that("Switching the same node is identity operation", {
  expect_that(result$newOrder, equals(oldOrder))
  expect_that(result$newLogOrderScores, equals(oldOrderScores))
})

test_that("Switching node orders stops on invalid input", {
  expect_that(flipNodesInOrder(oldOrder, oldOrderScores, placesToSwitch=integer(0), maxParents=3, mockLogLocalStructureScore), throws_error())
  expect_that(flipNodesInOrder(oldOrder, oldOrderScores, placesToSwitch=1, maxParents=3, mockLogLocalStructureScore), throws_error())
  expect_that(flipNodesInOrder(oldOrder, oldOrderScores, placesToSwitch=1:3, maxParents=3, mockLogLocalStructureScore), throws_error())
  expect_that(flipNodesInOrder(oldOrder, oldOrderScores, placesToSwitch=c(100,1000), maxParents=3, mockLogLocalStructureScore), throws_error())
})



# A more difficult case

maxParents <- 2
oldOrder <- c(1, 2, 3, 4, 5, 6)
oldOrderScores <- c(0, 0, 0, 0, 0, 0)

placesToSwitch <- c(2, 5)
newOrder <- c(1, 5, 3, 4, 2, 6)

mockLogLocalStructureScore <- function(node, vParents, vOrder) {
  #cat("Local score requested for node", node, "parents", vParents, fill=T)
  
  # As a switched node, must be computed from scratch
  if (node == 2) {
    if (length(vParents) == 0) return(-0.5)
    else if (all(vParents == 1)) return(-1)
    else if (all(vParents == 5)) return(-5)
    else if (all(vParents == 3)) return(-3)
    else if (all(vParents == 4)) return(-4)
    else if (all(vParents == c(1,5))) return(-1-5)
    else if (all(vParents == c(1,3))) return(-1-3)
    else if (all(vParents == c(1,4))) return(-1-4)
    else if (all(vParents == c(3,5))) return(-3-5)
    else if (all(vParents == c(4,5))) return(-4-5)
    else if (all(vParents == c(3,4))) return(-3-4)
  }
  
  # Terms corresponding to sets with old parent 2 must be subracted and 
  # terms with new parent 5 must be added
  if (node == 3) {
    if (all(vParents == 2)) return(-2) # for subtraction
    else if (all(vParents == c(1,2))) return(-1-2) # for subtraction
    else if (all(vParents == 5)) return(-5) # for adding
    else if (all(vParents == c(1,5))) return(-1-5)  # for adding
  }
  
  # Terms corresponding to sets with old parent 2 must be subtracted and 
  # terms with new parent 5 must be added
  if (node == 4) {
    if (all(vParents == 2)) return(-2)
    else if (all(vParents == c(1,2))) return(-1-2)
    else if (all(vParents == c(2,3))) return(-2-3)
    else if (all(vParents == 5)) return(-5)
    else if (all(vParents == c(1,5))) return(-1-5)
    else if (all(vParents == c(3,5))) return(-3-5)
  }
  
  # Must be computed from scratch
  if (node == 5) {
    if (length(vParents) == 0) return(-0.5)
    else if (all(vParents == 1)) return(-1)
  }
  
  stop("Illegal node or parent set")
}

newOrderScores <- c(0, 
                    log(exp(-0.5) + exp(-1) + exp(-5) + exp(-3) + exp(-4) + exp(-1-5) + exp(-1-3) + exp(-1-4) + exp(-3-5) + exp(-4-5) + exp(-3-4)), 
                    log(exp(0) - exp(-2) - exp(-1-2) + exp(-5) + exp(-1-5)), 
                    log(exp(0) - exp(-2) - exp(-1-2) - exp(-2-3) + exp(-5)+ exp(-1-5) + exp(-3-5)),
                    log(exp(-0.5) + exp(-1)),
                    0)

result <- flipNodesInOrder(oldOrder, oldOrderScores, placesToSwitch, maxParents, mockLogLocalStructureScore)
test_that("Switching node orders computes only the necessary changes. Case max parents = 2.", {
  expect_that(result$newOrder, equals(newOrder))
  expect_that(result$newLogOrderScores, equals(newOrderScores))
})


# Test the same case as above with an algorithm that computes terms of all affected nodes from scratch.
# This way the subraction exp(oldLogScore) - exp(logTermsRelatedToOldParent) + - exp(logTermsRelatedToNewParent)
# cannot produce a negative number and log(score) = NaN. The algorithm can also be faster.

mockLogLocalStructureScore <- function(node, vParents, vOrder) {
  #cat("Local score requested for node", node, "parents", vParents, fill=T)
  
  if (node == 2) {
    if (length(vParents) == 0) return(-0.5)
    else if (all(vParents == 1)) return(-1)
    else if (all(vParents == 5)) return(-5)
    else if (all(vParents == 3)) return(-3)
    else if (all(vParents == 4)) return(-4)
    else if (all(vParents == c(1,5))) return(-1-5)
    else if (all(vParents == c(1,3))) return(-1-3)
    else if (all(vParents == c(1,4))) return(-1-4)
    else if (all(vParents == c(3,5))) return(-3-5)
    else if (all(vParents == c(4,5))) return(-4-5)
    else if (all(vParents == c(3,4))) return(-3-4)
  }
  
  if (node == 3) {
    if (length(vParents) == 0) return(-0.5)
    else if (all(vParents == 1)) return(-1) 
    else if (all(vParents == 5)) return(-5)
    else if (all(vParents == c(1,5))) return(-1-5) 
  }
  
  if (node == 4) {
    if (length(vParents) == 0) return(-0.5)
    else if (all(vParents == 1)) return(-1)
    else if (all(vParents == 5)) return(-5)
    else if (all(vParents == 3)) return(-3)
    else if (all(vParents == c(1,5))) return(-1-5)
    else if (all(vParents == c(1,3))) return(-1-3)
    else if (all(vParents == c(3,5))) return(-3-5)
  }
  
  if (node == 5) {
    if (length(vParents) == 0) return(-0.5)
    else if (all(vParents == 1)) return(-1)
  }
  
  stop("Illegal node or parent set")
}

newOrderScores <- c(0, 
                    log(exp(-0.5) + exp(-1) + exp(-5) + exp(-3) + exp(-4) + exp(-1-5) + exp(-1-3) + exp(-1-4) + exp(-3-5) + exp(-4-5) + exp(-3-4)), 
                    log(exp(-0.5) + exp(-1) + exp(-5) + exp(-1-5)), 
                    log(exp(-0.5) + exp(-1) + exp(-5) + exp(-3) + exp(-1-5) + exp(-1-3) + exp(-3-5)),
                    log(exp(-0.5) + exp(-1)),
                    0)

result <- flipNodesInOrder2(oldOrder, oldOrderScores, placesToSwitch, maxParents, mockLogLocalStructureScore)
test_that("[Algorithm 2] Switching node orders computes only the necessary changes. Case max parents = 1.", {
  expect_that(result$newOrder, equals(newOrder))
  expect_that(result$newLogOrderScores, equals(newOrderScores))
})
test_that("[Algorithm 2] Switching node orders stops on invalid input", {
  expect_that(flipNodesInOrder2(oldOrder, oldOrderScores, placesToSwitch=integer(0), maxParents=3, mockLogLocalStructureScore), throws_error())
  expect_that(flipNodesInOrder2(oldOrder, oldOrderScores, placesToSwitch=1, maxParents=3, mockLogLocalStructureScore), throws_error())
  expect_that(flipNodesInOrder2(oldOrder, oldOrderScores, placesToSwitch=1:3, maxParents=3, mockLogLocalStructureScore), throws_error())
  expect_that(flipNodesInOrder2(oldOrder, oldOrderScores, placesToSwitch=c(100,1000), maxParents=3, mockLogLocalStructureScore), throws_error())
})




# A case where only scores of switched nodes need to be computed

maxParents <- 2
oldOrder <- c(1, 2, 3, 4, 5, 6)
oldOrderScores <- c(0, 0, 0, 0, 0, 0)

placesToSwitch <- c(1, 2)
newOrder <- c(2, 1, 3, 4, 5, 6)

mockLogLocalStructureScore <- function(node, vParents, vOrder) {
  # As a switched node, must be computed from scratch
  if (node == 1) {
    if (length(vParents) == 0) return(-0.5)
    else if (all(vParents == 2)) return(-2)
  }
  
  # As a switched node, must be computed from scratch
  if (node == 2) {
    if (length(vParents) == 0) return(-0.5)
  }
  
  stop("Illegal node or parent set")
}

newOrderScores <- c(log(exp(-0.5) + exp(-2)), 
                    log(exp(-0.5)), 
                    0, 
                    0,
                    0,
                    0)

result <- flipNodesInOrder(oldOrder, oldOrderScores, placesToSwitch, maxParents, mockLogLocalStructureScore)
test_that("Switching node orders computes only the necessary changes. Case only switched nodes affected.", {
  expect_that(result$newOrder, equals(newOrder))
  expect_that(result$newLogOrderScores, equals(newOrderScores))
})

result <- flipNodesInOrder2(oldOrder, oldOrderScores, placesToSwitch, maxParents, mockLogLocalStructureScore)
test_that("[Algorithm 2] Switching node orders computes only the necessary changes. Case only switched nodes affected.", {
  expect_that(result$newOrder, equals(newOrder))
  expect_that(result$newLogOrderScores, equals(newOrderScores))
})
