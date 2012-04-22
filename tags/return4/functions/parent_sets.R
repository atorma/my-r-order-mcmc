# Returns the number of possible parent sets of given node
# consistent with given order. 
# 
# node: node index
# vOrder: ordering of node indexes
# size: vector defining allowed parent set sizes (e.g. 0:3)
getNumParentSets <- function(node, vOrder, size) {
  nodePos <- which(node == vOrder)
  size <- size[size < nodePos] # can't select more parents than ahead of our node in order
  numSets <- (min(size) == 0)*1 + sum(choose(nodePos-1, size[size > 0]))
  return(numSets)
}


# Returns a list of all parent sets of a node consistent with an order and sizes. 
#
# Parent sets are sorted ascending by node number.
#
# TODO an iterator object that does not hold all parent combinations in memory
# would be more memory efficient.
getParentSets <- function(node, vOrder, size) {
  nodePos <- which(node == vOrder)
  if (length(nodePos) != 1) stop("Node not found in order")
  if (min(size) < 0) stop("Negative set size")
  
  
  possibleParents <- vOrder[1:(nodePos-1)]
  size <- size[size < nodePos] # can't select more parents than ahead of our node in order
  
  numSets <- (min(size) == 0)*1 + sum(choose(nodePos-1, size[size > 0]))
  parentSets <- vector("list", numSets)
  setIndex <- 1
  
  if (min(size) == 0) { 
    parentSets[1] <- list(integer(0)) # combinations function can't handle "n choose 0"
    size <- size[size > 0]
    setIndex <- setIndex + 1
  }
  
  for (s in size) { # now 0 < s < pos
    sets <- combinations(length(possibleParents), s, possibleParents)
    for (i in 1:nrow(sets)) {
      parentSets[[setIndex]] <- sets[i,]
      setIndex <- setIndex + 1
    }
  }
  
  return(parentSets)
}


# Returns a list of all parent sets of a node consistent with an order and sizes and 
# including given parents.
#
# Parent sets are sorted ascending by node number.
getParentSetsIncludingParent <- function(node, vOrder, size, parents) {
  if (max(size) < length(parents)) {
    stop("Allowed parent set sizes smaller than set of required parents")
  }
  if (node %in% parents) {
    stop("Node cannot be its own parent")
  }
  cParents <- length(parents)
  if (cParents == 0) {
    return(getParentSets(node, vOrder, size))
  }
  
  nodePos <- which(vOrder == node)
  
  # If the parents are not consistent with the order, there are no consistent families
  if (!all(parents %in% vOrder[1:(nodePos-1)])) {
    return(list())
  }
  
  # First all sets of requested_size - number_of_requested-parents, then add requested parents to each
  tempSize <- size - cParents
  tempSize <- tempSize[tempSize >= 0]
  parentSets <- getParentSets(node, setdiff(vOrder[1:nodePos], parents), tempSize)
  for (i in 1:length(parentSets)) {
    parentSets[[i]] <- sort(c(parentSets[[i]], parents))
  }
  return(parentSets)
}