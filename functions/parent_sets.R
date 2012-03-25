# Returns a list of all parent sets of a node consistent with an order and sizes. 
# Parent sets are sorted ascending by node number.
#
# TODO an iterator object that does not hold all parent combinations in memory
# would be more memory efficient.
getParentSets <- function(node, vOrder, size) {
  nodePos <- which(node == vOrder)
  if (length(nodePos) != 1) stop("Node not found in order")
  if (min(size) < 0) stop("Negative set size")
  
  
  size <- size[size < nodePos] # can't select more parents than ahead of our node in order
  
  numSets <- (min(size) == 0)*1 + sum(choose(nodePos-1, size[size > 0]))
  parentSets <- vector("list", numSets)
  
  
  setIndex <- 1
  
  if (min(size) == 0) { 
    parentSets[1] <- list(integer(0)) # combinations function can't handle "n choose 0"
    size <- size[size > 0]
    setIndex <- setIndex + 1
  }
  
  possibleParents <- vOrder[1:(nodePos-1)]
  for (s in size) { # now 0 < s < pos
    sets <- combinations(length(possibleParents), s, possibleParents)
    for (i in 1:nrow(sets)) {
      parentSets[[setIndex]] <- sets[i,]
      setIndex <- setIndex + 1
    }
  }
  
  return(parentSets)
}

createCachedParentSetsProvider <- function(numNodes, maxParents) {
  combsCache <- hash()
  
  getKey <- function(nPossibleParents, nActualParents) {
    paste(c(nPossibleParents, nActualParents), collapse=' ')
  }
  
  if (maxParents > 0) {
    for (nPossibleParents in 1:(numNodes-1)) {
      for (nActualParents in 1:min(nPossibleParents, maxParents)) {
        key <- getKey(nPossibleParents, nActualParents)
        combsCache[[key]] <- combinations(nPossibleParents, nActualParents, 1:nPossibleParents)
      }
    }
  }
  
  cacheAccessor <- function(node, vOrder, size) {
    nodePos <- which(node == vOrder)
    if (length(nodePos) != 1) stop("Node not found in order")
    if (min(size) < 0) stop("Negative set size")
    
    parentSets <- list()
    if (min(size) == 0) { 
      parentSets[1] <- list(integer(0)) # we didn't cache this
      size <- size[size > 0]
    }
    size <- size[size < nodePos] # can't select more parents than ahead of our node in order
    
    possibleParents <- vOrder[1:(nodePos-1)]
    for (s in size) { # now 0 < s < pos
      key <- getKey(length(possibleParents), s)
      positionsOfParents <- combsCache[[key]]
      for (i in 1:nrow(positionsOfParents)) {
        parentSets <- c(parentSets, list(vOrder[positionsOfParents[i,]]))
      }
    }
    
    return(parentSets)
  }
  return(cacheAccessor)
}

# Returns a list of all parent sets of a node consistent with an order and sizes and 
# including given parents.
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