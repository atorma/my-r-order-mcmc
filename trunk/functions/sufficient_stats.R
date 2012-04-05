# Returns a function f(node, parents) such that f returns the 
# sufficient statistics matrix for the given node with the given
# parents.
#
# Makes client code needing sufficient statistics cleaner and 
# allows improving the performance of computing the sufficent stats
# later.
createSufficientStatsProvider <- function(cardinalities, mObs) {
  # Updated to user R's cross-tabulation (table function)
  # after discovering it. This makes it about 32 % faster than the old 
  # version that used countStates and countSufficientStats when
  # computing stats from scratch. Another improvement is caching
  # which makes it very fast when the same (node, parents) combination
  # is requested again. This is expected when computing sufficient
  # stats over posterior order samples.
  
  cache <- hash()
  getKey <- function(node, parents) {
    return( paste(c(parents, node), collapse=' ') )
  }
  
  allFactors <- list(ncol(mObs))
  for (node in 1:ncol(mObs)) {
    allFactors[[node]] <- factor(mObs[,node], levels=1:cardinalities[node])
  }
  
  # TODO could this be even more efficient if mObs is a data.frame
  # and we use tapply() or aggregate() or aggregate.data.frame() method?
  computeSuffStats <- function(node, parents) {
    factors <- list(1 + length(parents))
    factors[[1]] <- allFactors[[node]]
    
    f <- 2
    for (parent in parents) {
      factors[[f]] <- allFactors[[parent]]
      f <- f + 1
    }
    
    xtab <- table(factors)
    suffStats <- matrix(xtab, nrow=prod(cardinalities[parents]), ncol=cardinalities[node], byrow=TRUE)
    
    return(suffStats)
  }
  
  return(function(node, parents, parentsSorted=FALSE) {    

    if (!parentsSorted) {
      parents <- sort(parents) # for consistency and cache keys.
    }
    
    key <- getKey(node, parents)
    suffStats <- cache[[key]]
    
    if (is.null(suffStats)) {
      suffStats <- computeSuffStats(node, parents)
      cache[[key]] <- suffStats
    }
    return(suffStats)
    
  })
}