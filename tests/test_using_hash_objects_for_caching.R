context("Using hash objects for caching")

createCache <- function() {
  cache <- hash()
  
  getKey <- function(node, parents) {
    return( paste(c(parents, node), collapse=' ') )
  }
  
  
  return(function(node, parents, value=NULL) {    
    
    parents <- sort(parents) # for consistency and cache keys.
    key <- getKey(node, parents)
    cachedValue <- cache[[key]]
    
    if (is.null(cachedValue)) {
      cachedValue <- value
      cache[[key]] <- cachedValue
    }
    
    return(cachedValue)
    
  })
  
}

cache1 <- createCache()
cache2 <- createCache()

foo1 <- cache1(1, c(2, 3), 14)
foo2 <- cache2(1, c(2, 3), 99)
cache1(1, c(2, 3))
cache2(1, c(2, 3))

test_that("Caches are independent", {
  expect_that(foo1, equals(cache1(1, c(2, 3))))
  expect_that(foo2, equals(cache2(1, c(2, 3))))
})