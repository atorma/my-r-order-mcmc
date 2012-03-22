context("Computation of exact edge probabilities")

numNodes <- 3
maxParents <- 1
mOrders <- matrix(c(1, 2, 3
                    2, 1, 3
                    3, 2, 1
                    1, 3, 2
                    2, 3, 1
                    3, 1, 2))

mockLogLocalStructureScore <- function(node, vParents, vOrder) {
  #cat("Node", node, "parents", vParents, "order", vOrder, fill=T)
  
  if (all(vOrder == mOrders[1,])) {
    
    if (node == 2) {
      if (length(vParents) == 0) return(-1)
      else if (vParents == 1) return(-2)    
      
    } else if (node == 3) {
      if (length(vParents) == 0) return(-1)
      else if (vParents == 1) return(-2)
      else if (vParents == 2) return(-3)
      
    }
  }
  
  if (all(vOrder == mOrders[2,])) {
    
    if (node == 1) {
      if (length(vParents) == 0) return(-1)
      else if (vParents == 1) return(-2)    
      
    } else if (node == 3) {
      if (length(vParents) == 0) return(-1)
      else if (vParents == 2) return(-2)
      else if (vParents == 1) return(-3)
      
    }
  }
  
  if (all(vOrder == mOrders[3,])) {
    
    if (node == 2) {
      if (length(vParents) == 0) return(-1)
      else if (vParents == 3) return(-2)    
      
    } else if (node == 1) {
      if (length(vParents) == 0) return(-1)
      else if (vParents == 3 return(-2)
      else if (vParents == 2) return(-3)
      
    }
  }
  
  if (all(vOrder == mOrders[4,])) {
    
    if (node == 3) {
      if (length(vParents) == 0) return(-1)
      else if (vParents == 1) return(-2)    
      
    } else if (node == 2) {
      if (length(vParents) == 0) return(-1)
      else if (vParents == 1) return(-2)
      else if (vParents == 3) return(-3)
      
    }
  }
  
  if (all(vOrder == mOrders[5,])) {
    
    if (node == 3) {
      if (length(vParents) == 0) return(-1)
      else if (vParents == 2) return(-2)    
      
    } else if (node == 1) {
      if (length(vParents) == 0) return(-1)
      else if (vParents == 2) return(-2)
      else if (vParents == 3) return(-3)
      
    }
  }
  
  if (all(vOrder == mOrders[6,])) {
    
    if (node == 1) {
      if (length(vParents) == 0) return(-1)
      else if (vParents == 3) return(-2)    
      
    } else if (node == 2) {
      if (length(vParents) == 0) return(-1)
      else if (vParents == 3) return(-2)
      else if (vParents == 1) return(-3)
      
    }
  }
  
  stop("Invalid node, order, or parent set")
}
