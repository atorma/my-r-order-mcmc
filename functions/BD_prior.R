# Returns the scalar hyperparameter value for a uniform Bayesian Dirichlet Equivalent prior.
getAlphaBDeu <- function(equivalentSampleSize, cardinalityOfVariable, numberOfParentStates) {
  alpha <- equivalentSampleSize/(cardinalityOfVariable*numberOfParentStates);
  alpha[alpha == Inf] <- 1;
  return(alpha)
}

# Returns a matrix of BDEu prior parameters.
# Each matrix row corresponds to a parent configuration, each column to a node state.
getBDeuParams <- function(node, parents, cardinalities, equivalentSampleSize) {
  if (length(cardinalities) == 1) {
    cardinalities <- rep(cardinalities, max(c(node, parents)) )
  }
  
  cardinalityOfVariable <- cardinalities[node]
  numberOfParentStates <- prod(cardinalities[parents])
  alphas <- matrix(getAlphaBDeu(equivalentSampleSize, cardinalityOfVariable, numberOfParentStates), nrow=numberOfParentStates, ncol=cardinalityOfVariable)
  return(alphas)
}

# Returns a function f(node, parents) that returns the BDeu Dirichlet prior parameter
# matrix for the node and parents.
# 
# cardinalities:  cardinalities[i] is the cardinality of domain of node i
# equivalentSampleSize: Equivalent sample size of the BDeu prior
createBDeuPriorParamsProvider <- function(cardinalities, equivalentSampleSize) {
  return(function(node, parents) {
    getBDeuParams(node, parents, cardinalities, equivalentSampleSize)
  }) 
}


# Returns a function f(node, parents) that returns the K2-type BD Dirichlet prior parameter
# matrix for the node and parents.
# 
# cardinalities:  cardinalities[i] is the cardinality of domain of node i
# alpha: Constant scalar value of all Dirichlet prior parameters
createBDK2PriorParamsProvider <- function(cardinalities, alpha) {
  return(function(node, parents) {
    matrix(alpha, nrow=prod(cardinalities[parents]), ncol=cardinalities[node])
  }) 
}
