# Returns the scalar hyperparameter value for a uniform Bayesian Dirichlet Equivalent prior.
getAlphaBDEu <- function(equivalentSampleSize, cardinalityOfVariable, numberOfParentStates) {
  alpha <- equivalentSampleSize/(cardinalityOfVariable*numberOfParentStates);
  alpha[alpha == Inf] <- 1;
  return(alpha)
}

# Returns a matrix of BDEu prior parameters.
# Each matrix row corresponds to a parent configuration, each column to a node state.
getBDEuParams <- function(node, parents, cardinalities, equivalentSampleSize) {
  if (length(cardinalities) == 1) {
    cardinalities <- rep(cardinalities, max(c(node, parents)) )
  }
  
  cardinalityOfVariable <- cardinalities[node]
  numberOfParentStates <- prod(cardinalities[parents])
  alphas <- matrix(getAlphaBDEu(equivalentSampleSize, cardinalityOfVariable, numberOfParentStates), nrow=numberOfParentStates, ncol=cardinalityOfVariable)
  return(alphas)
}