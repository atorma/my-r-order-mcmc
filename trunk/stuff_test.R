


#tests
cPossibleStates <- 3

mObs <- matrix(0, 5, 3)
mObs[1,] <- c(2, 2, 3)
mObs[2,] <- c(3, 2, 2)
mObs[3,] <- c(3, 1, 2)
mObs[4,] <- c(2, 2, 3)
mObs[5,] <- c(2, 1, 3)

# Local score for node 1 without parents
# Expected
lgamma(1) - lgamma(6) + lgamma(1/3+0) + lgamma(1/3+3) + lgamma(1/3+2) - 3*lgamma(1/3)
# Actual
getLogLocalStructureBelief(1, cPossibleStates, integer(0), mObs)

# Local score for node 1 with parent 2
# Expected
lgamma(3/9) - lgamma(3/9+2) + lgamma(1/9+0) + lgamma(1/9+1) + lgamma(1/9+1) - 3*lgamma(1/9)  +  lgamma(3/9) - lgamma(3/9+3) + lgamma(1/9+0) + lgamma(1/9+2) + lgamma(1/9+1) - 3*lgamma(1/9)  +  lgamma(3/9) - lgamma(3/9+0) + 3*lgamma(3/9+0) - 3*lgamma(3/9)
# Actual
getLogLocalStructureBelief(1, cPossibleStates, 2, mObs)




