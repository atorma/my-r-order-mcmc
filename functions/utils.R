# Computes log(sum(exp(x)) in a way where numerical stability is improved.
#
# With large number of observations log scores are large negative and cause underflow
# when computing exp(x). Transformation 
#
# log(sum(exp(x)) = m + log(sum(exp(x - m))),
#
# where m = min(x), improves numerical stability.
#
# exponentials: vector x in log(sum(exp(x))
getLogSumOfExponentials <- function(x) {
  m <- max(x)
  logSumExp <- m + log(sum(exp(x - m)))
  return(logSumExp)
  
  # Alternative is to use log(1+x) and knowing that exp(x[xmax] - x[xmax]) = exp(0) = 1.
  #xmax <- which.max(x)
  #return(log1p( sum(exp(x[-xmax] - x[xmax])) ) + x[xmax])
}

# Computes ROC curve. First column is false positive rate, 
# 2nd column is true positive rate 
getRocCurve <- function(mPredictedEdgeProbs, mActualAdjacencies) {
  thresholds <- seq(from=0, to=1, by=0.001)
  roc <- matrix(NA, length(thresholds), 2)
  for (i in 1:length(thresholds)) {
    predicted <- mPredictedEdgeProbs >= thresholds[i]
    TP <- sum(predicted*mActualAdjacencies)
    FP <- sum(predicted*!mActualAdjacencies)
    P <- sum(mActualAdjacencies)
    N <- sum(!mActualAdjacencies)
    roc[i,1] <- FP/N
    roc[i,2] <- TP/P
  }
  return(roc)
}

# Computes AUROC from the result of getRocCurve.
getAuroc <- function(roc) {
  # Order the points increasing by false positive rate
  ordering <- order(roc[,1])
  roc <- roc[ordering,]
  
  lastPoint <- nrow(roc)
  # Numeric integration using trapezoidial rule
  a <- roc[-lastPoint, 1]
  fa <- roc[-lastPoint, 2]
  b <- roc[-1, 1]
  fb <- roc[-1, 2]
  auroc <- sum((b - a)*(fa + fb)/2)
  return(auroc)
}

# Computes Kullback-Leibler divergence criterion.
#
# p: true probabilities
# q: approximated probabilities
getKLDivergence <- function(p, q) {
  if (length(p) != length(q)) stop("p and q must have same length")
  
  p <- p/sum(p)
  q <- q/sum(q)
  
  sum(p*log(p/q))
}