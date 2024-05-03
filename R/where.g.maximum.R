where.g.maximum <- function(m, k, n) {
  if (m == 0 && k == n) {
    return(NULL) # because best.p is non-unique
  } else if (m == 0) {
    return(0)
  } else if (k == n) {
    return(1)
  } else {
    a <- (lchoose(n-1, m-1) - lchoose(n-1, k)) / (k - m + 1)
    best.p <- 1 - 1/(1 + exp(a))
    return(best.p)
  }
}
