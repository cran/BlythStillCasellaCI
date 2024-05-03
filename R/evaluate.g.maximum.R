evaluate.g.maximum <- function(m, k, n) {
  if (m == 0 || k == n) {
    max.value <- 1
  } else {
    a <- (lchoose(n-1, m-1) - lchoose(n-1, k)) / (k - m + 1)
    best.p <- 1 - 1/(1 + exp(a))
    max.value <- g(m, k, n, best.p)
  }
  return(max.value)
}
