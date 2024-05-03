g <- function(m, k, n, p) {
  sum(dbinom(m:k, size = n, prob = p))
}

