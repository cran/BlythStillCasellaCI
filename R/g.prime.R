g.prime <- function(m, k, n, p) {
  # i.e., dg/dp
  #n*p^(m-1)*(1-p)^(n-m)*(choose(n-1,m-1)-choose(n-1,k)*(p/(1-p))^(k-m+1))
  q <- (1 - p)
  out <- exp( lchoose(n-1, m-1) + (m-1)*log(p) + (n-m)*log(q) )
  out <- out - exp( lchoose(n-1, k) + k*log(p) + (n-k-1)*log(q) )
  out <- n*out
  out
}
