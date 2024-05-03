calc.total.length <- function(L.vec, U.vec=NULL) {
  
  n <- length(L.vec)-1
  stopifnot(n>0)
  
  if(is.null(U.vec)) {
    U.vec <- rev(1-L.vec)
  }
  stopifnot(length(U.vec)==(n+1))
  
  width.vec <- (U.vec-L.vec)
  stopifnot(all(width.vec>=0))

  total.length <- sum(width.vec)
  return(total.length)
}
