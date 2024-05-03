evaluate.coverage <- function(L.vec, U.vec=NULL, digits=10) {
  n <- length(L.vec)-1
  stopifnot(n>0)
  
  L.vec <- round(L.vec, digits)
  if(is.null(U.vec)) {
    U.vec <- rev(1-L.vec)
  }
  U.vec <- round(U.vec, digits)
  stopifnot(length(U.vec)==(n+1))
  
  coverage.at.left.endpt.vec  <- rep(NA, n+1)
  coverage.at.right.endpt.vec <- rep(NA, n+1)
  
  for(i in 1:length(L.vec)) { # evaluate coverage at all left endpoints
    left.endpt <- L.vec[i]
    if(left.endpt==0) {
      coverage.at.left.endpt.vec[i] <- ifelse(L.vec[1]==0,1,0)
    } else {
      X.range <- range(which(L.vec<left.endpt & left.endpt<=U.vec))
      m <- X.range[1]-1
      k <- X.range[2]-1
      coverage.at.left.endpt.vec[i] <- g(m,k,n,p=left.endpt)
    }
  }
  
  for(i in 1:length(U.vec)) { # evaluate coverage at all right endpoints
    right.endpt <- U.vec[i]
    if(right.endpt==1) {
      coverage.at.right.endpt.vec[i] <- ifelse(U.vec[n+1]==1,1,0)
    } else {
      X.range <- range(which(U.vec>right.endpt & L.vec<=right.endpt))
      m <- X.range[1]-1
      k <- X.range[2]-1
      coverage.at.right.endpt.vec[i] <- g(m,k,n,p=right.endpt)
    }
  }
  output <- cbind(coverage.at.left.endpt.vec,coverage.at.right.endpt.vec)
  rownames(output) <- paste("X=", c(0:n))
  return(output)
}
