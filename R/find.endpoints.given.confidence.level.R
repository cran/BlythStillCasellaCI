find.endpoints.given.confidence.level <- function(m,k,n,alpha=0.05,digits=4) {
  left.endpt  <- NULL
  right.endpt <- NULL
  
  interval.length <- 10^(-digits)
  
  if(alpha==0) {
    if(m!=0 || k!=n) {
      return(NULL)
    } else {
      left.endpt  <- 0
      right.endpt <- 1
    }
  }
  
  determine.next.x <- function(x0,a,b) {
    # Newton-Raphson method
    g0       <- g(m,k,n,x0)
    g.prime0 <- g.prime(m,k,n,x0)
    if(is.na(g.prime0) || g.prime0==0 || is.infinite(g.prime0)) {
      return((a+b)/2)
    } else {
      f0 <- g0-(1-alpha) # search for root for g(x)-(1-alpha)
      x1 <- x0-f0/g.prime0
      if(x1>=a && x1<=b) {
        return(x1)
      } else {
        return((a+b)/2)
      }
    }
  }
  
  search.for.right.endpoint <- function() {
    while((b-a)>=interval.length) {
      x1 <- determine.next.x(x0,a,b)
      # update a and b
      g1 <- g(m,k,n,x1)
      if(g1<=(1-alpha)) {
        b <- x1
        if(a<(x1-interval.length*0.9)) { # avoid getting stuck with b <- x1 <- x1 <- ...
          a1 <- x1-interval.length*0.9
          ga1 <- g(m,k,n,a1)
          if(ga1>(1-alpha)) {
            a <- a1
          } else {
            b  <- a1
            x1 <- a1
          }
        }
      } else {
        a <- x1
        if(b>(x1+interval.length*0.9)) { # avoid getting stuck with a <- x1 <- x1 <- ...
          b1 <- x1+interval.length*0.9
          gb1 <- g(m,k,n,b1)
          if(gb1<=(1-alpha)) {
            b <- b1
          } else {
            a  <- b1
            x1 <- b1
          }
        }
      }
      x0 <- x1
    }
    
    bb <- floor(b*10^digits)/10^digits # bb<=b
    if(bb<=a) {
      right.endpt <<- bb
    } else {
      gbb <- g(m,k,n,bb)
      if(gbb<(1-alpha)) {
        aa <- floor(a*10^digits)/10^digits # aa<=a
        right.endpt <<- aa
      } else {
        right.endpt <<- bb
      }
    }
  }
  
  search.for.left.endpoint <- function() {
    while((b-a)>=interval.length) {
      x1 <- determine.next.x(x0,a,b)
      # update a and b
      g1 <- g(m,k,n,x1)
      if(g1<=(1-alpha)) {
        a <- x1
        if(b>(x1+interval.length*0.9)) { # avoid getting stuck with a <- x1 <- x1 <- ...
          b1 <- x1+interval.length*0.9
          gb1 <- g(m,k,n,b1)
          if(gb1>(1-alpha)) {
            b <- b1
          } else {
            a  <- b1
            x1 <- b1
          }
        }
      } else {
        b <- x1
        if(a<(x1-interval.length*0.9)) { # avoid getting stuck with b <- x1 <- x1 <- ...
          a1 <- x1-interval.length*0.9
          ga1 <- g(m,k,n,a1)
          if(ga1<=(1-alpha)) {
            a <- a1
          } else {
            b  <- a1
            x1 <- a1
          }
        }
      }
      x0 <- x1
    }
    
    aa <- ceiling(a*10^digits)/10^digits # aa>=a
    if(aa>=b) {
      left.endpt <<- aa
    } else {
      gaa <- g(m,k,n,aa)
      if(gaa<(1-alpha)) {
        bb <- ceiling(b*10^digits)/10^digits # bb>=b
        left.endpt <<- bb
      } else {
        left.endpt <<- aa
      }
    }
  }
  
  ## --
  if(m==0 && k==n) {
    left.endpt  <- 0
    right.endpt <- 1
  } else if(m==0) { # m==0 && k<n
    left.endpt  <- 0
    # search for right endpoint
    a <- 0
    b <- 1
    x0 <- 1/2
    g0 <- g(m,k,n,x0)
    if(g0<=(1-alpha)) {
      b <- x0
    } else {
      a <- x0
    }
    search.for.right.endpoint()
    
  } else if(k==n) { # m>0 && k==n
    right.endpt <- 1
    # search for left endpoint
    a <- 0
    b <- 1
    x0 <- 1/2
    g0 <- g(m,k,n,x0)
    if(g0<=(1-alpha)) {
      a <- x0
    } else {
      b <- x0
    }
    search.for.left.endpoint()
    
  } else { # m>0 and k<n
    max.coverage <- evaluate.g.maximum(m,k,n)
    if(max.coverage<(1-alpha)) {
      return(NULL)
    } else {
      g.max <- where.g.maximum(m,k,n)
      ## search for left endpoint
      a <- 0
      b <- g.max
      x0 <- (a+b)/2
      g0 <- g(m,k,n,x0)
      if(g0<=(1-alpha)) {
        a <- x0
      } else {
        b <- x0
      }
      search.for.left.endpoint()
      
      ## search for right endpoint
      a <- g.max
      b <- 1
      x0 <- (a+b)/2
      g0 <- g(m,k,n,x0)
      
      if(g0<=(1-alpha)) {
        b <- x0
      } else {
        a <- x0
      }
      search.for.right.endpoint()
      
      if(left.endpt>right.endpt) {
        return(NULL)
      }
    }
  }
  return(c(left.endpt,right.endpt))
}

