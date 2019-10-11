#' @importFrom stats dbinom qbeta

g <- function(m,k,n,p) {
  sum(dbinom(m:k,size=n,prob=p))
}

g.prime <- function(m,k,n,p) { # i.e., dg/dp
  #n*p^(m-1)*(1-p)^(n-m)*(choose(n-1,m-1)-choose(n-1,k)*(p/(1-p))^(k-m+1))
  q <- (1-p)
  out <- exp(lchoose(n-1,m-1)+(m-1)*log(p)+(n-m)*log(q))
  out <- out - exp(lchoose(n-1,k)+k*log(p)+(n-k-1)*log(q))
  out <- n*out
  out
}

evaluate.g.maximum <- function(m,k,n) {
  if(m==0 || k==n) {
    max.value <- 1
  } else {
    a <- (lchoose(n-1,m-1)-lchoose(n-1,k))/(k-m+1)
    best.p <- 1-1/(1+exp(a))
    max.value <- g(m,k,n,best.p)
  }
  return(max.value)
}

where.g.maximum <- function(m,k,n) {
  if(m==0 && k==n) {
    return(NULL) # because best.p is non-unique
  } else if(m==0) {
    return(0)
  } else if(k==n) {
    return(1)
  } else {
    a <- (lchoose(n-1,m-1)-lchoose(n-1,k))/(k-m+1)
    best.p <- 1-1/(1+exp(a))
    return(best.p)
  }
}

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

evaluate.coverage <- function(L.vec, U.vec=NULL, digits=10) {
  n <- length(L.vec)-1
  stopifnot(n>0)

  L.vec <- round(L.vec*10^digits)/10^digits
  if(is.null(U.vec)) {
    U.vec <- rev(1-L.vec)
  }
  U.vec <- round(U.vec*10^digits)/10^digits
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

# --
# main function
#            -- #
#' Blyth-Still-Casella exact binomial confidence intervals
#'
#' computes the Blyth-Still-Casella exact binomial confidence intervals based on a refining procedure proposed by George Casella (1986).
#'
#' @param n number of trials
#' @param X number of successes (optional)
#' @param alpha confidence level = 1 - alpha
#' @param digits number of significant digits after the decimal point
#' @param CIs.init initial confidence intervals from which the refinement procedure begins
#' (default starts from Clopper-Pearson confidence intervals)
#' @param additional.info additional information about the types of interval endpoints and their possible range is provided if TRUE (default = FALSE)
#' @return If \code{X} is specified, the corresponding confidence interval will be returned, otherwise a list of n + 1 confidence intervals will be returned.
#' @return If \code{additional.info = FALSE}, only a list of confidence interval(s) will be returned. For any conincidental endpoint, midpoint of its range will be displayed.
#' @return If \code{additional.info = TRUE}, the following lists will be returned:
#' @return \tabular{ll}{
#'    \code{CI}                          \tab a list of confidence intervals \cr
#'    \code{coinc.index}                 \tab indices of coincidental lower endpoints (L.Index) and their corresponding upper endpoints (U.index)\cr
#'    \code{endpoint.type}               \tab whether the endpoint is coincidental (C) or non-coincidental (NC)\cr
#'    \code{range}                       \tab range for each endpoint\cr
#' }
#' @examples
#' # to obtain 95% CIs for n = 30 and X = 0 to 30
#' blyth.still.casella(n = 30, alpha = 0.05, digits = 4)
#'
#' # to obtain 90% CIs, endpoint types, indices of coincidental enpoints (if any),
#' # and range of each endpoint for n = 30 and X = 23
#' blyth.still.casella(n = 30, X = 23, alpha = 0.05, digits = 4, additional.info = TRUE)
#'
#' # use initial confidence intervals defined by the user instead of Clopper-Pearson CIs
#' # CIs.input needs to be a (n + 1) x 2 matrix with sufficient coverage
#' CIs.input <- matrix(c(0,1), nrow = 11, ncol = 2, byrow = TRUE) # start with [0,1] intervals
#' blyth.still.casella(n = 10, alpha = 0.05, digits = 4, CIs.init = CIs.input, additional.info = TRUE)
#'
#' # use summary function to see the range for each endpoint
#' output <- blyth.still.casella(n = 5, alpha = 0.1, digits = 4, additional.info = TRUE)
#' summary(output)
#' @export
blyth.still.casella <- function(n,
                                X=NULL,
                                alpha=0.05,
                                digits=2,
                                CIs.init=NULL,
                                additional.info=FALSE) {

  stopifnot(alpha>=0 && alpha<1)
  stopifnot(n==floor(n) && n>0)
  if(!is.null(X)) {
    stopifnot(X==floor(X) && X>=0 && X<=n)
  }

  if(alpha==0) {
    CI.mat <- matrix(NA, nrow=n+1, ncol=2)
    CI.mat[,1] <- 0
    CI.mat[,2] <- 1
    rownames(CI.mat) <- paste0("X=",0:n)
    colnames(CI.mat) <- c("L","U")

    Endpoint.Type.mat <- matrix("NC",nrow=n+1,ncol=2)
    rownames(Endpoint.Type.mat) <- paste0("X=",0:n)
    colnames(Endpoint.Type.mat) <- c("L","U")

    Range.mat <- matrix(NA,nrow=n+1,ncol=4)
    Range.mat[,1:2] <- 0
    Range.mat[,3:4] <- 1
    rownames(Range.mat) <- paste0("X=",0:n)
    colnames(Range.mat) <- c("L.min","L.max","U.min","U.max")

    if(!is.null(X)) {
      CI.mat            <- CI.mat[X+1,,drop=FALSE]
      Endpoint.Type.mat <- Endpoint.Type.mat[X+1,,drop=FALSE]
      Range.mat         <- Range.mat[X+1,,drop=FALSE]
    }

    if(!additional.info) {
      result <- CI.mat
      class(result) <- "customclass"
      return(result)
    } else {
      result <- list(CI=CI.mat,Endpoint.Type=Endpoint.Type.mat,Range=Range.mat)
      class(result) <- "customclass"
      return(result)
    }
  }

  if(n==1) {
    CI.mat <- matrix(NA, nrow=2, ncol=2)
    CI.mat[1,] <- c(             0,max(1-alpha,0.5))
    CI.mat[2,] <- c(min(alpha,0.5),               1)
    CI.mat[1,2] <- ceiling(CI.mat[1,2]*10^digits)/10^digits
    CI.mat[2,1] <- floor(CI.mat[2,1]*10^digits)/10^digits
    #CI.mat[1,2] <- ceiling((1-alpha)*10^digits)/10^digits
    #CI.mat[2,1] <- floor(alpha*10^digits)/10^digits
    rownames(CI.mat) <- c("X=0","X=1")
    colnames(CI.mat) <- c("L","U")

    Endpoint.Type.mat <- matrix("NC",nrow=2,ncol=2)
    rownames(Endpoint.Type.mat) <- c("X=0","X=1")
    colnames(Endpoint.Type.mat) <- c("L","U")

    Range.mat <- cbind(CI.mat[,1],CI.mat[,1],CI.mat[,2],CI.mat[,2])
    rownames(Range.mat) <- c("X=0","X=1")
    colnames(Range.mat) <- c("L.min","L.max","U.min","U.max")

    if(!is.null(X)) {
      CI.mat            <- CI.mat[X+1,,drop=FALSE]
      Endpoint.Type.mat <- Endpoint.Type.mat[X+1,,drop=FALSE]
      Range.mat         <- Range.mat[X+1,,drop=FALSE]
    }

    if(!additional.info) {
      result <- CI.mat
      class(result) <- "customclass"
      return(result)
    } else {
      result <- list(CI=CI.mat,Endpoint.Type=Endpoint.Type.mat,Range=Range.mat)
      class(result) <- "customclass"
      return(result)
    }
  }

  digits.save <- digits

  if(!is.null(CIs.init)) {
    stopifnot(ncol(CIs.init)==2 && nrow(CIs.init)==(n+1))
    L.vec <-   floor(CIs.init[,1]*10^digits)/10^digits
    U.vec <- ceiling(CIs.init[,2]*10^digits)/10^digits
    stopifnot(all(diff(L.vec)>=0))
    stopifnot(all(diff(U.vec)>=0))
    min.coverage <- min(evaluate.coverage(L.vec,U.vec,digits)[[1]], evaluate.coverage(L.vec,U.vec,digits)[[2]])
    if(min.coverage<(1-alpha)) stop("Initial set of confidence intervals does not have sufficient coverage probability");
  } else {
    # Clopper-Pearson exact confidence intervals
    L.vec <- qbeta(alpha/2,0:n,(n+1):1)
    U.vec <- qbeta(1-alpha/2,1:(n+1),n:0)
    L.vec <-   floor(L.vec*10^digits)/10^digits
    U.vec <- ceiling(U.vec*10^digits)/10^digits
  }

  move.noncoincidental.endpoint <- function() { # move a noncoincidental lower endpoint (indexed by i) to the right
    L.at.i <- L.vec[i+1]
    j <- which(U.vec[0:(i-1)+1]>L.at.i)
    if(length(j)==0) {
      i <<- i-1
      return(invisible(NULL))
    }# if occurs, the endpoint is not movable
    j <- j[1]-1
    # determine first.touch
    # L.at.i can possibly touch one of the following endpoints:
    #   (1) L.at.(i+1), non-moving
    #   (2) U.at.i    , non-moving # pathological interval
    #   (3) U.at.j    , non-moving
    #   (4) 0.5       , if i==(n/2) or (i+j)==n
    first.touch.vers <- rep(Inf,4)
    if(i<n)                  { # L.at.(i+1) exists
      first.touch.vers[1] <- L.vec[i+2]
    }
    if(i!=(n/2))             { # U.at.i not moving
      first.touch.vers[2] <- U.vec[i+1]
    }
    if((i+j)!=n)             { # U.at.j not moving
      first.touch.vers[3] <- U.vec[j+1]
    }
    if(i==(n/2) || (i+j)==n) { # U.at.i moving or U.at.j moving
      first.touch.vers[4] <- 0.5
    }
    first.touch <- min(first.touch.vers)
    if(first.touch==Inf) {
      i <<- i-1
      return(invisible(NULL))
    }

    if(g(j,i-1,n,first.touch)<(1-alpha)) { # insufficient coverage
      endpts <- find.endpoints.given.confidence.level(j,i-1,n,alpha,digits)
      if(is.null(endpts)) {
        i <<- (i-1)
      } else {
        right.endpt <- endpts[2]
        if(L.at.i < right.endpt) {
          L.vec[i+1]   <<- right.endpt
          U.vec[n-i+1] <<- round((1-right.endpt)*10^digits)/10^digits
          coincidental.endpoint.idx[i+1] <<- NA
          i <<- (i-1)
        } else { # i.e., L.vec[i+1]>=right.endpt and cannot be moved
          coincidental.endpoint.idx[i+1] <<- NA
          i <<- (i-1)
        }
      }
    } else { # sufficient coverage
      if(L.at.i < first.touch) { # then move to first.touch
        L.vec[i+1]   <<- first.touch
        U.vec[n-i+1] <<- round((1-first.touch)*10^digits)/10^digits
        # reset endpoint type to be non-coincidental
        coincidental.endpoint.idx[i+1] <<- NA
      } else {
        i <<- i-1
      }
    }

  }

  move.coincidental.endpoint <- function() { # move a coincidental lower endpoint (indexed by i) to the right
    L.at.i <- L.vec[i+1]
    j.equal <- which(U.vec[0:(i-1)+1]==L.at.i)-1

    # to remove index that have been assigned as coincidental endpoints from previous iteration
    if(i<n) {
      j.not.taken <- setdiff(j.equal, coincidental.endpoint.idx[(i+2):(n+1)])
      if(!all(is.na(coincidental.endpoint.idx[(i+2):(n+1)]))) { # if not all elements are NA
        # force index to be smaller than already assigned index
        j.not.taken <- j.not.taken[j.not.taken < min(coincidental.endpoint.idx[(i+2):(n+1)], na.rm = T)]
      }
    }else {
      j.not.taken <- j.equal
    }
    if(length(j.not.taken)>0) {
      coincidental.endpoint.idx[i+1] <<- j.not.taken[length(j.not.taken)]
      j <- j.not.taken[length(j.not.taken)]
    }else{
      stop("No matching j is found")
    }

    if(j==(n-i)) {
      i <<- i-1
      return(invisible(NULL)) # although coincidental, not movable
    }

    if(g(j,i-1,n,L.at.i)<(1-alpha) || g(j+1,i,n,L.at.i)<(1-alpha)) {
      # j is not the index for coincidental point, since coverage is not met
      # to find the index whose coverage is sufficient
      new.j <- j
      k <- length(j.not.taken) # move to the second most large j.equal
      while (g(new.j,i-1,n,L.at.i)<(1-alpha) || g(new.j+1,i,n,L.at.i)<(1-alpha)) {
        if(k<1) {
          stop("Error: None of the endpoints satisfy the coverage")
        }
        else{
          k <- k-1
          new.j <- j.not.taken[k]
        }
      }
      coincidental.endpoint.idx[i+1] <<- new.j
      i <<- i-1
      return(invisible(NULL)) # although coincidental, not movable
    }

    # determine first.touch
    # L.at.i can possibly touch one of the following endpoints:
    #   (1) L.at.(i+1), non-moving
    #   (2) U.at.i    , non-moving # pathological interval
    #   (3) 0.5       , if i+j<n
    # U.at.j can possibly touch one of the following endpoints:
    #   (4) U.at.(j+1), non-moving
    #   (5) 0.5       , if i+j<n
    # when either case (3) or (5) occurs: i.e., L.at.i=L.at.(n-j)=U.at.j=U.at.(n-i)=0.5
    #
    # Case (6): l_i & u_j are actually separable

    first.touch.vers <- rep(Inf,6)
    if(i<n && (i+1)!=(n-j)) { # L.at.(i+1) not moving
      first.touch.vers[1] <- L.vec[i+2]
    }
    if(i!=(n-i)) { # U.at.i not moving
      first.touch.vers[2] <- U.vec[i+1]
    }
    if((i+j)<n) { # i<(n-j)
      first.touch.vers[3] <- 0.5
    }
    if((j+1)!=(n-i)) { # U.at.(j+1) not moving
      first.touch.vers[4] <- U.vec[j+2]
    }
    if((i+j)<n) { # j<(n-i) # this case is repetitive
      first.touch.vers[5] <- 0.5
    }
    # Case (6)
    if(j<(i-1)) {
      endpts.sep.set <- find.endpoints.given.confidence.level(j+1,i-1,n,alpha,digits) # defines the range in which l_i & u_j are separable
      if(!is.null(endpts.sep.set)) {
        if(endpts.sep.set[2]>L.at.i) {
          g.max <- where.g.maximum(j+1,i-1,n)
          first.choice  <- round(g.max*10^digits)/10^digits
          second.choice <- floor(endpts.sep.set[2]*10^digits)/10^digits
          if(first.choice>endpts.sep.set[2] || first.choice<=L.at.i) {
            first.choice <- NULL
          }
          if(second.choice<=L.at.i) {
            second.choice <- NULL
          }

          if(!is.null(first.choice)) {
            first.touch.vers[6] <- first.choice
          } else if(!is.null(second.choice)) {
            first.touch.vers[6] <- second.choice
          }
        }
      }
    }
    first.touch <- min(first.touch.vers)

    if(first.touch==Inf) {
      i <<- i-1
      return(invisible(NULL))
    }

    if(g(j,i-1,n,first.touch)<(1-alpha) || g(j+1,i,n,first.touch)<(1-alpha)) { # insufficient coverage
      endpts.set1 <- find.endpoints.given.confidence.level(j,i-1,n,alpha,digits)
      endpts.set2 <- find.endpoints.given.confidence.level(j+1,i,n,alpha,digits)

      if(is.null(endpts.set1) || is.null(endpts.set2)) {
        i <<- (i-1)
      } else {
        right.endpt  <- endpts.set1[2] # set 1's right endpoint
        if(right.endpt < endpts.set2[1] || right.endpt <= L.at.i) {
          i <<- (i-1)
        } else {
          L.vec[i+1]   <<- right.endpt
          U.vec[n-i+1] <<- round((1-right.endpt)*10^digits)/10^digits
          U.vec[j+1]   <<- right.endpt
          L.vec[n-j+1] <<- round((1-right.endpt)*10^digits)/10^digits
          i <<- (i-1)
        }
      }
    } else { # sufficient coverage
      if(L.at.i < first.touch) {
        L.vec[i+1]   <<- first.touch
        U.vec[n-i+1] <<- round((1-first.touch)*10^digits)/10^digits
        U.vec[j+1]   <<- first.touch
        L.vec[n-j+1] <<- round((1-first.touch)*10^digits)/10^digits
        # i stays at i
      } else {
        i <<- (i-1)
      }
    }

  }

  L.no.longer.changes <- FALSE
  iter <- 0
  while(!L.no.longer.changes) {
    iter <- iter + 1
    L.old <- L.vec
    i <- n
    coincidental.endpoint.idx <- rep(NA, n+1) # a matrix recording coincidental points

    while(i>0) {
      L.at.i <- L.vec[i+1]
      j.equal <- which(U.vec[0:(i-1)+1]==L.at.i)-1
      # to remove index that have been assigned as coincidental endpoints from previous iteration
      if(i<n){
        j.not.taken <- setdiff(j.equal, coincidental.endpoint.idx[(i+2):(n+1)])
        if(!all(is.na(coincidental.endpoint.idx[(i+2):(n+1)]))) { # if not all elements are NA
          # force index to be smaller than the smallest index that has been assigned
          j.not.taken <- j.not.taken[j.not.taken < min(coincidental.endpoint.idx[(i+2):(n+1)], na.rm = T)]
        }
      }else {
        j.not.taken <- j.equal
      }
      j <- j.not.taken

      if(length(j)>0) {
        j <- j[length(j)] # last element
        # determine if separable
        separable <- TRUE
        if(j==(i-1)) {
          separable <- FALSE
        } else {
          if(g(j+1,i-1,n,L.at.i)<(1-alpha)) {
            separable <- FALSE
          }
        }

        if(separable) {
          move.noncoincidental.endpoint()
        } else {
          move.coincidental.endpoint()
        }
      } else { # L.at.i does not equal to any U.at.(0..(i-1))
        move.noncoincidental.endpoint()
      }

    }

    L.new <- L.vec
    if(all(L.new==L.old)) L.no.longer.changes <- TRUE
  }

  # determine range for pairs of coincidental endpoints
  coincidental.endpt.mat <- cbind(c(0:n), coincidental.endpoint.idx)
  colnames(coincidental.endpt.mat) <- NULL
  # coincidental.endpt.mat.reduced: only include index for coincidental endpoints
  keep.ind <- which(apply(coincidental.endpt.mat,1,sum)!=n & !is.na(coincidental.endpt.mat[,2]))
  if(length(keep.ind)==0) {
    coincidental.endpt.mat.reduced <- c()
  } else {
    coincidental.endpt.mat.reduced <- coincidental.endpt.mat[keep.ind,,drop=FALSE]
  }
  if(!is.null(coincidental.endpt.mat.reduced)) {
    colnames(coincidental.endpt.mat.reduced) <- c("L.index", "U.index")
  }

  # find range for coinc. endpoints and store them in range.mat
  Range.mat <- cbind(L.vec,L.vec,U.vec,U.vec)
  if(!is.null(coincidental.endpt.mat.reduced)) {
    for(k in 1:nrow(coincidental.endpt.mat.reduced)) { # note that it is 'k', instead of 'i'
      i <- coincidental.endpt.mat.reduced[k,1]
      j <- coincidental.endpt.mat.reduced[k,2]

      endpts.set1 <- find.endpoints.given.confidence.level(j,i-1,n,alpha,digits)
      endpts.set2 <- find.endpoints.given.confidence.level(j+1,i,n,alpha,digits)

      if(is.null(endpts.set1) || is.null(endpts.set2)) {
        stop("Error in deciding allowable range for coincidental endpoints.
             Please consider increasing the argument 'digits'.")
      }

      allowable.range <- c(endpts.set2[1],endpts.set1[2])
      Range.mat[i+1,1:2]   <- allowable.range
      Range.mat[j+1,3:4]   <- allowable.range
      Range.mat[n-i+1,3:4] <- round(rev(1-allowable.range)*10^digits)/10^digits
      Range.mat[n-j+1,1:2] <- round(rev(1-allowable.range)*10^digits)/10^digits
      }
  }

  # enforce monotonicity of endpoints in Range.mat
  while(TRUE) {
    no.more.adjustment <- TRUE
    if(!is.null(coincidental.endpt.mat.reduced)) {
      for(k in 1:nrow(coincidental.endpt.mat.reduced)) { # note that it is 'k', instead of 'i'
        i <- coincidental.endpt.mat.reduced[k,1]
        j <- coincidental.endpt.mat.reduced[k,2]

        range.lower.limit <- Range.mat[i+1,1]
        range.upper.limit <- Range.mat[i+1,2]

        # aliases
        rll <- range.lower.limit
        rul <- range.upper.limit

        rll.candidates.vec <- c()
        rul.candidates.vec <- c()

        if(i+j<n) {
          rul.candidates.vec <- c(rul.candidates.vec, 0.5)
        }
        if(i+j>n) {
          rll.candidates.vec <- c(rll.candidates.vec, 0.5)
        }

        if(i>0) {
          rll.candidates.vec <- c(rll.candidates.vec, Range.mat[i,1])
        }
        if(j>0) {
          rll.candidates.vec <- c(rll.candidates.vec, Range.mat[j,3])
        }

        if(i<n) {
          rul.candidates.vec <- c(rul.candidates.vec, Range.mat[i+2,2])
        }
        if(j<n) {
          rul.candidates.vec <- c(rul.candidates.vec, Range.mat[j+2,4])
        }

        if(rll<max(rll.candidates.vec)) {
          no.more.adjustment <- FALSE
          rll <- max(rll.candidates.vec)
          Range.mat[  i+1,1] <- rll
          Range.mat[  j+1,3] <- rll
          Range.mat[n-i+1,4] <- 1-rll
          Range.mat[n-j+1,2] <- 1-rll
        }

        if(rul>min(rul.candidates.vec)) {
          no.more.adjustment <- FALSE
          rul <- min(rul.candidates.vec)
          Range.mat[  i+1,2] <- rul
          Range.mat[  j+1,4] <- rul
          Range.mat[n-i+1,3] <- 1-rul
          Range.mat[n-j+1,1] <- 1-rul
        }
      }
    }

    if(no.more.adjustment)
      break;
  }

  rownames(Range.mat) <- paste0("X=",0:n)
  colnames(Range.mat) <- c("L.min","L.max","U.min","U.max")

  # update L.vec & U.vec using mid-ranges
  L.vec <- round((Range.mat[,1]+Range.mat[,2])/2*10^digits)/10^digits
  U.vec <- rev(round((1-L.vec)*10^digits)/10^digits)
  if(!is.null(coincidental.endpt.mat.reduced)) { # realign coincidental endpoints after rounding
    for(k in 1:nrow(coincidental.endpt.mat.reduced)) {
      i <- coincidental.endpt.mat.reduced[k,1]
      j <- coincidental.endpt.mat.reduced[k,2]
      realigned.endpt <- max(L.vec[i+1],U.vec[j+1])
      L.vec[i+1] <- U.vec[j+1] <- realigned.endpt
    }
  }
  CI.mat <- cbind(L.vec,U.vec)
  rownames(CI.mat) <- paste0("X=",0:n)
  colnames(CI.mat) <- c("L","U")

  Endpoint.Type.mat <- matrix("NC",nrow=n+1,ncol=2)
  half.range.mat <- matrix(NA, nrow = n + 1, ncol = 2)
  if(!is.null(coincidental.endpt.mat.reduced)) {
    coincidental.Ls <- coincidental.endpt.mat.reduced[,1]
    coincidental.Us <- coincidental.endpt.mat.reduced[,2]
    Endpoint.Type.mat[coincidental.Ls+1,1] <- "C"
    Endpoint.Type.mat[coincidental.Us+1,2] <- "C"
    half.range.mat[coincidental.Ls+1,1] <- (Range.mat[coincidental.Ls+1,1] + Range.mat[coincidental.Ls+1,2])/2
    half.range.mat[coincidental.Us+1,2] <- (Range.mat[coincidental.Us+1,1] + Range.mat[coincidental.Us+1,2])/2
  }
  rownames(Endpoint.Type.mat) <- paste0("X=",0:n)
  colnames(Endpoint.Type.mat) <- c("L","U")
  rownames(half.range.mat) <- paste0("X=",0:n)
  colnames(half.range.mat) <- c("L","U")

  if(digits.save!=digits) { # digits.save < digits (which =4)
    digits    <- digits.save   # restore original precision
    CI.mat[,1]    <-   floor(CI.mat[,1]*10^digits)/10^digits
    CI.mat[,2]    <- ceiling(CI.mat[,2]*10^digits)/10^digits

    Range.mat[,1] <- pmin(ceiling(Range.mat[,1]*10^digits)/10^digits, CI.mat[,1])
    Range.mat[,2] <- pmax(  floor(Range.mat[,2]*10^digits)/10^digits, CI.mat[,1])
    Range.mat[,3] <- pmin(ceiling(Range.mat[,3]*10^digits)/10^digits, CI.mat[,2])
    Range.mat[,4] <- pmax(  floor(Range.mat[,4]*10^digits)/10^digits, CI.mat[,2])
  }

  # if user specified the value of X
  if(!is.null(X)) {
    CI.mat            <- CI.mat[X+1,,drop=FALSE]
    if(length(c(which(coincidental.endpt.mat.reduced[,1] == X),
                which(coincidental.endpt.mat.reduced[,2] == X))) == 0){
      coincidental.endpt.mat.reduced <- NULL
    }else{
      coincidental.endpt.mat.reduced <- coincidental.endpt.mat.reduced[c(which(coincidental.endpt.mat.reduced[,1] == X),
                                                                         which(coincidental.endpt.mat.reduced[,2] == X)),,drop=FALSE]
    }
    Endpoint.Type.mat <- Endpoint.Type.mat[X+1,,drop=FALSE]
    Range.mat         <- Range.mat[X+1,,drop=FALSE]
  }

  if(!additional.info) {
    result <- CI.mat
    class(result) <- "customclass"
    return(result)
  } else {
    result <- list(CI=CI.mat, coinc.index = coincidental.endpt.mat.reduced, endpoint.type=Endpoint.Type.mat, range=Range.mat)
    class(result) <- "customclass"
    return(result)
  }
}

# to customize the summary function
add.sign.to.number <- function(x, digits) {
  x.c <- lapply(x, function(x) {if(!is.na(x)) {x <- formatC(x, format = 'f', digits = digits)} else {x<- ""}})
  x.plus.sign <- lapply(x.c, function(x) {if(x != "") x <- paste(c(bquote("\U00B1"), x), collapse = " ")})
  x.plus.sign[sapply(x.plus.sign, is.null)] <- ""
  x.plus.sign <- unlist(x.plus.sign)
  return(x.plus.sign)
}

#' @export
summary.customclass <- function(object, ...) {
  arguments <- list(...)
  if (!is.null(names(object))) { # check if x has Range matrix
    ls <- object$CI[,1]
    us <- object$CI[,2]
    # to figure out the number of digits used for calculation
    lus <- c(ls, us)
    if(all(lus %in% c(0,1))){ # if alpha = 0
      l.for.digits <- 0
    } else { # if alpha != 0
      l.for.digits <- max(lus[!lus %in% c(0,1)])
    }
    digits <- nchar(sub("^.+[.]","",l.for.digits))
    # calculate ranges
    L.range <- floor(apply(object$range, 1, function(x) if(x[1] != x[2]) {1/2*(x[2] - x[1])} else NA)*10^digits)/10^digits
    L.range.plus.sign <- add.sign.to.number(L.range, digits)
    U.range <- floor(apply(object$range, 1, function(x) if(x[3] != x[4]) {1/2*(x[4] - x[3])} else NA)*10^digits)/10^digits
    U.range.plus.sign <- add.sign.to.number(U.range, digits)
    output <- data.frame(L = ls, L.range = L.range.plus.sign,
                         U = us, U.range = U.range.plus.sign)
    return(output)
  } else {
    output <- object
    class(output) <- "matrix"
    return(summary(output))
  }
}

