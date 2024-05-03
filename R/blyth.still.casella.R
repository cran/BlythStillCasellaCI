#' Blyth-Still-Casella Exact Binomial Confidence Intervals
#'
#' `blyth.still.casella()` computes Blyth-Still-Casella exact binomial confidence intervals based on a refining procedure proposed by George Casella (1986).
#'
#' @importFrom stats dbinom qbeta
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
#'    \code{ConfidenceInterval}          \tab a list of confidence intervals \cr
#'    \code{CoincidenceEndpoint}         \tab indices of coincidental lower endpoints (L.Index) and their corresponding upper endpoints (U.index)\cr
#'    \code{Range}                       \tab range for each endpoint\cr
#' }
#' @examples
#' # to obtain 95% CIs for n = 30 and X = 0 to 30
#' blyth.still.casella(n = 30, alpha = 0.05, digits = 4)
#'
#' # to obtain 90% CIs, endpoint types, indices of coincidental endpoints (if any),
#' # and range of each endpoint for n = 30 and X = 23
#' blyth.still.casella(n = 30, X = 23, alpha = 0.05, digits = 4, additional.info = TRUE)
#'
#' # use initial confidence intervals defined by the user instead of Clopper-Pearson CIs
#' # CIs.input needs to be a (n + 1) x 2 matrix with sufficient coverage
#' CIs.input <- matrix(c(0,1), nrow = 11, ncol = 2, byrow = TRUE) # start with [0,1] intervals
#' blyth.still.casella(n = 10, alpha = 0.05, digits = 4, CIs.init = CIs.input, additional.info = TRUE)
#'
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

  # warning message
  if(additional.info && digits < 4) {
    warn.msg <-
      paste0("For a more accurate determination of pairs of coincidental ",
             "endpoints, please consider\n",
             "    increasing the number of significant digits after the ",
             "decimal point (digits) to at\n",
             "    least 4.\n")
    warning(warn.msg)
  }

  if(alpha==0) {
    CI.mat <- matrix(NA, nrow=n+1, ncol=2)
    CI.mat[,1] <- 0
    CI.mat[,2] <- 1
    rownames(CI.mat) <- paste0("X=",0:n)
    colnames(CI.mat) <- c("L","U")

    coincidental.endpt.mat.reduced <- NULL

    Range.mat <- matrix(NA,nrow=n+1,ncol=4)
    Range.mat[,1:2] <- 0
    Range.mat[,3:4] <- 1
    rownames(Range.mat) <- paste0("X=",0:n)
    colnames(Range.mat) <- c("L.min","L.max","U.min","U.max")

    if(!is.null(X)) {
      CI.mat            <- CI.mat[X+1,,drop=FALSE]
      # coincidental.endpt.mat.reduced <- NULL
      Range.mat         <- Range.mat[X+1,,drop=FALSE]
    }

    if(!additional.info) {
      result <- CI.mat
      return(result)
    } else {
      result <- list(ConfidenceInterval=CI.mat,
                     CoincidentalEndpoint=coincidental.endpt.mat.reduced,
                     Range=Range.mat)
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

    coincidental.endpt.mat.reduced <- NULL

    Range.mat <- cbind(CI.mat[,1],CI.mat[,1],CI.mat[,2],CI.mat[,2])
    rownames(Range.mat) <- c("X=0","X=1")
    colnames(Range.mat) <- c("L.min","L.max","U.min","U.max")

    if(!is.null(X)) {
      CI.mat            <- CI.mat[X+1,,drop=FALSE]
      # coincidental.endpt.mat.reduced <- NULL
      Range.mat         <- Range.mat[X+1,,drop=FALSE]
    }

    if(!additional.info) {
      result <- CI.mat
      return(result)
    } else {
      result <- list(ConfidenceInterval=CI.mat,
                     CoincidentalEndpoint=coincidental.endpt.mat.reduced,
                     Range=Range.mat)
      return(result)
    }
  }

  if(!is.null(CIs.init)) {
    stopifnot(ncol(CIs.init)==2 && nrow(CIs.init)==(n+1))
    L.vec <-   floor(CIs.init[,1]*10^digits)/10^digits
    U.vec <- ceiling(CIs.init[,2]*10^digits)/10^digits
    stopifnot(all(diff(L.vec)>=0))
    stopifnot(all(diff(U.vec)>=0))
    min.coverage <- min(evaluate.coverage(L.vec,U.vec,digits)[,1],
                        evaluate.coverage(L.vec,U.vec,digits)[,2])
    if(min.coverage<(1-alpha)) {
      err.msg <- paste0("Initial set of confidence intervals does not have ",
                        "sufficient coverage probability")
      stop(err.msg)
    }
  } else {
    # Clopper-Pearson exact confidence intervals
    L.vec <- qbeta(alpha/2,0:n,(n+1):1)
    U.vec <- qbeta(1-alpha/2,1:(n+1),n:0)
    L.vec <-   floor(L.vec*10^digits)/10^digits
    U.vec <- ceiling(U.vec*10^digits)/10^digits
  }

  move.noncoincidental.endpoint <- function() { # move a noncoincidental
                                                # lower endpoint (indexed by i)
                                                # to the right
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
    #   (2) U.at.j    , non-moving
    #   (3) 0.5       , if (i+j)==n
    first.touch.vers <- rep(Inf,3)
    if(i<n)                  { # L.at.(i+1) exists
      first.touch.vers[1] <- L.vec[i+2]
    }
    first.touch.vers[2] <- U.vec[j+1]
    if((i+j)==n)             {
      first.touch.vers[3] <- 0.5
    }
    first.touch <- min(first.touch.vers)

    if(g(j,i-1,n,first.touch)<(1-alpha)) { # insufficient coverage
      endpts <- find.endpoints.given.confidence.level(j,i-1,n,alpha,digits)
      if(is.null(endpts)) {
        i <<- (i-1)
      } else {
        right.endpt <- endpts[2]
        if(L.at.i < right.endpt) {
          L.vec[i+1]   <<- right.endpt
          U.vec[n-i+1] <<- round(1-right.endpt, digits)
          i <<- (i-1)
        } else { # i.e., L.vec[i+1]>=right.endpt and cannot be moved
          i <<- (i-1)
        }
      }
    } else { # sufficient coverage
      if(L.at.i < first.touch) { # then move to first.touch
        L.vec[i+1]   <<- first.touch
        U.vec[n-i+1] <<- round(1-first.touch, digits)
      } else {
        i <<- i-1
      }
    }

  }

  move.coincidental.endpoint <- function() { # move a coincidental lower
                                             # endpoint (indexed by i) to the
                                             # right
    L.at.i <- L.vec[i+1]
    j.equal <- which(U.vec[0:(i-1)+1]==L.at.i)-1
    j1 <- min(j.equal)
    j2 <- max(j.equal)

    if(j1<j2) {
      if(g(j2,i-1,n,L.at.i)<(1-alpha)) { # then cannot be moved
        i <<- i-1
        return(invisible(NULL))
      }
    }
    j <- j2 # move L.at.i and U.at.j to the right together

    if(j==(n-i)) {
      i <<- i-1
      return(invisible(NULL)) # although coincidental, not movable
    }

    # determine first.touch
    # L.at.i/U.at.j can possibly touch one of the following endpoints:
    #   (1) L.at.(i+1), non-moving
    #   (2) U.at.(j+1), non-moving
    #   (3) 0.5       , if i+j=n-1
    #
    # Case (4): l_i & u_j are actually separable

    first.touch.vers <- rep(Inf,4)
    if(i<n) {
      first.touch.vers[1] <- L.vec[i+2] # L.at.(i+1)
    }
    first.touch.vers[2] <- U.vec[j+2] # U.at.(j+1)
    if((i+j)==(n-1)) {
      first.touch.vers[3] <- 0.5
    }
    # Case (4)
    if(j<(i-1)) {
      endpts.sep.set <-
        find.endpoints.given.confidence.level(j+1,i-1,n,alpha,digits)
        # defines the range in which l_i & u_j are separable
      if(!is.null(endpts.sep.set)) {
        if(endpts.sep.set[1]>L.at.i) {
          first.touch.vers[4] <- endpts.sep.set[1]
        }
      }
    }
    first.touch <- min(first.touch.vers)

    if(first.touch==Inf) {
      i <<- i-1
      return(invisible(NULL))
    }

    if(g(j,i-1,n,first.touch)<(1-alpha) ||
       g(j+1,i,n,first.touch)<(1-alpha)) { # insufficient coverage
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
          U.vec[n-i+1] <<- round(1-right.endpt, digits)
          U.vec[j+1]   <<- right.endpt
          L.vec[n-j+1] <<- round(1-right.endpt, digits)
          i <<- (i-1)
        }
      }
    } else { # sufficient coverage
      if(L.at.i < first.touch) {
        L.vec[i+1]   <<- first.touch
        U.vec[n-i+1] <<- round(1-first.touch, digits)
        U.vec[j+1]   <<- first.touch
        L.vec[n-j+1] <<- round(1-first.touch, digits)
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

    while(i>0) {
      L.at.i <- L.vec[i+1]
      j.equal <- which(U.vec[0:(i-1)+1]==L.at.i)-1

      if(length(j.equal)>0) {
        j <- j.equal[length(j.equal)] # last element
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

  # determine pairs of coincidental endpoints
  coincidental.endpt.mat <- cbind(c(0:n), rep(NA,n+1))
  j.set <- 0:(n-1)
  for(i in n:1) {
    j.set <- intersect(j.set, 0:(i-1))
    if(length(j.set)==0) break;

    if(!is.na(coincidental.endpt.mat[i+1,2])) next;

    L.at.i <- L.vec[i+1]
    j.equal <- j.set[U.vec[j.set+1]==L.at.i]

    if(length(j.equal)>0) {
      j <- j.equal[length(j.equal)] # last element
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
        # move.noncoincidental.endpoint()
      } else {
        # move.coincidental.endpoint()
        if((i+j)!=n) {
          coincidental.endpt.mat[  i+1,2] <- j
          coincidental.endpt.mat[n-j+1,2] <- (n-i)
          j.set <- setdiff(j.set, c(j,n-i))
        }
      }
    }

  }

  coincidental.endpt.mat.reduced <-
    coincidental.endpt.mat[!is.na(coincidental.endpt.mat[,2]),,drop=FALSE]

  if(nrow(coincidental.endpt.mat.reduced)==0) {
    coincidental.endpt.mat.reduced <- NULL
  } else {
    colnames(coincidental.endpt.mat.reduced) <- c("L.index", "U.index")
  }

  range.err.msg.flag <- FALSE
  # find range for coincidental endpoints and store them in Range.mat
  Range.mat <- cbind(L.vec,L.vec,U.vec,U.vec)
  if(!is.null(coincidental.endpt.mat.reduced)) {
    for(k in 1:nrow(coincidental.endpt.mat.reduced)) { # note that it is 'k',
                                                       #   instead of 'i'
      i <- coincidental.endpt.mat.reduced[k,1]
      j <- coincidental.endpt.mat.reduced[k,2]

      endpts.set1 <- find.endpoints.given.confidence.level(j,i-1,n,alpha,digits)
      endpts.set2 <- find.endpoints.given.confidence.level(j+1,i,n,alpha,digits)

      err.msg.flag <- FALSE                            # reset flag
      if(is.null(endpts.set1) || is.null(endpts.set2)) {
        err.msg <-
          paste0("Failure in deciding allowable range for coincidental endpoints. ",
                 "Please\n",
                 "    consider increasing the argument 'digits'.\n")
        warning(err.msg)
        err.msg.flag <- TRUE
        range.err.msg.flag <- TRUE

        # stop("Error in deciding "allowable range for coincidental endpoints.
        #      Please consider increasing the argument 'digits'.")
      }

      if(err.msg.flag) {
        allowable.range <- c(L.vec[i+1],L.vec[i+1])
      } else {
        if(endpts.set2[1]<=endpts.set1[2]) {
          allowable.range <- c(endpts.set2[1],endpts.set1[2])
        } else {
          warning(paste0("Numerical error in evaluating allowable range for ",
                         "coincidental endpoints.\n"))
          allowable.range <- c(L.vec[i+1],L.vec[i+1])
        }
      }
      Range.mat[i+1,1:2]   <- allowable.range
      Range.mat[j+1,3:4]   <- allowable.range
      Range.mat[n-i+1,3:4] <- round(rev(1-allowable.range),digits)
      Range.mat[n-j+1,1:2] <- round(rev(1-allowable.range),digits)
      }
  }

  # enforce monotonicity of endpoints in Range.mat
  while(TRUE) {
    no.more.adjustment <- TRUE
    if(!is.null(coincidental.endpt.mat.reduced)) {
      for(k in 1:nrow(coincidental.endpt.mat.reduced)) { # note that it is 'k',
                                                         #   instead of 'i'
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
          Range.mat[n-i+1,4] <- round(1-rll,digits)
          Range.mat[n-j+1,2] <- round(1-rll,digits)
        }

        if(rul>min(rul.candidates.vec)) {
          no.more.adjustment <- FALSE
          rul <- min(rul.candidates.vec)
          Range.mat[  i+1,2] <- rul
          Range.mat[  j+1,4] <- rul
          Range.mat[n-i+1,3] <- round(1-rul,digits)
          Range.mat[n-j+1,1] <- round(1-rul,digits)
        }
      }
    }

    if(no.more.adjustment)
      break;
  }

  rownames(Range.mat) <- paste0("X=",0:n)
  colnames(Range.mat) <- c("L.min","L.max","U.min","U.max")

  # update L.vec & U.vec using mid-ranges
  L.vec <- round((Range.mat[,1]+Range.mat[,2])/2, digits)
  U.vec <- round(rev(1-L.vec),digits)
  # realign coincidental endpoints after rounding
  if(!is.null(coincidental.endpt.mat.reduced)) {
    for(k in 1:nrow(coincidental.endpt.mat.reduced)) {
      i <- coincidental.endpt.mat.reduced[k,1]
      j <- coincidental.endpt.mat.reduced[k,2]
      # U.vec[j+1] <- L.vec[i+1] # should already equal to each other
      L.vec[n-j+1] <- U.vec[n-i+1] <- round(1-L.vec[i+1],digits)
    }
  }

  CI.mat <- cbind(L.vec,U.vec)
  rownames(CI.mat) <- paste0("X=",0:n)
  colnames(CI.mat) <- c("L","U")

  # if user specified the value of X
  if(!is.null(X)) {
    CI.mat            <- CI.mat[X+1,,drop=FALSE]
    if(length(c(which(coincidental.endpt.mat.reduced[,1] == X),
                which(coincidental.endpt.mat.reduced[,2] == X))) == 0){
      coincidental.endpt.mat.reduced <- NULL
    }else{
      coincidental.endpt.mat.reduced <-
        coincidental.endpt.mat.reduced[
          c(which(coincidental.endpt.mat.reduced[,1] == X),
            which(coincidental.endpt.mat.reduced[,2] == X)),,drop=FALSE]
    }
    Range.mat         <- Range.mat[X+1,,drop=FALSE]
  }

  if(!additional.info) {
    result <- CI.mat
    return(result)
  } else {
    result <- list(ConfidenceInterval=CI.mat,
                   CoincidenceEndpoint=coincidental.endpt.mat.reduced,
                   Range=Range.mat)
    return(result)
  }

  # return(list(ConfidenceInterval=CI.mat,
  #             CoincidentalEndpoint=coincidental.endpt.mat.reduced,
  #             Range=Range.mat))
}

