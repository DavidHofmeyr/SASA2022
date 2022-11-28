if(!('FKSUM'%in%installed.packages())) install.packages('FKSUM')
library(FKSUM)

# Objective function for evaluating negative entropy,
# alternatively pseudo-likelihood
# arguments:
# v = projection vector
# X = (whitened) data matrix
# h = bandwidth for density estimates

obj_ica <- function(v, X, h){
  nv <- sqrt(sum(v^2))
  p <- X%*%v/nv
  n <- length(p)
  fhat <- fk_sum(p, rep(1/n/h, n), h, nbin = 2000)
  sum(log(fhat))/n
}

# Gradient of negative entropy objective
# arguments are the same as for the objective

grad_ica <- function(v, X, h){
  nv <- sqrt(sum(v^2))
  p <- X%*%v/nv
  n <- length(p)
  fhat <- fk_sum(p, rep(1/n/h, n), h, nbin = 2000)
  d1 <- fk_sum(p, rep(1/n^2/h^2, n), h, nbin = 2000, type = 'dksum')/fhat
  d2 <- fk_sum(p, 1/n^2/h^2/fhat, h, nbin = 2000, type = 'dksum')
  -c(d1 + d2)%*%(X-p%*%t(v)/nv)/nv
}

# Function to optimise the projection matrix
# in order to minimise the negative entropy
# subject to orthonormality
# The only argument is the (un-whitened) data matrix

myica <- function(X, ncomp = NULL){
  n <- nrow(X)
  d <- ncol(X)
  if(is.null(ncomp)) ncomp = d
  E <- eigen(cov(X))
  K <- t(t(E$vectors[,1:ncomp])/sqrt(E$values[1:ncomp]))
  Y <- (X-matrix(colMeans(X), n, d, byrow = T))%*%K
  Y0 <- Y
  V <- matrix(0, ncomp, ncomp)
  h <- .75/n^.2
  for(j in 1:ncomp){
    v <- rnorm(ncomp-j+1)
    for(iter in 1:20){
      fval <- obj_ica(v, Y, h)
      dir <- c(grad_ica(v, Y, h))
      dir <- dir/sqrt(sum(dir^2))
      stp <- .5
      repeat{
        fnew <- obj_ica(v+stp*dir, Y, h)
        if((fnew-fval) > .0001*fval) break
        else if(stp < 1e-9) break
        else stp <- stp/3
      }
      if(stp > 1e-9) v <- v + stp*dir
      else break
    }
    if(j > 1){
      v <- Null(V[,1:(j-1)])%*%v
    }
    V[,j] <- v/sqrt(sum(v^2))
    Y <- Y0%*%Null(V[,1:j])
  }
  list(S = Y0%*%V, K = K, W = V)
}

# The following function is useful
# only for the experiment

plotimage = function(x) image(matrix(x, 481, 321, byrow = TRUE), col = rgb(1:100/100, 1:100/100, 1:100/100), xaxt = 'n', yaxt = 'n')