fCmat <- function(n, i, ii) 
{
  ## get constranits
  A <- matrix(0, n - 2, n - 1)
  A[ , ii] <- 1
  if (ii < i){
    for (k in 1:(n - 1)) {
      if (k < ii) A[k, k] <- -1
      if (k > ii) A[k - 1, k] <- -1
    }
  } else {
    for (k in 1:(n - 1)) {
      if (k < ii) A[k, k] <- 1
      if (k > ii) A[k - 1, k] <- 1
    }
  }
  A
} 
 
fRHS <- function(X, i, ii) 
{
  ## get RHS values in the constraints
  n <- dim(X)[1]
  if (ii < i) {
    if (ii > 1) {
      b <- 1 - c(X[1:(ii - 1), ii], X[ii,(ii + 1):n][-(i - ii)])
    } else {
      b <- 1 - c(X[ii,(ii + 1):n][-(i - ii)])
    }
  } else {
    if(ii < n - 1) {
      b <- 1 + c(X[1:ii, (ii + 1)][-i], X[(ii + 1), (ii + 2):n])
    } else {
      b <- 1 + c(X[1:ii, (ii + 1)][-i])
    }
  }
  b
}
 
objfeval <- function(X, rho, K) {
  ## objective function evaluation
  indices <- upper.tri(rho)
  sum((rho[indices] - K) * X[indices])
}
 
LauGreenAlg0 <- function(K, rho)
{
  ##############################################
  # Bayesian Model Based Clustering Procedures #
  # John W. Lau and Peter J. Greeny            #
  # June 10, 2006                              #
  # Algorithm 1, p. 19                         #
  ##############################################
  #require(lpSolve)  ## to solve the binary integer programming problem
  n   <- nrow(rho)            ## n subjects
  X   <- matrix(0.0, n, n)    ## initial partition 
  val <- objfeval(X, rho, K)  ## objective function value
  repeat {
    for (i in 1:n) {
      ## coeffs in \sum X_{j,k}(rho_{j,k} - K)
      f.obj <- rep(NA, n)
      for (k in 1:n) {
        if (k < i)  ## (1,i), (2,i), ..., (i-1,i)
          f.obj[k] <- rho[k, i] - K
        if (k > i)  ## (i,i+1), (i,i+2), ..., (i,n)
          f.obj[k] <- rho[i, k] - K
      }
      f.obj <- f.obj[-i]
      ## constraints coeffs
      ## X_{ij} + X_{ik} - X_{jk} \leq 1 for \{i\neq j\neq k\neq i:i,j,k\in {1,\ldots,n}\}
      f.con <- do.call(rbind, lapply(as.matrix(1:(n - 1)), function(ii) fCmat(n, i, ii)))
      ## constraints inequlities
      f.dir <- rep("<=", nrow(f.con))
      ## RHS of the constraints
      f.rhs <- c(t(do.call(rbind, lapply(as.matrix(1:(n - 1)), function(ii) fRHS(X, i, ii)))))    
      ## solve the binary integer programming problem
      LP   <- lp(direction = "max", objective.in = f.obj, const.mat = f.con, 
                 const.dir = f.dir, const.rhs = f.rhs, all.bin = TRUE)
      xnew <- LP$solution
      ## update results into X
      for (k in 1:(n - 1)) {
        if (k < i) {  ## (1,i), (2,i), ..., (i-1,i)
          X[k, i] <- xnew[k] 
        } else {      ## (i,i+1), (i,i+2), ..., (i,n)
          X[i, k + 1] <- xnew[k]
        }
      }
    }
    valnew <- objfeval(X, rho, K)
    if (valnew > val) {  
      ## val < valnew  
      val <- valnew 
    } else { 
      break 
    }
  }
  ## symmetrize X
  for (i in 2:n) {
    for(ii in 1:(i - 1)) {
      X[i, ii] <- X[ii, i]
    }
  }
  ## return
  list(partition = X[lower.tri(X)],
       val       = valnew)
}
 
LGalg_R <- function (Kgrid, rho) 
{
  ## estimates the incidence matrix
  ngrid <- length(Kgrid)
  I     <- nrow(rho)
  ovals <- rep(NA, ngrid)
  parts <- matrix(NA, ngrid, I*(I - 1)/2) 
  for (k in 1:ngrid) {
    tmp <- LauGreenAlg0(Kgrid[k], rho)
    ovals[k]   <- tmp$val
    parts[k, ] <- tmp$partition 
    cat("Kgrid = ", k, " ... \n", sep = "")
  }
  as.matrix(parts[which.max(ovals)[1], ])
}