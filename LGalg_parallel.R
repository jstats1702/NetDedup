LGalg_parallel_R <- function (Kgrid, rho, path.code) 
{
  ## estimates the incidence matrix
  require(doParallel)
  require(foreach)
  rho   <- as.matrix(rho)
  ngrid <- length(Kgrid)
  ## set cores and cluster
  dc <- detectCores()            ## detect number of cores
  nc <- min(dc, ngrid)           ## set working number of cores
  cl <- makeCluster(nc)          ## create cluster with desired number of cores
  registerDoParallel(cl)         ## register cluster
  parts <- foreach (k = 1:ngrid, .combine = rbind, .packages = "lpSolve") %dopar% {
    source(paste0(path.code, "LG_functions.R"))
    tmp <- LauGreenAlg0(Kgrid[k], rho)
    c(tmp$val, tmp$partition)  # objective value and partition vec by cols
  }
  stopCluster(cl)
  ## return
  ovals <- parts[ ,  1]
  parts <- parts[ , -1]
  parts[which.max(ovals)[1], ]
}