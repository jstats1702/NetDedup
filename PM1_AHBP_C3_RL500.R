rm(list = ls (all.names = TRUE))

####################
# general settings #
####################

#loc     <- "C:/RL1_paper/"
loc     <- "/soe/jsosamar/RL1_paper/"

model   <- "PM1"

prior   <- "AHBP" 

case    <- "C3"

dataset <- "RL500"


###### paths
path.code <- paste0(loc, "code/")
path.data <- paste0(loc, "data/")
path.outs <- paste0(loc, "outputs/", model, "_", prior, "_", case, "_", dataset, "/")
if (!dir.exists(path.outs)) dir.create(path.outs)


###### data
load(paste0(path.data, "data_", dataset, ".RData"))


##### libs and source
library(Rcpp)
source   (paste0(loc, "code/matching_R.R"))
sourceCpp(paste0(path.code, model, "_", prior, "_Cfunctions.cpp"))


##### mcmc
E     <- 1/100
a_psi <- 1
b_psi <- a_psi * (1/E - 1)

nburn <- 5e+05
nsams <- 1e+05
nskip <- 1
ndisp <- floor(0.01 * (nburn + nskip * nsams))

Q2    <- floor(I/2)
Er1.0 <- floor(I/2)
E     <- ( I - Er1.0 )/( 2 * Q2 )
CV    <- 0.5
a_vt  <- ( 1/E - 1 - CV^2 )/( (1 + 1/E - 1) * CV^2 )
b_vt  <- a_vt * (1/E - 1)

ptm <- proc.time()
set.seed(1)
#MCMC (a_psi, b_psi, Xi, a_vt, b_vt, Mell, Pe, nburn, nsams, nskip, ndisp, path.outs)
proc.time() - ptm


##### inference
ptm <- proc.time()
cstats_R (model, prior, case, dataset, loc)
proc.time() - ptm

ptm <- proc.time()
matching0_R(Xi, model, prior, case, dataset, loc)
proc.time() - ptm