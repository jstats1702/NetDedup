rm(list = ls (all.names = TRUE))

####################
# general settings #
####################

#loc     <- "C:/RL1_paper/"
loc     <- "/soe/jsosamar/RL1_paper/"

model   <- "PM1"

prior   <- "EPP" 

case    <- "C1"

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

Er1.0 <- 2 * N - I
f     <- function (x, I, Er1) 2*x*log((x + I)/x) - (I + Er1)
E     <- uniroot(f = f, lower = 1e-06, upper = 1e+06, I = I, Er1 = Er1.0)$root
CV    <- 0.05
a_t   <- 1/CV^2
b_t   <- 1/(E * CV^2)

ptm <- proc.time()
set.seed(1)
MCMC (a_psi, b_psi, Xi, a_t, b_t, Mell, Pe, nburn, nsams, nskip, ndisp, path.outs)
proc.time() - ptm


##### inference
ptm <- proc.time()
cstats_R (model, prior, case, dataset, loc)
proc.time() - ptm

ptm <- proc.time()
matching0_R(Xi, model, prior, case, dataset, loc)
proc.time() - ptm