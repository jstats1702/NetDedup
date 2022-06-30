rm(list = ls (all.names = TRUE))

####################
# general settings #
####################

#loc     <- "C:/RL1_paper/"
loc     <- "/soe/jsosamar/RL1_paper/"

model   <- "NM1"

prior   <- "AHBP" 

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
nburn <- 1e+06
nsams <- 50000
nskip <- 20
ndisp <- floor(0.01 * (nburn + nskip * nsams))

Q2    <- floor(I/2)
Er1.0 <- 2 * N - I    
E     <- ( I - Er1.0 )/( 2 * Q2 )
CV    <- 0.05
a_vt  <- ( 1/E - 1 - CV^2 )/( (1 + 1/E - 1) * CV^2 )
b_vt  <- a_vt * (1/E - 1)

ptm <- proc.time()
set.seed(1)
MCMC (Xi, a_vt, b_vt, K, Y, nburn, nsams, nskip, ndisp, path.outs)
proc.time() - ptm


##### log-likelihood plot
ll.data <- as.matrix(read.table(paste0(path.outs, "ll_chain.txt")))

pdf(paste0(path.outs, "loglik.pdf"), height = 5, width = 10, pointsize = 12)
par(mfrow = c(1, 1), mar = c(4, 4, 3, 2) - 0.0, mgp = c(2, 1, 0), oma = c(0, 0, 0, 0))
plot(ll.data, type = "l", ylab = "Log-likelihood", main = " ", col = "gray")
dev.off()


##### inference
ptm <- proc.time()
cstats_R (model, prior, case, dataset, loc)
proc.time() - ptm

ptm <- proc.time()
matching0_R(Xi, model, prior, case, dataset, method = "PM", loc)
proc.time() - ptm

ptm <- proc.time()
matching0_R(Xi, model, prior, case, dataset, method = "LG", loc)
proc.time() - ptm

ptm <- proc.time()
matching0_R(Xi, model, prior, case, dataset, method = "MPMMS", loc)
proc.time() - ptm