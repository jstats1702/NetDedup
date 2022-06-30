rm(list = ls (all.names = TRUE))

####################
# general settings #
####################

loc <- "C:/RL1_paper/"
#loc <- "/soe/jsosamar/RL1_paper/"

dataset <- "RL500"


###### paths
path.code <- paste0(loc, "code/")
path.outs <- paste0(loc, "outputs/")


##### libs and source
library(Rcpp)
#library(gtools)

source   (paste0(path.code, "prior_Rfunctions.R"))

sourceCpp(paste0(path.code, "prior_Cfunctions.cpp"))

## 1  : AUP
## 2  : AHUP
## 3  : AHBP C1
## 4  : AHBP C2
## 5  : AHBP C3
## 6  : EPP  C1
## 7  : EPP  C2
## 8  : EPP  C3
## 9  : NBNBP
## 10 : NBDP


##### sim pars
nsim <- 5000
M    <- 250
I    <- 500
N    <- 450

p    <- 10
Q2   <- floor(I/2)
S.data  <- matrix(NA, nsim, p)
S.distr <- matrix(NA, I+1,  p)


##### 1 : AUP
j <- 1

q  <- rep(NA, Q2+1)
for (i in 0:Q2) q[i+1] <- logpr2_unif_R(I, i)
q <- exp(q - max(q))/sum(exp(q - max(q)))

for(i in 1:nsim) {
     set.seed(i)
     r2 <- sample(x = 0:Q2, size = 1, replace = F, prob = q)
     Xi <- gen_Xi_allelic_R(I, r2)
     S.data[i,j] <- sum(table(Xi) == 1)
     rm(r2, Xi)
     if (i%%(0.1*nsim) == 0) cat(round(100*i/nsim, 2), "% done ...", "\n", sep = "") 
}

#for (i in 0:I) S.distr[i+1,1] <- sum(S.data[,1] == i)/nsim
#Er1.AUP <- sum( 0:I * S.distr[,1] )
Er1.AUP <- I - 2 * sum( q * (0:Q2) )  ## E[No. of singletons] under unif
Er1.0   <- 2 * N - I                  ## true number of singletons


##### 2 : AHUP
j <- 2

for(i in 1:nsim) {
     r2 <- floor(runif(n = 1, min = 0, max = Q2 + 1))
     Xi <- gen_Xi_allelic_R(I, r2)
     S.data[i,j] <- sum(table(Xi) == 1)
     rm(Xi, r2)
     if (i%%(0.1*nsim) == 0) cat(round(100*i/nsim, 2), "% done ...", "\n", sep = "") 
}


##### 3 : AHBP C1
j <- 3

E    <- ( I - Er1.0 )/( 2 * Q2 )
CV   <- 0.05
a_vt <- ( 1/E - 1 - CV^2 )/( (1 + 1/E - 1) * CV^2 )
b_vt <- a_vt * (1/E - 1)

for(i in 1:nsim) {
     vt <- rbeta(n = 1, shape1 = a_vt, shape2 = b_vt)
     r2 <- rbinom(n = 1, size = Q2, prob = vt)
     Xi <- gen_Xi_allelic_R(I, r2)
     S.data[i,j] <- sum(table(Xi) == 1)
     rm(Xi, r2, vt)
     if (i%%(0.1*nsim) == 0) cat(round(100*i/nsim, 2), "% done ...", "\n", sep = "") 
}


##### 4 : AHBP C2
j <- 4

E    <- ( I - Er1.0 )/( 2 * Q2 )
CV   <- 0.5
a_vt <- ( 1/E - 1 - CV^2 )/( (1 + 1/E - 1) * CV^2 )
b_vt <- a_vt * (1/E - 1)

for(i in 1:nsim) {
     vt <- rbeta(n = 1, shape1 = a_vt, shape2 = b_vt)
     r2 <- rbinom(n = 1, size = Q2, prob = vt)
     Xi <- gen_Xi_allelic_R(I, r2)
     S.data[i,j] <- sum(table(Xi) == 1)
     rm(Xi, r2, vt)
     if (i%%(0.1*nsim) == 0) cat(round(100*i/nsim, 2), "% done ...", "\n", sep = "") 
}


##### 5 : AHBP C3
j <- 5

E    <- 0.5
CV   <- 0.5
a_vt <- ( 1/E - 1 - CV^2 )/( (1 + 1/E - 1) * CV^2 )
b_vt <- a_vt * (1/E - 1)

for(i in 1:nsim) {
     vt <- rbeta(n = 1, shape1 = a_vt, shape2 = b_vt)
     r2 <- rbinom(n = 1, size = Q2, prob = vt)
     Xi <- gen_Xi_allelic_R(I, r2)
     S.data[i,j] <- sum(table(Xi) == 1)
     rm(Xi, r2, vt)
     if (i%%(0.1*nsim) == 0) cat(round(100*i/nsim, 2), "% done ...", "\n", sep = "") 
}


##### 6 : EPP C1
j <- 6

f   <- function (x, I, Er1) 2*x*log((x + I)/x) - (I + Er1)
E   <- uniroot(f = f, lower = 1e-06, upper = 1e+06, I = I, Er1 = Er1.0)$root
CV  <- 0.05
a_t <- 1/CV^2
b_t <- 1/(E * CV^2)

for(i in 1:nsim) {
     theta <- rgamma(n = 1, shape = a_t, rate = b_t)
     Xi <- gen_Xi_EPP_R(I, theta)
     S.data[i,j] <- sum(table(Xi) == 1)
     rm(Xi, theta)
     if (i%%(0.1*nsim) == 0) cat(round(100*i/nsim, 2), "% done ...", "\n", sep = "") 
}


##### 7 : EPP C2
j <- 7

f   <- function (x, I, Er1) 2*x*log((x + I)/x) - (I + Er1)
E   <- uniroot(f = f, lower = 1e-06, upper = 1e+06, I = I, Er1 = Er1.0)$root
CV  <- 0.5
a_t <- 1/CV^2
b_t <- 1/(E * CV^2)

for(i in 1:nsim) {
     theta <- rgamma(n = 1, shape = a_t, rate = b_t)
     Xi <- gen_Xi_EPP_R(I, theta)
     S.data[i,j] <- sum(table(Xi) == 1)
     rm(Xi, theta)
     if (i%%(0.1*nsim) == 0) cat(round(100*i/nsim, 2), "% done ...", "\n", sep = "") 
}


##### 8 : EPP C2
j <- 8

f   <- function (x, I, Er1) 2*x*log((x + I)/x) - (I + Er1)
E   <- uniroot(f = f, lower = 1e-06, upper = 1e+06, I = I, Er1 = Q2)$root
CV  <- 0.5
a_t <- 1/CV^2
b_t <- 1/(E * CV^2)

for(i in 1:nsim) {
     theta <- rgamma(n = 1, shape = a_t, rate = b_t)
     Xi <- gen_Xi_EPP_R(I, theta)
     S.data[i,j] <- sum(table(Xi) == 1)
     rm(Xi, theta)
     if (i%%(0.1*nsim) == 0) cat(round(100*i/nsim, 2), "% done ...", "\n", sep = "") 
}


##### 9 : NBNBP
j <- 9

E     <- I/2
a     <- 1
q     <- 1 - 1/E
a_eta <- 1
b_eta <- 1
a_t   <- 2
b_t   <- 2

for(i in 1:nsim) {
     eta   <- rgamma(n = 1, shape = a_eta, rate = b_eta) 
     theta <- rbeta(n = 1, shape1 = a_t, shape2 = b_t)
     Xi <- gen_Xi_NBNBP_C(a, q, eta, theta, M, I)
     S.data[i,j] <- sum(table(Xi) == 1)
     rm(Xi, eta, theta)
     if (i%%(0.1*nsim) == 0) cat(round(100*i/nsim, 2), "% done ...", "\n", sep = "") 
}


##### 10 : NBDP
j <- 10

E     <- I/2
a     <- 1
q     <- 1 - 1/E
al    <- 1
Mu0   <- dgeom(x = 1:I, prob = 0.5)/(1 - dgeom(x = 0, prob = 0.5))
Mu0   <- c( Mu0, 1 - sum(Mu0) )

for(i in 1:nsim) {
     Xi <- gen_Xi_NBDP_C(a, q, al, Mu0, M, I) 
     S.data[i,j] <- sum(table(Xi) == 1)
     rm(Xi)
     if (i%%(0.1*nsim) == 0) cat(round(100*i/nsim, 2), "% done ...", "\n", sep = "") 
}


save(S.data, file = paste0(path.outs, "prior_sim_", dataset, ".Rdata"))
load(paste0(path.outs, "prior_sim_", dataset, ".Rdata"))


for (i in 0:I) for (j in 1:p) S.distr[i+1,j] <- sum(S.data[,j] == i)/nsim

ES <- rep(NA, p)
for (j in 1:p) ES[j] <- sum( 0:I * S.distr[,j] )

prior.labs <- c("Uniform","AHUP","AHBP - C1", "AHBP - C2", "AHBP - C3", "EPP - C1", "EPP - C2", "EPP - C3", "NBNBP","NBDP")

##### histograms
pdf(paste0(path.outs, "prior_sim_", dataset, ".pdf"), width = 25, height = 10, pointsize = 20)
#windows(width = 20, height = 10)
par(mfrow = c(2, 5), mar = c(4, 3, 3, 2) + 1.0, mgp = c(2, 1, 0), oma = c(0, 0, 0, 0))
for (j in 1:p) {
     hist(S.data[,j], freq = F, breaks = 25, xlim=c(0,I), ylim = range(S.distr, na.rm = T), cex.axis = 0.5, las = 1, 
          xlab = "Number of singletons", ylab = "Probability", main = paste0(prior.labs[j], "\n Mean = ", round(ES[j], 1)))
}
dev.off()
