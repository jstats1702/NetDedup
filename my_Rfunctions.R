dyx_R <- function (i, ii, I)
{ 
     ## returns (vectorized) linear index from a lower triangular matrix position
     I*(I - 1)/2 - (I - ii + 1)*(I - ii)/2 + i - ii
}

procrus_R <- function (Z, Z0) 
{
     ## Procrustes transform; gives rotation,reflection,trranslation of Z closest to Z0
     ## Author: Peter Hoff
     K <- ncol(Z)
     for (i in 1:K) Z[ , i] <- Z[ , i] - mean(Z[ , i]) + mean(Z0[ , i])  ## translation
     A  <- t(Z) %*% Z0 %*% t(Z0) %*% Z
     eA <- eigen(A, symmetric = T)
     Ahalf <- eA$vec[ , 1:K] %*% diag(sqrt(eA$val[1:K])) %*% t(eA$vec[ , 1:K])
     Z %*% solve(Ahalf) %*% t(Z) %*% Z0
     #t(t(Z0) %*% Z %*% solve(Ahalf) %*% t(Z)) 
}

vec_to_symmat_R <- function (x)
{
     ## returns a symmetric matrix given the vec lower triangular part of it 
     ## ones in the main diagonal 
     N <- length(x)
     n <- (1 + sqrt(8*N + 1))/2
     out <- matrix(0, n, n)
     out[lower.tri(out)] <- c(x)
     out <- out + t(out)
     diag(out) <- 1
     out
}

relabel_R <- function (Xi) 
{
     ## re-labels latent individuals from 1 to N
     ## where N is the number of unique elements in Xi
     Xi0 <- Xi  # Xi copy
     uXi0 <- sort(unique(Xi0))
     N <- length(uXi0)  # n latent inidviduals 
     Xi <- matrix(NA, length(Xi), 1)
     for (n in 1:N) 
          Xi[Xi0 == uXi0[n]] <- n
     Xi
}

gen_Xi_R <- function (I, N)
{
     actors <- 1:I
     counts <- rep(0, I)
     Xi     <- matrix(NA, I, 1)
     for (n in 1:N) {
          i     <- sample(x = actors[is.na(Xi)],  size = 1, replace = F)
          Xi[i] <- sample(x = actors[counts < 1], size = 1, replace = F)
          counts[Xi[i]] = counts[Xi[i]] + 1 
     }
     for (n in 1:(I-N)) {
          if (length(actors[is.na(Xi)]) > 1)
               i <- sample(x = actors[is.na(Xi)], size = 1, replace = F)
          else
               i <- actors[is.na(Xi)]
          Xi[i] <- sample(x = actors[counts == 1], size = 1, replace = F)
          counts[Xi[i]] = counts[Xi[i]] + 1 
     }
     relabel_R(Xi)
}

get_Ustar_R <- function (U, Xi) {
     ## creates U* (N x K), where N is the number of latent individuals
     Xi <- relabel_R(Xi)
     N <- max(Xi)
     K <- ncol(U)
     Ustar <- matrix(NA, N, K)
     for (n in 1:N) {
          id <- which(Xi == n)
          Ustar[n, ] <- U[id[1], ]
     }
     Ustar
}

get_U_R <- function (Ustar, Xi) {
     ## creates U (I x K), where I is the number of actors
     uXi <- unique(Xi)
     N <- length(uXi)
     U <- matrix(NA, length(Xi), ncol(Ustar))
     for (n in 1:N) {
          id <- which(Xi == n)
          for (kk in 1:length(id)) 
               U[id[kk], ] <- Ustar[n, ]
     }
     U
}

get_matches_R <- function (Xi)
{
     I <- length(Xi)
     N <- length(unique(Xi))
     nmatch  <- I - N
     matches <- matrix(NA, nmatch, 2)
     cont <- 1
     for (n in 1:N) {
          id <- which(Xi == n)
          if (length(id) > 1) {
               matches[cont, ] <- c(id)
               cont <- cont + 1
          }
     }
     matches[order(matches[ , 1]), ]
}

ll_iter_net_R <- function (I, beta, U, Y) 
{
     ## computes the log-likelihood in a given iteration of the MCMC
     out <- 0
     for (i in 2:I) {  ## lower triangular indices
          for (ii in 1:(i-1)) {
               eta <- beta - sqrt(sum( (U[i, ] - U[ii, ])^2 ))
               out <- out + eta * Y[dyx_R(i, ii, I)] - log(1 + exp(eta))
          }
     }
     out
}

ll_iter_pro_R <- function(Mell, Vartheta, W, Pe) 
{
     ## computes the log-likelihood in a given iteration of the MCMC
     ## profile part of log-likelihood
     I <- nrow(Pe)
     L <- ncol(Pe)
     out <- 0.0
     for (l in 1:L) {
          lo <- sum(Mell[1:l]) - Mell[l]
          for (i in 1:I) {
               if (W[i, l] == 1) {
                    out = out + log(Vartheta[lo + Pe[i, l]])
               }
          }
     }
     out
}

lower_tri_keep_R <- function (I, Samip) 
{
     out <- rep(FALSE, I*(I-1)/2)
     for (i in 1:I) {
          if (i > 2) {
               if (i %in% Samip) {
                    for (ii in 1:(i-1))
                         out[dyx_R(i, ii, I)] <- TRUE
               }
          }
          if (i < I) {
               if (i %in% Samip) {
                    for (ii in (i+1):I)
                         out[dyx_R(ii, i, I)] <- TRUE
               }
          }
     }
     out
}

accuracy_R <- function(truth, predicted)
{
     ## Compute the Matthews correlation coefficient (MCC) score
     ## Jeff Hebert 9/1/2016
     ## truth = vector of true outcomes, 1 = Positive, 0 = Negative
     ## predicted = vector of predicted outcomes, 1 = Positive, 0 = Negative
     ## function returns MCC
     ## F_{2} weighs recall higher than precision 
     ## F_{0.5} weighs recall lower than precision
     TP <- sum(truth == 1 & predicted == 1)
     TN <- sum(truth == 0 & predicted == 0)
     FP <- sum(truth == 0 & predicted == 1)
     FN <- sum(truth == 1 & predicted == 0)
     ## measures
     ##MCC <- ( TP*TN - FP*FN ) / exp( 0.5*(log(TP + FP) + log(TP + FN) + log(TN + FP) + log(TN + FN)) )
     NN  <- TN + TP + FN + FP
     SS  <- (TP + FN)/NN
     PP  <- (TP + FP)/NN
     MCC <- ( TP/NN - SS * PP ) / exp( 0.5 * (log(PP) + log(SS) + log(1-SS) + log(1-PP)) ) 
     P   <- TP / (TP + FP)  
     R   <- TP / (TP + FN)
     FDR <- FP / (TP + FP)
     FNR <- FN / (TP + FN)
     Fh  <- ( 1 + 0.5^2 ) * TP / ( (1 + 0.5^2) * TP + 0.5^2 * FN + FP )
     F1  <- ( 1 + 1.0^2 ) * TP / ( (1 + 1.0^2) * TP + 1.0^2 * FN + FP )  
     F2  <- ( 1 + 2.0^2 ) * TP / ( (1 + 2.0^2) * TP + 2.0^2 * FN + FP )  
     list (MCC = MCC,
           P   = P,
           R   = R,
           FDR = FDR,
           FNR = FNR,
           Fh  = Fh,
           F1  = F1,
           F2  = F2)
}

load_matching_results0 <- function (model, prior, case, dataset, method, loc)
{
     path.outs <- paste0(loc, "outputs/", model, "_", prior, "_", case, "_", dataset, "/")
     load(file = paste0(path.outs, model, "_", prior, "_", case, "_", dataset, "_matching_results_", method, ".RData"))
     out
}

load_matching_results1 <- function (model, dataset, method, loc)
{
     rbind(
          load_matching_results0(model, prior = "AUP",   case = "UC", dataset, method, loc),
          #load_matching_results0(model, prior = "AHUP",  case = "UC", dataset, method, loc),
          #load_matching_results0(model, prior = "AHBP",  case = "C1", dataset, method, loc),
          load_matching_results0(model, prior = "AHBP",  case = "C2", dataset, method, loc),
          load_matching_results0(model, prior = "AHBP",  case = "C3", dataset, method, loc),
          load_matching_results0(model, prior = "AHBP3", case = "C2", dataset, method, loc),
          #load_matching_results0(model, prior = "EPP",   case = "C1", dataset, method, loc),
          load_matching_results0(model, prior = "EPP",   case = "C2", dataset, method, loc),
          load_matching_results0(model, prior = "EPP",   case = "C3", dataset, method, loc),
          load_matching_results0(model, prior = "NBDP",  case = "UC", dataset, method, loc),
          load_matching_results0(model, prior = "NBNBP", case = "UC", dataset, method, loc) 
     )
}

load_cstats_results0 <- function (model, prior, case, dataset, loc)
{
     path.outs <- paste0(loc, "outputs/", model, "_", prior, "_", case, "_", dataset, "/")
     load(file = paste0(path.outs, model, "_", prior, "_", case, "_", dataset, "_cstats_results.RData"))
     out
}

load_cstats_results1 <- function (model, dataset, loc)
{
     rbind(
          load_cstats_results0(model, prior = "AUP",   case = "UC", dataset, loc),
          #load_cstats_results0(model, prior = "AHUP",  case = "UC", dataset, loc),
          #load_cstats_results0(model, prior = "AHBP",  case = "C1", dataset, loc),
          load_cstats_results0(model, prior = "AHBP",  case = "C2", dataset, loc),
          load_cstats_results0(model, prior = "AHBP",  case = "C3", dataset, loc),
          load_cstats_results0(model, prior = "AHBP3",  case = "C2", dataset, loc),
          #load_cstats_results0(model, prior = "EPP",   case = "C1", dataset, loc),
          load_cstats_results0(model, prior = "EPP",   case = "C2", dataset, loc),
          load_cstats_results0(model, prior = "EPP",   case = "C3", dataset, loc),
          load_cstats_results0(model, prior = "NBNBP", case = "UC", dataset, loc),
          load_cstats_results0(model, prior = "NBDP",  case = "UC", dataset, loc)
     )
}