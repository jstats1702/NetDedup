gen_Vartheta_R <- function(Mell, Alpha)
{
   require(gtools)
   L <- length(Mell)
   M <- sum(Mell)
   Vartheta <- matrix(NA, M, 1)
   for (l in 1:L) {
      idx <- (sum(Mell[1:l]) - Mell[l] + 1):sum(Mell[1:l])
      Vartheta[idx] <- rdirichlet(1, Alpha[idx])
   }
   Vartheta
}

gen_W_R <- function(I, Psi)
{
   L <- length(Psi)
   W <- matrix(NA, I, L)
   for (l in 1:L) {
      for (i in 1:I) {
         W[i, l] <- rbinom(n = 1, size = 1, prob = Psi[l])
      }
   }
   W
}

gen_Pi_R <- function(N, L, Mell, Vartheta)
{
   Pi <- matrix(NA, N, L)
   for (l in 1:L) {
      vt <- Vartheta[(sum(Mell[1:l]) - Mell[l] + 1):sum(Mell[1:l])]
      for(n in 1:N)
         Pi[n, l] <- which(rmultinom(n = 1, size = 1, prob = vt) == 1)
   }
   Pi
}

datagen_PNM1_R <- function(I, N, K, sigsq, beta, Mell, Psi, Alpha) 
{
   ## NETWORK parameters
   Xi    <- gen_Xi_R(I, N)
   Ustar <- matrix(rnorm(N * K, 0.0, sqrt(sigsq)), N, K)
   for (k in 1:K) Ustar[ , k] <- Ustar[ , k] - mean(Ustar[ , k])  ## center values
   U     <- get_U_R(Ustar, Xi)
   ## PROFILE parameters
   L        <- length(Mell)
   Vartheta <- gen_Vartheta_R(Mell, Alpha)
   W        <- gen_W_R(I, Psi)
   Pi       <- gen_Pi_R(N, L, Mell, Vartheta)
   ## PROFILE data
   Pe <- matrix(data = NA, nrow = I, ncol = L)
   for (l in 1:L) {
        vt <- Vartheta[(sum(Mell[1:l]) - Mell[l] + 1):sum(Mell[1:l])]
        for (i in 1:I) {
             if (W[i, l] == 0) {
                  Pe[i, l] <- Pi[Xi[i], l]  ## NO distortion!
             } else {  # w_{i,l,j} = 1
                  Pe[i, l] <- which(rmultinom(n = 1, size = 1, prob = vt) == 1)
             }
        }
        rm(vt)
   }
   ## NETWORK data  
   NN    <- I*(I - 1)/2  ## n dyads
   Theta <- matrix(NA, NN, 1)
   Y     <- matrix(NA, NN, 1)
   for (i in 2:I) {  ## lower triangular indices
      for (ii in 1:(i-1)) {
         m        <- dyx_R(i, ii, I)
         eta      <- beta - sqrt(sum((Ustar[Xi[i], ] - Ustar[Xi[ii], ])^2))
         Theta[m] <- 1/(1 + exp(-eta)) 
         Y[m]     <- rbinom(1, 1, Theta[m])
      }
   }
   ## log-likelihood
   ll <- ll_iter_pro_R(Mell, Vartheta, W, Pe) + ll_iter_net_R(I, beta, U, Y)
   ## return
   list (sigsq    = sigsq,
         beta     = beta,
         Ustar    = Ustar,
         Xi       = Xi,
         Theta    = Theta,
         Y        = Y,
         Alpha    = Alpha,
         Psi      = Psi,
         Vartheta = Vartheta,
         W        = W,
         Pi       = Pi,
         Pe       = Pe,
         ll       = ll,
         I        = I,
         N        = N,
         K        = K,
         Mell     = Mell,
         L        = L)
}