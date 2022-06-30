datagen_NM1_R <- function (I, N, K, sigsq, beta, Xi) 
{
     ## NETWORK parameters
     #Xi    <- gen_Xi_R(I, N)
     Ustar <- matrix(rnorm(N * K, 0.0, sqrt(sigsq)), N, K)
     for (k in 1:K) Ustar[ , k] <- Ustar[ , k] - mean(Ustar[ , k])  ## center values
     U <- get_U_R(Ustar, Xi)
     ## NETWORK data  
     NN    <- I*(I - 1)/2  ## n dyads
     Theta <- matrix(NA, NN, 1)
     Y     <- matrix(NA, NN, 1)
     for (i in 2:I) {  ## lower triangular indices
          for (ii in 1:(i-1)) {
               m        <- dyx_R(i, ii, I)
               eta      <- beta - sqrt(sum( (U[i,] - U[ii,])^2 ))
               Theta[m] <- 1/(1 + exp(-eta)) 
               Y[m]     <- rbinom(1, 1, Theta[m])
          }
     }
     ## log-likelihood
     ll <- ll_iter_net_R(I, beta, U, Y) 
     ## return
     list (sigsq = sigsq,
           beta  = beta,
           Ustar = Ustar,
           Xi    = Xi,
           Y     = Y,
           Theta = Theta,
           ll    = ll,
           I     = I,
           N     = N,
           K     = K)
}