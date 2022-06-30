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

logpr2_unif_R <- function(I, r2) {
     if (r2 == 0) {
          out <- 0
     } else {
          if (r2 == floor(I/2)) {
               out <- sum(log(1:I)) - r2*log(2) - sum(log(1:r2))
          } else {
               out <- sum(log(1:I)) - r2*log(2) - sum(log(1:r2)) - sum(log(1:(I-2*r2)))
          }
     }
     out
}

gen_Xi_allelic_R <- function (I, r2)
{
     r2counts <- rep(0, 2*r2)
     r2s      <- sample(x = 1:I, size = 2*r2, replace = F)
     Xi       <- rep(NA, I)
     if (r2 > 0) {
          for (i in 1:r2) {
               tmp <- sample(x = r2s[r2counts < 1], size = 2, replace = F)
               Xi[tmp] <- i
               r2counts[r2s %in% tmp] <- r2counts[r2s %in% tmp] + 1
          }
     }
     Xi[is.na(Xi)] <- (r2+1):(I-r2)
     c(relabel_R(Xi))
}

gen_Xi_EPP_R <- function (I, theta)
{
     K  <- 1
     nk <- c(1)
     Xi <- c(1)
     for (j in 2:I) {
          pk <- nk/(j-1+theta)
          pp <- theta/(j-1+theta)
          s  <- sample(x = c(1:K, K+1), size = 1, replace = F, prob = c(pk, pp))
          Xi[j] <- s
          if (s == (K+1)) {
               K  <- K+1
               nk <- c(nk, 1)
          } else {
               nk[s] <- nk[s] + 1
          }
     }
     c(relabel_R(Xi))
}

gen_Xi_NBNBP_R <- function (a, q, r, p,  M, I) 
{
     ## Sampler for NBNB prior
     ## Produces one sample from NBNB
     ## Needs M gibbs samples bc sampling is not exact
     w <- (q * (1-p)^r) / (1 - (1-p)^r)
     cluster <- 1:I
     for (j in 1:M) {
          clusterA <- cluster
          for (i in 1:I) {
               clusterA[i]  <- 0
               clusterA[-i] <- as.numeric(factor(cluster[-i], labels = seq(1, length(unique(cluster[-i])))))
               nk <- table(clusterA[-i])
               clusterA[i] <- sample(x = 1:(max(clusterA) + 1), size = 1, prob = c(r + nk, (max(clusterA) + a) * w * r))
               cluster <- clusterA
          }
     }
     return( cluster )
}