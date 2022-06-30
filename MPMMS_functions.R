MPMMSAlg_R <- function (Xi.data, threshold)
{
     S <- dim(Xi.data)[1]
     I <- dim(Xi.data)[2]
     ## Most Probable MMS (MPMMS)
     MPMMS <- vector("list", I)
     for (i in 1:I) {
          #maXimal matching set (MMS)
          MMS <- vector("list", S)
          for (s in 1:S) MMS[[s]] <- which(Xi.data[s, ] == Xi.data[s,i])
          uMMS <- unique(MMS)
          nMMS <- length(uMMS)
          freq <- rep(NA, nMMS)
          for (k in 1:nMMS) {
               cont <- 0
               for (l in 1:S) if(identical(uMMS[k], MMS[l])) cont <- cont + 1
               freq[k] <- cont
          }
          #maxMSSidx  <- which.max(freq)
          #MPMMS[[i]] <- uMMS[[which.max(unlist(lapply(uMMS[maxMSSidx], FUN = length)))]]
          maxMSSidx <- which(freq/S > threshold)
          if (length(maxMSSidx) > 0) {
               MPMMS[[i]] <- uMMS[[which.max(unlist(lapply(uMMS[maxMSSidx], FUN = length)))]]
          } else {
               MPMMS[[i]] <- as.integer(c(i))               
          }
          if (i%%floor(0.1*I) == 0) cat(round(100*i/I,1), "% done ... \n", sep = "")
     }
     ## linkage structure
     link.str <- rep(NA, I)
     n <- 1
     for (i in 1:I) {
          if (is.na(link.str[i])) {
               for (ii in 1:I) if (identical(MPMMS[[i]], MPMMS[[ii]])) link.str[ii] <- n 
               n <- n + 1
          }
     }
     ## clusts
     clusts <- matrix(0, I, I)
     for (i in 2:I) {
          for (j in 1:(i - 1)) {
               if (link.str[i] == link.str[j]) clusts[i,j] <- 1
          }
     }
     clusts <- clusts + t(clusts)
     clusts[lower.tri(clusts)]
}