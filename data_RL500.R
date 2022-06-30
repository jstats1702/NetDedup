rm(list = ls (all.names = TRUE))

####################
# general settings #
####################

loc     <- "C:/RL1_paper/"

model   <- "PNM1"

dataset <- "RL500"


##### paths
path.code <- paste0(loc, "code/")
path.data <- paste0(loc, "data/")


##### libs and source
library(RecordLinkage)

source(paste0(path.code, "my_plots.R"))
source(paste0(path.code, "my_Rfunctions.R"))
source(paste0(path.code, "NM1_datagen.R"))


###################
# data generation #
###################

data(RLdata500)

#write.table(x = RLdata500[ , c(1:7)], file = paste0(path.data, dataset, "_data0.txt"), quote = F, row.names = F)
#NA fields included as zeros
RLdata500 <- read.table(file = paste0(path.data, dataset, "_data0.txt"), header = T)
RLdata500 <- RLdata500[ , c(5, 6, 7)]

I  <- dim(RLdata500)[1]
L  <- dim(RLdata500)[2]
Xi <- as.matrix(identity.RLdata500 + 1)
#Xi <- relabel_R(Xi)
N  <- length(unique(Xi))

Mell <- matrix(NA, L, 1)
for (l in 1:L) Mell[l] <- length(unique(RLdata500[ , l]))


##### TRUE LINKAGE STRUCTURE
matches <- get_matches_R(Xi)
nmatch  <- nrow(matches)


##### PROFILE DATA (CODIFIED)
Pe0 <- RLdata500
Pe  <- matrix(NA, I, L)
for(l in 1:L) {
     nivs <- unique(Pe0[ , l])
     for (m in 1:length(nivs)) {
          cont <- which(nivs == nivs[m])
          for (i in 1:I) {
               if (Pe0[i, l] == nivs[m]) Pe[i, l] <- cont
          }
     }
}
rm(Pe0, nivs, cont)


##### NETWORK DATA
K      <- 2
sigsq  <- 178
beta   <- 10

set.seed(1)
mydata <- datagen_NM1_R(I, N, K, sigsq, beta, Xi)
Ustar  <- mydata$Ustar
Theta  <- mydata$Theta
Y      <- mydata$Y
U      <- get_U_R(Ustar, Xi)

## connect with a single link unconnected records
Yadj  <- matrix(0, I, I)
Yadj [lower.tri(Yadj)]  <- c(Y)
Yadj  <- Yadj  + t(Yadj)

zro <- (1:I)[ Yadj %*% rep(1, I) == 0 ]
if (length(zro) > 0) {
     for (i in zro) {
          D  <- rep(0, I)
          for (ii in 1:I) D[ii] <- sqrt(sum( (U[i, ] - U[ii, ])^2 ))
          idx <- which((D > 0) & (D == min(D[-i])))[1]
          Yadj[i, idx] <- 1
          Yadj[idx, i] <- 1
     }
}
rm(zro, i, D, idx)

isSymmetric(Yadj)

Y <- Yadj[lower.tri(Yadj)]

save(Xi, matches, nmatch, Pe, beta, U, Ustar, Y, Theta, Mell, I, N, L, K, file = paste0(path.data, "data_", dataset, ".RData"))


##### NETWORK plots
Ttrue <- matrix(0, I, I)
Ttrue[lower.tri(Ttrue)] <- c(Theta)
Ttrue <- Ttrue + t(Ttrue) 

pdf(paste0(path.data, dataset, "_adj_mat_K_", K, ".pdf"), height = 10, width = 10, pointsize = 12)
par(mfrow = c(1, 1), mar = c(4, 4, 3, 2) - 0.0, mgp = c(2, 1, 0), oma = c(0, 0, 0, 0))
adjacency.plot(mat = Yadj, labs = NA, show.grid = F, tick = 0, main = paste0("Adjacency Matrix", " RL500"))
dev.off()

pdf(paste0(path.data, dataset, "_true_inter_probs_K_", K, ".pdf"), height = 10, width = 10, pointsize = 12)
par(mfrow = c(1, 1), mar = c(4, 4, 3, 2) - 0.0, mgp = c(2, 1, 0), oma = c(0, 0, 0, 0))
heat.plot0(mat = Ttrue, labs = NA, show.grid = F, tick = 0, main = paste0("True interaction probabilities", " RL500"))
dev.off()

pdf(paste0(path.data, dataset, "_true_lat_positions_K_", K, ".pdf"), height = 10, width = 10, pointsize = 12)
par(mfrow = c(1, 1), mar = c(4, 4, 3, 2) - 0.0, mgp = c(2, 1, 0), oma = c(0, 0, 0, 0))
urange <- max(abs(range(Ustar))) * c(-1,1)
plot(NA, xlab = "", ylab = "", type = "n", xlim = urange, ylim = urange, main = paste0("True Latent Positions", " RL500"))
abline(h = 0, v = 0, col = "lightgray")
for (n in 1:N) text(x = Ustar[n, 1], y = Ustar[n, 2], labels = n, col = length(which(Xi == n)), cex = 0.5)
dev.off()




require(igraph)

g <- graph.adjacency(adjmatrix = Yadj, mode = "undirected", diag = F) 

transitivity(graph = g, type = "global")

assortativity.degree(graph = g, directed = F)

sum(Y)/length(Y)
