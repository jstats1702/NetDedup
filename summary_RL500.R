rm(list = ls (all.names = TRUE))

####################
# general settings #
####################

library(xtable)

loc     <- "C:/RL1_paper/"

dataset <- "RL500"

source(paste0(loc, "code/my_Rfunctions.R"))

rn <- 2

out <- rbind( load_matching_results1(model = "PM1",  dataset, method = "LG", loc),
              load_matching_results1(model = "PNM1", dataset, method = "LG", loc) )

cstats  <- rbind( load_cstats_results1(model = "PM1",  dataset, loc),
                  load_cstats_results1(model = "PNM1", dataset, loc) )

tab.out1 <- out
tab.out1[,c(5,6,7,10)] <- round(matrix(as.numeric(as.matrix(out[,c(5,6,7,10)])), 16, 4), rn)
tab.out1 <- tab.out1[,-c(1,3,4,5,8,9)]

tab.stats1 <- cstats
tab.stats1[,c(4,5)] <- round(matrix(as.numeric(as.matrix(cstats[,c(4,5)])), 16, 2), rn)
tab.stats1 <- tab.stats1[,-c(1,2,3)]


# PM1
xtable(cbind(tab.out1[1:8,], tab.stats1[1:8,]), digits = c(rep(0,2), rep(rn, 5)), align = rep("c", 6 + 1) )


# PNM1
xtable(cbind(tab.out1[9:16,], tab.stats1[9:16,]), digits = c(rep(0,2), rep(rn, 5)), align = rep("c", 6 + 1) )

#############
# Network 2 #
#############

loc     <- "C:/RL1_paper_Net_2/"
out     <- load_matching_results1(model = "PNM1", dataset, method = "LG", loc)
cstats  <- load_cstats_results1(model = "PNM1", dataset, loc)

tab.out1 <- out
tab.out1[,c(5,6,7,10)] <- round(matrix(as.numeric(as.matrix(out[,c(5,6,7,10)])), 16, 4), rn)
tab.out1 <- tab.out1[,-c(1,3,4,5,8,9)]

tab.stats1 <- cstats
tab.stats1[,c(4,5)] <- round(matrix(as.numeric(as.matrix(cstats[,c(4,5)])), 16, 2), rn)
tab.stats1 <- tab.stats1[,-c(1,2,3)]

# PNM1
xtable(cbind(tab.out1, tab.stats1), digits = c(rep(0,2), rep(rn, 5)), align = rep("c", 6 + 1) )
