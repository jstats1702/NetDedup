cstats_R <- function (model, prior, case, dataset, loc)
{
     rn <- 2
     path.outs <- paste0(loc, "outputs/", model, "_", prior, "_", case, "_", dataset, "/")
     M   <- as.matrix(read.table(paste0(path.outs, "stats_out.txt")))
     out <- c( mean( M[,1] ), sd( M[,1] ) )
     out <- c( model, prior, case, round( out, rn ) )
     out <- as.data.frame(t(as.matrix( out )))
     colnames(out) <- c( "model" , "prior" , "case" , "N.pm", "N.sd")
     save(out, file = paste0(path.outs, model, "_", prior, "_", case, "_", dataset, "_cstats_results.RData"))
}

matching0_R <- function (Xi, model, prior, case, dataset, loc)
{
     rn     <- 7
     method <- "LG"
     cat("Working on model ", model, ", prior ", prior, ", case ", case, ", data ", dataset, ", METHOD ", method, ", ...\n", sep = "" )
     ##### paths
     path.code <- paste0(loc, "code/")
     path.outs <- paste0(loc, "outputs/", model, "_", prior, "_", case, "_", dataset, "/")
     ##### libs and source
     library(Rcpp)
     sourceCpp(paste0(path.code, "my_Cfunctions.cpp"))
     source   (paste0(path.code, "my_Rfunctions.R"))
     source   (paste0(path.code, "LG_functions.R"))
     source   (paste0(path.code, "LGalg_parallel.R"))
     ##### data
     inmat.hat <- c( as.matrix(read.table(paste0(path.outs, "inmat_out.txt"))) )
     truth     <- c( get_inmat(Xi) )  ## true incidence matrix
     ###### Xi inference, clustering (lower triangular)
     LGmin     <- 0.5
     LGmax     <- 0.99
     nLGgrid   <- 5
     predicted <- LGalg_parallel_R(Kgrid = seq(LGmin, LGmax, length = nLGgrid), rho = vec_to_symmat_R(inmat.hat), path.code)
     save(predicted, file = paste0(path.outs, model, "_", prior, "_", case, "_", dataset, "_matching_data_", method, ".RData"))
     AC  <- accuracy_R(truth, predicted)
     out <- c( model, prior, case, method, round( c( AC$MCC, AC$R, AC$P, AC$FNR, AC$FDR, AC$F1 ), rn ) )
     ##### return
     out <- as.data.frame(t(as.matrix( out )))
     colnames(out) <- c( "model" , "prior" , "case" , "method", "MCC" , "Recall" , "Precision" , "FNR" , "FDR" , "F1-Score" )
     save(out, file = paste0(path.outs, model, "_", prior, "_", case, "_", dataset, "_matching_results_", method, ".RData"))
}