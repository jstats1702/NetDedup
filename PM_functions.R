PMAlg_R <- function (Xi.data, truth) 
{
     S <- dim(Xi.data)[1]
     AC.data <- matrix(NA, S, 6)
     for(s in 1:S) {
          inmat     <- vec_to_symmat_R(get_inmat( Xi.data[s, ] ))
          predicted <- inmat[lower.tri(inmat)]
          AC        <- accuracy_R(truth, predicted)
          AC.data[s, ] <- c( AC$MCC, AC$R, AC$P, AC$FNR, AC$FDR, AC$F1 )
          if (s%%floor(0.1*S) == 0) cat(round(100*s/S,1), "% done ... \n", sep = "")
     }
     AC.data
}