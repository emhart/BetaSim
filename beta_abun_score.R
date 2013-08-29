#' Takes a matrix of species and abundance classes and return log-odds predation ratio
#' @param mat a matrix of sp x abun generation by betasim()
#' @param reverse True or False.  If True the predator will prefer rare species, otherwise it will prefer common species
#' @details Log odds are set between pred.alpha and -pred.alpha.  Odds are based on relative abundance across all species within a sit.
#' 
#' @return log.odds a matrix of log odds ratios for predation scores.

beta_abun_score <- function(mat,reverse = FALSE){
  comm <- as.vector(mat)
  abun.score <- seq(-1,1,length=length(unique(comm)))
  #abun.score <- seq(-1,1,length=sum(mat))
  Fn <- ecdf(comm)
  quants <- Fn(comm)
  if(reverse){
    quants <- 1-Fn(comm)
  }
  log.odds <- quantile(abun.score,quants)
  dim(log.odds) <- dim(mat)
  
  return(log.odds)  
}


