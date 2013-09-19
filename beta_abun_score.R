#' Takes a matrix of species and abundance classes and return log-odds predation ratio
#' @param mat a matrix of sp x abun generation by betasim()
#' @details Log odds are set between pred.alpha and -pred.alpha.  Odds are based on relative abundance across all species within a sit.
#' 
#' @return log.odds a matrix of log odds ratios for predation scores.

beta_abun_score <- function(mat,alpha=1){
  comm <- as.vector(mat)
  abun.score <- seq(1,-1,length=length(unique(comm)))
  
  #abun.score <- seq(-1,1,length=sum(mat))
  Fn <- ecdf(comm)
  quants <- Fn(comm)
  log.odds <- quantile(abun.score,quants)
  log.odds <- log.odds*alpha
  if(alpha==0){
   log.odds <- rep(1,length(log.odds)) 
  }
  dim(log.odds) <- dim(mat)

  return(log.odds)  
}


