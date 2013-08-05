
calcbeta <- function(m,method="jaccard",trap.errors=TRUE,
                     distances=c("centroid","pairwise")) {
  ## calculate beta distances, either
  distances <- match.arg(distances)
  nsites <- nrow(m)
  ## HACK: use binary with jaccard, otherwise not
  vv <- vegdist(m,method=method,binary=(method=="jaccard"))
  nvals <- if (distances=="centroid") nsites else (nsites*(nsites-1)/2)
  if (all(vv==0)) return(rep(0,nvals)) ## all identical sites
  ##  if (any(vv==0)) NA else
  if (distances=="centroid") {
    if (!trap.errors) {
      retval <- betadisper(vv,rep("blank",nsites))$distances
    } else {
      tt <- try(betadisper(vv,rep("blank",nsites)))
      if (inherits(tt,"try-error")) {
        retval <- rep(NA,nsites)
      } else retval <- tt$distances
    } 
  } else {
    retval <- c(unclass(vv)) 
  }
  retval
}
