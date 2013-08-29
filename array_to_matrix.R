#' Converts a dataframe output from betasim() or beta_predation() to a matrix
#' based on earlier code by B. Bolker
#' 
#' @param d a multidimensional array from betasim() or beta_predation()
#' @return a species x site matrix
#' 

array_to_matrix <- function(d){
  d2 <- as.data.frame.table(d)

  ## rearrange data frame to species matrix
  d3 <- transform(d2,sp=paste(species,abund,sep="")) ## species names
  d4 <- dcast(d3,site~sp,value.var="Freq")[,-1] ## recast & drop site column
  dnames <- list(dimnames(d)[["site"]],colnames(d4))
  d5 <- as.matrix(d4)
  dimnames(d5) <- dnames  ## restore row/column names

  return(d5)
}