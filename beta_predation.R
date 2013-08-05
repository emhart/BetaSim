
beta_predation <-function(d,
                          pred.spec=c("none","focal","abund"),
                          pred.rand=c("binom","none"),
                          pred.gamma=-Inf,      ## predation intensity (-Inf==none)
                          pred.alpha=0,         ## predator preference parameter
                          ## if pred.spec=="abund" 0=neutral, -1=pref. rare, +1=pref. common
                          ## if pred.spec=="focal" 0=neutral, -1=pref non-endemic, +1=pref endemic
                          pred.gamma.sd=0       ## patch-to-patch variability in predation intensity (0=none)) 
                          )
                          {
  if (any(floor(d)!=d)) {
    warning("rounding to do predator selection")
    d <- round(d)
  }  
  ## rules for predator impact:
  ## if non-specific:
  ##   if pred.gamma.sd==0:
  ##      equivalent to rarefying with prob=plogis(pred.gamma)
  ##   if pred.gamma.sd>0:
  ##      rarefy with prob(site)=plogis(rnorm(nsite,pred.gamma,pred.gamma.sd))
  ## if focal:
  ##      as above, but add +pred.alpha to endemic species, -pred.alpha to other species
  ## if abund:
  ##      as above, but add +pred.alpha*(scaled p.abund)
  pred.logis <- d             ## copy shape of species table
  pred.logis[] <- pred.gamma  ## replace values with baseline predation value
  n.abund <- dim(d[,,1])[2]
  
  if (pred.gamma.sd>0) {
    sitepred <- rnorm(n.site,sd=pred.gamma.sd)
    ## dim 3 == site
    pred.logis <- pred.logis + sitepred[slice.index(pred.logis,3)]
  }
  
  ### Abundance based predation needs to be fixed to handle duplicate abundance values
  
  if (pred.spec=="focal") {
    ## dim 3 == site, dim 1 == origin
    origsite <- (slice.index(pred.logis,1)==slice.index(pred.logis,3))
    pred.logis[] <- pred.logis[] + ifelse(origsite,pred.alpha,-pred.alpha)
  } else if (pred.spec=="cabund") {
   ### Abundance score is calculated based on relative abundance of each taxa.
    
    abundscore <- seq(1,-1,length=n.abund)[slice.index(pred.logis,2)]
    pred.logis[] <- pred.logis[] + pred.alpha*abundscore
    for(i in 1:dim(d)[3]){
      index <- t(apply(d[,,i],1,order,decreasing=T))
      for(j in 1:dim(index)[1]){
        pred.logis[j,,i] <- pred.logis[j,index[j,],i]
      }
    }
    
    
  } else if (pred.spec=="rabund") {
    
    abundscore <- seq(-1,1,length=n.abund)[slice.index(pred.logis,2)]
    pred.logis[] <- pred.logis[] + pred.alpha*abundscore
    ### Now make sure that reductions follow the actual values
    for(i in 1:dim(d)[3]){
      index <- t(apply(d[,,i],1,order,decreasing=T))
      for(j in 1:dim(index)[1]){
        pred.logis[j,,i] <- pred.logis[j,index[j,],i]
      }
    }
    
  }
  ## survival probability is COMPLEMENTARY plogis() ...
  survprob <- plogis(pred.logis,lower.tail=FALSE)
  d[] <- if (pred.rand=="binom") {
    rbinom(length(d),size=d,prob=survprob)
  } else d*survprob
  return(d)
}