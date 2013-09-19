####Examples...
require(abind)
require(reshape2)
require(vegan)
require(ggplot2)

tmp <- betasim(n.site=2,p.mix=c(0,1),n.abund=2,spcat = 5 ,rand="poisson")
calcbeta(df_to_matrix(betasim(n.site=10,p.mix=.5,rand="poisson")))

#### Mix of sites and probabilities of mixing
pmix <- seq(0,1,length=10)
pabun <- c(.02,.05,.2,.5)
pop_sizes <- seq(1,4,length=50)
no_pred_df <- matrix(0,nrow=0,ncol=3)
for(i in pabun){
  for(j in pop_sizes){
    z <- mean(replicate(10,median(calcbeta(df_to_matrix(betasim(n.site=10, n.indiv.site=10^j,p.mix=i,rand="poisson"))))))
    no_pred_df <- rbind(no_pred_df,c(i,j,z))           
  }
}

no_pred_df <- as.data.frame(no_pred_df)
colnames(no_pred_df) <- c("pmix","population","beta")



ggplot(no_pred_df,aes(x=population,y=beta,group=pmix,colour=as.factor(pmix)))+geom_point()+geom_path()

#### Predator simulations, normal abundance reduction

pred_df <- matrix(0,nrow=0,ncol=3)
for(i in pabun){
  for(j in pop_sizes){
    z <- mean(replicate(50,median(calcbeta(df_to_matrix(beta_predation(betasim(n.site=10, n.indiv.site=10^j,p.mix=i,rand="poisson"),pred.gamma = 2, pred.spec = "none", pred.rand = "binom")  )))))
    pred_df <- rbind(pred_df,c(i,j,z))           
  }
}

pred_df <- as.data.frame(pred_df)
colnames(pred_df) <- c("pmix","population","beta")



ggplot(pred_df,aes(x=population,y=beta,group=pmix,colour=as.factor(pmix)))+geom_point()+geom_path()



#### Predator simulations, normal abundance reduction

pred_df_focal <- matrix(0,nrow=0,ncol=3)
for(i in pabun){
  for(j in pop_sizes){
    z <- mean(replicate(50,median(calcbeta(df_to_matrix(beta_predation(betasim(n.site=10, n.indiv.site=10^j,p.mix=i,rand="poisson"),pred.gamma = 2, pred.spec = "focal",pred.alpha = 1, pred.rand = "binom")  )))))
    pred_df_focal <- rbind(pred_df_focal,c(i,j,z))           
  }
}

pred_df_focal<- as.data.frame(pred_df_focal)
colnames(pred_df_focal) <- c("pmix","population","beta")



ggplot(pred_df_focal,aes(x=population,y=beta,group=pmix,colour=as.factor(pmix)))+geom_point()+geom_path()


###

pred_df_abund <- matrix(0,nrow=0,ncol=3)
for(i in pabun){
  for(j in pop_sizes){
    z <- mean(replicate(50,median(calcbeta(df_to_matrix(beta_predation(betasim(n.site=10, n.indiv.site=10^j,p.mix=i,rand="poisson"),pred.gamma = 2, pred.spec = "abund",pred.alpha = 1, pred.rand = "binom")  )))))
    pred_df_abund <- rbind(pred_df_abund,c(i,j,z))           
  }
}

pred_df_abund<- as.data.frame(pred_df_abund)
colnames(pred_df_abund) <- c("pmix","population","beta")



ggplot(pred_df_abund,aes(x=population,y=beta,group=pmix,colour=as.factor(pmix)))+geom_point()+geom_path()

#### Bind frames together....

types <- c(rep("No Predation",200),rep("Equal reduction",200), rep("Focal species",200),rep("Most abund pref.",200))

f.df <- cbind(rbind(no_pred_df,pred_df,pred_df_focal,pred_df_abund),types)


ggplot(f.df,aes(x=population,y=beta,group=pmix,colour=as.factor(pmix)))+geom_point()+geom_path()+facet_wrap(~types)+theme_bw()+scale_colour_manual("Mixing probability",values=c("Black","Red","Blue","purple")) + xlab(expression("Population size in Log"[10])) + ylab("Beta (mean Jaccard's dist)")
