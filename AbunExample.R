###Abundance examples.
require(abind)
require(reshape2)
require(vegan)
require(ggplot2)

###generalist predator
## Change values of pred.gamma


#### Mix of sites and probabilities of mixing
pop_sizes <- seq(1.5,4,length=10)
generalist <- matrix(0,nrow=0,ncol=3)
gamma.vals <- round(seq(-2,2,length=10),2)
for(i in pop_sizes){
  for(j in gamma.vals){
    z <- mean(replicate(10,median(calcbeta(array_to_matrix(beta_predation(betasim(n.site=10, n.indiv.site=10^i,p.mix=1,rand="poisson",spcat=9,n.abund=1),pred.gamma = j, pred.spec = "none", pred.rand = "binom"))))))
    generalist <- rbind(generalist,c(i,j,z))           
  }
}

generalist <- as.data.frame(generalist)
colnames(generalist) <- c("population","predgamma","beta")
ggplot(generalist,aes(x=population,y=beta,group=predgamma,colour=as.factor(predgamma)))+geom_point()+geom_path()+theme_bw()

test <- generalist

### Now lets look at alpha on preference for common species


#### Mix of sites and probabilities of mixing
pop_sizes <- seq(1.5,4,length=30)
common <- matrix(0,nrow=0,ncol=3)
alpha.vals <- round(seq(0,2,length=10),2)
for(i in pop_sizes){
  for(j in alpha.vals){
    z <- mean(replicate(10,median(calcbeta(array_to_matrix(beta_predation(betasim(n.site=10, n.indiv.site=10^i,p.mix=.3,rand="poisson"),pred.gamma = 2,pred.alpha=j, pred.spec = "cabund", pred.rand = "binom"))))))
    common <- rbind(common,c(i,j,z))           
  }
}

common <- as.data.frame(common)
colnames(common) <-  c("population","alpha","beta")
ggplot(common,aes(x=population,y=beta,group=alpha,colour=as.factor(alpha)))+geom_point()+geom_path()+theme_bw()
ggplot(c2,aes(x=population,y=beta,group=alpha,colour=as.factor(alpha)))+geom_point()+geom_path()+theme_bw()


#### Mix of sites and probabilities of mixing
pop_sizes <- seq(1.5,4,length=30)
rare <- matrix(0,nrow=0,ncol=3)
alpha.vals <- round(seq(0,5,length=10),2)
for(i in pop_sizes){
  for(j in alpha.vals){
    z <- mean(replicate(10,median(calcbeta(array_to_matrix(beta_predation(betasim(n.site=10, n.indiv.site=10^i,p.mix=.3,rand="poisson"),pred.gamma = 2,pred.alpha=j, pred.spec = "rabund", pred.rand = "binom"))))))
    rare <- rbind(rare,c(i,j,z))           
  }
}

rare <- as.data.frame(rare)
colnames(rare) <- c("population","alpha","beta")
ggplot(rare,aes(x=population,y=beta,group=alpha,colour=as.factor(alpha)))+geom_point()+geom_path()+theme_bw()



### add a no predator scenario...

#### Mix of sites and probabilities of mixing
pop_sizes <- seq(1.5,4,length=30)
no_pred<- matrix(0,nrow=0,ncol=3)

for(i in pop_sizes){
 
    z <- mean(replicate(10,median(calcbeta(array_to_matrix(betasim(n.site=10, n.indiv.site=10^i,p.mix=.3,rand="poisson"))))))
    no_pred <- rbind(no_pred,c(i,NA,z))           
}
no_pred <- as.data.frame(no_pred)
colnames(no_pred) <- c("population","predgamma","beta")
