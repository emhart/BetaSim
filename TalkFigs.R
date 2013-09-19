#### Set-up heat plots

### Set up abundance classes
### Set up species categories
nabun <- c(60,30,10,5,3,2,1)
nsp <-  sort(nabun)



mix <- seq(0,1,length=11)
gamma <- seq(0,4, length=11)

# Define matrix of mixing probabilities

mixcols <- rbind(c(1,2),c(0,1))
mixcol.names <- c("Source: Rare","Source: Common")

### Set mean population size
### Population is log 10
pop <- 2.8

### Set site size
n.site = 10

### Set alpha

alpha <- c(-1,0,1)

output <- matrix(ncol=4,nrow=0)
text.vec <- vector()

for(i in 1:length(mix)){
  
  for(j in 1:length(gamma)){
    
     for(k in 1:length(alpha)){
       
       for(l in 1:2){
         mixmat <- c(0,0)
         mixmat[l] <- mix[i]  
    z <- mean(replicate(10,median(calcbeta(array_to_matrix(beta_predation(betasim(n.site=n.site, n.abund=2,spcat=10,n.indiv.site=10^pop,p.mix=mixmat,rand="poisson"),pred.gamma = gamma[j],pred.alpha=alpha[k], pred.spec = "abund", pred.rand = "binom"))))))
    
    output <- rbind(output,c(mix[i],gamma[j],alpha[k],z) )
    text.vec <- c(text.vec, mixcol.names[l])
       }
     }
  }
}

output <- as.data.frame(output)
colnames(output) <- c("Mixing","gamma","alpha","beta")
output$bsource <- text.vec
output$alpha <- as.character(output$alpha)
output$alpha[output$alpha=="-1"] <- "Rare Pred. Pref"
output$alpha[output$alpha=="0"] <- "No Pred. Pref"
output$alpha[output$alpha=="1"] <- "Common Pred. Pref"





ggplot(output,aes(y=Mixing*100, x=gamma,fill=beta)) + geom_raster()+scale_fill_continuous("Beta",low="yellow",high="red")+theme_bw()+facet_grid(bsource~alpha)+xlab("Strength of Predation")+ylab("Percent mixing of non-source class")+theme(strip.background = element_rect(fill="white"))
  


map <- map_data('world')
ggplot(map, aes(long, lat, group=group,fill=region)) + geom_path() 



test <- map
test <- subset(test,test$long > -170 & test$long < -30)
test <- subset(test,test$lat > 19 & test$lat < 83)

ggplot(data = test, aes(long, lat, group=group,fill=region)) + geom_path() + geom_polygon()




