#### SCRATCH CODE WORKSPACE######

library(plyr)


tmp <- betasim(n.site=3, n.indiv.site=100,n.abun=1,spcat = 9,p.mix=1,rand="poisson")
plot(radfit(df_to_matrix(tmp)))

tmp.p <-beta_predation(tmp,pred.spec="none",pred.rand="binom",pred.gamma=.1)


beta_predation(betasim(n.site=10, n.indiv.site=10^i,p.mix=.3,rand="poisson"),pred.gamma = 2,pred.alpha=j, pred.spec = "cabund", pred.rand = "binom")


z <- df_to_matrix(beta_predation(betasim(n.site=10, n.indiv.site=10^i,p.mix=.3,rand="poisson"),pred.gamma = 2,pred.alpha=j, pred.spec = "cabund", pred.rand = "binom"))
calcbeta(z)


mat1 <- c(10,10,10,10,10,10,10,10,10)
dim(mat1) <- c(3,3)
mat2 <- c(1000,500,2,1)
dim(mat2) <- c(1,4)
beta_abun_score(mat1)
beta_abun_score(mat2)