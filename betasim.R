betasim <- function(n.abund=5,        ## number of abundance categories
                    p.abund=0.5,      ## relative ranking of subsequent abundance categories
                    diff.abund=NULL,  ## alternative parameterization: difference between most and least abundant category
                    spcat=2,          ## number of species per abundance category PER SITE (overdetermined)
                    totsp=30,         ## total species pool size (FIXME)
                    n.site=3,         ## number of sites
                    n.indiv.site=100,      ## individuals per site
                    n.indiv.tot=NULL,      ## total number of individuals
                    p.mix=rep(0,n.abund),  ## default: no mixing
                    ## Probability of across-site mixing for each abundance category:
                    ##   0=completely endemic, 1=completely scrambled
                    seed=NULL,
                    rand=c("none","multinom","poisson"),
                    ## randomization type: 'none' to keep numbers==expected numbers
                    ##                     'multinom' keeps total per site fixed
                    ##                     'poisson' allows Poisson number per site
                    rarefy=1

) {
  ## species category/site combinations
  ## was: LETTER.spcat where LETTER was associated with site,
  ## e.g.
  if ((n.mix <- length(p.mix))>1 && n.mix!=n.abund) {
    stop("mismatch between mixing proportions and number of abundance classes")
  }
  if (spcat<1 && !isTRUE(all.equal(round(n.site*spcat),n.site*spcat)))
    stop("if spcat<1, n.site*spcat must be an integer")
  if (!is.null(n.indiv.tot)) {
    if (!missing(n.indiv.site)) stop("must specify at most one of n.indiv.site and n.indiv.tot")
    n.indiv.site <- n.indiv.tot/n.site
    if (floor(n.indiv.site) != n.indiv.site) {
      warning("rounding number of individuals per site")
      n.indiv.site <- round(n.indiv.site)
    }
  }
  p.mix <- rep(p.mix,length.out=n.abund)  ## replicate p.mix (if necessary)
  rand <- match.arg(rand)  ## check/expand 'rand' argument
  if (!is.null(seed)) set.seed(seed)
  ## Generate proportions in each category (on each reef):
  if (missing(p.abund) && !is.null(diff.abund))
    p.abund <- exp(log(diff.abund)/n)  ## n^th root of diff.abund
  avec <- p.abund^(0:(n.abund-1))
  avec <- avec/sum(avec)
  ## FIXME: adjust warning
  ## if (totsp %% (n.site * n.abund) != 0) 
  ##        warning("totsp is not an even multiple of n.site*n.abund")
  ## FIXME: we could round() if necessary ...
  if (missing(spcat)) {
    spcat <- totsp/(n.site*n.abund)  ## number of endemic species per abundance class per site
  } else if (!missing(n.abund) && !missing(n.site) && !missing(totsp))
    stop("must specify at most three of spcat, n.abund, n.site, totsp")
  ## FIXME: allow filling-in accordingly
  
  ## Set up three-way array:
  ##  dim 1: species class and number within class
  ##         [lower-case letters + lower case roman numerals]: see spcat
  ##         if spcat >= 1, then lower-case letters correspond to the original endemic site
  ##         (a=1, b=2, etc.); otherwise, lower-case letters may be found at a range of sites
  ##  dim 2: abundance class [0 to $n-1$]
  ##  dim 3: site [numbers -- WAS upper case letters]
  ##
  totspecies <- round(n.site*spcat)
  ## n.site <- 20; spcat <- 0.25; totspecies <- 5  ## 5 LETTERS, 1 roman code
  ## n.site <- 5; spcat <- 1; totspecies <- 5   ## 5 LETTERS, 1 roman code
  ## n.site <- 5; spcat <- 2; totspecies <- 10  ## 5 LETTERS, 2 roman codes
  ## ugh.  there should be a better way to get this pattern ...
  betw <- function(x,a,b) {
    y <- x > a & x <= b
    storage.mode(y) <- "numeric"
    y
  }
  
  m <- matrix(rep(seq(n.site),totspecies),ncol=totspecies)
  m2 <- t(betw(ceiling(m-(col(m)-1)/spcat),0,ceiling(1/spcat)))
  sp1 <- rep(letters[seq(totspecies/ceiling(spcat))],each=ceiling(spcat))
  sp2 <- rep(tolower(as.roman(seq(ceiling(spcat)))),length.out=totspecies)
  dimnames(m2) <- list(species=paste(sp1,sp2,sep="."),site=1:n.site)
  a1 <- abind(replicate(n.abund, m2, simplify=FALSE), along=3)
  a2 <- aperm(a1,c(1,3,2))
  dimnames(a2)[[2]] <- 0:(n.abund-1)
  names(dimnames(a2)) <- c("species","abund","site")
  sp_array <- d <- a2
  ## Assign proportions to endemic sites:
  d <- sweep(d,2,avec,"*")
  ## Normalize:
  d <- sweep(d,3,apply(d,3,sum),"/")
  if (FALSE) {
    ## examine:
    as.data.frame.table(d["a.i",,"1",drop=FALSE])  ## endemic species
    as.data.frame.table(d["b.i",,"1",drop=FALSE])  ## non-endemic species
  }
  ## mix among sites
  for (j in 1:n.abund) {
    pvec0 <- rep(p.mix[j]/n.site,n.site) ## redistributed individuals
    for (i in 1:totspecies) {
      n.exp <- sum(d[i,j,])
      endem <- (d[i,j,]>0)
      pvec <- pvec0+endem*((1-p.mix[j])/sum(endem))
      d[i,j,] <- n.exp*pvec
    }
  }     
  if (rand=="multinom") {  ## fixed number of individuals per site
    for (i in 1:n.site) {
      d[,,i] <- rmultinom(1,size=n.indiv.site,prob=d[,,i])
    }
  } else if (rand=="poisson") {  ## allow variation
    d[] <- rpois(length(d),lambda=n.indiv.site*d)
  }  else if (rand=="none") {
    d[] <- n.indiv.site*d
  }
  return(d)
}

