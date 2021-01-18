
# Start  Inf.Dorf.calc1() function
###############################################################################
# Brianna Hitt - 12-06-19
# This function is the same as Inf.Dorf.OTC1(), but no longer finds the 
#   optimal testing configuration. It only calculates operating
#   characteristics for a specific testing configuration.

Inf.Dorf.calc1 <- function(p, Se, Sp, group.sz, pool.szs, 
                           alpha=2, a, trace=TRUE, print.time=TRUE, ...){
  
  start.time<-proc.time()
  
  N <- group.sz
  
  # build a vector of probabilities for a heterogeneous population
  if(length(p)==1){
    p.vec <- expectOrderBeta(p=p, alpha=alpha, grp.sz=N, ...)
  } else if(length(p)>1){
    p.vec <- sort(p)
    alpha <- NA
  }
  
  # calculate descriptive measures for informative Dorfman testing, 
  #   given a configuration/set of pool sizes
  save.info <- inf.dorf.measures(prob=p.vec, se=Se, sp=Sp, N=N, 
                                 pool.sizes=pool.szs)
  
  # extract accuracy measures for each individual
  ET <- save.info$e
  all.ind.testerror <- save.info$summary[,-1]
  ind.testerror <- get.unique.index(all.ind.testerror[a, ], 
                                    which(colnames(all.ind.testerror)=="PSp"), 
                                    rowlabel = a)[,-1]
  colnames(ind.testerror) <- c("PSe", "PSP", "PPPV", "PNPV", "individuals")
  PSe.vec <- save.info$summary[,3]
  PSp.vec <- save.info$summary[,4]
  
  # calculate overall accuracy measures
  PSe <- sum(p.vec*PSe.vec)/sum(p.vec)
  PSp <- sum((1-p.vec)*(PSp.vec))/sum(1-p.vec)
  PPPV <- sum(p.vec*PSe.vec)/sum(p.vec*PSe.vec + (1-p.vec)*(1-PSp.vec))
  PNPV <- sum((1-p.vec)*PSp.vec)/sum((1-p.vec)*PSp.vec + p.vec*(1-PSe.vec))
  
  save.it <- c(alpha, N, ET, ET/N, PSe, PSp, PPPV, PNPV)
  
  # put accuracy measures in a matrix for easier display of results
  acc.ET <- matrix(data=save.it[5:8], nrow=1, ncol=4, 
                   dimnames=list(NULL, c("PSe", "PSp", "PPPV", "PNPV")))
  
  # create input accuracy value matrices for output display
  Se.display <- matrix(data = Se, nrow = 1, ncol = 2, 
                       dimnames = list(NULL, "Stage"=1:2))
  Sp.display <- matrix(data = Sp, nrow = 1, ncol = 2, 
                       dimnames = list(NULL, "Stage"=1:2))
  
  # print time elapsed, if print.time == TRUE
  if(print.time){
    time.it(start.time)
  }
  
  list("algorithm"="Informative two-stage hierarchical testing",
       "prob"=list(p), "alpha"=alpha, "Se"=Se.display, "Sp"=Sp.display,
       "Config"=list("Block.sz"=save.it[2], "pool.szs"=pool.szs), 
       "p.vec"=p.vec, "ET"=save.it[3], "value"=save.it[4], 
       "Accuracy"=list("Individual"=ind.testerror, "Overall"=acc.ET))
}

###############################################################################
