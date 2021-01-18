
# Start  Inf.D3.calc1() function
###############################################################################
# Brianna Hitt - 12-06-19
# This function is the same as Inf.D3.OTC1(), but no longer finds the 
#   optimal testing configuration. It only calculates operating
#   characteristics for a specific testing configuration.

Inf.D3.calc1 <- function(p, Se, Sp, group.sz, pool.szs, 
                         alpha=2, a, trace=TRUE, print.time=TRUE, ...){
  
  start.time<-proc.time()
  
  I <- group.sz
  
  # build a vector of probabilities for a heterogeneous population
  if(length(p)==1){
    p.vec <- expectOrderBeta(p=p, alpha=alpha, grp.sz=I, ...)
  } else if(length(p)>1){
    p.vec <- sort(p)
    alpha <- NA
  }
  order.for.p <- 1:I
  
  # call hierarchical.desc2() for the configuration
  save.info <- hierarchical.desc2(p=p.vec[order.for.p], se=Se, sp=Sp, 
                                  I2=pool.szs, order.p=FALSE)
  
  # extract accuracy measures for each individual
  ET <- save.info$ET
  ind.testerror <- get.unique.index(save.info$individual.testerror[a,], 
                                    which(colnames(save.info$individual.testerror)=="psp.vec"), 
                                    rowlabel = a)[,-1]
  colnames(ind.testerror) <- c("PSe", "PSP", "PPPV", "PNPV", "individuals")
  PSe.vec <- save.info$individual.testerror$pse.vec
  PSp.vec <- save.info$individual.testerror$psp.vec
  
  # calculate overall accuracy measures
  group.testerror <- save.info$group.testerror
  names(group.testerror) <- NULL
  PSe <- group.testerror[1]
  PSp <- group.testerror[2]
  PPPV <- group.testerror[3]
  PNPV <- group.testerror[4]
  
  save.it <- c(alpha, I, ET, ET/I, PSe, PSp, PPPV, PNPV)
  
  # put accuracy measures in a matrix for easier display of results
  acc.ET <- matrix(data=save.it[5:8], nrow=1, ncol=4, 
                   dimnames=list(NULL, c("PSe", "PSp", "PPPV", "PNPV")))
  
  # create input accuracy value matrices for output display
  Se.display <- matrix(data = Se, nrow = 1, ncol = 3, 
                       dimnames = list(NULL, "Stage"=1:3))
  Sp.display <- matrix(data = Sp, nrow = 1, ncol = 3, 
                       dimnames = list(NULL, "Stage"=1:3))
  
  # print time elapsed, if print.time == TRUE
  if(print.time){
    time.it(start.time)
  }
  
  list("algorithm"="Informative three-stage hierarchical testing",
       "prob"=list(p), "alpha"=alpha, "Se"=Se.display, "Sp"=Sp.display, 
       "Config"=list("Stage1"=save.it[2], "Stage2"=pool.szs), 
       "p.vec"=p.vec, "ET"=save.it[3], "value"=save.it[4], 
       "Accuracy"=list("Individual"=ind.testerror, "Overall"=acc.ET))
}

###############################################################################

