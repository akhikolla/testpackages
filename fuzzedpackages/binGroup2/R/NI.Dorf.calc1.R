
# Start  NI.Dorf.calc1() function
###################################################################
# Brianna Hitt - 12-06-19
# This function is the same as NI.Dorf(), but no longer finds the 
#   optimal testing configuration. It only calculates operating
#   characteristics for a specific testing configuration.

NI.Dorf.calc1 <- function(p, Se, Sp, group.sz, a, 
                          trace=TRUE, print.time=TRUE, ...){
  
  start.time<-proc.time()
  
  I <- group.sz

  # generate a probability vector for homogeneous population
  p.vec <- rep(x=p[1], times=I)
  order.for.p <- 1:I
  
  # calculate descriptive measures for two-stage hierarchical testing
  save.info <- hierarchical.desc2(p=p.vec[order.for.p], 
                                  se=Se, sp=Sp, I2=NULL, order.p=FALSE)
  
  # extract ET, PSe, PSp and calculate the MAR function
  ET <- save.info$ET

  # for non-informative Dorfman (two-stage hierarchical) testing, 
  #   all individuals have the same testing accuracy measures
  check <- check.all.equal(save.info$individual.testerror, 
                           which(colnames(save.info$individual.testerror)=="psp.vec"))
  if (is.null(check)) {
    ind.testerror <- get.unique.index(save.info$individual.testerror[a,], 
                                      which(colnames(save.info$individual.testerror)=="psp.vec"), 
                                      rowlabel = a)[,-1]
  } else {
    ind.testerror <- check[,-1]
  }
  colnames(ind.testerror) <- c("PSe", "PSP", "PPPV", "PNPV", "individuals")
  
  group.testerror <- save.info$group.testerror
  names(group.testerror) <- NULL
  PSe <- group.testerror[1]
  PSp <- group.testerror[2]
  PPPV <- group.testerror[3]
  PNPV <- group.testerror[4]

  save.it <- c(p[1], I, ET, ET/I, PSe, PSp, PPPV, PNPV)

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
  
  list("algorithm"="Non-informative two-stage hierarchical testing",
       "prob"=p[1], "Se"=Se.display, "Sp"=Sp.display, 
       "Config"=list("Stage1"=save.it[2]), "p.vec"=p.vec, 
       "ET"=save.it[3], "value"=save.it[4], 
       "Accuracy"=list("Individual"=ind.testerror, "Overall"=acc.ET))
}

###################################################################
