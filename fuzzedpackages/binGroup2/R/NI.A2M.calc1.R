
# Start  NI.A2M.calc1() functions
###################################################################
# Brianna Hitt - 12-06-19
# This function is the same as NI.A2M(), but no longer finds the 
#   optimal testing configuration. It only calculates operating
#   characteristics for a specific testing configuration.
# Brianna Hitt - 03-18-2020
# Revised - Array.Measures() and MasterPool.Array.Measures() now allow 
#   Se/Sp to vary across stages of testing. Array.Measures() needs to 
#   be called with a subset of the Se/Sp vector (no master pool testing). 
# Changes to Array.Measures() can be viewed in arrayFunctions.R and 
#   changes to MasterPool.Array.Measures() can be viewed in arrayMasterFunctions.R.

NI.A2M.calc1 <- function(p, Se, Sp, group.sz, a, 
                         trace=TRUE, print.time=TRUE, ...){
  
  start.time<-proc.time()
  
  I <- group.sz
  N <- I^2
  
  # build a matrix of probabilities
  p.mat <- matrix(data=p[1], nrow=I, ncol=I)
  
  # calculate the measures for array testing, with and without master pooling
  # based on equations provided by Kim et al. (2007) where n represents the 
  #   row/column size and n^2 represents the array size
  # Se/Sp are (master pool, row/col test, individual test)
  #   do not need master pool Se/Sp for Array.Measures()
  save.info.A2 <- Array.Measures(p=p.mat, se=Se[-1], sp=Sp[-1])
  #   do need all Se/Sp for MasterPool.Array.Measures()
  save.info.A2M <- MasterPool.Array.Measures(results=save.info.A2, n=I, 
                                             pmat=p.mat, Se=Se, Sp=Sp)
  
  # extract accuracy measures for each individual, with and without 
  #   master pooling
  ET <- save.info.A2M$ET
  # use the transpose of each matrix to read individuals across each row, 
  #   so that individuals 1 through J come from the first row, individuals
  #   J+1 through 2*J come from the second row, and so on.
  all.ind.testerror <- data.frame("p"=as.numeric(t(p.mat)), 
                                  "pse.vec"=as.numeric(t(save.info.A2M$PSe)), 
                                  "psp.vec"=as.numeric(t(save.info.A2M$PSp)), 
                                  "pppv.vec"=as.numeric(t(save.info.A2M$PPV)), 
                                  "pnpv.vec"=as.numeric(t(save.info.A2M$NPV)))
  check <- check.all.equal(all.ind.testerror, 
                           which(colnames(all.ind.testerror)=="psp.vec"))
  if (is.null(check)) {
    ind.testerror <- get.unique.index(all.ind.testerror[a, ], 
                                      which(colnames(all.ind.testerror)=="psp.vec"), 
                                      rowlabel = a)[,-1]
  } else {
    ind.testerror <- check[, -1]
  }
  colnames(ind.testerror) <- c("PSe", "PSP", "PPPV", "PNPV", "individuals")
  PSe.mat <- save.info.A2M$PSe
  PSp.mat <- save.info.A2M$PSp

  # calculate overall accuracy measures
  PSe <- sum(p.mat*PSe.mat)/sum(p.mat)
  PSp <- sum((1-p.mat)*(PSp.mat))/sum(1-p.mat)
  PPPV <- sum(p.mat*PSe.mat)/sum(p.mat*PSe.mat + (1-p.mat)*(1-PSp.mat))
  PNPV <- sum((1-p.mat)*PSp.mat)/sum((1-p.mat)*PSp.mat + p.mat*(1-PSe.mat))
  
  save.it <- c(p[1], I, N, ET, ET/N, PSe, PSp, PPPV, PNPV)
  
  # put accuracy measures in a matrix for easier display of results
  acc.ET <- matrix(data=save.it[6:9], nrow=1, ncol=4, 
                   dimnames=list(NULL, c("PSe", "PSp", "PPPV", "PNPV")))

  # create input accuracy value matrices for output display
  Se.display <- matrix(data = Se, nrow = 1, ncol = 3, 
                       dimnames = list(NULL, "Test"=c("Master Pool", "Row/Column", "Individual")))
  Sp.display <- matrix(data = Sp, nrow = 1, ncol = 3, 
                       dimnames = list(NULL, "Test"=c("Master Pool", "Row/Column", "Individual")))
  # use below if Se/Sp for row and column testing are allowed to differ
  # Se.display <- matrix(data = Se, nrow = 1, ncol = 4, 
  #                      dimnames = list(NULL, "Test"=c("Master Pool", "Row", "Column", "Individual")))
  # Sp.display <- matrix(data = Sp, nrow = 1, ncol = 4, 
  #                      dimnames = list(NULL, "Test"=c("Master Pool", "Row", "Column", "Individual")))
  
  # print time elapsed, if print.time == TRUE
  if(print.time){
    time.it(start.time)
  }
  
  list("algorithm"="Non-informative array testing with master pooling",
       "prob"=p[1], "Se"=Se.display, "Sp"=Sp.display, 
       "Config"=list("Array.dim"=save.it[2], "Array.sz"=save.it[3]), 
       "p.mat"=p.mat, "ET"=save.it[4], "value"=save.it[5], 
       "Accuracy"=list("Individual"=ind.testerror, "Overall"=acc.ET))
}

###################################################################
