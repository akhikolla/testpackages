
# Start  NI.Array.calc1() function
###################################################################
# Brianna Hitt - 12-06-19
# This function is the same as NI.Array(), but no longer finds the 
#   optimal testing configuration. It only calculates operating
#   characteristics for a specific testing configuration.

NI.Array.calc1 <- function(p, Se, Sp, group.sz, a, 
                           trace=TRUE, print.time=TRUE, ...){
  
  start.time<-proc.time()
  
  I <- group.sz
  N <- I^2
  
  # build a matrix of probabilities
  # this is the same for an overall probability p and for a vector p
  p.mat <- matrix(data=p[1], nrow=I, ncol=I)
  
  # call Array.Measures to calculate descriptive measures for the given array size
  save.info <- Array.Measures(p=p.mat, se=Se, sp=Sp)
  
  # extract accuracy measures for each individual
  ET <- save.info$ET
  # use the transpose of each matrix to read individuals across each row, 
  #   so that individuals 1 through J come from the first row, individuals
  #   J+1 through 2*J come from the second row, and so on.
  all.ind.testerror <- data.frame("p"=as.numeric(t(p.mat)), 
                                  "pse.vec"=as.numeric(t(save.info$PSe)), 
                                  "psp.vec"=as.numeric(t(save.info$PSp)), 
                                  "pppv.vec"=as.numeric(t(save.info$PPV)), 
                                  "pnpv.vec"=as.numeric(t(save.info$NPV)))
  check <- check.all.equal(all.ind.testerror, 
                           which(colnames(all.ind.testerror)=="psp.vec"))
  if (is.null(check)) {
    ind.testerror <- get.unique.index(all.ind.testerror[a,], 
                                      which(colnames(all.ind.testerror)=="psp.vec"), 
                                      rowlabel = a)[,-1]
  } else {
    ind.testerror <- check[, -1]
  }
  colnames(ind.testerror) <- c("PSe", "PSP", "PPPV", "PNPV", "individuals")
  
  PSe.mat <- save.info$PSe
  PSp.mat <- save.info$PSp

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
  Se.display <- matrix(data = Se, nrow = 1, ncol = 2, 
                       dimnames = list(NULL, "Test"=c("Row/Column", "Individual")))
  Sp.display <- matrix(data = Sp, nrow = 1, ncol = 2, 
                       dimnames = list(NULL, "Test"=c("Row/Column", "Individual")))
  # use below if Se/Sp for row and column testing is allowed to differ
  # Se.display <- matrix(data = Se, nrow = 1, ncol = 3, 
  #                      dimnames = list(NULL, "Test"=c("Row", "Column", "Individual")))
  # Sp.display <- matrix(data = Sp, nrow = 1, ncol = 3, 
  #                      dimnames = list(NULL, "Test"=c("Row", "Column", "Individual")))
  
  # print time elapsed, if print.time == TRUE
  if(print.time){
    time.it(start.time)
  }
  
  list("algorithm"="Non-informative array testing without master pooling",
       "prob"=p[1], "Se"=Se.display, "Sp"=Sp.display, 
       "Config"=list("Array.dim"=save.it[2], "Array.sz"=save.it[3]), 
       "p.mat"=p.mat, "ET"=save.it[4], "value"=save.it[5], 
       "Accuracy"=list("Individual"=ind.testerror, "Overall"=acc.ET))
}

###################################################################
