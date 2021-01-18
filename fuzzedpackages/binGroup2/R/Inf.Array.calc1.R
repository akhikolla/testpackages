
# Start  Inf.Array.calc1() function
###############################################################################
# Brianna Hitt - 12-06-19
# This function is the same as Inf.Array.OTC1(), but no longer finds the 
#   optimal testing configuration. It only calculates operating
#   characteristics for a specific testing configuration.

Inf.Array.calc1 <- function(p, Se, Sp, group.sz, alpha=2, a, 
                            trace=TRUE, print.time=TRUE, ...){
  
  start.time<-proc.time()
  
  I <- group.sz
  N <- I^2
    
  # build a vector of probabilities for a heterogeneous population
  if(length(p)==1){
    p.vec <- expectOrderBeta(p=p, alpha=alpha, grp.sz=N, ...)
  } else if(length(p)>1){
    p.vec <- sort(p)
    alpha <- NA
  }
  
  # build a matrix of probabilities using the gradient design
  p.ga <- informativeArrayProb(prob.vec=p.vec, nr=I, nc=I, method="gd")
  
  # call Array.Measures() to calculate descriptive measures for the given 
  #   array size
  save.info <- Array.Measures(p=p.ga, se=Se, sp=Sp)
  
  # extract accuracy measures for each individual
  ET <- save.info$ET
  # use the transpose of each matrix to read individuals across each row, 
  #   so that individuals 1 through J come from the first row, individuals
  #   J+1 through 2*J come from the second row, and so on.
  all.ind.testerror <- data.frame("p"=as.numeric(t(p.ga)), 
                                  "pse.vec"=as.numeric(t(save.info$PSe)), 
                                  "psp.vec"=as.numeric(t(save.info$PSp)), 
                                  "pppv.vec"=as.numeric(t(save.info$PPV)), 
                                  "pnpv.vec"=as.numeric(t(save.info$NPV)))
  ind.testerror <- get.unique.index(all.ind.testerror[a,], 
                                    which(colnames(all.ind.testerror)=="psp.vec"), 
                                    rowlabel = a)[,-1]
  colnames(ind.testerror) <- c("PSe", "PSP", "PPPV", "PNPV", "individuals")
  
  PSe.mat <- save.info$PSe
  PSp.mat <- save.info$PSp
  
  # calculate overall accuracy measures
  PSe <- sum(p.ga*PSe.mat)/sum(p.ga)
  PSp <- sum((1-p.ga)*(PSp.mat))/sum(1-p.ga)
  PPPV <- sum(p.ga*PSe.mat)/sum(p.ga*PSe.mat + (1-p.ga)*(1-PSp.mat))
  PNPV <- sum((1-p.ga)*PSp.mat)/sum((1-p.ga)*PSp.mat + p.ga*(1-PSe.mat))
  
  save.it <- c(alpha, I, N, ET, ET/N, PSe, PSp, PPPV, PNPV)
  
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
  
  list("algorithm"="Informative array testing without master pooling",
       "prob"=list(p), "alpha"=alpha, "Se"=Se.display, "Sp"=Sp.display, 
       "Config"=list("Array.dim"=save.it[2], "Array.sz"=save.it[3]), 
       "p.mat"=p.ga, "ET"=save.it[4], "value"=save.it[5], 
       "Accuracy"=list("Individual"=ind.testerror, "Overall"=acc.ET))
}

###################################################################
