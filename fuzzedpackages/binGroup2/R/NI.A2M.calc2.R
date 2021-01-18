# Start  NI.A2M.calc2() function
###################################################################
# Brianna Hitt - 12-31-19
# This function is the same as TwoDisease.NI.A2M(), but no longer 
#   finds the optimal testing configuration. It only calculates 
#   operating characteristics for a specific testing configuration.

###############################################################################
# Uses Rcpp functions corresponding to the manuscript entitled "Array Testing 
#   for Multiplex Assays" (Hou et al. 2018)
# Note: some supporting functions are written in standalone C++ file for 
#   faster calculation. Latest update of R (3.4.2), Rtools (Rtools34), R 
#   packages "Rcpp" and "RcppArmadillo" are required.
###############################################################################

# Results for array testing with master pooling
NI.A2M.calc2 <- function(p.vec, Se, Sp, group.sz, 
                         trace=TRUE, print.time=TRUE, ...){
  
  start.time<-proc.time()
  
  I <- group.sz
  N <- I^2
  
  # calculate results for a single array size
  results <- ARRAY_master(p=p.vec, SE=Se, SP=Sp, n=I)
  
  # save the results for each array size in a row of the results matrix
  results.ATM <- c(I, N, results$ET_ATM*N, results$ET_ATM, 
                   results$PSE1_ATM, results$PSE2_ATM,
                   results$PSP1_ATM, results$PSP2_ATM, results$PPV1_ATM,
                   results$PPV2_ATM, results$NPV1_ATM, results$NPV2_ATM)
  
  acc.ATM <- matrix(data=results.ATM[5:12], nrow=2, ncol=4, byrow=FALSE, 
                    dimnames=list("Disease" = 1:2, 
                                  c("PSe", "PSp", "PPPV", "PNPV")))
  
  # Brianna 03.04.2020 
  # Added individual accuracy measures for consistency with other algorithms
  # Will need to edit this when the ARRAY function calculates 
  #   individual accuracy measures
  ind.acc1 <- data.frame(PSe = acc.ATM[1,1], PSp = acc.ATM[1,2], 
                         PPPV = acc.ATM[1,3], PNPV = acc.ATM[1,4], 
                         individuals = "All")
  ind.acc2 <- data.frame(PSe = acc.ATM[2,1], PSp = acc.ATM[2,2], 
                         PPPV = acc.ATM[2,3], PNPV = acc.ATM[2,4], 
                         individuals = "All") 
  
  # create input accuracy value matrices for output display
  Se.display <- matrix(data = Se, nrow = 2, ncol = 3, 
                       dimnames = list("Disease" = 1:2, 
                                       "Test" = c("Master Pool", "Row/Column", 
                                                  "Individual")))
  Sp.display <- matrix(data = Sp, nrow = 2, ncol = 3, 
                       dimnames = list("Disease" = 1:2, 
                                       "Test" = c("Master Pool", "Row/Column", 
                                                  "Individual")))
  # use below if Se/Sp for row and column testing is allowed to differ
  # Se.display <- matrix(data = Se, nrow = 2, ncol = 4, 
  #                      dimnames = list("Disease" = c("1", "2"), 
  #                                      c("Master Pool Testing", "Row Testing", 
  #                                        "Column Testing", "Individual Testing")))
  # Sp.display <- matrix(data = Sp, nrow = 2, ncol = 4, 
  #                      dimnames = list("Disease" = c("1", "2"), 
  #                                      c("Master Pool Testing", "Row Testing", 
  #                                        "Column Testing", "Individual Testing")))
  
  joint.p <- matrix(data=rep(p.vec, each=results.ATM[2]), nrow=4, ncol=results.ATM[2], 
                    byrow=TRUE, dimnames=list(c("00", "10", "01", "11"), 
                                              as.character(1:results.ATM[2])))
  
  # print time elapsed, if print.time == TRUE
  if(print.time){
    time.it(start.time)
  }
  
  # reorganize results
  list("algorithm"="Non-informative array testing with master pooling",
       "prob.vec"=p.vec, "Se"=Se.display, "Sp"=Sp.display,
       "Config"=list("Array.dim"=results.ATM[1], "Array.sz"=results.ATM[2]),
       "p.mat" = joint.p, "ET"=results.ATM[3], "value"=results.ATM[4],
       "Accuracy"=list("Disease 1 Individual" = ind.acc1, 
                       "Disease 2 Individual" = ind.acc2, 
                       "Overall" = acc.ATM))
}

#
