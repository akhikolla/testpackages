# Start  NI.Array.calc2() function
###################################################################
# Brianna Hitt - 12-31-19
# This function is the same as TwoDisease.NI.Array(), but no longer 
#   finds the optimal testing configuration. It only calculates 
#   operating characteristics for a specific testing configuration.

###############################################################################
# Uses Rcpp functions corresponding to the manuscript entitled "Array Testing 
#   for Multiplex Assays" (Hou et al. 2018)
# Note: some supporting functions are written in standalone C++ file for 
#   faster calculation. Latest update of R (3.4.2), Rtools (Rtools34), R 
#   packages "Rcpp" and "RcppArmadillo" are required.
###############################################################################

# Results for array testing without master pooling
NI.Array.calc2 <- function(p.vec, Se, Sp, group.sz, 
                           trace=TRUE, print.time=TRUE, ...){
  
  start.time<-proc.time()
  
  I <- group.sz
  N <- I^2
  
  # calculate results for a single array size
  # results <- ARRAY(p=p.vec, SE=Se, SP=Sp, n=i)
  # USE BELOW UNTIL ARRAY ALLOWS FOR SE, SP TO VARY ACROSS STAGES
  # results <- ARRAY(p=p.vec, SE=c(Se[1,1], Se[2,1]), SP=c(Sp[1,1], Sp[2,1]), n=I)
  results <- ARRAY_nomaster(p=p.vec, SE=Se, SP=Sp, n=I)
  
  # save the results for each array size in a row of the results matrix
  results.AT <- c(I, N, results$ET_AT*N, results$ET_AT, 
                  results$PSE1_AT, results$PSE2_AT,
                  results$PSP1_AT, results$PSP2_AT, results$PPV1_AT,
                  results$PPV2_AT, results$NPV1_AT, results$NPV2_AT)
  
  acc.AT <- matrix(data=results.AT[5:12], nrow=2, ncol=4, byrow=FALSE, 
                   dimnames=list("Disease" = 1:2, 
                                 c("PSe", "PSp", "PPPV", "PNPV")))
  
  # Brianna 03.04.2020 
  # Added individual accuracy measures for consistency with other algorithms
  # Will need to edit this when the ARRAY function calculates 
  #   individual accuracy measures
  ind.acc1 <- data.frame(PSe = acc.AT[1,1], PSp = acc.AT[1,2], 
                         PPPV = acc.AT[1,3], PNPV = acc.AT[1,4], 
                         individuals = "All")
  ind.acc2 <- data.frame(PSe = acc.AT[2,1], PSp = acc.AT[2,2], 
                         PPPV = acc.AT[2,3], PNPV = acc.AT[2,4], 
                         individuals = "All") 
  
  # create input accuracy value matrices for output display
  Se.display <- matrix(data = Se, nrow = 2, ncol = 2, 
                       dimnames = list("Disease" = 1:2, 
                                       "Test" = c("Row/Column", "Individual")))
  Sp.display <- matrix(data = Sp, nrow = 2, ncol = 2, 
                       dimnames = list("Disease" = 1:2, 
                                       "Test" = c("Row/Column", "Individual")))
  # use below if Se/Sp for row and column testing is allowed to differ
  # Se.display <- matrix(data = Se, nrow = 2, ncol = 3, 
  #                      dimnames = list("Disease" = c("1", "2"), 
  #                                      c("Row Testing", "Column Testing", 
  #                                        "Individual Testing")))
  # Sp.display <- matrix(data = Sp, nrow = 2, ncol = 3, 
  #                      dimnames = list("Disease" = c("1", "2"), 
  #                                      c("Row Testing", "Column Testing", 
  #                                        "Individual Testing")))
  
  joint.p <- matrix(data=rep(p.vec, each=results.AT[2]), nrow=4, ncol=results.AT[2], 
                    byrow=TRUE, dimnames=list(c("00", "10", "01", "11"), 
                                              as.character(1:results.AT[2])))
  
  # print time elapsed, if print.time == TRUE
  if(print.time){
    time.it(start.time)
  }
  
  # reorganize results
  list("algorithm"="Non-informative array testing without master pooling",
       "prob.vec"=p.vec, "Se"=Se.display, "Sp"=Sp.display,
       "Config"=list("Array.dim"=results.AT[1], "Array.sz"=results.AT[2]),
       "p.mat" = joint.p, "ET"=results.AT[3], "value"=results.AT[4],
       "Accuracy"=list("Disease 1 Individual" = ind.acc1, 
                       "Disease 2 Individual" = ind.acc2, 
                       "Overall" = acc.AT))
}

#