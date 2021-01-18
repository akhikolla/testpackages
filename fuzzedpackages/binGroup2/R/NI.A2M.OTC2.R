# Start  NI.A2M.OTC2() function
###################################################################
# Brianna Hitt

# Brianna Hitt - 04.02.2020
# Changed cat() to message()

###############################################################################
# Uses Rcpp functions corresponding to the manuscript entitled "Array Testing 
#   for Multiplex Assays" (Hou et al. 2018)
# Note: some supporting functions are written in standalone C++ file for 
#   faster calculation. Latest update of R (3.4.2), Rtools (Rtools34), R 
#   packages "Rcpp" and "RcppArmadillo" are required.
###############################################################################
# Finds the OTC for array testing with master pooling using an assay that 
#   tests for two diseases
# Allows the user to specify a range of array sizes instead of a only a 
#   maximum array size (like in the original ARRAY function from Hou et al.)
# Reorganizes results from the format of the ARRAY function
###############################################################################

# Results for array testing with master pooling

# Brianna Hitt - 04.02.2020
# Changed cat() to message

NI.A2M.OTC2 <- function(p.vec, Se, Sp, group.sz, updateProgress=NULL, 
                        trace=TRUE, print.time=TRUE, ...){
  
  start.time <- proc.time()
  
  set.of.I <- group.sz
  
  results.ATM <- matrix(data=NA, nrow=length(set.of.I), ncol=12)
  count <- 1
  
  for(I in set.of.I){
    N <- I^2
    
    # calculate results for a single array size
    results <- ARRAY_master(p=p.vec, SE=Se, SP=Sp, n=I)
    
    # save the results for each array size in a row of the results matrix
    results.ATM[count,] <- c(I, N, results$ET_ATM*N,
                             results$ET_ATM, results$PSE1_ATM, results$PSE2_ATM,
                             results$PSP1_ATM, results$PSP2_ATM, results$PPV1_ATM,
                             results$PPV2_ATM, results$NPV1_ATM, results$NPV2_ATM)

    if(is.function(updateProgress)){
      updateText <- paste0("Row/Column Size = ", I, ", Array Size = ", N)
      updateProgress(value = count/(length(set.of.I)+1), detail=updateText)
    }
    
    # print progress, if trace == TRUE
    if(trace){
      cat("Row/Column Size = ", I, ", Array Size = ", N, "\n", sep="")
    }
    count <- count + 1
  }
  
  # save the results for each initial array size
  if(length(set.of.I)==1){
    configs.ATM <- NA
  } else{
    configs.ATM <- results.ATM
    colnames(configs.ATM) <- c("I", "N", "ET", "value", "PSe1", "PSe2", "PSp1", 
                               "PSp2", "PPPV1", "PPPV2", "PNPV1", "PNPV2")
    configs.ATM <- convert.config(algorithm="A2M", results=configs.ATM, diseases=2)
  }
  
  opt.ATM <- results.ATM[results.ATM[,4]==min(results.ATM[,4])]
  
  acc.ATM <- matrix(data=opt.ATM[5:12], nrow=2, ncol=4, byrow=FALSE, 
                    dimnames=list("Disease"=c("1", "2"), c("PSe", "PSp", 
                                                           "PPPV", "PNPV")))
  
  # create input accuracy value matrices for output display
  Se.display <- matrix(data = Se, nrow = 2, ncol = 3, 
                       dimnames = list("Disease" = 1:2, 
                                       "Test" = c("Master Pool", "Row/Column", "Individual")))
  Sp.display <- matrix(data = Sp, nrow = 2, ncol = 3, 
                       dimnames = list("Disease" = 1:2, 
                                       "Test" = c("Master Pool", "Row/Column", "Individual")))
  # use below if Se/Sp for row and column testing is allowed to differ
  # Se.display <- matrix(data = Se, nrow = 2, ncol = 4, 
  #                      dimnames = list("Disease" = c("1", "2"), 
  #                                      c("Master Pool Testing", "Row Testing", 
  #                                        "Column Testing", "Individual Testing")))
  # Sp.display <- matrix(data = Sp, nrow = 2, ncol = 4, 
  #                      dimnames = list("Disease" = c("1", "2"), 
  #                                      c("Master Pool Testing", "Row Testing", 
  #                                        "Column Testing", "Individual Testing")))
  
  joint.p <- matrix(data=rep(p.vec, each=opt.ATM[2]), nrow=4, ncol=opt.ATM[2], 
                    byrow=TRUE, dimnames=list(c("00", "10", "01", "11"), 
                                              as.character(1:opt.ATM[2])))
  
  # print time elapsed, if print.time == TRUE
  if(print.time){
    time.it(start.time)
  }
  
  # reorganize results
  list("algorithm"="Non-informative array testing with master pooling",
       "prob.vec"=p.vec, "Se"=Se.display, "Sp"=Sp.display,
       "opt.ET"=list("OTC"=list("Array.dim"=opt.ATM[1], "Array.sz"=opt.ATM[2]),
                     "p.mat" = joint.p, "ET"=opt.ATM[3], "value"=opt.ATM[4],
                     "Accuracy"=acc.ATM), 
       "Configs"=configs.ATM)
  
}

#
