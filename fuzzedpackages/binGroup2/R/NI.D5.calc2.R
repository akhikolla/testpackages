# Start  NI.D3.calc2() function
###################################################################
# Brianna Hitt - 02-18-2020
# This function calculates operating characteristics for a 
#   specific testing configuration using non-informative 
#   five-stage hierarchical testing.

NI.D5.calc2 <- function(p.vec, Se, Sp, ordering, group.mem, a, 
                        trace=TRUE, print.time=TRUE, ...){
  
  start.time <- proc.time()
  
  I <- ncol(group.mem)
  
  # generate a joint probability vector for each individual
  joint.p <- matrix(data=rep(p.vec, times=I), nrow=4, ncol=I, byrow=FALSE, 
                    dimnames=list(joint.p.row.labels(ordering), 
                                  as.character(1:I)))
  
  save.info <- TwoDisease.Hierarchical(joint.p=joint.p, group.mem=group.mem, 
                                       Se=Se, Sp=Sp, ordering=ordering, 
                                       a=1:I, accuracy=TRUE)
  
  ind.acc1.all <- save.info$Disease1$Individual
  rownames(ind.acc1.all) <- ind.acc1.all[, "Individual"]
  check1 <- check.all.equal(ind.acc1.all, 
                            which(colnames(ind.acc1.all)=="PSp")) 
  if (is.null(check1)) {
    ind.acc1 <- get.unique.index(ind.acc1.all[a,], 
                                 which(colnames(ind.acc1.all)=="PSp"), 
                                 rowlabel = a)[,-1]
  } else {
    ind.acc1 <- check1[,-1]
  }
  
  ind.acc2.all <- save.info$Disease2$Individual
  rownames(ind.acc2.all) <- ind.acc2.all[, "Individual"]
  check2 <- check.all.equal(ind.acc2.all, 
                            which(colnames(ind.acc2.all)=="PSp")) 
  if (is.null(check2)) {
    ind.acc2 <- get.unique.index(ind.acc2.all[a,], 
                                 which(colnames(ind.acc2.all)=="PSp"), 
                                 rowlabel = a)[,-1]
  } else {
    ind.acc2 <- check2[,-1]
  }
  
  acc.dis1 <- save.info$Disease1$Overall
  acc.dis2 <- save.info$Disease2$Overall
  
  save.it <- c(I, save.info$Expt, save.info$Expt/I, acc.dis1, acc.dis2)
  
  accuracy.mat <- matrix(data=save.it[4:11], nrow=2, ncol=4, byrow=TRUE, 
                         dimnames=list("Disease" = 1:2, 
                                       c("PSe", "PSp", "PPPV", "PNPV")))
  
  Se.display <- matrix(data = Se, nrow = 2, ncol = 5, 
                       dimnames = list("Disease" = 1:2, "Stage" = 1:5))
  Sp.display <- matrix(data = Sp, nrow = 2, ncol = 5, 
                       dimnames = list("Disease" = 1:2, "Stage" = 1:5))

  # print time elapsed, if print.time == TRUE
  if(print.time){
    time.it(start.time)
  }
  
  list("algorithm"="Non-informative five-stage hierarchical testing",
       "prob.vec"=p.vec, "Se"=Se.display, "Sp"=Sp.display, 
       "Config"=list("Stage1"=save.it[1], "Stage2"=get.pools(group.mem[2,]), 
                     "Stage3"=get.pools(group.mem[3,]), "Stage4"=get.pools(group.mem[4,])), 
       "p.mat"=joint.p, "ET"=save.it[2], "value"=save.it[3], 
       "Accuracy"=list("Disease 1 Individual"=ind.acc1, 
                       "Disease 2 Individual"=ind.acc2, 
                       "Overall"=accuracy.mat))
}
