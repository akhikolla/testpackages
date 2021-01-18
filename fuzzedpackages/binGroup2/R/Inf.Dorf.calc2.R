# Start  Inf.Dorf.calc2() function
###################################################################
# Brianna Hitt

Inf.Dorf.calc2 <- function(alpha=NULL, probabilities=NULL, 
                           Se, Sp, ordering, group.mem, a, 
                           trace=TRUE, print.time=TRUE, ...){
  start.time <- proc.time()
  
  N <- ncol(group.mem)

  # generate a joint probability vector for each individual
  if(!is.null(probabilities)){
    joint.p <- probabilities
  } else{
    joint.p <- joint.p.create(alpha=alpha, I=N)
    row.names(joint.p) <- joint.p.row.labels(ordering)
    colnames(joint.p) <- as.character(1:N)
  }
  
  save.info <- TwoDisease.Hierarchical(joint.p=joint.p, group.mem=group.mem, 
                                       Se=Se, Sp=Sp, ordering=ordering, 
                                       a=1:N, accuracy=TRUE)
  ind.acc1.all <- save.info$Disease1$Individual
  rownames(ind.acc1.all) <- ind.acc1.all[, "Individual"]
  ind.acc1 <- get.unique.index(ind.acc1.all[a,], 
                               which(colnames(ind.acc1.all)=="PSp"), 
                               rowlabel = a)[,-1]
  ind.acc2.all <- save.info$Disease2$Individual
  rownames(ind.acc2.all) <- ind.acc2.all[, "Individual"]
  ind.acc2 <- get.unique.index(ind.acc2.all[a,], 
                               which(colnames(ind.acc2.all)=="PSp"), 
                               rowlabel = a)[,-1]
  acc.dis1 <- save.info$Disease1$Overall
  acc.dis2 <- save.info$Disease2$Overall

  save.it <- c(N, save.info$Expt, save.info$Expt/N, acc.dis1, acc.dis2)
  
  accuracy.mat <- matrix(data=save.it[4:11], nrow=2, ncol=4, byrow=TRUE, 
                         dimnames=list("Disease" = 1:2, 
                                       c("PSe", "PSp", "PPPV", "PNPV")))

  Se.display <- matrix(data = Se, nrow = 2, ncol = 2, 
                       dimnames = list("Disease" = 1:2, "Stage" = 1:2))
  Sp.display <- matrix(data = Sp, nrow = 2, ncol = 2, 
                       dimnames = list("Disease" = 1:2, "Stage" = 1:2))
  
  # print time elapsed, if print.time == TRUE
  if(print.time){
    time.it(start.time)
  }
  
  list("algorithm"="Informative two-stage hierarchical testing",
       "joint.p"=probabilities, "alpha.vec"=alpha, 
       "Se"=Se.display, "Sp"=Sp.display, 
       "Config"=list("Block.sz"=save.it[1], "pool.szs"=get.pools(group.mem[1,])), 
       "p.mat"=joint.p, "ET"=save.it[2], "value"=save.it[3], 
       "Accuracy"=list("Disease 1 Individual"=ind.acc1, 
                       "Disease 2 Individual"=ind.acc2, 
                       "Overall"=accuracy.mat))
}

#
