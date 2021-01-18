# Start  NI.Dorf.OTC2() function
###################################################################
# Brianna Hitt

# Brianna Hitt - 04.02.2020
# Changed cat() to message()

NI.Dorf.OTC2 <- function(p.vec, Se, Sp, ordering, group.sz, 
                         updateProgress=NULL, trace=TRUE, print.time=TRUE, ...){
  start.time <- proc.time()
  
  set.of.I <- sort(group.sz)
  save.it <- matrix(data=NA, nrow=length(set.of.I), ncol=11)
  count <- 1
  
  for(I in set.of.I){
    # generate a joint probability vector for each individual
    joint.p <- matrix(data=rep(p.vec, times=I), nrow=4, ncol=I, byrow=FALSE, 
                      dimnames=list(joint.p.row.labels(ordering), 
                                    as.character(1:I)))
    
    # calculate descriptive measures for two-stage hierarchical testing
    group.mem <- matrix(data = c(rep(x = 1, times = I), 1:I), nrow = 2, 
                        ncol = I, byrow = TRUE, 
                        dimnames = list(Stage = 1:2, Individual = 1:I))
    save.info <- TwoDisease.Hierarchical(joint.p=joint.p, group.mem=group.mem, 
                                         Se=Se, Sp=Sp, ordering=ordering, a=1, 
                                         accuracy=TRUE)
    acc.dis1 <- save.info$Disease1$Overall
    acc.dis2 <- save.info$Disease2$Overall
    
    save.it[count,] <- c(I, save.info$Expt, save.info$Expt/I, 
                         acc.dis1, acc.dis2)
    
    if(is.function(updateProgress)){
      updateText <- paste0("Initial Group Size = ", I)
      updateProgress(value = count/(length(set.of.I)+1), detail=updateText)
    }
    
    # print progress, if trace == TRUE
    if(trace){
      cat("Initial Group Size =", I, "\n")
    }
    
    count <- count + 1
  }
  
  # save the results for each initial group size
  if(length(set.of.I)==1){
    configs <- NA
  } else{
    configs <- as.matrix(save.it[,])[order(save.it[,3]),]
    colnames(configs) <- c("I", "ET", "value", "PSe1", "PSp1", 
                           "PPPV1", "PNPV1", "PSe2", "PSp2", "PPPV2", "PNPV2")
    configs <- convert.config(algorithm="D2", results=configs, diseases=2)
  }
  
  # find the optimal configuration, with the minimum value
  results <- save.it[save.it[,3]==min(save.it[,3]),]
  joint.p <- matrix(data=rep(p.vec, each=results[1]), nrow=4, ncol=results[1], 
                    byrow=TRUE, 
                    dimnames=list(joint.p.row.labels(ordering), 
                                  as.character(1:results[1])))
  accuracy.mat <- matrix(data=results[4:11], nrow=2, ncol=4, byrow=TRUE, 
                         dimnames=list("Disease"=c("1", "2"), 
                                       c("PSe", "PSp", "PPPV", "PNPV")))
  opt.ET <- list("OTC"=list("Stage1"=results[1]), "p.mat"=joint.p, 
                 "ET"=results[2], "value"=results[3], 
                 "Accuracy"=accuracy.mat)
  Se.display <- matrix(data = Se, nrow = 2, ncol = 2, 
                       dimnames = list("Disease" = c("1", "2"), 
                                       c("Stage 1", "Stage 2")))
  Sp.display <- matrix(data = Sp, nrow = 2, ncol = 2, 
                       dimnames = list("Disease" = c("1", "2"), 
                                      c("Stage 1", "Stage 2")))
  
  # print time elapsed, if print.time == TRUE
  if(print.time){
    time.it(start.time)
  }
  
  list("algorithm"="Non-informative two-stage hierarchical testing",
       "prob.vec"=p.vec, "Se"=Se.display, "Sp"=Sp.display, "opt.ET"=opt.ET, 
       "Configs"=configs)
}

#