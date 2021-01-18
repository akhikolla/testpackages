# Start  Inf.Dorf.OTC2() function
###################################################################
# Brianna Hitt

# Brianna Hitt - 04.02.2020
# Changed cat() to message()

Inf.Dorf.OTC2 <- function(alpha=NULL, probabilities=NULL, Se, Sp, 
                          ordering, group.sz, updateProgress=NULL, 
                          trace=TRUE, print.time=TRUE, ...){
  start.time <- proc.time()
  
  set.of.blocks <- sort(group.sz)
  save.ET <- matrix(data=NA, nrow=length(set.of.blocks), 
                    ncol=5*max(set.of.blocks)+11)
  count <- 1
  
  for(N in set.of.blocks){
    # generate a joint probability vector for each individual
    if(!is.null(probabilities)){
      joint.p <- probabilities
    } else{
      joint.p <- joint.p.create(alpha=alpha, I=N)
      row.names(joint.p) <- joint.p.row.labels(ordering)
      colnames(joint.p) <- as.character(1:N)
    }
    
    # generate a matrix of all possible configurations/sets of pool sizes
    # the parts() function requires loading the partitions library
    # do not include the first column, which would test the initial group twice
    possible.groups <- parts(n=N)[,-1]
    
    save.it <- matrix(data=NA, nrow=ncol(possible.groups), 
                      ncol=5*max(set.of.blocks)+11)
    counter <- 1
    
    for(c in 1:ncol(possible.groups)){
      #develop the group membership matrix for each configuration
      gp.sizes <- possible.groups[,c]
      # remove the first row of the group membership matrix since the initial 
      #   pool is not tested in informative two-stage hierarchical testing
      group.mem <- group.membership.mat(gp.sizes)[-1,]
      rownames(group.mem) <- c("1", "2")
      # stage1.sizes <- group.membership.mat(gp.sizes)[2:3,]
      # group.mem <- matrix(data = c(stage1.sizes, 1:I), nrow = 2, ncol = I, 
      #                     byrow = TRUE, 
      #                     dimnames = list(Stage = 1:2, Individual = 1:I))
      
      save.info <- TwoDisease.Hierarchical(joint.p=joint.p, group.mem=group.mem, 
                                           Se=Se, Sp=Sp, ordering=ordering, a=1, 
                                           accuracy=TRUE)
      acc.dis1 <- save.info$Disease1$Overall
      acc.dis2 <- save.info$Disease2$Overall
      # acc.dis1 <- rep(NA, 4)
      # acc.dis2 <- rep(NA, 4)
      
      save.it[counter,] <- c(joint.p, rep(NA, max(0, 4*max(set.of.blocks)-length(joint.p))),
                             N, save.info$Expt, save.info$Expt/N, acc.dis1, acc.dis2, 
                             gp.sizes, rep(0, max(set.of.blocks)-length(gp.sizes)))
      counter <- counter + 1
    }
    
    # save the top configurations for each initial group size
    num.top <- 10
    if(N==set.of.blocks[1]){
      top.configs <- as.matrix(save.it[order(save.it[,(4*max(set.of.blocks)+3)]),])[1:min(nrow(save.it), num.top),]
    } else if(N>set.of.blocks[1]){
      top.configs <- rbind(top.configs, as.matrix(save.it[order(save.it[,(4*max(set.of.blocks)+3)]),])[1:min(nrow(save.it), num.top),])
    }
    
    # find the best configuration for each initial block size N, out of all possible configurations
    save.ET[count,] <- save.it[save.it[,(4*max(set.of.blocks)+3)]==min(save.it[,(4*max(set.of.blocks)+3)]),]
    
    if(is.function(updateProgress)){
      updateText <- paste0("Block Size = ", N)
      updateProgress(value = count/(length(set.of.blocks)+1), detail=updateText)
    }
    
    # print progress, if trace == TRUE
    if(trace){
      cat("Block Size =", N, "\n")
    }
    count <- count + 1
  }
  
  # reorder matrix of top configurations by E(T)/N
  top.configs <- as.matrix(top.configs)[order(top.configs[,(4*max(set.of.blocks)+3)]),]
  colnames(top.configs) <- c(rep(x="joint.p", times=4*max(set.of.blocks)),
                             "N", "ET", "value", "PSe1", "PSp1",  "PPPV1", 
                             "PNPV1", "PSe2", "PSp2", "PPPV2", "PNPV2", 
                             rep(x="pool.sz", times=max(set.of.blocks)))
  top.configs <- convert.config(algorithm="ID2", results=top.configs, diseases=2)
  
  # save the results for each initial group size
  if(length(set.of.blocks)==1){
    configs <- NA
  } else{
    configs <- as.matrix(save.ET)[order(save.ET[,(4*max(set.of.blocks)+3)]),]
    colnames(configs) <- c(rep(x="joint.p", times=4*max(set.of.blocks)),
                           "N", "ET", "value", "PSe1", "PSp1",  "PPPV1", 
                           "PNPV1", "PSe2", "PSp2", "PPPV2", "PNPV2", 
                           rep(x="pool.sz", times=max(set.of.blocks)))
    configs <- convert.config(algorithm="ID2", results=configs, diseases=2)
  }
  
  # find the optimal configuration, with the minimum value
  result.ET <- save.ET[save.ET[,(4*max(set.of.blocks)+3)]==min(save.ET[,(4*max(set.of.blocks)+3)]),]
  joint.p <- matrix(data=result.ET[1:(4*max(set.of.blocks))][!is.na(result.ET[1:(4*max(set.of.blocks))])], 
                    nrow=4, ncol=result.ET[4*max(set.of.blocks)+1])
  row.names(joint.p) <- joint.p.row.labels(ordering)
  colnames(joint.p) <- as.character(1:result.ET[4*max(set.of.blocks)+1])
  accuracy.mat <- matrix(data=result.ET[4*max(set.of.blocks)+4:11], nrow=2, ncol=4, byrow=TRUE, 
                         dimnames=list("Disease"=c("1", "2"), c("PSe", "PSp", "PPPV", "PNPV")))
  opt.ET <- list("OTC"=list("Block.sz"=result.ET[4*max(set.of.blocks)+1], 
                            "pool.szs"=(result.ET[(4*max(set.of.blocks)+12):length(result.ET)])[result.ET[(4*max(set.of.blocks)+12):length(result.ET)]!=0]), 
                 "p.mat"=joint.p, 
                 "ET"=result.ET[4*max(set.of.blocks)+2], "value"=result.ET[4*max(set.of.blocks)+3], 
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
  
  list("algorithm"="Informative two-stage hierarchical testing",
       "joint.p"=probabilities, "alpha.vec"=alpha, "Se"=Se.display, 
       "Sp"=Sp.display, "opt.ET"=opt.ET, "Configs"=configs, 
       "Top.Configs"=top.configs)
}
