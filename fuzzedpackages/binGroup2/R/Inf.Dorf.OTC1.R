
# Start  inf.dorf.measures() function
###############################################################################
#    Brianna Hitt - 4-18-17
#    Purpose: calculates the testing expenditure and accuracy measures for 
#             informative Dorfman testing, given a testing configuration; 
#             Informative Dorfman testing is the same as pool-specific optimal 
#             Dorfman testing, described by McMahan et al. (2012), except for 
#             that it attempts to find the optimal testing configuration by 
#             examining all possible configurations rather than using the 
#             greedy algorithm proposed by McMahan et al. (2012)
#      calls: characteristics.pool() - calculates the testing expenditure for a 
#               given configuration
#             accuracy.dorf() - calculates measures of testing accuracy for 
#               informative Dorfman testing
#      inputs: p = probability that an individual is infected, can be an 
#                overall probability of disease or a vector of individual 
#                probabilities
#              Se = sensitivity of the diagnostic test
#              Sp = specificity of the diagnostic test
#              N = block size/initial group  size, up to 50
#              pool.sizes = a configuration/set of pool sizes from a matrix of 
#                all possible configurations, generated using the three-stage 
#                hierarchical setup
#      outputs: list of the testing expenditure (e), testing variance (v), and
#               measures of testing accuracy (summary)
#      notes: much of the code for this function, including the 
#               characteristics.pool() and accuracy.dorf() functions are from 
#               Chris McMahan's programs, provided with "Informative Dorfman 
#               screening" by McMahan, Tebbs, and Bilder (2012)

inf.dorf.measures <- function(prob, se, sp, N, pool.sizes){

  # Saves original ordering
  ind.order<-order(prob)

  # Orders subjects, required under all Informative algorithms
  prob<-sort(prob)

  # Determines the number of subjects being screened and sets up vectors for 
  #   storing summary measures
  N<-length(prob)
  pool.id<-rep(-1000, N)
  PSe<-rep(-100, N)
  PSp<-rep(-100, N)
  PPV<-rep(-100, N)
  NPV<-rep(-100, N)

  # pool.sizes is a single configuration/set of pool sizes from the matrix of 
  #   all possible configurations from the three-stage hierarchical testing 
  #   setup
  psz <- pool.sizes
  J <- length(psz)

  # Finds measures pool by pool
  psz<-c(psz, 0)
  lower<-1
  upper<-psz[1]
  vec.e<-rep(-1000, J)
  vec.v<-rep(-1000, J)
  for(i in 1:J){
    p.pool<-prob[lower:upper]
    pool.id[lower:upper]<-rep(i, length(p.pool))

    # calculates the testing expenditure and variance per individual
    res<-characteristics.pool(p=p.pool, se=se, sp=sp)
    vec.e[i]<-res$e
    vec.v[i]<-res$v

    # calculates measures of testing accuracy per individual
    res.acc<-accuracy.dorf(p=p.pool, se=se, sp=sp)
    PSe[lower:upper]<-res.acc$PSe
    PSp[lower:upper]<-res.acc$PSp
    PPV[lower:upper]<-res.acc$PPV
    NPV[lower:upper]<-res.acc$NPV

    lower<-1+upper
    upper<-upper+psz[i+1]
  }

  # Finds total expectation and variation
  res.e<-sum(vec.e)
  res.v<-sum(vec.v)

  # Returns all subjects to original ordering, along with their corresponding 
  #   measures
  prob<-prob[order(ind.order)]
  pool.id<-pool.id[order(ind.order)]
  PSe<-PSe[order(ind.order)]
  PSp<-PSp[order(ind.order)]
  PPV<-PPV[order(ind.order)]
  NPV<-NPV[order(ind.order)]

  res.mat<-matrix(c(pool.id, prob, PSe, PSp, PPV, NPV), nrow=N, ncol=6, 
                  byrow=FALSE, dimnames=list(as.character(1:N) , 
                                             c("pool", "probability", 
                                               "PSe", "PSp", "PPV", "NPV")))

  prob<-prob[ind.order]

  list("e"=res.e, "v"=res.v, "summary"=res.mat)
}

###############################################################################




# Start  Inf.Dorf.OTC1() function
###############################################################################

# Brianna Hitt - 4-18-17
# Updated: Brianna Hitt - 06-20-18

# Brianna Hitt - 04.02.2020
# Changed cat() to message()

Inf.Dorf.OTC1 <- function(p, Se, Sp, group.sz, obj.fn, weights=NULL, alpha=2, 
                          updateProgress=NULL, trace=TRUE, print.time=TRUE, ...){
  
  start.time<-proc.time()

  set.of.blocks <- group.sz

  save.ET <- matrix(data=NA, nrow=length(set.of.blocks), 
                    ncol=2*max(set.of.blocks)+8)
  save.MAR <- matrix(data=NA, nrow=length(set.of.blocks), 
                     ncol=2*max(set.of.blocks)+8)
  save.GR1 <- matrix(data=NA, nrow=length(set.of.blocks), 
                     ncol=2*max(set.of.blocks)+8)
  save.GR2 <- matrix(data=NA, nrow=length(set.of.blocks), 
                     ncol=2*max(set.of.blocks)+8)
  save.GR3 <- matrix(data=NA, nrow=length(set.of.blocks), 
                     ncol=2*max(set.of.blocks)+8)
  save.GR4 <- matrix(data=NA, nrow=length(set.of.blocks), 
                     ncol=2*max(set.of.blocks)+8)
  save.GR5 <- matrix(data=NA, nrow=length(set.of.blocks), 
                     ncol=2*max(set.of.blocks)+8)
  save.GR6 <- matrix(data=NA, nrow=length(set.of.blocks), 
                     ncol=2*max(set.of.blocks)+8)

  count <- 1

  for(N in set.of.blocks){
    # build a vector of probabilities for a heterogeneous population
    if(length(p)==1){
      p.vec <- expectOrderBeta(p=p, alpha=alpha, grp.sz=N, ...)
    } else if(length(p)>1){
      p.vec <- sort(p)
      alpha <- NA
    }

    # generate a matrix of all possible configurations/sets of pool sizes
    # the parts() function requires loading the partitions library
    # do not include the first column, which would test the initial group twice
    possible.groups <- parts(n=N)[,-1]

    save.it <- matrix(data=NA, nrow=ncol(possible.groups), 
                      ncol=2*max(set.of.blocks)+15)

    counter <- 1
    for(c in 1:ncol(possible.groups)){
      pool.sizes <- possible.groups[,c]
      # calculate descriptive measures for informative Dorfman testing, 
      #   given a configuration/set of pool sizes
      save.info <- inf.dorf.measures(prob=p.vec, se=Se, sp=Sp, N=N, 
                                     pool.sizes=pool.sizes[pool.sizes!=0])

      # extract the configuration/pool sizes
      row.names(save.info$summary)=NULL
      pool.sz <- table(save.info$summary[,1])
      row.names(pool.sz)=NULL

      # extract accuracy measures for each individual
      ET <- save.info$e
      PSe.vec <- save.info$summary[,3]
      PSp.vec <- save.info$summary[,4]
      if("MAR" %in% obj.fn){
        MAR <- MAR.func(ET=ET, p.vec=p.vec, PSe.vec=PSe.vec, PSp.vec=PSp.vec)
      } else{MAR <- NA}

      # calculate overall accuracy measures
      PSe <- sum(p.vec*PSe.vec)/sum(p.vec)
      PSp <- sum((1-p.vec)*(PSp.vec))/sum(1-p.vec)
      PPPV <- sum(p.vec*PSe.vec)/sum(p.vec*PSe.vec + (1-p.vec)*(1-PSp.vec))
      PNPV <- sum((1-p.vec)*PSp.vec)/sum((1-p.vec)*PSp.vec + p.vec*(1-PSe.vec))

      # for each row in the matrix of weights, calculate the GR function
      if(is.null(dim(weights))){
        GR1 <- NA
        GR2 <- NA
        GR3 <- NA
        GR4 <- NA
        GR5 <- NA
        GR6 <- NA
      } else{
        GR1 <- GR.func(ET=ET, p.vec=p.vec, PSe.vec=PSe.vec, PSp.vec=PSp.vec, 
                       D1=weights[1,1], D2=weights[1,2])
        if(dim(weights)[1]>=2){
          GR2 <- GR.func(ET=ET, p.vec=p.vec, PSe.vec=PSe.vec, PSp.vec=PSp.vec, 
                         D1=weights[2,1], D2=weights[2,2])
        } else{GR2 <- NA}
        if(dim(weights)[1]>=3){
          GR3 <- GR.func(ET=ET, p.vec=p.vec, PSe.vec=PSe.vec, PSp.vec=PSp.vec, 
                         D1=weights[3,1], D2=weights[3,2])
        } else{GR3 <- NA}
        if(dim(weights)[1]>=4){
          GR4 <- GR.func(ET=ET, p.vec=p.vec, PSe.vec=PSe.vec, PSp.vec=PSp.vec, 
                         D1=weights[4,1], D2=weights[4,2])
        } else{GR4 <- NA}
        if(dim(weights)[1]>=5){
          GR5 <- GR.func(ET=ET, p.vec=p.vec, PSe.vec=PSe.vec, PSp.vec=PSp.vec, 
                         D1=weights[5,1], D2=weights[5,2])
        } else{GR5 <- NA}
        if(dim(weights)[1]>=6){
          GR6 <- GR.func(ET=ET, p.vec=p.vec, PSe.vec=PSe.vec, PSp.vec=PSp.vec, 
                         D1=weights[6,1], D2=weights[6,2])
        } else{GR6 <- NA}
      }

      save.it[counter,] <- c(p.vec, rep(NA, max(0, max(set.of.blocks)-length(p.vec))), 
                             alpha, N, ET, ET/N, MAR, GR1/N, GR2/N, GR3/N, 
                             GR4/N, GR5/N, GR6/N, PSe, PSp, PPPV, PNPV, pool.sz, 
                             rep(0, max(0, max(set.of.blocks)-length(pool.sz))))
      counter <- counter + 1

    }

    # save the top configurations from each initial group/block size
    num.top <- 10
    if(N==set.of.blocks[1]){
      if(obj.fn[1]=="ET"){
        top.configs <- as.matrix(save.it[order(save.it[,(max(set.of.blocks)+4)]),])[1:min(nrow(save.it), num.top), c(1:(max(set.of.blocks)+3),(max(set.of.blocks)+4),(max(set.of.blocks)+12):ncol(save.it))]
      } else if(obj.fn[1]=="MAR"){
        top.configs <- as.matrix(save.it[order(save.it[,(max(set.of.blocks)+5)]),])[1:min(nrow(save.it), num.top), c(1:(max(set.of.blocks)+3),(max(set.of.blocks)+5),(max(set.of.blocks)+12):ncol(save.it))]
      } else if(obj.fn[1]=="GR"){
        top.configs <- as.matrix(save.it[order(save.it[,(max(set.of.blocks)+6)]),])[1:min(nrow(save.it), num.top), c(1:(max(set.of.blocks)+3),(max(set.of.blocks)+6),(max(set.of.blocks)+12):ncol(save.it))]
      }
    } else if(N>set.of.blocks[1]){
      if(obj.fn[1]=="ET"){
        top.configs <- rbind(top.configs, as.matrix(save.it[order(save.it[,(max(set.of.blocks)+4)]),])[1:min(nrow(save.it), num.top), c(1:(max(set.of.blocks)+3),(max(set.of.blocks)+4),(max(set.of.blocks)+12):ncol(save.it))])
      } else if(obj.fn[1]=="MAR"){
        top.configs <- rbind(top.configs, as.matrix(save.it[order(save.it[,(max(set.of.blocks)+5)]),])[1:min(nrow(save.it), num.top), c(1:(max(set.of.blocks)+3),(max(set.of.blocks)+5),(max(set.of.blocks)+12):ncol(save.it))])
      } else if(obj.fn[1]=="GR"){
        top.configs <- rbind(top.configs, as.matrix(save.it[order(save.it[,(max(set.of.blocks)+6)]),])[1:min(nrow(save.it), num.top), c(1:(max(set.of.blocks)+3),(max(set.of.blocks)+6),(max(set.of.blocks)+12):ncol(save.it))])
      }
    }

    # find the best configuration for each block size N, out of all possible configurations
    save.ET[count,] <- save.it[save.it[,(max(set.of.blocks)+4)]==min(save.it[,(max(set.of.blocks)+4)]), c(1:(max(set.of.blocks)+3),(max(set.of.blocks)+4),(max(set.of.blocks)+12):ncol(save.it))]
    if(class(try(save.MAR[count,] <- save.it[save.it[,(max(set.of.blocks)+5)]==min(save.it[,(max(set.of.blocks)+5)]), c(1:(max(set.of.blocks)+3),(max(set.of.blocks)+5),(max(set.of.blocks)+12):ncol(save.it))],silent=T))!="try-error"){
      save.MAR[count,] <- save.it[save.it[,(max(set.of.blocks)+5)]==min(save.it[,(max(set.of.blocks)+5)]), c(1:(max(set.of.blocks)+3),(max(set.of.blocks)+5),(max(set.of.blocks)+12):ncol(save.it))]
    }
    if(class(try(save.GR1[count,] <- save.it[save.it[,(max(set.of.blocks)+6)]==min(save.it[,(max(set.of.blocks)+6)]), c(1:(max(set.of.blocks)+3),(max(set.of.blocks)+6),(max(set.of.blocks)+12):ncol(save.it))],silent=T))!="try-error"){
      save.GR1[count,] <- save.it[save.it[,(max(set.of.blocks)+6)]==min(save.it[,(max(set.of.blocks)+6)]), c(1:(max(set.of.blocks)+3),(max(set.of.blocks)+6),(max(set.of.blocks)+12):ncol(save.it))]
    }
    if(class(try(save.GR2[count,] <- save.it[save.it[,(max(set.of.blocks)+7)]==min(save.it[,(max(set.of.blocks)+7)]), c(1:(max(set.of.blocks)+3),(max(set.of.blocks)+7),(max(set.of.blocks)+12):ncol(save.it))],silent=T))!="try-error"){
      save.GR2[count,] <- save.it[save.it[,(max(set.of.blocks)+7)]==min(save.it[,(max(set.of.blocks)+7)]), c(1:(max(set.of.blocks)+3),(max(set.of.blocks)+7),(max(set.of.blocks)+12):ncol(save.it))]
    }
    if(class(try(save.GR3[count,] <- save.it[save.it[,(max(set.of.blocks)+8)]==min(save.it[,(max(set.of.blocks)+8)]), c(1:(max(set.of.blocks)+3),(max(set.of.blocks)+8),(max(set.of.blocks)+12):ncol(save.it))],silent=T))!="try-error"){
      save.GR3[count,] <- save.it[save.it[,(max(set.of.blocks)+8)]==min(save.it[,(max(set.of.blocks)+8)]), c(1:(max(set.of.blocks)+3),(max(set.of.blocks)+8),(max(set.of.blocks)+12):ncol(save.it))]
    }
    if(class(try(save.GR4[count,] <- save.it[save.it[,(max(set.of.blocks)+9)]==min(save.it[,(max(set.of.blocks)+9)]), c(1:(max(set.of.blocks)+3),(max(set.of.blocks)+9),(max(set.of.blocks)+12):ncol(save.it))],silent=T))!="try-error"){
      save.GR4[count,] <- save.it[save.it[,(max(set.of.blocks)+9)]==min(save.it[,(max(set.of.blocks)+9)]), c(1:(max(set.of.blocks)+3),(max(set.of.blocks)+9),(max(set.of.blocks)+12):ncol(save.it))]
    }
    if(class(try( save.GR5[count,] <- save.it[save.it[,(max(set.of.blocks)+10)]==min(save.it[,(max(set.of.blocks)+10)]), c(1:(max(set.of.blocks)+3),(max(set.of.blocks)+10),(max(set.of.blocks)+12):ncol(save.it))],silent=T))!="try-error"){
      save.GR5[count,] <- save.it[save.it[,(max(set.of.blocks)+10)]==min(save.it[,(max(set.of.blocks)+10)]), c(1:(max(set.of.blocks)+3),(max(set.of.blocks)+10),(max(set.of.blocks)+12):ncol(save.it))]
    }
    if(class(try(save.GR6[count,] <- save.it[save.it[,(max(set.of.blocks)+11)]==min(save.it[,(max(set.of.blocks)+11)]), c(1:(max(set.of.blocks)+3),(max(set.of.blocks)+11),(max(set.of.blocks)+12):ncol(save.it))],silent=T))!="try-error"){
      save.GR6[count,] <- save.it[save.it[,(max(set.of.blocks)+11)]==min(save.it[,(max(set.of.blocks)+11)]), c(1:(max(set.of.blocks)+3),(max(set.of.blocks)+11),(max(set.of.blocks)+12):ncol(save.it))]
    }

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

  # reorder matrix of top configurations by E(T)/I
  top.configs <- top.configs[order(top.configs[,(max(set.of.blocks)+4)]),]
  colnames(top.configs) <- c(rep(x = "p", times = max(set.of.blocks)), 
                             "alpha", "N", "ET", "value", "PSe", "PSp", 
                             "PPPV", "PNPV", rep(x = "pool.sz", times=max(set.of.blocks)))
  top.configs <- convert.config(algorithm="ID2", results=top.configs)

  # save the best configuration for each initial group/block size
  if(length(set.of.blocks)==1){
    configs <- NA
  } else{
    if(obj.fn[1]=="ET"){
      configs <- as.matrix(save.ET)[order(save.ET[,(max(set.of.blocks)+4)]),]
    } else if(obj.fn[1]=="MAR"){
      configs <- as.matrix(save.MAR)[order(save.MAR[,(max(set.of.blocks)+4)]),]
    } else if(obj.fn[1]=="GR"){
      configs <- as.matrix(save.GR1)[order(save.GR1[,(max(set.of.blocks)+4)]),]
    }

    colnames(configs) <- c(rep(x = "p", times = max(set.of.blocks)), 
                           "alpha", "N", "ET", "value", "PSe", "PSp", 
                           "PPPV", "PNPV", rep(x = "pool.sz", times=max(set.of.blocks)))
    configs <- convert.config(algorithm="ID2", results=configs)
  }

  # find the optimal configuration over all block sizes considered
  result.ET <- save.ET[save.ET[,(max(set.of.blocks)+4)]==min(save.ET[,(max(set.of.blocks)+4)]),]
  result.MAR <- save.MAR[save.MAR[,(max(set.of.blocks)+4)]==min(save.MAR[,(max(set.of.blocks)+4)]),]
  result.GR1 <- save.GR1[save.GR1[,(max(set.of.blocks)+4)]==min(save.GR1[,(max(set.of.blocks)+4)]),]
  result.GR2 <- save.GR2[save.GR2[,(max(set.of.blocks)+4)]==min(save.GR2[,(max(set.of.blocks)+4)]),]
  result.GR3 <- save.GR3[save.GR3[,(max(set.of.blocks)+4)]==min(save.GR3[,(max(set.of.blocks)+4)]),]
  result.GR4 <- save.GR4[save.GR4[,(max(set.of.blocks)+4)]==min(save.GR4[,(max(set.of.blocks)+4)]),]
  result.GR5 <- save.GR5[save.GR5[,(max(set.of.blocks)+4)]==min(save.GR5[,(max(set.of.blocks)+4)]),]
  result.GR6 <- save.GR6[save.GR6[,(max(set.of.blocks)+4)]==min(save.GR6[,(max(set.of.blocks)+4)]),]

  p.vec.ET <- (result.ET[1:max(set.of.blocks)])[!is.na(result.ET[1:max(set.of.blocks)])]
  if("MAR" %in% obj.fn){
    p.vec.MAR <- (result.MAR[1:max(set.of.blocks)])[!is.na(result.MAR[1:max(set.of.blocks)])]
  } else{p.vec.MAR <- NA}
  if(is.null(dim(weights))){
    p.vec.GR1 <- NA
    p.vec.GR2 <- NA
    p.vec.GR3 <- NA
    p.vec.GR4 <- NA
    p.vec.GR5 <- NA
    p.vec.GR6 <- NA
  } else{
    p.vec.GR1 <- (result.GR1[1:max(set.of.blocks)])[!is.na(result.GR1[1:max(set.of.blocks)])]
    if(dim(weights)[1]>=2){
      p.vec.GR2 <- (result.GR2[1:max(set.of.blocks)])[!is.na(result.GR2[1:max(set.of.blocks)])]
    } else{p.vec.GR2 <- NA}
    if(dim(weights)[1]>=3){
      p.vec.GR3 <- (result.GR3[1:max(set.of.blocks)])[!is.na(result.GR3[1:max(set.of.blocks)])]
    } else{p.vec.GR3 <- NA}
    if(dim(weights)[1]>=4){
      p.vec.GR4 <- (result.GR4[1:max(set.of.blocks)])[!is.na(result.GR4[1:max(set.of.blocks)])]
    } else{p.vec.GR4 <- NA}
    if(dim(weights)[1]>=5){
      p.vec.GR5 <- (result.GR5[1:max(set.of.blocks)])[!is.na(result.GR5[1:max(set.of.blocks)])]
    } else{p.vec.GR5 <- NA}
    if(dim(weights)[1]>=6){
      p.vec.GR6 <- (result.GR6[1:max(set.of.blocks)])[!is.na(result.GR6[1:max(set.of.blocks)])]
    } else{p.vec.GR6 <- NA}
  }
  
  # put accuracy measures in a matrix for easier display of results
  acc.ET <- matrix(data=result.ET[(max(set.of.blocks)+5:8)], nrow=1, ncol=4, 
                   dimnames=list(NULL, c("PSe", "PSp", "PPPV", "PNPV")))
  acc.MAR <- matrix(data=result.MAR[(max(set.of.blocks)+5:8)], nrow=1, ncol=4, 
                    dimnames=list(NULL, c("PSe", "PSp", "PPPV", "PNPV")))
  acc.GR1 <- matrix(data=result.GR1[(max(set.of.blocks)+5:8)], nrow=1, ncol=4, 
                    dimnames=list(NULL, c("PSe", "PSp", "PPPV", "PNPV")))
  acc.GR2 <- matrix(data=result.GR2[(max(set.of.blocks)+5:8)], nrow=1, ncol=4, 
                    dimnames=list(NULL, c("PSe", "PSp", "PPPV", "PNPV")))
  acc.GR3 <- matrix(data=result.GR3[(max(set.of.blocks)+5:8)], nrow=1, ncol=4, 
                    dimnames=list(NULL, c("PSe", "PSp", "PPPV", "PNPV")))
  acc.GR4 <- matrix(data=result.GR4[(max(set.of.blocks)+5:8)], nrow=1, ncol=4, 
                    dimnames=list(NULL, c("PSe", "PSp", "PPPV", "PNPV")))
  acc.GR5 <- matrix(data=result.GR5[(max(set.of.blocks)+5:8)], nrow=1, ncol=4, 
                    dimnames=list(NULL, c("PSe", "PSp", "PPPV", "PNPV")))
  acc.GR6 <- matrix(data=result.GR6[(max(set.of.blocks)+5:8)], nrow=1, ncol=4, 
                    dimnames=list(NULL, c("PSe", "PSp", "PPPV", "PNPV")))

  # create a list of results for each objective function
  opt.ET <- list("OTC"=list("Block.sz"=result.ET[(max(set.of.blocks)+2)],
                            "pool.szs"=(result.ET[(max(set.of.blocks)+9):length(result.ET)])[result.ET[(max(set.of.blocks)+9):length(result.ET)]!=0]), 
                 "p.vec"=p.vec.ET, "ET"=result.ET[(max(set.of.blocks)+3)], 
                 "value"=result.ET[(max(set.of.blocks)+4)], "Accuracy"=acc.ET)
  opt.MAR <- list("OTC"=list("Block.sz"=result.MAR[(max(set.of.blocks)+2)], 
                             "pool.szs"=(result.MAR[(max(set.of.blocks)+9):length(result.MAR)])[result.MAR[(max(set.of.blocks)+9):length(result.MAR)]!=0]), 
                  "p.vec"=p.vec.MAR, "ET"=result.MAR[(max(set.of.blocks)+3)], 
                  "value"=result.MAR[(max(set.of.blocks)+4)], "Accuracy"=acc.MAR)
  opt.GR1 <- list("OTC"=list("Block.sz"=result.GR1[(max(set.of.blocks)+2)], 
                             "pool.szs"=(result.GR1[(max(set.of.blocks)+9):length(result.GR1)])[result.GR1[(max(set.of.blocks)+9):length(result.GR1)]!=0]), 
                  "p.vec"=p.vec.GR1, "ET"=result.GR1[(max(set.of.blocks)+3)], 
                  "value"=result.GR1[(max(set.of.blocks)+4)], "Accuracy"=acc.GR1)
  opt.GR2 <- list("OTC"=list("Block.sz"=result.GR2[(max(set.of.blocks)+2)], 
                             "pool.szs"=(result.GR2[(max(set.of.blocks)+9):length(result.GR2)])[result.GR2[(max(set.of.blocks)+9):length(result.GR2)]!=0]), 
                  "p.vec"=p.vec.GR2, "ET"=result.GR2[(max(set.of.blocks)+3)], 
                  "value"=result.GR2[(max(set.of.blocks)+4)], "Accuracy"=acc.GR2)
  opt.GR3 <- list("OTC"=list("Block.sz"=result.GR3[(max(set.of.blocks)+2)], 
                             "pool.szs"=(result.GR3[(max(set.of.blocks)+9):length(result.GR3)])[result.GR3[(max(set.of.blocks)+9):length(result.GR3)]!=0]), 
                  "p.vec"=p.vec.GR3, "ET"=result.GR3[(max(set.of.blocks)+3)], 
                  "value"=result.GR3[(max(set.of.blocks)+4)], "Accuracy"=acc.GR3)
  opt.GR4 <- list("OTC"=list("Block.sz"=result.GR4[(max(set.of.blocks)+2)], 
                             "pool.szs"=(result.GR4[(max(set.of.blocks)+9):length(result.GR4)])[result.GR4[(max(set.of.blocks)+9):length(result.GR4)]!=0]), 
                  "p.vec"=p.vec.GR4, "ET"=result.GR4[(max(set.of.blocks)+3)], 
                  "value"=result.GR4[(max(set.of.blocks)+4)], "Accuracy"=acc.GR4)
  opt.GR5 <- list("OTC"=list("Block.sz"=result.GR5[(max(set.of.blocks)+2)], 
                             "pool.szs"=(result.GR5[(max(set.of.blocks)+9):length(result.GR5)])[result.GR5[(max(set.of.blocks)+9):length(result.GR5)]!=0]), 
                  "p.vec"=p.vec.GR5, "ET"=result.GR5[(max(set.of.blocks)+3)], 
                  "value"=result.GR5[(max(set.of.blocks)+4)], "Accuracy"=acc.GR5)
  opt.GR6 <- list("OTC"=list("Block.sz"=result.GR6[(max(set.of.blocks)+2)], 
                             "pool.szs"=(result.GR6[(max(set.of.blocks)+9):length(result.GR6)])[result.GR6[(max(set.of.blocks)+9):length(result.GR6)]!=0]), 
                  "p.vec"=p.vec.GR6, "ET"=result.GR6[(max(set.of.blocks)+3)], 
                  "value"=result.GR6[(max(set.of.blocks)+4)], "Accuracy"=acc.GR6)

  # create input accuracy value matrices for output display
  Se.display <- matrix(data = Se, nrow = 1, ncol = 2, 
                       dimnames = list(NULL, "Stage"=1:2))
  Sp.display <- matrix(data = Sp, nrow = 1, ncol = 2, 
                       dimnames = list(NULL, "Stage"=1:2))
  
  # create a list of results, including all objective functions
  opt.all <- list("opt.ET"=opt.ET, "opt.MAR"=opt.MAR, "opt.GR1"=opt.GR1, 
                  "opt.GR2"=opt.GR2,  "opt.GR3"=opt.GR3, "opt.GR4"=opt.GR4, 
                  "opt.GR5"=opt.GR5, "opt.GR6"=opt.GR6)
  # remove any objective functions not requested by the user
  opt.req <- Filter(function(x) !is.na(x$ET), opt.all)
  
  # print time elapsed, if print.time == TRUE
  if(print.time){
    time.it(start.time)
  }

  inputs <- list("algorithm"="Informative two-stage hierarchical testing",
                 "prob"=list(p), "alpha"=alpha, "Se"=Se.display, "Sp"=Sp.display)
  res <- c(inputs, opt.req)
  res[["Configs"]] <- configs
  res[["Top.Configs"]] <- top.configs
  res
}

###############################################################################
