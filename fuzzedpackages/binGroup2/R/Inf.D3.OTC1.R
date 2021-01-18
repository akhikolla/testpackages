# Start  Inf.D3.OTC1() function
###############################################################################
# Brianna Hitt - 05-01-17
# Updated: Brianna Hitt - 06-20-18

# Brianna Hitt - 04.02.2020
# Changed cat() to message()

Inf.D3.OTC1 <- function(p, Se, Sp, group.sz, obj.fn, weights=NULL, alpha=2, 
                        updateProgress=NULL, trace=TRUE, print.time=TRUE, ...){
  
  start.time<-proc.time()

  set.of.I <- group.sz

  save.ET <- matrix(data=NA, nrow=length(set.of.I), ncol=2*max(set.of.I)+8)
  save.MAR <- matrix(data=NA, nrow=length(set.of.I), ncol=2*max(set.of.I)+8)
  save.GR1 <- matrix(data=NA, nrow=length(set.of.I), ncol=2*max(set.of.I)+8)
  save.GR2 <- matrix(data=NA, nrow=length(set.of.I), ncol=2*max(set.of.I)+8)
  save.GR3 <- matrix(data=NA, nrow=length(set.of.I), ncol=2*max(set.of.I)+8)
  save.GR4 <- matrix(data=NA, nrow=length(set.of.I), ncol=2*max(set.of.I)+8)
  save.GR5 <- matrix(data=NA, nrow=length(set.of.I), ncol=2*max(set.of.I)+8)
  save.GR6 <- matrix(data=NA, nrow=length(set.of.I), ncol=2*max(set.of.I)+8)

  count <- 1

  for(I in set.of.I){
    # build a vector of probabilities for a heterogeneous population
    if(length(p)==1){
      p.vec <- expectOrderBeta(p=p, alpha=alpha, grp.sz=I, ...)
    } else if(length(p)>1){
      p.vec <- sort(p)
      alpha <- NA
    }

    # generate a matrix of all possible configurations/sets of pool sizes
    # the parts() function requires loading the partitions library
    # do not include the first column, which would test the initial group twice
    possible.groups <- parts(n=I)[,-1]

    save.it <- matrix(data=NA, nrow=ncol(possible.groups), 
                      ncol=2*max(set.of.I)+15)

    counter <- 1
    for(c in 1:ncol(possible.groups)){
      # extract the configuration, ordering, and group sizes for each column
      config <- c(possible.groups[,c], 1:I)
      order.for.p <- config[(1+I):(2*I)]
      gp.sizes <- config[1:I]

      # call hierarchical.desc2() for the configuration
      save.info <- hierarchical.desc2(p=p.vec[order.for.p], se=Se, sp=Sp, 
                                      I2=gp.sizes[gp.sizes!=0], order.p=FALSE)

      # extract accuracy measures for each individual
      ET <- save.info$ET
      PSe.vec <- save.info$individual.testerror$pse.vec
      PSp.vec <- save.info$individual.testerror$psp.vec
      if("MAR" %in% obj.fn){
        MAR <- MAR.func(ET=ET, p.vec=p.vec, PSe.vec=PSe.vec, PSp.vec=PSp.vec)
      } else{MAR <- NA}

      # calculate overall accuracy measures
      group.testerror <- save.info$group.testerror
      names(group.testerror) <- NULL
      PSe <- group.testerror[1]
      PSp <- group.testerror[2]
      PPPV <- group.testerror[3]
      PNPV <- group.testerror[4]

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

      save.it[counter,] <- c(p.vec, rep(NA, max(0, max(set.of.I)-length(p.vec))), 
                             alpha, I, ET, ET/I, MAR, GR1/I, GR2/I, GR3/I, 
                             GR4/I, GR5/I, GR6/I, PSe, PSp, PPPV, PNPV, gp.sizes, 
                             rep(0, max(0, max(set.of.I)-length(gp.sizes))))
      counter <- counter + 1

    }

    # save the top configurations for each initial group size
    num.top <- 10
    if(I == set.of.I[1]){
      if(obj.fn[1]=="ET"){
        top.configs <- as.matrix(save.it[order(save.it[,(max(set.of.I)+4)]),])[1:min(nrow(save.it), num.top), c(1:(max(set.of.I)+3),(max(set.of.I)+4),(max(set.of.I)+12):ncol(save.it))]
      } else if(obj.fn[1]=="MAR"){
        top.configs <- as.matrix(save.it[order(save.it[,(max(set.of.I)+5)]),])[1:min(nrow(save.it), num.top), c(1:(max(set.of.I)+3),(max(set.of.I)+5),(max(set.of.I)+12):ncol(save.it))]
      } else if(obj.fn[1]=="GR"){
        top.configs <- as.matrix(save.it[order(save.it[,(max(set.of.I)+6)]),])[1:min(nrow(save.it), num.top), c(1:(max(set.of.I)+3),(max(set.of.I)+6),(max(set.of.I)+12):ncol(save.it))]
      }
    } else if(I > set.of.I[1]){
      if(obj.fn[1]=="ET"){
        top.configs <- rbind(top.configs, as.matrix(save.it[order(save.it[,(max(set.of.I)+4)]),])[1:min(nrow(save.it), num.top), c(1:(max(set.of.I)+3),(max(set.of.I)+4),(max(set.of.I)+12):ncol(save.it))])
      } else if(obj.fn[1]=="MAR"){
        top.configs <- rbind(top.configs, as.matrix(save.it[order(save.it[,(max(set.of.I)+5)]),])[1:min(nrow(save.it), num.top), c(1:(max(set.of.I)+3),(max(set.of.I)+5),(max(set.of.I)+12):ncol(save.it))])
      } else if(obj.fn[1]=="GR"){
        top.configs <- rbind(top.configs, as.matrix(save.it[order(save.it[,(max(set.of.I)+6)]),])[1:min(nrow(save.it), num.top), c(1:(max(set.of.I)+3),(max(set.of.I)+6),(max(set.of.I)+12):ncol(save.it))])
      }
    }

    # find the best configuration for each initial group size I, out of all possible configurations
    save.ET[count,] <- save.it[save.it[,(max(set.of.I)+4)]==min(save.it[,(max(set.of.I)+4)]), c(1:(max(set.of.I)+3),(max(set.of.I)+4),(max(set.of.I)+12):ncol(save.it))]
    if(class(try(save.MAR[count,] <- save.it[save.it[,(max(set.of.I)+5)]==min(save.it[,(max(set.of.I)+5)]), c(1:(max(set.of.I)+3),(max(set.of.I)+5),(max(set.of.I)+12):ncol(save.it))],silent=T))!="try-error"){
      save.MAR[count,] <- save.it[save.it[,(max(set.of.I)+5)]==min(save.it[,(max(set.of.I)+5)]), c(1:(max(set.of.I)+3),(max(set.of.I)+5),(max(set.of.I)+12):ncol(save.it))]
    }
    if(class(try(save.GR1[count,] <- save.it[save.it[,(max(set.of.I)+6)]==min(save.it[,(max(set.of.I)+6)]), c(1:(max(set.of.I)+3),(max(set.of.I)+6),(max(set.of.I)+12):ncol(save.it))],silent=T))!="try-error"){
      save.GR1[count,] <- save.it[save.it[,(max(set.of.I)+6)]==min(save.it[,(max(set.of.I)+6)]), c(1:(max(set.of.I)+3),(max(set.of.I)+6),(max(set.of.I)+12):ncol(save.it))]
    }
    if(class(try(save.GR2[count,] <- save.it[save.it[,(max(set.of.I)+7)]==min(save.it[,(max(set.of.I)+7)]), c(1:(max(set.of.I)+3),(max(set.of.I)+7),(max(set.of.I)+12):ncol(save.it))],silent=T))!="try-error"){
      save.GR2[count,] <- save.it[save.it[,(max(set.of.I)+7)]==min(save.it[,(max(set.of.I)+7)]), c(1:(max(set.of.I)+3),(max(set.of.I)+7),(max(set.of.I)+12):ncol(save.it))]
    }
    if(class(try(save.GR3[count,] <- save.it[save.it[,(max(set.of.I)+8)]==min(save.it[,(max(set.of.I)+8)]), c(1:(max(set.of.I)+3),(max(set.of.I)+8),(max(set.of.I)+12):ncol(save.it))],silent=T))!="try-error"){
      save.GR3[count,] <- save.it[save.it[,(max(set.of.I)+8)]==min(save.it[,(max(set.of.I)+8)]), c(1:(max(set.of.I)+3),(max(set.of.I)+8),(max(set.of.I)+12):ncol(save.it))]
    }
    if(class(try(save.GR4[count,] <- save.it[save.it[,(max(set.of.I)+9)]==min(save.it[,(max(set.of.I)+9)]), c(1:(max(set.of.I)+3),(max(set.of.I)+9),(max(set.of.I)+12):ncol(save.it))],silent=T))!="try-error"){
      save.GR4[count,] <- save.it[save.it[,(max(set.of.I)+9)]==min(save.it[,(max(set.of.I)+9)]), c(1:(max(set.of.I)+3),(max(set.of.I)+9),(max(set.of.I)+12):ncol(save.it))]
    }
    if(class(try( save.GR5[count,] <- save.it[save.it[,(max(set.of.I)+10)]==min(save.it[,(max(set.of.I)+10)]), c(1:(max(set.of.I)+3),(max(set.of.I)+10),(max(set.of.I)+12):ncol(save.it))],silent=T))!="try-error"){
      save.GR5[count,] <- save.it[save.it[,(max(set.of.I)+10)]==min(save.it[,(max(set.of.I)+10)]), c(1:(max(set.of.I)+3),(max(set.of.I)+10),(max(set.of.I)+12):ncol(save.it))]
    }
    if(class(try(save.GR6[count,] <- save.it[save.it[,(max(set.of.I)+11)]==min(save.it[,(max(set.of.I)+11)]), c(1:(max(set.of.I)+3),(max(set.of.I)+11),(max(set.of.I)+12):ncol(save.it))],silent=T))!="try-error"){
      save.GR6[count,] <- save.it[save.it[,(max(set.of.I)+11)]==min(save.it[,(max(set.of.I)+11)]), c(1:(max(set.of.I)+3),(max(set.of.I)+11),(max(set.of.I)+12):ncol(save.it))]
    }

    if(is.function(updateProgress)){
      updateText <- paste0("Initial Group Size = ", I)
      updateProgress(value = count/(length(set.of.I)+1), detail=updateText)
    }
    
    # print the progress, if trace == TRUE
    if(trace){
      cat("Initial Group Size =", I, "\n")
    }
    
    count <- count + 1
  }

  # reorder matrix of top configurations by E(T)/I
  top.configs <- top.configs[order(top.configs[,(max(set.of.I)+4)]),]
  colnames(top.configs) <- c(rep(x = "p", times = max(set.of.I)), 
                             "alpha", "I", "ET", "value", "PSe", "PSp", 
                             "PPPV", "PNPV", rep(x = "pool.sz", times=max(set.of.I)))
  top.configs <- convert.config(algorithm="ID3", results=top.configs)

  # save the best configuration for each initial group size
  if(length(set.of.I)==1){
    configs <- NA
  } else{
    if(obj.fn[1]=="ET"){
      configs <- as.matrix(save.ET)[order(save.ET[,(max(set.of.I)+4)]),]
    } else if(obj.fn[1]=="MAR"){
      configs <- as.matrix(save.MAR)[order(save.MAR[,(max(set.of.I)+4)]),]
    } else if(obj.fn[1]=="GR"){
      configs <- as.matrix(save.GR1)[order(save.GR1[,(max(set.of.I)+4)]),]
    }

    colnames(configs) <- c(rep(x = "p", times = max(set.of.I)), 
                           "alpha", "I", "ET", "value", "PSe", "PSp", 
                           "PPPV", "PNPV", rep(x = "pool.sz", times=max(set.of.I)))
    configs <- convert.config(algorithm="ID3", results=configs)
  }

  # find the optimal testing configuration, over all initial group sizes considered
  result.ET <- save.ET[save.ET[,(max(set.of.I)+4)]==min(save.ET[,(max(set.of.I)+4)]),]
  result.MAR <- save.MAR[save.MAR[,(max(set.of.I)+4)]==min(save.MAR[,(max(set.of.I)+4)]),]
  result.GR1 <- save.GR1[save.GR1[,(max(set.of.I)+4)]==min(save.GR1[,(max(set.of.I)+4)]),]
  result.GR2 <- save.GR2[save.GR2[,(max(set.of.I)+4)]==min(save.GR2[,(max(set.of.I)+4)]),]
  result.GR3 <- save.GR3[save.GR3[,(max(set.of.I)+4)]==min(save.GR3[,(max(set.of.I)+4)]),]
  result.GR4 <- save.GR4[save.GR4[,(max(set.of.I)+4)]==min(save.GR4[,(max(set.of.I)+4)]),]
  result.GR5 <- save.GR5[save.GR5[,(max(set.of.I)+4)]==min(save.GR5[,(max(set.of.I)+4)]),]
  result.GR6 <- save.GR6[save.GR6[,(max(set.of.I)+4)]==min(save.GR6[,(max(set.of.I)+4)]),]

  p.vec.ET <- (result.ET[1:max(set.of.I)])[!is.na(result.ET[1:max(set.of.I)])]
  if("MAR" %in% obj.fn){
    p.vec.MAR <- (result.MAR[1:max(set.of.I)])[!is.na(result.MAR[1:max(set.of.I)])]
  } else{p.vec.MAR <- NA}
  if(is.null(dim(weights))){
    p.vec.GR1 <- NA
    p.vec.GR2 <- NA
    p.vec.GR3 <- NA
    p.vec.GR4 <- NA
    p.vec.GR5 <- NA
    p.vec.GR6 <- NA
  } else{
    p.vec.GR1 <- (result.GR1[1:max(set.of.I)])[!is.na(result.GR1[1:max(set.of.I)])]
    if(dim(weights)[1]>=2){
      p.vec.GR2 <- (result.GR2[1:max(set.of.I)])[!is.na(result.GR2[1:max(set.of.I)])]
    } else{p.vec.GR2 <- NA}
    if(dim(weights)[1]>=3){
      p.vec.GR3 <- (result.GR3[1:max(set.of.I)])[!is.na(result.GR3[1:max(set.of.I)])]
    } else{p.vec.GR3 <- NA}
    if(dim(weights)[1]>=4){
      p.vec.GR4 <- (result.GR4[1:max(set.of.I)])[!is.na(result.GR4[1:max(set.of.I)])]
    } else{p.vec.GR4 <- NA}
    if(dim(weights)[1]>=5){
      p.vec.GR5 <- (result.GR5[1:max(set.of.I)])[!is.na(result.GR5[1:max(set.of.I)])]
    } else{p.vec.GR5 <- NA}
    if(dim(weights)[1]>=6){
      p.vec.GR6 <- (result.GR6[1:max(set.of.I)])[!is.na(result.GR6[1:max(set.of.I)])]
    } else{p.vec.GR6 <- NA}
  }

  # put accuracy measures in a matrix for easier display of results
  acc.ET <- matrix(data=result.ET[(max(set.of.I)+5:8)], nrow=1, ncol=4, 
                   dimnames=list(NULL, c("PSe", "PSp", "PPPV", "PNPV")))
  acc.MAR <- matrix(data=result.MAR[(max(set.of.I)+5:8)], nrow=1, ncol=4, 
                    dimnames=list(NULL, c("PSe", "PSp", "PPPV", "PNPV")))
  acc.GR1 <- matrix(data=result.GR1[(max(set.of.I)+5:8)], nrow=1, ncol=4, 
                    dimnames=list(NULL, c("PSe", "PSp", "PPPV", "PNPV")))
  acc.GR2 <- matrix(data=result.GR2[(max(set.of.I)+5:8)], nrow=1, ncol=4, 
                    dimnames=list(NULL, c("PSe", "PSp", "PPPV", "PNPV")))
  acc.GR3 <- matrix(data=result.GR3[(max(set.of.I)+5:8)], nrow=1, ncol=4, 
                    dimnames=list(NULL, c("PSe", "PSp", "PPPV", "PNPV")))
  acc.GR4 <- matrix(data=result.GR4[(max(set.of.I)+5:8)], nrow=1, ncol=4, 
                    dimnames=list(NULL, c("PSe", "PSp", "PPPV", "PNPV")))
  acc.GR5 <- matrix(data=result.GR5[(max(set.of.I)+5:8)], nrow=1, ncol=4, 
                    dimnames=list(NULL, c("PSe", "PSp", "PPPV", "PNPV")))
  acc.GR6 <- matrix(data=result.GR6[(max(set.of.I)+5:8)], nrow=1, ncol=4, 
                    dimnames=list(NULL, c("PSe", "PSp", "PPPV", "PNPV")))
  
  # create a list of results for each objective function
  opt.ET <- list("OTC"=list("Stage1"=result.ET[(max(set.of.I)+2)], 
                            "Stage2"=(result.ET[(max(set.of.I)+9):length(result.ET)])[result.ET[(max(set.of.I)+9):length(result.ET)]!=0]), 
                 "p.vec"=p.vec.ET, "ET"=result.ET[(max(set.of.I)+3)], 
                 "value"=result.ET[(max(set.of.I)+4)], "Accuracy"=acc.ET)
  opt.MAR <- list("OTC"=list("Stage1"=result.MAR[(max(set.of.I)+2)], 
                             "Stage2"=(result.MAR[(max(set.of.I)+9):length(result.MAR)])[result.MAR[(max(set.of.I)+9):length(result.MAR)]!=0]), 
                  "p.vec"=p.vec.MAR, "ET"=result.MAR[(max(set.of.I)+3)], 
                  "value"=result.MAR[(max(set.of.I)+4)], "Accuracy"=acc.MAR)
  opt.GR1 <- list("OTC"=list("Stage1"=result.GR1[(max(set.of.I)+2)], 
                             "Stage2"=(result.GR1[(max(set.of.I)+9):length(result.GR1)])[result.GR1[(max(set.of.I)+9):length(result.GR1)]!=0]), 
                  "p.vec"=p.vec.GR1, "ET"=result.GR1[(max(set.of.I)+3)], 
                  "value"=result.GR1[(max(set.of.I)+4)], "Accuracy"=acc.GR1)
  opt.GR2 <- list("OTC"=list("Stage1"=result.GR2[(max(set.of.I)+2)], 
                             "Stage2"=(result.GR2[(max(set.of.I)+9):length(result.GR2)])[result.GR2[(max(set.of.I)+9):length(result.GR2)]!=0]), 
                  "p.vec"=p.vec.GR2, "ET"=result.GR2[(max(set.of.I)+3)], 
                  "value"=result.GR2[(max(set.of.I)+4)], "Accuracy"=acc.GR2)
  opt.GR3 <- list("OTC"=list("Stage1"=result.GR3[(max(set.of.I)+2)], 
                             "Stage2"=(result.GR3[(max(set.of.I)+9):length(result.GR3)])[result.GR3[(max(set.of.I)+9):length(result.GR3)]!=0]), 
                  "p.vec"=p.vec.GR3, "ET"=result.GR3[(max(set.of.I)+3)], 
                  "value"=result.GR3[(max(set.of.I)+4)], "Accuracy"=acc.GR3)
  opt.GR4 <- list("OTC"=list("Stage1"=result.GR4[(max(set.of.I)+2)], 
                             "Stage2"=(result.GR4[(max(set.of.I)+9):length(result.GR4)])[result.GR4[(max(set.of.I)+9):length(result.GR4)]!=0]), 
                  "p.vec"=p.vec.GR4,  "ET"=result.GR4[(max(set.of.I)+3)], 
                  "value"=result.GR4[(max(set.of.I)+4)], "Accuracy"=acc.GR4)
  opt.GR5 <- list("OTC"=list("Stage1"=result.GR5[(max(set.of.I)+2)], 
                             "Stage2"=(result.GR5[(max(set.of.I)+9):length(result.GR5)])[result.GR5[(max(set.of.I)+9):length(result.GR5)]!=0]), 
                  "p.vec"=p.vec.GR5, "ET"=result.GR5[(max(set.of.I)+3)], 
                  "value"=result.GR5[(max(set.of.I)+4)], "Accuracy"=acc.GR5)
  opt.GR6 <- list("OTC"=list("Stage1"=result.GR6[(max(set.of.I)+2)], 
                             "Stage2"=(result.GR6[(max(set.of.I)+9):length(result.GR6)])[result.GR6[(max(set.of.I)+9):length(result.GR6)]!=0]), 
                  "p.vec"=p.vec.GR6, "ET"=result.GR6[(max(set.of.I)+3)], 
                  "value"=result.GR6[(max(set.of.I)+4)], "Accuracy"=acc.GR6)

  # create input accuracy value matrices for output display
  Se.display <- matrix(data = Se, nrow = 1, ncol = 3, 
                       dimnames = list(NULL, "Stage"=1:3))
  Sp.display <- matrix(data = Sp, nrow = 1, ncol = 3, 
                       dimnames = list(NULL, "Stage"=1:3))
  
  # create a list of results, including all objective functions
  opt.all <- list("opt.ET"=opt.ET, "opt.MAR"=opt.MAR, "opt.GR1"=opt.GR1, 
                  "opt.GR2"=opt.GR2, "opt.GR3"=opt.GR3, "opt.GR4"=opt.GR4, 
                  "opt.GR5"=opt.GR5, "opt.GR6"=opt.GR6)
  # remove any objective functions not requested by the user
  opt.req <- Filter(function(x) !is.na(x$ET), opt.all)
  
  # print the time elapsed, if print.time == TRUE
  if(print.time){
    time.it(start.time)
  }

  inputs <- list("algorithm"="Informative three-stage hierarchical testing",
                 "prob"=list(p), "alpha"=alpha, "Se"=Se.display, "Sp"=Sp.display)
  res <- c(inputs, opt.req)
  res[["Configs"]] <- configs
  res[["Top.Configs"]] <- top.configs
  res
}

###############################################################################

