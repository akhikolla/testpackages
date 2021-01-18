
# Start  NI.Dorf.OTC1() function
###################################################################

# Brianna Hitt - 4-17-17
# Updated: Brianna Hitt - 6-20-18

# Brianna Hitt - 04.02.2020
# Changed cat() to message()

NI.Dorf.OTC1 <- function(p, Se, Sp, group.sz, obj.fn, weights=NULL, 
                         updateProgress=NULL, trace=TRUE, print.time=TRUE){
  
  start.time<-proc.time()

  set.of.I <- group.sz
  save.it <- matrix(data=NA, nrow=length(set.of.I), ncol=15)
  count <- 1

  for(I in set.of.I){
    # generate a probability vector for homogeneous population
    p.vec <- rep(x=p[1], times=I)

    # calculate descriptive measures for two-stage hierarchical testing
    save.info <- hierarchical.desc2(p=p.vec, se=Se, sp=Sp, I2=NULL, order.p=FALSE)

    # extract ET, PSe, PSp and calculate the MAR function
    ET <- save.info$ET
    PSe.vec <- save.info$individual.testerror$pse.vec
    PSp.vec <- save.info$individual.testerror$psp.vec
    if("MAR" %in% obj.fn){
      MAR <- MAR.func(ET=ET, p.vec=p.vec, PSe.vec=PSe.vec, PSp.vec=PSp.vec)
    } else{MAR <- NA}

    # for non-informative Dorfman (two-stage hierarchical) testing, all individuals have the same testing accuracy measures
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
      GR1 <- GR.func(ET=ET, p.vec=p.vec, PSe.vec=PSe.vec, PSp.vec=PSp.vec, D1=weights[1,1], D2=weights[1,2])
      if(dim(weights)[1]>=2){
        GR2 <- GR.func(ET=ET, p.vec=p.vec, PSe.vec=PSe.vec, PSp.vec=PSp.vec, D1=weights[2,1], D2=weights[2,2])
      } else{GR2 <- NA}
      if(dim(weights)[1]>=3){
        GR3 <- GR.func(ET=ET, p.vec=p.vec, PSe.vec=PSe.vec, PSp.vec=PSp.vec, D1=weights[3,1], D2=weights[3,2])
      } else{GR3 <- NA}
      if(dim(weights)[1]>=4){
        GR4 <- GR.func(ET=ET, p.vec=p.vec, PSe.vec=PSe.vec, PSp.vec=PSp.vec, D1=weights[4,1], D2=weights[4,2])
      } else{GR4 <- NA}
      if(dim(weights)[1]>=5){
        GR5 <- GR.func(ET=ET, p.vec=p.vec, PSe.vec=PSe.vec, PSp.vec=PSp.vec, D1=weights[5,1], D2=weights[5,2])
      } else{GR5 <- NA}
      if(dim(weights)[1]>=6){
        GR6 <- GR.func(ET=ET, p.vec=p.vec, PSe.vec=PSe.vec, PSp.vec=PSp.vec, D1=weights[6,1], D2=weights[6,2])
      } else{GR6 <- NA}
    }

    save.it[count,] <- c(p[1], I, ET, ET/I, MAR, GR1/I, GR2/I, GR3/I, GR4/I, GR5/I, GR6/I, PSe, PSp, PPPV, PNPV)
    
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
    if(obj.fn[1]=="ET"){
      configs <- as.matrix(save.it[, c(1:3,4,12:ncol(save.it))])[order(save.it[,4]),]
    } else if(obj.fn[1]=="MAR"){
      configs <- as.matrix(save.it[, c(1:3,5,12:ncol(save.it))])[order(save.it[,5]),]
    } else if(obj.fn[1]=="GR"){
      configs <- as.matrix(save.it[, c(1:3,6,12:ncol(save.it))])[order(save.it[,6]),]
    }

    colnames(configs) <- c("p", "I", "ET", "value", "PSe", "PSp", "PPPV", "PNPV")
    configs <- convert.config(algorithm="D2", results=configs)
  }

  # find the testing configuration with the minimum value, for each objective function
  result.ET <- save.it[save.it[,4]==min(save.it[,4]), c(1:3,4,12:ncol(save.it))]
  result.MAR <- save.it[save.it[,5]==min(save.it[,5]), c(1:3,5,12:ncol(save.it))]
  result.GR1 <- save.it[save.it[,6]==min(save.it[,6]), c(1:3,6,12:ncol(save.it))]
  result.GR2 <- save.it[save.it[,7]==min(save.it[,7]), c(1:3,7,12:ncol(save.it))]
  result.GR3 <- save.it[save.it[,8]==min(save.it[,8]), c(1:3,8,12:ncol(save.it))]
  result.GR4 <- save.it[save.it[,9]==min(save.it[,9]), c(1:3,9,12:ncol(save.it))]
  result.GR5 <- save.it[save.it[,10]==min(save.it[,10]), c(1:3,10,12:ncol(save.it))]
  result.GR6 <- save.it[save.it[,11]==min(save.it[,11]), c(1:3,11,12:ncol(save.it))]

  p.vec.ET <- rep(x=result.ET[1], times=result.ET[2])
  if("MAR" %in% obj.fn){
    p.vec.MAR <- rep(x=result.MAR[1], times=result.MAR[2])
  } else{p.vec.MAR <- NA}
  if(is.null(dim(weights))){
    p.vec.GR1 <- NA
    p.vec.GR2 <- NA
    p.vec.GR3 <- NA
    p.vec.GR4 <- NA
    p.vec.GR5 <- NA
    p.vec.GR6 <- NA
  } else{
    p.vec.GR1 <- rep(x=result.GR1[1], times=result.GR1[2])
    if(dim(weights)[1]>=2){
      p.vec.GR2 <- rep(x=result.GR2[1], times=result.GR2[2])
    } else{p.vec.GR2 <- NA}
    if(dim(weights)[1]>=3){
      p.vec.GR3 <- rep(x=result.GR3[1], times=result.GR3[2])
    } else{p.vec.GR3 <- NA}
    if(dim(weights)[1]>=4){
      p.vec.GR4 <- rep(x=result.GR4[1], times=result.GR4[2])
    } else{p.vec.GR4 <- NA}
    if(dim(weights)[1]>=5){
      p.vec.GR5 <- rep(x=result.GR5[1], times=result.GR5[2])
    } else{p.vec.GR5 <- NA}
    if(dim(weights)[1]>=6){
      p.vec.GR6 <- rep(x=result.GR6[1], times=result.GR6[2])
    } else{p.vec.GR6 <- NA}
  }
  
  # put accuracy measures in a matrix for easier display of results
  acc.ET <- matrix(data=result.ET[5:8], nrow=1, ncol=4, 
                   dimnames=list(NULL, c("PSe", "PSp", "PPPV", "PNPV")))
  acc.MAR <- matrix(data=result.MAR[5:8], nrow=1, ncol=4, 
                   dimnames=list(NULL, c("PSe", "PSp", "PPPV", "PNPV")))
  acc.GR1 <- matrix(data=result.GR1[5:8], nrow=1, ncol=4, 
                   dimnames=list(NULL, c("PSe", "PSp", "PPPV", "PNPV")))
  acc.GR2 <- matrix(data=result.GR2[5:8], nrow=1, ncol=4, 
                   dimnames=list(NULL, c("PSe", "PSp", "PPPV", "PNPV")))
  acc.GR3 <- matrix(data=result.GR3[5:8], nrow=1, ncol=4, 
                   dimnames=list(NULL, c("PSe", "PSp", "PPPV", "PNPV")))
  acc.GR4 <- matrix(data=result.GR4[5:8], nrow=1, ncol=4, 
                   dimnames=list(NULL, c("PSe", "PSp", "PPPV", "PNPV")))
  acc.GR5 <- matrix(data=result.GR5[5:8], nrow=1, ncol=4, 
                   dimnames=list(NULL, c("PSe", "PSp", "PPPV", "PNPV")))
  acc.GR6 <- matrix(data=result.GR6[5:8], nrow=1, ncol=4, 
                   dimnames=list(NULL, c("PSe", "PSp", "PPPV", "PNPV")))
  
  
  # create a list of results for each objective function
  opt.ET <- list("OTC"=list("Stage1"=result.ET[2]), "p.vec"=p.vec.ET, "ET"=result.ET[3], "value"=result.ET[4], "Accuracy"=acc.ET)
  opt.MAR <- list("OTC"=list("Stage1"=result.MAR[2]), "p.vec"=p.vec.MAR, "ET"=result.MAR[3], "value"=result.MAR[4], "Accuracy"=acc.MAR)
  opt.GR1 <- list("OTC"=list("Stage1"=result.GR1[2]), "p.vec"=p.vec.GR1, "ET"=result.GR1[3], "value"=result.GR1[4], "Accuracy"=acc.GR1)
  opt.GR2 <- list("OTC"=list("Stage1"=result.GR2[2]), "p.vec"=p.vec.GR2, "ET"=result.GR2[3], "value"=result.GR2[4], "Accuracy"=acc.GR2)
  opt.GR3 <- list("OTC"=list("Stage1"=result.GR3[2]), "p.vec"=p.vec.GR3, "ET"=result.GR3[3], "value"=result.GR3[4], "Accuracy"=acc.GR3)
  opt.GR4 <- list("OTC"=list("Stage1"=result.GR4[2]), "p.vec"=p.vec.GR4, "ET"=result.GR4[3], "value"=result.GR4[4], "Accuracy"=acc.GR4)
  opt.GR5 <- list("OTC"=list("Stage1"=result.GR5[2]), "p.vec"=p.vec.GR5, "ET"=result.GR5[3], "value"=result.GR5[4], "Accuracy"=acc.GR5)
  opt.GR6 <- list("OTC"=list("Stage1"=result.GR6[2]), "p.vec"=p.vec.GR6, "ET"=result.GR6[3], "value"=result.GR6[4], "Accuracy"=acc.GR6)

  # create input accuracy value matrices for output display
  Se.display <- matrix(data = Se, nrow = 1, ncol = 2, 
                       dimnames = list(NULL, "Stage"=1:2))
  Sp.display <- matrix(data = Sp, nrow = 1, ncol = 2, 
                       dimnames = list(NULL, "Stage"=1:2))
  
  # create a list of results, including all objective functions
  opt.all <- list("opt.ET"=opt.ET, "opt.MAR"=opt.MAR, "opt.GR1"=opt.GR1, "opt.GR2"=opt.GR2,
                  "opt.GR3"=opt.GR3, "opt.GR4"=opt.GR4, "opt.GR5"=opt.GR5, "opt.GR6"=opt.GR6)
  # remove any objective functions not requested by the user
  opt.req <- Filter(function(x) !is.na(x$ET), opt.all)
  
  # print time elapsed, if print.time == TRUE
  if(print.time){
    time.it(start.time)
  }

  inputs <- list("algorithm"="Non-informative two-stage hierarchical testing",
                 "prob"=p, "Se"=Se.display, "Sp"=Sp.display)
  res <- c(inputs, opt.req)
  res[["Configs"]] <- configs
  res
}

###################################################################
