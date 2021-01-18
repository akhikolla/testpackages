
# Start  NI.Array.OTC1() function
###################################################################

# Brianna Hitt - 05-01-17
# Updated: Brianna Hitt - 06-20-18

# Brianna Hitt - 04.02.2020
# Changed cat() to message()

NI.Array.OTC1 <- function(p, Se, Sp, group.sz, obj.fn, weights=NULL, 
                          updateProgress=NULL, 
                          trace=TRUE, print.time=TRUE, ...){

  start.time<-proc.time()

  set.of.I <- group.sz

  save.it <- matrix(data=NA, nrow=length(set.of.I), ncol=16)
  count <- 1

  for(I in set.of.I){
    N <- I^2

    # build a matrix of probabilities
    # this is the same for an overall probability p and for a vector p
    p.mat <- matrix(data=p[1], nrow=I, ncol=I)

    # call Array.Measures to calculate descriptive measures for the given array size
    save.info <- Array.Measures(p=p.mat, se=Se, sp=Sp)

    # extract accuracy measures for each individual
    ET <- save.info$ET
    PSe.mat <- save.info$PSe
    PSp.mat <- save.info$PSp
    if("MAR" %in% obj.fn){
      MAR <- MAR.func(ET=ET, p.vec=p.mat, PSe.vec=PSe.mat, PSp.vec=PSp.mat)
    } else{MAR <- NA}

    # calculate overall accuracy measures
    PSe <- sum(p.mat*PSe.mat)/sum(p.mat)
    PSp <- sum((1-p.mat)*(PSp.mat))/sum(1-p.mat)
    PPPV <- sum(p.mat*PSe.mat)/sum(p.mat*PSe.mat + (1-p.mat)*(1-PSp.mat))
    PNPV <- sum((1-p.mat)*PSp.mat)/sum((1-p.mat)*PSp.mat + p.mat*(1-PSe.mat))

    # for each row in the matrix of weights, calculate the GR function
    if(is.null(dim(weights))){
      GR1 <- NA
      GR2 <- NA
      GR3 <- NA
      GR4 <- NA
      GR5 <- NA
      GR6 <- NA
    } else{
      GR1 <- GR.func(ET=ET, p.vec=p.mat, PSe.vec=PSe.mat, PSp.vec=PSp.mat, D1=weights[1,1], D2=weights[1,2])
      if(dim(weights)[1]>=2){
        GR2 <- GR.func(ET=ET, p.vec=p.mat, PSe.vec=PSe.mat, PSp.vec=PSp.mat, D1=weights[2,1], D2=weights[2,2])
      } else{GR2 <- NA}
      if(dim(weights)[1]>=3){
        GR3 <- GR.func(ET=ET, p.vec=p.mat, PSe.vec=PSe.mat, PSp.vec=PSp.mat, D1=weights[3,1], D2=weights[3,2])
      } else{GR3 <- NA}
      if(dim(weights)[1]>=4){
        GR4 <- GR.func(ET=ET, p.vec=p.mat, PSe.vec=PSe.mat, PSp.vec=PSp.mat, D1=weights[4,1], D2=weights[4,2])
      } else{GR4 <- NA}
      if(dim(weights)[1]>=5){
        GR5 <- GR.func(ET=ET, p.vec=p.mat, PSe.vec=PSe.mat, PSp.vec=PSp.mat, D1=weights[5,1], D2=weights[5,2])
      } else{GR5 <- NA}
      if(dim(weights)[1]>=6){
        GR6 <- GR.func(ET=ET, p.vec=p.mat, PSe.vec=PSe.mat, PSp.vec=PSp.mat, D1=weights[6,1], D2=weights[6,2])
      } else{GR6 <- NA}
    }

    save.it[count,] <- c(p[1], I, N, ET, ET/N, MAR, GR1/N, GR2/N, GR3/N, GR4/N, GR5/N, GR6/N, PSe, PSp, PPPV, PNPV)

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
    configs <- NA
  } else{
    if(obj.fn[1]=="ET"){
      configs <- as.matrix(save.it[, c(1:4,5,13:ncol(save.it))])[order(save.it[,5]),]
    } else if(obj.fn[1]=="MAR"){
      configs <- as.matrix(save.it[, c(1:4,6,13:ncol(save.it))])[order(save.it[,6]),]
    } else if(obj.fn[1]=="GR"){
      configs <- as.matrix(save.it[, c(1:4,7,13:ncol(save.it))])[order(save.it[,7]),]
    }

    colnames(configs) <- c("p", "I", "N", "ET", "value", "PSe", "PSp", "PPPV", "PNPV")
    configs <- convert.config(algorithm="A2", results=configs)
  }

  # find the optimal testing configuration, over all array sizes considered
  result.ET <- save.it[save.it[,5]==min(save.it[,5]),c(1:4,5,13:ncol(save.it))]
  result.MAR <- save.it[save.it[,6]==min(save.it[,6]),c(1:4,6,13:ncol(save.it))]
  result.GR1 <- save.it[save.it[,7]==min(save.it[,7]),c(1:4,7,13:ncol(save.it))]
  result.GR2 <- save.it[save.it[,8]==min(save.it[,8]),c(1:4,8,13:ncol(save.it))]
  result.GR3 <- save.it[save.it[,9]==min(save.it[,9]),c(1:4,9,13:ncol(save.it))]
  result.GR4 <- save.it[save.it[,10]==min(save.it[,10]),c(1:4,10,13:ncol(save.it))]
  result.GR5 <- save.it[save.it[,11]==min(save.it[,11]),c(1:4,11,13:ncol(save.it))]
  result.GR6 <- save.it[save.it[,12]==min(save.it[,12]),c(1:4,12,13:ncol(save.it))]

  p.mat.ET <- matrix(data=result.ET[1], nrow=result.ET[2], ncol=result.ET[2])
  if("MAR" %in% obj.fn){
    p.mat.MAR <- matrix(data=result.MAR[1], nrow=result.MAR[2], ncol=result.MAR[2])
  } else{p.mat.MAR <- NA}
  if(is.null(dim(weights))){
    p.mat.GR1 <- NA
    p.mat.GR2 <- NA
    p.mat.GR3 <- NA
    p.mat.GR4 <- NA
    p.mat.GR5 <- NA
    p.mat.GR6 <- NA
  } else{
    p.mat.GR1 <- matrix(data=result.GR1[1], nrow=result.GR1[2], ncol=result.GR1[2])
    if(dim(weights)[1]>=2){
      p.mat.GR2 <- matrix(data=result.GR2[1], nrow=result.GR2[2], ncol=result.GR2[2])
    } else{p.mat.GR2 <- NA}
    if(dim(weights)[1]>=3){
      p.mat.GR3 <- matrix(data=result.GR3[1], nrow=result.GR3[2], ncol=result.GR3[2])
    } else{p.mat.GR3 <- NA}
    if(dim(weights)[1]>=4){
      p.mat.GR4 <- matrix(data=result.GR4[1], nrow=result.GR4[2], ncol=result.GR4[2])
    } else{p.mat.GR4 <- NA}
    if(dim(weights)[1]>=5){
      p.mat.GR5 <- matrix(data=result.GR5[1], nrow=result.GR5[2], ncol=result.GR5[2])
    } else{p.mat.GR5 <- NA}
    if(dim(weights)[1]>=6){
      p.mat.GR6 <- matrix(data=result.GR6[1], nrow=result.GR6[2], ncol=result.GR6[2])
    } else{p.mat.GR6 <- NA}
  }
  
  # put accuracy measures in a matrix for easier display of results
  acc.ET <- matrix(data=result.ET[6:9], nrow=1, ncol=4, 
                   dimnames=list(NULL, c("PSe", "PSp", "PPPV", "PNPV")))
  acc.MAR <- matrix(data=result.MAR[6:9], nrow=1, ncol=4, 
                    dimnames=list(NULL, c("PSe", "PSp", "PPPV", "PNPV")))
  acc.GR1 <- matrix(data=result.GR1[6:9], nrow=1, ncol=4, 
                    dimnames=list(NULL, c("PSe", "PSp", "PPPV", "PNPV")))
  acc.GR2 <- matrix(data=result.GR2[6:9], nrow=1, ncol=4, 
                    dimnames=list(NULL, c("PSe", "PSp", "PPPV", "PNPV")))
  acc.GR3 <- matrix(data=result.GR3[6:9], nrow=1, ncol=4, 
                    dimnames=list(NULL, c("PSe", "PSp", "PPPV", "PNPV")))
  acc.GR4 <- matrix(data=result.GR4[6:9], nrow=1, ncol=4, 
                    dimnames=list(NULL, c("PSe", "PSp", "PPPV", "PNPV")))
  acc.GR5 <- matrix(data=result.GR5[6:9], nrow=1, ncol=4, 
                    dimnames=list(NULL, c("PSe", "PSp", "PPPV", "PNPV")))
  acc.GR6 <- matrix(data=result.GR6[6:9], nrow=1, ncol=4, 
                    dimnames=list(NULL, c("PSe", "PSp", "PPPV", "PNPV")))

  # create a list of results for each objective function
  opt.ET <- list("OTC"=list("Array.dim"=result.ET[2], "Array.sz"=result.ET[3]), 
                 "p.mat"=p.mat.ET, "ET"=result.ET[4], "value"=result.ET[5], 
                 "Accuracy"=acc.ET)
  opt.MAR <- list("OTC"=list("Array.dim"=result.MAR[2], "Array.sz"=result.MAR[3]), 
                  "p.mat"=p.mat.MAR, "ET"=result.MAR[4], "value"=result.MAR[5], 
                  "Accuracy"=acc.MAR)
  opt.GR1 <- list("OTC"=list("Array.dim"=result.GR1[2], "Array.sz"=result.GR1[3]), 
                  "p.mat"=p.mat.GR1, "ET"=result.GR1[4], "value"=result.GR1[5], 
                  "Accuracy"=acc.GR1)
  opt.GR2 <- list("OTC"=list("Array.dim"=result.GR2[2], "Array.sz"=result.GR2[3]), 
                  "p.mat"=p.mat.GR2, "ET"=result.GR2[4], "value"=result.GR2[5], 
                  "Accuracy"=acc.GR2)
  opt.GR3 <- list("OTC"=list("Array.dim"=result.GR3[2], "Array.sz"=result.GR3[3]), 
                  "p.mat"=p.mat.GR3, "ET"=result.GR3[4], "value"=result.GR3[5], 
                  "Accuracy"=acc.GR3)
  opt.GR4 <- list("OTC"=list("Array.dim"=result.GR4[2], "Array.sz"=result.GR4[3]), 
                  "p.mat"=p.mat.GR4, "ET"=result.GR4[4], "value"=result.GR4[5], 
                  "Accuracy"=acc.GR4)
  opt.GR5 <- list("OTC"=list("Array.dim"=result.GR5[2], "Array.sz"=result.GR5[3]), 
                  "p.mat"=p.mat.GR5, "ET"=result.GR5[4], "value"=result.GR5[5], 
                  "Accuracy"=acc.GR5)
  opt.GR6 <- list("OTC"=list("Array.dim"=result.GR6[2], "Array.sz"=result.GR6[3]), 
                  "p.mat"=p.mat.GR6, "ET"=result.GR6[4], "value"=result.GR6[5], 
                  "Accuracy"=acc.GR6)

  # create input accuracy value matrices for output display
  Se.display <- matrix(data = Se, nrow = 1, ncol = 2, 
                       dimnames = list(NULL, "Test"=c("Row/Column", "Individual")))
  Sp.display <- matrix(data = Sp, nrow = 1, ncol = 2, 
                       dimnames = list(NULL, "Test"=c("Row/Column", "Individual")))
  # use below if Se/Sp for row and column testing is allowed to differ
  # Se.display <- matrix(data = Se, nrow = 1, ncol = 3, 
  #                      dimnames = list(NULL, "Test"=c("Row", "Column", "Individual")))
  # Sp.display <- matrix(data = Sp, nrow = 1, ncol = 3, 
  #                      dimnames = list(NULL, "Test"=c("Row", "Column", "Individual")))
  
  # create a list of results, including all objective functions
  opt.all <- list("opt.ET"=opt.ET, "opt.MAR"=opt.MAR, "opt.GR1"=opt.GR1, "opt.GR2"=opt.GR2,
                  "opt.GR3"=opt.GR3, "opt.GR4"=opt.GR4, "opt.GR5"=opt.GR5, "opt.GR6"=opt.GR6)
  # remove any objective functions not requested by the user
  opt.req <- Filter(function(x) !is.na(x$ET), opt.all)
  
  # print time elapsed, if print.time == TRUE
  if(print.time){
    time.it(start.time)
  }

  inputs <- list("algorithm"="Non-informative array testing without master pooling",
                 "prob"=p, "Se"=Se.display, "Sp"=Sp.display)
  res <- c(inputs, opt.req)
  res[["Configs"]] <- configs
  res
}

###################################################################
