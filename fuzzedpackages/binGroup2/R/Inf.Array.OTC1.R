# Start  Inf.Array.OTC1() function
###############################################################################
# Brianna Hitt - 05-01-17
# Updated: Brianna Hitt - 06-20-18

# Brianna Hitt - 04.02.2020
# Changed cat() to message()

Inf.Array.OTC1 <- function(p, Se, Sp, group.sz, obj.fn, weights=NULL, alpha=2, 
                           updateProgress=NULL, trace=TRUE, print.time=TRUE, ...){
  
  start.time<-proc.time()

  set.of.I <- group.sz

  save.it <- matrix(data=NA, nrow=length(set.of.I), ncol=max(set.of.I)^2+16)
  count <- 1

  for(I in set.of.I){
    N <- I^2

    # build a vector of probabilities for a heterogeneous population
    if(length(p)==1){
      p.vec <- expectOrderBeta(p=p, alpha=alpha, grp.sz=N, ...)
    } else if(length(p)>1){
      p.vec <- sort(p)
      alpha <- NA
    }

    # build a matrix of probabilities using the gradient design
    p.ga <- informativeArrayProb(prob.vec=p.vec, nr=I, nc=I, method="gd")

    # call Array.Measures() to calculate descriptive measures for the given 
    #   array size
    save.info <- Array.Measures(p=p.ga, se=Se, sp=Sp)

    # extract accuracy measures for each individual
    ET <- save.info$ET
    PSe.mat <- save.info$PSe
    PSp.mat <- save.info$PSp
    if("MAR" %in% obj.fn){
      MAR <- MAR.func(ET=ET, p.vec=p.ga, PSe.vec=PSe.mat, PSp.vec=PSp.mat)
    } else{MAR <- NA}

    # calculate overall accuracy measures
    PSe <- sum(p.ga*PSe.mat)/sum(p.ga)
    PSp <- sum((1-p.ga)*(PSp.mat))/sum(1-p.ga)
    PPPV <- sum(p.ga*PSe.mat)/sum(p.ga*PSe.mat + (1-p.ga)*(1-PSp.mat))
    PNPV <- sum((1-p.ga)*PSp.mat)/sum((1-p.ga)*PSp.mat + p.ga*(1-PSe.mat))

    # for each row in the matrix of weights, calculate the GR function
    if(is.null(dim(weights))){
      GR1 <- NA
      GR2 <- NA
      GR3 <- NA
      GR4 <- NA
      GR5 <- NA
      GR6 <- NA
    } else{
      GR1 <- GR.func(ET=ET, p.vec=p.ga, PSe.vec=PSe.mat, PSp.vec=PSp.mat, 
                     D1=weights[1,1], D2=weights[1,2])
      if(dim(weights)[1]>=2){
        GR2 <- GR.func(ET=ET, p.vec=p.ga, PSe.vec=PSe.mat, PSp.vec=PSp.mat, 
                       D1=weights[2,1], D2=weights[2,2])
      } else{GR2 <- NA}
      if(dim(weights)[1]>=3){
        GR3 <- GR.func(ET=ET, p.vec=p.ga, PSe.vec=PSe.mat, PSp.vec=PSp.mat, 
                       D1=weights[3,1], D2=weights[3,2])
      } else{GR3 <- NA}
      if(dim(weights)[1]>=4){
        GR4 <- GR.func(ET=ET, p.vec=p.ga, PSe.vec=PSe.mat, PSp.vec=PSp.mat, 
                       D1=weights[4,1], D2=weights[4,2])
      } else{GR4 <- NA}
      if(dim(weights)[1]>=5){
        GR5 <- GR.func(ET=ET, p.vec=p.ga, PSe.vec=PSe.mat, PSp.vec=PSp.mat, 
                       D1=weights[5,1], D2=weights[5,2])
      } else{GR5 <- NA}
      if(dim(weights)[1]>=6){
        GR6 <- GR.func(ET=ET, p.vec=p.ga, PSe.vec=PSe.mat, PSp.vec=PSp.mat, 
                       D1=weights[6,1], D2=weights[6,2])
      } else{GR6 <- NA}
    }

    save.it[count,] <- c(p.vec, rep(NA, max(0, max(set.of.I)^2-length(p.vec))), 
                         alpha, I, N, ET, ET/N, MAR, GR1/N, GR2/N, GR3/N, 
                         GR4/N, GR5/N, GR6/N, PSe, PSp, PPPV, PNPV)
    
    if(is.function(updateProgress)){
      updateText <- paste0("Row/Column Size=", I, ", Array Size=", N)
      updateProgress(value = count/(length(set.of.I)+1), detail=updateText)
    }
    
    # print the progress, if trace == TRUE
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
      configs <- (save.it[, c(1:(max(set.of.I)^2+4),(max(set.of.I)^2+5),
                              (max(set.of.I)^2+13):ncol(save.it))])[order(save.it[,(max(set.of.I)^2+5)]),]
    } else if(obj.fn[1]=="MAR"){
      configs <- (save.it[, c(1:(max(set.of.I)^2+4),(max(set.of.I)^2+6),
                              (max(set.of.I)^2+13):ncol(save.it))])[order(save.it[,(max(set.of.I)^2+6)]),]
    } else if(obj.fn[1]=="GR"){
      configs <- (save.it[, c(1:(max(set.of.I)^2+4),(max(set.of.I)^2+7),
                              (max(set.of.I)^2+13):ncol(save.it))])[order(save.it[,(max(set.of.I)^2+7)]),]
    }

    colnames(configs) <- c(rep(x = "p", times = max(set.of.I)^2), 
                           "alpha", "I", "N", "ET", "value", "PSe", "PSp", 
                           "PPPV", "PNPV")
    configs <- convert.config(algorithm="IA2", results=configs)
  }

  # find the optimal testing configuration, over all array sizes considered
  result.ET <- save.it[save.it[,(max(set.of.I)^2+5)]==min(save.it[,(max(set.of.I)^2+5)]), 
                       c(1:(max(set.of.I)^2+4),(max(set.of.I)^2+5),(max(set.of.I)^2+13):ncol(save.it))]
  result.MAR <- save.it[save.it[,(max(set.of.I)^2+6)]==min(save.it[,(max(set.of.I)^2+6)]), 
                        c(1:(max(set.of.I)^2+4),(max(set.of.I)^2+6),(max(set.of.I)^2+13):ncol(save.it))]
  result.GR1 <- save.it[save.it[,(max(set.of.I)^2+7)]==min(save.it[,(max(set.of.I)^2+7)]), 
                        c(1:(max(set.of.I)^2+4),(max(set.of.I)^2+7),(max(set.of.I)^2+13):ncol(save.it))]
  result.GR2 <- save.it[save.it[,(max(set.of.I)^2+8)]==min(save.it[,(max(set.of.I)^2+8)]), 
                        c(1:(max(set.of.I)^2+4),(max(set.of.I)^2+8),(max(set.of.I)^2+13):ncol(save.it))]
  result.GR3 <- save.it[save.it[,(max(set.of.I)^2+9)]==min(save.it[,(max(set.of.I)^2+9)]), 
                        c(1:(max(set.of.I)^2+4),(max(set.of.I)^2+9),(max(set.of.I)^2+13):ncol(save.it))]
  result.GR4 <- save.it[save.it[,(max(set.of.I)^2+10)]==min(save.it[,(max(set.of.I)^2+10)]), 
                        c(1:(max(set.of.I)^2+4),(max(set.of.I)^2+10),(max(set.of.I)^2+13):ncol(save.it))]
  result.GR5 <- save.it[save.it[,(max(set.of.I)^2+11)]==min(save.it[,(max(set.of.I)^2+11)]), 
                        c(1:(max(set.of.I)^2+4),(max(set.of.I)^2+11),(max(set.of.I)^2+13):ncol(save.it))]
  result.GR6 <- save.it[save.it[,(max(set.of.I)^2+12)]==min(save.it[,(max(set.of.I)^2+12)]), 
                        c(1:(max(set.of.I)^2+4),(max(set.of.I)^2+12),(max(set.of.I)^2+13):ncol(save.it))]

  p.ga.ET <- informativeArrayProb(prob.vec=(result.ET[1:max(set.of.I)^2])[!is.na(result.ET[1:max(set.of.I)^2])], 
                                    nr=result.ET[max(set.of.I)^2+2], nc=result.ET[max(set.of.I)^2+2], method="gd")
  if("MAR" %in% obj.fn){
    p.ga.MAR <- informativeArrayProb(prob.vec=(result.MAR[1:max(set.of.I)^2])[!is.na(result.MAR[1:max(set.of.I)^2])], 
                                       nr=result.MAR[max(set.of.I)^2+2], nc=result.MAR[max(set.of.I)^2+2], method="gd")
  } else{p.ga.MAR <- NA}
  if(is.null(dim(weights))){
    p.ga.GR1 <- NA
    p.ga.GR2 <- NA
    p.ga.GR3 <- NA
    p.ga.GR4 <- NA
    p.ga.GR5 <- NA
    p.ga.GR6 <- NA
  }  else{
    p.ga.GR1 <- informativeArrayProb(prob.vec=(result.GR1[1:max(set.of.I)^2])[!is.na(result.GR1[1:max(set.of.I)^2])], 
                                       nr=result.GR1[max(set.of.I)^2+2], nc=result.GR1[max(set.of.I)^2+2], method="gd")
    if(dim(weights)[1]>=2){
      p.ga.GR2 <- informativeArrayProb(prob.vec=(result.GR2[1:max(set.of.I)^2])[!is.na(result.GR2[1:max(set.of.I)^2])], 
                                         nr=result.GR2[max(set.of.I)^2+2], nc=result.GR2[max(set.of.I)^2+2], method="gd")
    } else{p.ga.GR2 <- NA}
    if(dim(weights)[1]>=3){
      p.ga.GR3 <- informativeArrayProb(prob.vec=(result.GR3[1:max(set.of.I)^2])[!is.na(result.GR3[1:max(set.of.I)^2])], 
                                         nr=result.GR3[max(set.of.I)^2+2], nc=result.GR3[max(set.of.I)^2+2], method="gd")
    } else{p.ga.GR3 <- NA}
    if(dim(weights)[1]>=4){
      p.ga.GR4 <- informativeArrayProb(prob.vec=(result.GR4[1:max(set.of.I)^2])[!is.na(result.GR4[1:max(set.of.I)^2])], 
                                         nr=result.GR4[max(set.of.I)^2+2], nc=result.GR4[max(set.of.I)^2+2], method="gd")
    } else{p.ga.GR4 <- NA}
    if(dim(weights)[1]>=5){
      p.ga.GR5 <- informativeArrayProb(prob.vec=(result.GR5[1:max(set.of.I)^2])[!is.na(result.GR5[1:max(set.of.I)^2])], 
                                         nr=result.GR5[max(set.of.I)^2+2], nc=result.GR5[max(set.of.I)^2+2], method="gd")
    } else{p.ga.GR5 <- NA}
    if(dim(weights)[1]>=6){
      p.ga.GR6 <- informativeArrayProb(prob.vec=(result.GR6[1:max(set.of.I)^2])[!is.na(result.GR6[1:max(set.of.I)^2])], 
                                         nr=result.GR6[max(set.of.I)^2+2], nc=result.GR6[max(set.of.I)^2+2], method="gd")
    } else{p.ga.GR6 <- NA}
  }

  # put accuracy measures in a matrix for easier display of results
  acc.ET <- matrix(data=result.ET[(max(set.of.I)^2+6:9)], nrow=1, ncol=4, 
                   dimnames=list(NULL, c("PSe", "PSp", "PPPV", "PNPV")))
  acc.MAR <- matrix(data=result.MAR[(max(set.of.I)^2+6:9)], nrow=1, ncol=4, 
                    dimnames=list(NULL, c("PSe", "PSp", "PPPV", "PNPV")))
  acc.GR1 <- matrix(data=result.GR1[(max(set.of.I)^2+6:9)], nrow=1, ncol=4, 
                    dimnames=list(NULL, c("PSe", "PSp", "PPPV", "PNPV")))
  acc.GR2 <- matrix(data=result.GR2[(max(set.of.I)^2+6:9)], nrow=1, ncol=4, 
                    dimnames=list(NULL, c("PSe", "PSp", "PPPV", "PNPV")))
  acc.GR3 <- matrix(data=result.GR3[(max(set.of.I)^2+6:9)], nrow=1, ncol=4, 
                    dimnames=list(NULL, c("PSe", "PSp", "PPPV", "PNPV")))
  acc.GR4 <- matrix(data=result.GR4[(max(set.of.I)^2+6:9)], nrow=1, ncol=4, 
                    dimnames=list(NULL, c("PSe", "PSp", "PPPV", "PNPV")))
  acc.GR5 <- matrix(data=result.GR5[(max(set.of.I)^2+6:9)], nrow=1, ncol=4, 
                    dimnames=list(NULL, c("PSe", "PSp", "PPPV", "PNPV")))
  acc.GR6 <- matrix(data=result.GR6[(max(set.of.I)^2+6:9)], nrow=1, ncol=4, 
                    dimnames=list(NULL, c("PSe", "PSp", "PPPV", "PNPV")))
  
  # create a list of results for each objective function
  opt.ET <- list("OTC"=list("Array.dim"=result.ET[(max(set.of.I)^2+2)], 
                            "Array.sz"=result.ET[(max(set.of.I)^2+3)]), 
                 "p.mat"=p.ga.ET, "ET"=result.ET[(max(set.of.I)^2+4)], 
                 "value"=result.ET[(max(set.of.I)^2+5)], "Accuracy"=acc.ET)
  opt.MAR <- list("OTC"=list("Array.dim"=result.MAR[(max(set.of.I)^2+2)], 
                             "Array.sz"=result.MAR[(max(set.of.I)^2+3)]), 
                  "p.mat"=p.ga.MAR, "ET"=result.MAR[(max(set.of.I)^2+4)], 
                  "value"=result.MAR[(max(set.of.I)^2+5)], "Accuracy"=acc.MAR)
  opt.GR1 <- list("OTC"=list("Array.dim"=result.GR1[(max(set.of.I)^2+2)], 
                             "Array.sz"=result.GR1[(max(set.of.I)^2+3)]), 
                  "p.mat"=p.ga.GR1, "ET"=result.GR1[(max(set.of.I)^2+4)], 
                  "value"=result.GR1[(max(set.of.I)^2+5)], "Accuracy"=acc.GR1)
  opt.GR2 <- list("OTC"=list("Array.dim"=result.GR2[(max(set.of.I)^2+2)], 
                             "Array.sz"=result.GR2[(max(set.of.I)^2+3)]), 
                  "p.mat"=p.ga.GR2, "ET"=result.GR2[(max(set.of.I)^2+4)], 
                  "value"=result.GR2[(max(set.of.I)^2+5)], "Accuracy"=acc.GR2)
  opt.GR3 <- list("OTC"=list("Array.dim"=result.GR3[(max(set.of.I)^2+2)], 
                             "Array.sz"=result.GR3[(max(set.of.I)^2+3)]), 
                  "p.mat"=p.ga.GR3, "ET"=result.GR3[(max(set.of.I)^2+4)], 
                  "value"=result.GR3[(max(set.of.I)^2+5)], "Accuracy"=acc.GR3)
  opt.GR4 <- list("OTC"=list("Array.dim"=result.GR4[(max(set.of.I)^2+2)], 
                             "Array.sz"=result.GR4[(max(set.of.I)^2+3)]), 
                  "p.mat"=p.ga.GR4, "ET"=result.GR4[(max(set.of.I)^2+4)], 
                  "value"=result.GR4[(max(set.of.I)^2+5)], "Accuracy"=acc.GR4)
  opt.GR5 <- list("OTC"=list("Array.dim"=result.GR5[(max(set.of.I)^2+2)], 
                             "Array.sz"=result.GR5[(max(set.of.I)^2+3)]), 
                  "p.mat"=p.ga.GR5, "ET"=result.GR5[(max(set.of.I)^2+4)], 
                  "value"=result.GR5[(max(set.of.I)^2+5)], "Accuracy"=acc.GR5)
  opt.GR6 <- list("OTC"=list("Array.dim"=result.GR6[(max(set.of.I)^2+2)], 
                             "Array.sz"=result.GR6[(max(set.of.I)^2+3)]), 
                  "p.mat"=p.ga.GR6, "ET"=result.GR6[(max(set.of.I)^2+4)], 
                  "value"=result.GR6[(max(set.of.I)^2+5)], "Accuracy"=acc.GR6)

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
  opt.all <- list("opt.ET"=opt.ET, "opt.MAR"=opt.MAR, "opt.GR1"=opt.GR1, 
                  "opt.GR2"=opt.GR2, "opt.GR3"=opt.GR3, "opt.GR4"=opt.GR4, 
                  "opt.GR5"=opt.GR5, "opt.GR6"=opt.GR6)
  # remove any objective functions not requested by the user
  opt.req <- Filter(function(x) !is.na(x$ET), opt.all)
  
  # print the time elapsed, if print.time == TRUE
  if(print.time){
    time.it(start.time)
  }

  inputs <- list("algorithm"="Informative array testing without master pooling", 
                 "prob"=list(p), "alpha"=alpha, "Se"=Se.display, "Sp"=Sp.display)
  res <- c(inputs, opt.req)
  res[["Configs"]] <- configs
  res
}

###################################################################
