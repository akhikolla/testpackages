bias <- function(workingDir) {
  setwd(workingDir)
  iter <- NULL
  n <- 0
  
  grid <- expand.grid(ML=FALSE,
                      #iter=1:20000,
                      iter=1:2,
                      n = c(10,100,1000),
                      #rho = c(-0.99,seq(-0.95,0.95,by=0.05), 0.99),
                      rho = c(-0.99,seq(-0.95,0.95,by=0.2), 0.99),
                      fast=TRUE)
  grid$reset <- TRUE
  grid <- subset(grid, ! ( (iter > 100) & (n == 1000) | ( (iter > 1000) & (n==100) ) ) )
  bias <- wCorrSim(n=grid$n, rho=grid$rho, ML=grid$ML, fast=grid$fast, reset=grid$reset, usew=FALSE)
  
  save(bias, file="bias.RData")
  
  
  bias$rmse <- sqrt( (bias$est - bias$rho)^2 )
  bias$bias <- bias$est - bias$rho

  aggbiasA <- aggregate(bias ~ n + rho + type, data=bias, FUN=mean, na.rm=TRUE)
  aggbiasB <- aggregate(rmse ~ n + rho + type, data=bias, FUN=mean, na.rm=TRUE)
  aggbias <- merge(aggbiasA, aggbiasB, by=c("n", "rho", "type"), all=TRUE)
  names(aggbias)[names(aggbias) == "bias"] <- "bias.mean"
  names(aggbias)[names(aggbias) == "rmse"] <- "rmse.mean"
  aggbias2 <- aggregate(rmse  ~ n+type, data=bias, FUN=mean, na.rm=TRUE)
  names(aggbias2)[names(aggbias2) == "rmse"] <- "rmse.mean"
  
  save(aggbias, aggbias2, file="aggbias.RData")
}

fast <- function(workingDir){
  setwd(workingDir)

  grid <- expand.grid(fast=c(TRUE,FALSE),
                      iter=1:10,
                      n = c(10,100,1000),
                      rho = c(-0.99,seq(-0.95,0.95,by=0.05), 0.99),
                      ML=FALSE)
  grid$reset <- grid$fast
  fast <- wCorrSim(n=grid$n, rho=grid$rho, ML=grid$ML, fast=grid$fast, reset=grid$reset, usew=FALSE, outstr="fast")
  
  save(fast, file="fast.RData")
  
  if(FALSE) {
    x <- rnorm(10)
    y <- rnorm(10)+x/2
    afast <- weightedCorr(x,y, method="Pearson", weights=rep(1,10), fast=TRUE, ML=FALSE)
    aslow <- weightedCorr(x,y, method="Pearson", weights=rep(1,10), fast=FALSE, ML=FALSE)
    c(afast-aslow,afast,aslow)
    afast2 <- weightedCorr(x,y, method="Pearson", weights=rep(1,10), fast=TRUE, ML=FALSE)
    c(afast-afast2,afast,afast2)
    afast <- weightedCorr(x,y, method="Spearman", weights=rep(1,10), fast=TRUE, ML=FALSE)
    aslow <- weightedCorr(x,y, method="Spearman", weights=rep(1,10), fast=FALSE, ML=FALSE)
    c( (afast-aslow)*10^16,afast,aslow)
    M <- 0 + (x > 0)
    P <- 9 + (y > 0)
    afast <- weightedCorr(x,y, method="polychoric", weights=rep(1,10), fast=TRUE, ML=FALSE)
    aslow <- weightedCorr(x,y, method="polychoric", weights=rep(1,10), fast=FALSE, ML=FALSE)
    c( (afast-aslow),afast,aslow)
  }
  
  fast$i <- rep(1:(nrow(fast)/2),each=2)
  mfast <- merge(subset(fast,fast),
                 subset(fast,!fast, c("i", "est")),
                 by="i",
                 suffixes=c(".fast",".slow"))
  mfast$fast <- NULL
  mfast$absdrho <- pmax(abs(mfast$est.fast - mfast$est.slow), 1E-16)
  aggfast <- aggregate(absdrho ~ n + rho + type, data=mfast, FUN=mean, na.rm=TRUE)
  names(aggfast)[names(aggfast) == "absdrho"] <- "absdrho.mean"
  save(aggfast, file="aggfast.RData")
  
}

ML <- function(workingDir) {
  type <- NULL
  grid <- expand.grid(ML=c(TRUE,FALSE),
                      iter=1:500,
                      n = c(10,100,1000),
                      rho = c(-0.99,seq(-0.95,0.95,by=0.05), 0.99),
                      fast=TRUE)
  grid$reset <- grid$ML
  ML <- wCorrSim(n=grid$n, rho=grid$rho, ML=grid$ML, fast=grid$fast, reset=grid$reset, usew=FALSE)
  
  save(ML, file="ML.RData")
  

  ml <- subset(ML, type %in% c("Polychoric", "Polyserial"))
  ml$rmse <- (ml$est - ml$rho)^2
  
  aggml <- aggregate(rmse ~ n + rho + type + ML, data=ml, FUN=mean, na.rm=TRUE)
  names(aggml)[names(aggml) == "rmse"] <- "rmse.mean"

  aggml$rmse.mean <- sqrt(aggml$rmse.mean)
  aggml$ml <- ifelse(aggml$ML==TRUE, "ML=TRUE", "ML=FALSE")
  aggml$nt <- factor(paste("n=",aggml$n))
  
  
  ml$i <- rep(1:(nrow(ml)/2),each=2)
  mml <- merge(subset(ml,ML),
               subset(ml,!ML, c("i", "est")),
               by="i",
               suffixes=c(".ml",".nonml"))
  mml$absd <- abs(mml$est.ml - mml$est.nonml)
  aggt1_0 <- aggregate(absd ~ type + n + ML, data=subset(mml, type=="Polychoric"), FUN=mean, na.rm=TRUE)
  names(aggt1_0)[names(aggt1_0) == "absd"] <- "absd.mean"
  aggt1_0$ML <- NULL
  
  aggt1 <- aggregate(rmse ~ type + n + ML, data=subset(ml, type=="Polychoric"), FUN=mean, na.rm=TRUE)
  names(aggt1)[names(aggt1) == "rmse"] <- "rmse.mean"
  
  
  
  aggt2_0 <- aggregate(absd ~ type + n + ML, data=subset(mml, type=="Polyserial"), FUN=mean, na.rm=TRUE)
  names(aggt2_0)[names(aggt2_0) == "absd"] <- "absd.mean"
  aggt2_0$ML <- NULL
  
  aggt2 <- aggregate(rmse ~ type + n + ML, data=subset(ml, type=="Polyserial"), FUN=mean, na.rm=TRUE)
  names(aggt2)[names(aggt2) == "rmse"] <- "rmse.mean"
  aggt2$rmse.mean <- sqrt(aggt2$rmse.mean)
  
  save(aggml, aggt1_0, aggt1, aggt2_0, aggt2, file="aggML.RData")
  
}

ntime <- function(workingDir) {
  setwd(workingDir)
  ntime1 <- ntime2 <- NULL
  
  grid <- expand.grid(ML=FALSE,
                      iter=1:5,
                      n = round(10^seq(1,6,by=0.25)),
                      rho = c(-0.99,seq(-0.95,0.95,by=0.05), 0.99),
                      fast=TRUE)
  ntime <- wCorrSim(n=grid$n, rho=grid$rho, ML=grid$ML, fast=grid$fast, reset=TRUE, usew=FALSE, outstr="ntime")
  
  save(ntime, file="ntime.RData")
  
  ######
  
  setwd(workingDir)
  grid <- expand.grid(ML=FALSE,
                      iter=1:5,
                      n = round(10^seq(6.5,7,by=0.5)),
                      rho = c(-0.99,seq(-0.95,0.95,by=0.05), 0.99),
                      fast=TRUE)
  ntime2 <- wCorrSim(n=grid$n, rho=grid$rho, ML=grid$ML, fast=grid$fast, reset=TRUE, usew=FALSE, outstr="ntime2")
  save(ntime2, file="ntime2.RData")
  
  ###### 
  setwd(workingDir)
  load("ntime1.RData")
  load("ntime2.RData")
  ntime <- rbind(ntime1, ntime2)
  save(ntime, file="ntime.RData")
  
  aggTime <- aggregate(t ~ n + type, data=ntime, FUN=mean, na.rm=TRUE)
  names(aggTime)[names(aggTime) == "t"] <- "t.mean"
  aggTime$t.mean <- ifelse(aggTime$t.mean==0, 0.001,aggTime$t.mean)
  save(aggTime, file="aggTime.RData")
  
}

speed <- function(workingDir) {
  setwd(workingDir)

  grid1 <- expand.grid(fast=c(TRUE,FALSE),
                       ML=c(TRUE,FALSE),
                       iter=80,
                       n = round(10^seq(1,4.75,by=0.25)),
                       rho = c(-0.99,seq(-0.95,0.95,by=0.05), 0.99))
  
  grid2 <- expand.grid(fast=c(TRUE,FALSE),
                       ML=c(TRUE,FALSE),
                       iter=20,
                       n = round(10^seq(5,7,by=0.5)),
                       rho = c(-0.99,seq(-0.95,0.95,by=0.05), 0.99))
  grid <- rbind(grid1, grid2)
  
  grid$reset <- (grid$ML) & (grid$fast)
  speed1 <- wCorrSim(n=grid$n, rho=grid$rho, ML=grid$ML, fast=grid$fast, reset=grid$reset, usew=FALSE,outstr="speed")
  
  save(speed, file="speed.RData")
  
  
  ##############
  setwd(workingDir)

  grid1 <- expand.grid(fast=c(TRUE,FALSE),
                       ML=c(TRUE,FALSE),
                       iter=80,
                       n = round(10^seq(1,4.75,by=0.25)),
                       rho = c(-0.99,seq(-0.95,0.95,by=0.05), 0.99))
  grid <- rbind(grid1)
  
  grid$reset <- (grid$ML) & (grid$fast)
  speed1 <- wCorrSim(n=grid$n, rho=grid$rho, ML=grid$ML, fast=grid$fast, reset=grid$reset, usew=FALSE,outstr="speed1")
  
  save(speed1, file="speed1.RData")
  
  ###################
  require(wCorr)
  setwd(workingDir)

  grid2 <- expand.grid(fast=c(TRUE,FALSE),
                       ML=c(TRUE,FALSE),
                       iter=20,
                       n = round(10^seq(5,6,by=0.5)),
                       rho = c(-0.99,seq(-0.95,0.95,by=0.05), 0.99))
  grid <- rbind(grid2)
  
  grid$reset <- (grid$ML) & (grid$fast)
  speed2 <- wCorrSim(n=grid$n, rho=grid$rho, ML=grid$ML, fast=grid$fast, reset=grid$reset, usew=FALSE,outstr="speed2")
  save(speed2, file="speed2.RData")
  
  
  ###################
  require(wCorr)
  setwd(workingDir)

  grid2 <- expand.grid(fast=c(TRUE,FALSE),
                       ML=c(TRUE,FALSE),
                       iter=20,
                       n = round(10^seq(6.5,6.5,by=0.5)),
                       rho = c(-0.99,seq(-0.95,0.95,by=0.05), 0.99))
  grid <- rbind(grid2)
  
  grid$reset <- (grid$ML) & (grid$fast)
  speed3 <- wCorrSim(n=grid$n, rho=grid$rho, ML=grid$ML, fast=grid$fast, reset=grid$reset, usew=FALSE,outstr="speed3")
  save(speed3, file="speed3.RData")
  
  
  ###################
  require(wCorr)
  setwd(workingDir)

  grid2 <- expand.grid(fast=c(TRUE,FALSE),
                       ML=c(TRUE,FALSE),
                       iter=10,
                       n = round(10^seq(7,7,by=0.5)),
                       rho = c(-0.99,seq(-0.95,0.95,by=0.05), 0.99))
  grid <- rbind(grid2)
  
  grid$reset <- (grid$ML) & (grid$fast)
  speed4 <- wCorrSim(n=grid$n, rho=grid$rho, ML=grid$ML, fast=grid$fast, reset=grid$reset, usew=FALSE,outstr="speed4")
  save(speed4, file="speed4.RData")
  
  ####
  setwd(workingDir)
  load("speed1.RData")
  load("speed2.RData")
  speed <- rbind(speed1, speed2)
  load("speed3.RData")
  speed <- rbind(speed, speed3)
  load("speed4.RData")
  speed <- rbind(speed, speed4)
  save(speed, file="speed.RData")
  
  speed$class <- ifelse(speed$ML, "ML=T,", "ML=F,")
  speed$class <- paste0(speed$class, ifelse(speed$fast, "fast=T", "fast=F"))
  speed$t <- pmax(speed$t, 0.001)
  aggSpeed <- aggregate(t ~ n + type + class, data=speed, FUN=mean, na.rm=TRUE)
  names(aggSpeed)[names(aggSpeed) == "t"] <- "t.mean"
  save(aggSpeed, file="aggSpeed.RData")
  
}

spearmanSpeed <- function(workingDir) {
  ###################
  require(wCorr)
  setwd(workingDir)

  grid <- expand.grid(usew=c(FALSE,TRUE),
                      iter=1:2,
                      n = c(10,100,1000, 10000),
                      rho = c(-0.99,seq(-0.95,0,by=0.05)))
  
  grid$reset <- !grid$usew
  spear <- spearmanSim(n=grid$n, rho=grid$rho, usew=grid$usew, outstr="spear")
  save(spear, file="spear.RData")
  
  load("../wgtvn.RData")
  
  spear$rho <- spear$SpearmanOrig
  spear$SpearmanOrig <- NULL
  spear$N <- NULL
  wgtvn <- wgtvn[wgtvn$type!= "Spearman",]
  
  wgt <- rbind(wgtvn, spear)
  wgt$mserho <- (wgt$est - wgt$rho)^2
  
  aggWgtvn <- aggregate(mserho ~ n + usew + type, data=wgt, FUN=mean, na.rm=TRUE)
  names(aggWgtvn)[names(aggWgtvn) == "mserho"] <- "mserho.mean"
  aggWgtvn$rmserho <- sqrt(aggWgtvn$mserho)
  aggWgtvn$weight <- ifelse(aggWgtvn$usew, "Weighted", "Unweighted")
  
  save(aggWgtvn, file="aggWgtvn.RData")
}
wgtvrho <- function(workingDir) {
  setwd(workingDir)
  grid <- expand.grid(usew=c(FALSE,TRUE),
                      iter=1:100,
                      n = c(10,100,1000),
                      rho = c(-0.99,seq(-0.95,0.95,by=0.05), 0.99))
  grid$reset <- !grid$usew
  wgtvrho <- wCorrSim(n=grid$n, rho=grid$rho, ML=FALSE, fast=TRUE, reset=TRUE, usew=grid$usew)
  
  save(wgtvrho, file="wgtvrho.RData")
  
  wgt <- wgtvrho
  wgt$absdrho <- abs(wgt$est - wgt$rho)
  
  aggWgtvrho <- aggregate(absdrho ~ rho + usew + type, data=wgt, FUN=mean, na.rm=TRUE)
  names(aggWgtvrho)[names(aggWgtvrho) == "absdrho"] <- "absdrho.mean"
  aggWgtvrho$weight <- ifelse(aggWgtvrho$usew, "Weighted", "Unweighted")
  save(aggWgtvrho, file="aggWgtvrho.RData")
}


wgtvn <- function(workingDir) {
  setwd(workingDir)

  grid <- expand.grid(usew=c(FALSE,TRUE),
                      iter=1:20,
                      n = c(10,100,1000,10000),
                      rho = c(-0.99,seq(-0.95,0,by=0.05)))
  grid$reset <- !grid$usew
  
  wgtvn <- wCorrSim(n=grid$n, rho=grid$rho, ML=FALSE, fast=TRUE, reset=TRUE, usew=grid$usew)
  save(wgtvn, file="wgtvn.RData")
}

# createRDA <- function() {
#   devtools::use_data(aggfast,aggml, aggt1, aggt1_0, aggt2, aggt2_0, aggSpeed, aggbias, aggbias2,
#                      aggTime, aggWgtvrho, aggWgtvn,speed, internal = TRUE, overwrite = TRUE)
# }
