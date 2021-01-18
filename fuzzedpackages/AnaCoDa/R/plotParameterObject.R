#' Plot Parameter 
#' 
#' @param x A parameter object
#' 
#' @param what Which aspect of the parameter to plot. Default value is
#' "Mutation".
#' 
#' @param samples Number of samples to plot using the posterior mean. Default
#' value is 100.
#'
#' @param mixture.name a vector with names/descriptions of the mixture distributions in the parameter object
#'
#' @param with.ci Plot with or without confidence intervals. Default value
#' is TRUE
#' 
#' @param ... Arguments to be passed to methods, such as graphical parameters.
#' 
#' @return This function has no return value.
#' 
#' @description \code{plot} graphs the mutation or selection parameter for a ROC or FONSE
#' parameter object for each mixture element.
#' 
#' @details Graphs are based off the last # samples for the posterior mean.
#' 
plot.Rcpp_ROCParameter <- function(x, what = "Mutation", samples = 100, mixture.name = NULL, with.ci = TRUE, ...)
{
  plotParameterObject(x, what = what, samples= samples, mixture.name=mixture.name, with.ci=with.ci, ...)
}


#' Plot Parameter 
#' 
#' @param x A parameter object
#' 
#' @param what Which aspect of the parameter to plot. Default value is
#' "Mutation".
#' 
#' @param samples Number of samples to plot using the posterior mean. Default
#' value is 100.
#'
#' @param mixture.name a vector with names/descriptions of the mixture distributions in the parameter object
#'
#' @param with.ci Plot with or without confidence intervals. Default value
#' is TRUE
#' 
#' @param ... Arguments to be passed to methods, such as graphical parameters.
#' 
#' @return This function has no return value.
#' 
#' @description \code{plot} graphs the mutation or selection parameter for a ROC or FONSE
#' parameter object for each mixture element.
#' 
#' @details Graphs are based off the last # samples for the posterior mean.
#' 
plot.Rcpp_FONSEParameter <- function(x, what = "Mutation", samples = 100, mixture.name = NULL, with.ci = TRUE, ...)
{
  plotParameterObject(x, what = what, samples=samples,mixture.name = mixture.name, with.ci=with.ci, ...)
}

plot.Rcpp_PAParameter <- function(x, what = "Mutation", samples = 100, mixture.name = NULL, with.ci = TRUE, ...)
{
  #to_plot = ifelse(what == "Alpha", "Mutation", "Selection")
  plotParameterObject(x, what = what, samples=samples,mixture.name = mixture.name, with.ci=with.ci, ...)
  #plotPA(x)
}

### NOT EXPOSED
plotParameterObject <- function(x, what = "Mutation", samples = 100, mixture.name = NULL, with.ci = TRUE, ...){
  numMixtures <- x$numMixtures
  means <- data.frame(matrix(0,ncol=numMixtures,nrow=40))
  sd.values <- data.frame(matrix(0,ncol=numMixtures*2,nrow=40))
  names.aa <- aminoAcids()
  paramType <- ifelse(what == "Mutation", 0, 1)
  #cat("ParamType: ", paramType, "\n")
  
  for (mixture in 1:numMixtures) {
    # get codon specific parameter
    count <- 1
    for (aa in names.aa) {
      if (aa == "M" || aa == "W" || aa == "X") next
      codons <- AAToCodon(aa, T)
      for (i in 1:length(codons))
      {
       means[count,mixture] <- x$getCodonSpecificPosteriorMean(mixture, samples,codons[i], paramType, TRUE, log_scale=FALSE)
       tmp <- x$getCodonSpecificQuantile(mixture,samples, codons[i], paramType, c(0.025, 0.975), TRUE, log_scale=FALSE)
        
        ## This approach to storing the quantiles may seem unconventional, but I actually found it to be the most straight forward approach
        ## for plotting later.
        sd.values[count,mixture] <- tmp[1]
        sd.values[count,mixture+numMixtures] <- tmp[2]
        count <- count + 1
      }
    }
  }
  ## Begin graphing
  mat <- matrix(rep(0,numMixtures*numMixtures),
                nrow = numMixtures, ncol = numMixtures, byrow = TRUE)
  count <- 1
  for(i in 1:numMixtures){
    for(j in 1:numMixtures){
      if(i<=j){
        mat[i,j] <-count
        count <- count + 1
      }
    }
  }
  nf <- layout(mat,widths=c(rep(5,numMixtures)),heights=c(rep(5,numMixtures)),respect=FALSE)
  par(mar=c(1,1,1,1))
  for(i in 1:numMixtures){
    for(j in 1:numMixtures){
      if(i==j)
      {
        plot(NULL, xlim=c(0,1), ylim=c(0,1), ylab="", xlab="",xaxt='n',yaxt='n',ann=FALSE)
        if(is.null(mixture.name)){
          text(x = 0.5, y = 0.5, paste0("Mixture\nElement",i), 
               cex = 1.6, col = "black")
        }else{
          text(x = 0.5, y = 0.5, mixture.name[i], 
               cex = 1.6, col = "black")
        }
      }
      else if (i < j){
        if(with.ci){
          plot(means[,j],means[,i],ann=FALSE,xlim=range(cbind(sd.values[,j],sd.values[,j+numMixtures])),ylim=range(cbind(sd.values[,i],sd.values[,i+numMixtures])))
          upper.panel.plot(means[,j],means[,i],sd.x=cbind(sd.values[,j],sd.values[,j+numMixtures]),sd.y=cbind(sd.values[,i],sd.values[,i+numMixtures]))
          #confidenceInterval.plot(x = means[,j],y = mean[,i], sd.x=sd.values[,j],sd.y=sd.values[,i])
        } else{
          plot(means[,j],means[,i],ann=FALSE,xlim=range(means[,j]),ylim=range(means[,i]))
          upper.panel.plot(means[,j],means[,i])
        }
      }
    }
  }
}




#TODO: should PA's ploting be here as well?

upper.panel.plot <- function(x, y, sd.x=NULL, sd.y=NULL, ...){
  abline(0, 1, col = "blue", lty = 2)
  points(x, y, ...)
  if(!is.null(sd.y)){
    y.up <- sd.y[,2]
    y.low <- sd.y[,1]
    epsilon <- range(x, na.rm = T) * 0.1
    segments(x, y.low, x, y.up, ...)
  }
  if(!is.null(sd.x)){
    x.up <- sd.x[,2]
    x.low <- sd.x[,1]
    epsilon <- range(y, na.rm = T) * 0.1
    segments(x.low, y, x.up, y, ...)
  }  
  
  lm.line <- lm(y~x, na.action = "na.exclude")
  abline(lm.line, col="blue", lwd = 2)
  
  R2 <- summary(lm.line)$r.squared
  
  b <- lm.line$coef[2]
  rho <- ifelse(b > 0, sqrt(R2), -sqrt(R2)) #make sure rho has correct sign
  

  if(!is.null(sd.x))
  {
    xlim <- range(sd.x, na.rm = T)
  }else{
    xlim <- range(x,na.rm = T)
  }
  if(!is.null(sd.y))
  {
    ylim <- range(sd.y, na.rm = T)
  }else{
    ylim <- range(y,na.rm=T)
  }
  
  width <- xlim[2] - xlim[1]
  height <- ylim[2] - ylim[1]
  
  std.error <- summary(lm.line)$coefficients[4]
  slope <- round(summary(lm.line)$coefficients[2], 3)
  intercept <- round(summary(lm.line)$coefficients[1], 3)
  t <- (slope - 1)/std.error
  
  if((t > qt(1-(0.05/2), lm.line$df.residual - 1))||(t < qt((0.05/2),lm.line$df.residual-1))){
    eq <- paste0("y = ", sprintf("%.3f", intercept), " + ", sprintf("%.3f", slope), "x *")
    text(xlim[1] + width * 0.01, ylim[2] - height * 0.2, eq, pos = 4, cex = 1.5)
  }else{
    eq <- paste0("y = ", sprintf("%.3f", intercept), " + ", sprintf("%.3f", slope), "x")
    text(xlim[1] + width * 0.01, ylim[2] - height * 0.2, eq, pos = 4, cex = 1.5)
  } 
  if(b > 0){
    text(xlim[2] - width * 0.04, ylim[1] + height * 0.05,
         parse(text = paste0("rho == ", sprintf("%.4f", rho))),
         pos = 2, cex = 1.5, font = 2)
  }else{
    text(xlim[2] - width * 0.04, ylim[2] - height * 0.05,
         parse(text = paste0("rho == ", sprintf("%.4f", rho))),
         pos = 2, cex = 1.5, font = 2)
  }
}


lower.panel.plot <- function(x, y, ...)
{
  
}


confidenceInterval.plot <- function(x, y, sd.x=NULL, sd.y=NULL, ...){
  points(x, y, ...)
  if(!is.null(sd.y)){
    y.up <- sd.y[,2]
    y.low <- sd.y[,1]
    epsilon <- range(x, na.rm = T) * 0.1
    segments(x, y.low, x, y.up, ...)
  }
  if(!is.null(sd.x)){
    x.up <- sd.x[,2]
    x.low <- sd.x[,1]
    epsilon <- range(y, na.rm = T) * 0.1
    segments(x.low, y, x.up, y, ...)
  }  

  lm.line <- lm(y~x, na.action = "na.exclude")

  b <- lm.line$coef[2]

  xlim <- range(x, na.rm = T)
  ylim <- range(y, na.rm = T)

  width <- xlim[2] - xlim[1]
  height <- ylim[2] - ylim[1]

  std.error <- summary(lm.line)$coefficients[4]
  slope <- round(summary(lm.line)$coefficients[2], 3)
  intercept <- round(summary(lm.line)$coefficients[1], 3)
  t <- (slope - 1)/std.error
}


plotPA <- function(parameter,genome,samples=100,mixture=1){
    #cat("hello")
    cat <- mixture
trace <- parameter$getTraceObject()
proposal <- FALSE
alphaList <- numeric (61)
lambdaPrimeList <- numeric (61)
waitingTimes <- numeric(61)
alpha.ci <- matrix(0, ncol=2, nrow=61)
lambdaPrime.ci <- matrix(0, ncol=2, nrow=61)
psiList <- numeric(length(genome))
ids <- numeric(length(genome))
codonList <- codons()
for (i in 1:61)
{
  codon <- codonList[i]
  alphaList[i] <- parameter$getCodonSpecificPosteriorMean(cat, samples * 0.5, codon, 0, FALSE,log_scale=F)
  alphaTrace <- trace$getCodonSpecificParameterTraceByMixtureElementForCodon(1, codon, 0, FALSE)
  alpha.ci[i,] <- quantile(alphaTrace[(samples * 0.5):samples], probs = c(0.025,0.975))
  
  
  lambdaPrimeList[i] <- parameter$getCodonSpecificPosteriorMean(cat, samples * 0.5, codon, 1, FALSE, log_scale=F)
  lambdaPrimeTrace <- trace$getCodonSpecificParameterTraceByMixtureElementForCodon(1, codon, 1, FALSE)
  lambdaPrime.ci[i,] <- quantile(lambdaPrimeTrace[(samples * 0.5):samples], probs = c(0.025,0.975))
  waitingTimes[i] <- alphaList[i] * lambdaPrimeList[i]
}

waitRates <- numeric(61)
for (i in 1:61) {
  waitRates[i] <- (1.0/waitingTimes[i])
}


for (geneIndex in 1:length(genome)) {
  psiList[geneIndex] <- parameter$getSynthesisRatePosteriorMeanForGene(samples * 0.5, geneIndex, 1)
}

for (i in 1:length(genome))
{
  g <- genome$getGeneByIndex(i, FALSE)
  ids[i] <- g$id
}

#Plot confidence intervals for alpha and lambda prime
plot(NULL, NULL, xlim=range(1:61, na.rm = T), ylim=range(alpha.ci), 
     main = "Confidence Intervals for Alpha Parameter", xlab = "Codons", 
     ylab = "Estimated values", axes=F) 
confidenceInterval.plot(x = 1:61, y = alphaList, sd.y = alpha.ci)
axis(2)
axis(1, tck = 0.02, labels = codonList[1:61], at=1:61, las=2, cex.axis=.6)

plot(NULL, NULL, xlim=range(1:61, na.rm = T), ylim=range(lambdaPrime.ci), 
     main = "Confidence Intervals for LambdaPrime Parameter", xlab = "Codons", 
     ylab = "Estimated values", axes=F) 
confidenceInterval.plot(x = 1:61, y = lambdaPrimeList, sd.y = lambdaPrime.ci)
axis(2)
axis(1, tck = 0.02, labels = codonList[1:61], at=1:61, las=2, cex.axis=.6)
#dev.off()
}

#' Plots ACF for codon specific parameter traces
#' 
#' @param parameter object of class Parameter
#' @param csp indicates which parameter to calculate the autocorrelation. Must be Mutation (the default, ROC, FONSE), Selection (ROC, FONSE), Alpha (PA, PANSE), LambdaPrime (PA, PANSE), NSERate (PA, PANSE)"
#' @param numMixtures indicates the number of CSP mixtures used
#' @param samples number of samples at the end of the trace used to calculate the acf
#' @param lag.max Maximum amount of lag to calculate acf. Default is 10*log10(N), where N i the number of observations.
#' @param plot logical. If TRUE (default) a plot of the acf is created
#' 
#' @description The function calculates and by defaults plots the acf and estimates the autocorrelation in the trace 
#' 
#' 
#' @seealso \code{\link{acfMCMC}}

acfCSP <- function(parameter, csp = "Mutation", numMixtures = 1, samples = NULL, lag.max = 40, plot=TRUE)
{
  if (csp == "Mutation" || csp == "Alpha")
  {
    paramType <- 0
  } else if (csp == "Selection" || csp == "LambdaPrime" || csp == "Lambda"){
    paramType <- 1
  } else if (csp == "NSERate"){
    paramType <- 2
  } else{
    stop("csp must take one of the following values: Mutation (ROC, FONSE), Selection (ROC, FONSE), Alpha (PA, PANSE), LambdaPrime (PA, PANSE), NSERate (PA, PANSE)")
  }
 
  
  acf.list <- list()
  names.aa <- aminoAcids()
  trace <- parameter$getTraceObject()
  if(is.null(samples))
  { 
    samples <- round(10*log10(length(trace))) 
  }

  ref.codon <- ifelse(csp %in% c("Selection","Mutation"),TRUE,FALSE)

  for (aa in names.aa)
  {
    if (aa == "X")
      next
     if ((aa == "M" || aa == "W") && ref.codon) ## If ROC or FONSE, skip amino acids without synonyms
      next
    codons <- AAToCodon(aa, ref.codon) ## If ROC or FONSE, skip reference codon
    codon.list <- list()
    for (i in 1:length(codons))
    {
      mix.list <- list()
      for (j in 1:numMixtures)
      {
        csp.trace <- trace$getCodonSpecificParameterTraceByMixtureElementForCodon(j, codons[i], paramType, ref.codon)
        csp.trace <- csp.trace[(length(csp.trace)-samples):length(csp.trace)]
        csp.acf <- acf(x = csp.trace, lag.max = lag.max, plot = FALSE)
        mix.list[[j]] <- csp.acf
        if (plot)
        {
          header <- paste(csp, aa, codons[i], "Mixture:", j, sep = " ")
          plot(x = csp.acf, xlab = "Lag time", ylab = "Autocorrelation", main = header)
        }
      }
      codon.list[[codons[i]]] <-mix.list
    }
    acf.list[[aa]] <- codon.list
  }
  return(acf.list)
}


# NOT EXPOSED
# plots to data.frames / matrices with CSP estimates.
plotCSPdf <- function(df1, df2, xlab = "", ylab = "", main = "")
{
  fill.in <- data.frame(Codon=as.character(codons()))
  df1 <- merge(x = df1, y = fill.in, by = "Codon", all = T)
  df2 <- merge(x = df2, y = fill.in, by = "Codon", all = T)
  
  for(i in 1:64){
    df1$AA[i] <- codonToAA(df1$Codon[i])
    df2$AA[i] <- codonToAA(df2$Codon[i])
  }
  df1[is.na(df1)] <- 0
  df1 <- df1[order(df1$AA), ]
  df2[is.na(df2)] <- 0
  df2 <- df2[order(df2$AA), ]
  
    
  aas <- aminoAcids()
  n.aa <- length(aas)
  kt <- rep(0, n.aa)
  dn <- rep(0, n.aa)
  for(j in 1:n.aa)
  {
    aa <- aas[j]
    if(aa == "W" || aa == "M" || aa == "X") next
    aa.pos <- which(df1$AA == aa)
    df1[aa.pos,4] <- df1[aa.pos,4] - mean(df1[aa.pos,3])
    df2[aa.pos,4] <- df2[aa.pos,4] - mean(df2[aa.pos,3])
    
    df1[aa.pos,5] <- df1[aa.pos,5] - mean(df1[aa.pos,3])
    df2[aa.pos,5] <- df2[aa.pos,5] - mean(df2[aa.pos,3])

    df1[aa.pos,3] <- df1[aa.pos,3] - mean(df1[aa.pos,3])
    df2[aa.pos,3] <- df2[aa.pos,3] - mean(df2[aa.pos,3])
  }
  
  xlim <- range(df1[, 4:5])
  ylim <- range(df2[, 4:5])
  
  plot(NULL, NULL, axes=F, xlab = "", ylab = "", xlim = xlim, ylim = ylim)
  points(df1[, 3], df2[,3], pch = 19, col = rgb(0,0,0,0.7))
  
  segments(x0 = df1[,3], y0 = df2[,4], x1 = df1[,3], y1 = df2[,5], lwd=2, col=adjustcolor("black", 0.5))
  segments(x0 = df1[,4], y0 = df2[,3], x1 = df1[,5], y1 = df2[,3], lwd=2, col=adjustcolor("black", 0.5))
  
  type2.reg <- lmodel2::lmodel2(df2[, 3] ~ df1[,3], nperm = 10, range.y = "interval", range.x = "interval")
  intercept <- type2.reg$regression.results[4, 2]
  slope <- type2.reg$regression.results[4, 3]
  r <- type2.reg$r
  
  if(slope > 0){
    text(x = xlim[1] + abs(xlim[1])*0.1, y = ylim[2] - abs(ylim[2])*0.2, adj = 0,
         labels = bquote("y = " ~.(round(intercept, 2)) ~ " + " ~.(round(slope, 2)) ~ "x"), cex = 0.75)
    text(x = xlim[1] + abs(xlim[1])*0.1, y = (ylim[2] - abs(ylim[2])*0.2) - 0.4, adj = 0,
         labels = bquote(rho ~ " = " ~.(round(r, 2))), cex = 0.75)
  }else{
    text(x = xlim[1] + abs(xlim[1])*0.2, y = ylim[1] + abs(ylim[1])*0.2 + 0.4, adj = 0,
         labels = bquote("y = " ~.(round(intercept, 2)) ~ " + " ~.(round(slope, 2)) ~ "x"), cex = 0.75)
    text(x = xlim[1] + abs(xlim[1])*0.2, y = (ylim[1] + abs(ylim[1])*0.2), adj = 0,
         labels = bquote(rho ~ " = " ~.(round(r, 2))), cex = 0.75)  
  }
  
  abline(a = intercept, b = slope, lty = 1, lwd = 2, col = adjustcolor("red", 0.7))
  abline(a = 0, b = 1, lty = 2, lwd = 1)
  abline(h = 0, lty = 3, lwd = 1)
  abline(v = 0, lty = 3, lwd = 1)
  
  title(main = main)
  
  mtext(text = xlab, side = 1, line = 2.5, font = 2, cex=1.5)
  mtext(text = ylab, side = 2, line = 2.5, font = 2, cex=1.5)
  
  axis(side = 1, tick = T, font.axis = 2, lwd=2, las=1)
  axis(side = 2, tick = T, font.axis = 2, lwd=2, las=1)
}
