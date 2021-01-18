#' Plot Model Object
#' 
#' @param x An Rcpp model object initialized with \code{initializeModelObject}.
#' @param genome An Rcpp genome object initialized with \code{initializeGenomeObject}.
#' @param samples The number of samples in the trace
#' @param mixture The mixture for which to graph values.
#' @param simulated A boolean value that determines whether to use the simulated genome.
#' @param ... Optional, additional arguments.
#' For this function, a possible title for the plot in the form of a list if set with "main".
#'  
#' @return This function has no return value.
#'  
#' @description Plots traces from the model object such as synthesis rates for each gene.
#' Will work regardless of whether or not expression/synthesis rate levels are being
#' estimated. If you wish to plot observed/empirical values, these values MUST be set
#' using the initial.expression.values parameter found in initializeParameterObject.
#' Otherwise, the expression values plotted will just be SCUO values estimated upon
#' initialization of the Parameter object.

plot.Rcpp_ROCModel <- function(x, genome = NULL, samples = 100, mixture = 1, 
                               simulated = FALSE, ...)
{
  model <- x
  opar <- par(no.readonly = T) 
  
  input_list <- as.list(list(...))
  
  if("main" %in% names(input_list)){
    main <- input_list$main
    input_list$main <- NULL
  }else{
    main <- ""
  }
  
  mat <- matrix(c(rep(1, 4), 2:21, rep(22, 4)),
                nrow = 7, ncol = 4, byrow = TRUE)
  mat <- cbind(rep(23, 7), mat, rep(24, 7))
  nf <- layout(mat, c(3, rep(8, 4), 2), c(3, 8, 8, 8, 8, 8, 3), respect = FALSE)
  ### Plot title.
  par(mar = c(0, 0, 0, 0))
  plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE)
  text(0.5, 0.6, main)
  text(0.5, 0.4, date(), cex = 0.6)
  
  num.genes <- length(genome)
  parameter <- model$getParameter()

  mixtureAssignment <- unlist(lapply(1:num.genes,  function(geneIndex){parameter$getEstimatedMixtureAssignmentForGene(samples, geneIndex)}))
  genes.in.mixture <- which(mixtureAssignment == mixture)
  expressionCategory <- parameter$getSynthesisRateCategoryForMixture(mixture)
  
  # need expression values to know range
  num.genes <- length(genes.in.mixture)
  expressionValues <- unlist(lapply(genes.in.mixture, function(geneIndex){
    parameter$getSynthesisRatePosteriorMeanForGene(samples, geneIndex, FALSE)
    }))  
  expressionValues <- log10(expressionValues)
  genome <- genome$getGenomeForGeneIndices(genes.in.mixture, simulated)
  
  names.aa <- aminoAcids()
  for(aa in names.aa)
  {
    if(aa == "M" || aa == "W" || aa == "X") next
    codon.probability <- calculateProbabilityVector(parameter,model,expressionValues,mixture,samples,aa,model.type="ROC")
    xlimit <- plotSinglePanel(parameter, model, genome, expressionValues, samples, mixture, aa,codon.probability = codon.probability)
    box()
    main.aa <- aa #TODO map to three letter code
    text(mean(xlimit), 1, main.aa, cex = 1.5)
    if(aa %in% c("A", "F", "K", "Q", "V")){
      axis(2, las=1)
    }
    if(aa %in% c("T", "V", "Y", "Z")){
      axis(1)
    }
    if(aa %in% c("A", "C", "D", "E")){
      axis(3)
    }
    if(aa %in% c("E", "I", "P", "T")){
      axis(4, las=1)
    }
    axis(1, tck = 0.02, labels = FALSE)
    axis(2, tck = 0.02, labels = FALSE)
    axis(3, tck = 0.02, labels = FALSE)
    axis(4, tck = 0.02, labels = FALSE)    
  }
  
  ## adding a histogram of phi values to plot
  hist.values <- hist(expressionValues, plot=FALSE, nclass=30)
  plot(hist.values, axes = FALSE, main = "", xlab = "", ylab = "")
  axis(1)
  axis(4, las=1)
  
  ### Plot xlab.
  plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE)
  text(0.5, 0.2, expression("log"[10]~"(Protein Synthesis Rate"~phi~")"))  
  #text(0.5, 0.5, "Production Rate (log10)")
  
  ### Plot ylab.
  plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE)
  text(0.5, 0.5, "Propotion", srt = 90)
  
  par(opar)
}

#' Plot Model Object
#' 
#' @param x An Rcpp model object initialized with \code{initializeModelObject}.
#'
#' @param genome An Rcpp genome object initialized with \code{initializeGenomeObject}.
#'
#' @param samples The number of samples in the trace
#'
#' @param mixture The mixture for which to graph values.
#'
#' @param simulated A boolean value that determines whether to use the simulated genome.
#'
#' @param codon.window A boolean value that determines the codon window to use for calculating codon frequencies. If NULL (the default), use complete sequences.
#'
#' @param ... Optional, additional arguments.
#' For this function, a possible title for the plot in the form of a list if set with "main".
#'  
#' @return This function has no return value.
#'
#' @description Plots traces from the model object such as synthesis rates for each gene.
#' Will work regardless of whether or not expression/synthesis rate levels are being
#' estimated. If you wish to plot observed/empirical values, these values MUST be set
#' using the initial.expression.values parameter found in initializeParameterObject.
#' Otherwise, the expression values plotted will just be SCUO values estimated upon
#' initialization of the Parameter object.
#'
plot.Rcpp_FONSEModel <- function(x, genome, samples = 100, mixture = 1, 
                               simulated = FALSE, codon.window = NULL,...)
{
  model <- x
  opar <- par(no.readonly = T) 
  
  input_list <- as.list(list(...))
  
  if("main" %in% names(input_list)){
    main <- input_list$main
    input_list$main <- NULL
  }else{
    main <- ""
  }
  
  mat <- matrix(c(rep(1, 4), 2:21, rep(22, 4)),
                nrow = 7, ncol = 4, byrow = TRUE)
  mat <- cbind(rep(23, 7), mat, rep(24, 7))
  nf <- layout(mat, c(3, rep(8, 4), 2), c(3, 8, 8, 8, 8, 8, 3), respect = FALSE)
  ### Plot title.
  par(mar = c(0, 0, 0, 0))
  plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE)
  text(0.5, 0.6, main)
  text(0.5, 0.4, date(), cex = 0.6)
  
  num.genes <- length(genome)
  parameter <- model$getParameter()
  
  mixtureAssignment <- unlist(lapply(1:num.genes,  function(geneIndex){parameter$getEstimatedMixtureAssignmentForGene(samples, geneIndex)}))
  genes.in.mixture <- which(mixtureAssignment == mixture)
  expressionCategory <- parameter$getSynthesisRateCategoryForMixture(mixture)

  
  # need expression values to know range
  num.genes <- length(genes.in.mixture)
  expressionValues <- unlist(lapply(genes.in.mixture, function(geneIndex){
      parameter$getSynthesisRatePosteriorMeanForGene(samples, geneIndex, FALSE)
    }))  

  expressionValues <- log10(expressionValues)
  genome <- genome$getGenomeForGeneIndices(genes.in.mixture, simulated)
  genes <- genome$getGenes(simulated)
  genome$clear()
  if (is.null(codon.window))
  {
    codon.window <- seq(1,100000)
  } else if (length(codon.window) == 2)
  {
    codon.window <- seq(codon.window[1],codon.window[2])
  }
  for (i in 1:length(genes))
  {
    dna <- genes[[i]]$seq
    start <- seq(1, nchar(dna), 3)
    stop  <- pmin(start + 2, nchar(dna))
    codons <- substring(dna,start,stop)
    codons <- codons[codon.window]
    codons <- codons[which(is.na(codons) == F)]
    dna <- paste(codons,collapse='')
    genes[[i]]$seq <- dna
    genome$addGene(genes[[i]],simulated)
  }
  names.aa <- aminoAcids()
  for(aa in names.aa)
  {
    if(aa == "M" || aa == "W" || aa == "X") next
    codon.probability <- calculateProbabilityVector(parameter,model,expressionValues,mixture,samples,aa,model.type="FONSE",codon.window = codon.window)
    xlimit <- plotSinglePanel(parameter, model, genome, expressionValues, samples, mixture, aa,codon.probability = codon.probability)
    box()
    main.aa <- aa #TODO map to three letter code
    text(mean(xlimit), 1, main.aa, cex = 1.5)
    if(aa %in% c("A", "F", "K", "Q", "V")){
      axis(2, las=1)
    }
    if(aa %in% c("T", "V", "Y", "Z")){
      axis(1)
    }
    if(aa %in% c("A", "C", "D", "E")){
      axis(3)
    }
    if(aa %in% c("E", "I", "P", "T")){
      axis(4, las=1)
    }
    axis(1, tck = 0.02, labels = FALSE)
    axis(2, tck = 0.02, labels = FALSE)
    axis(3, tck = 0.02, labels = FALSE)
    axis(4, tck = 0.02, labels = FALSE)    
  }
  
  ## adding a histogram of phi values to plot
  hist.values <- hist(expressionValues, plot=FALSE, nclass=30)
  plot(hist.values, axes = FALSE, main = "", xlab = "", ylab = "")
  axis(1)
  axis(4, las=1)
  
  ### Plot xlab.
  plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE)
  text(0.5, 0.2, expression("log"[10]~"(Protein Synthesis Rate"~phi~")"))  
  #text(0.5, 0.5, "Production Rate (log10)")
  
  ### Plot ylab.
  plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE)
  text(0.5, 0.5, "Propotion", srt = 90)
  
  par(opar)
}

calculateProbabilityVector <- function(parameter,model,expressionValues,mixture,samples,aa,model.type="ROC",codon.window = c(1,300))
{
  codons <- AAToCodon(aa, T)
  
  # get codon specific parameter
  selection <- vector("numeric", length(codons))
  mutation <- vector("numeric", length(codons))
  for (i in 1:length(codons))
  {
    selection[i] <- parameter$getCodonSpecificPosteriorMean(mixture, samples, codons[i], 1, T,log_scale = F)
    mutation[i] <- parameter$getCodonSpecificPosteriorMean(mixture, samples, codons[i], 0, T,log_scale = F)
  }
  
  # calculate codon probabilities with respect to phi
  expression.range <- range(expressionValues)
  phis <- seq(from = expression.range[1], to = expression.range[2], by = 0.01)
  if (model.type == "ROC")
  {
    codonProbability <- lapply(10^phis,  
                               function(phi){
                                 model$CalculateProbabilitiesForCodons(mutation, selection, phi)
                               })
    codonProbability <- do.call("rbind",codonProbability)
  } else if (model.type == "FONSE")
  {
    codonProbability <- vector(mode = "list", length = length(codon.window))
    for (i in 1:length(codon.window))
    {
      codonProbability.tmp <- lapply(10^phis,
                               function(phi)
                                 {
                                 model$CalculateProbabilitiesForCodons(mutation, selection, phi,codon.window[i])
                               })
      codonProbability[[i]] <- do.call("rbind",codonProbability.tmp)
    }
    codonProbability <- Reduce("+",codonProbability)/length(codon.window)
  }
  return(codonProbability)
}



# NOT EXPOSED
plotSinglePanel <- function(parameter, model, genome, expressionValues, samples, mixture, aa,codon.probability)
{
  codons <- AAToCodon(aa, T)
  # 
  # get codon specific parameter
  selection <- vector("numeric", length(codons))
  mutation <- vector("numeric", length(codons))
  for (i in 1:length(codons))
  {
    selection[i] <- parameter$getCodonSpecificPosteriorMean(mixture, samples, codons[i], 1, T, log_scale = F)
    mutation[i] <- parameter$getCodonSpecificPosteriorMean(mixture, samples, codons[i], 0, T, log_scale = F)
  }
  # 
  expression.range <- range(expressionValues)
  phis <- seq(from = expression.range[1], to = expression.range[2], by = 0.01)
  # codonProbability <- lapply(10^phis,  
  #                            function(phi){
  #                              model$CalculateProbabilitiesForCodons(mutation, selection, phi)
  #                            })
  codonProbability <- codon.probability
  #get codon counts
  codons <- AAToCodon(aa, F)
  codonCounts <- vector("list", length(codons))
  for(i in 1:length(codons))
  {
    codonCounts[[i]] <- genome$getCodonCountsPerGene(codons[i])
  }
  codonCounts <- do.call("cbind", codonCounts)
  # codon proportions
  codonCounts <- codonCounts / rowSums(codonCounts)
  codonCounts[is.nan(codonCounts)] <- NA # necessary if AA does not appear in gene
  
  # make empty plot
  xlimit <- range(expressionValues, na.rm = T)
  plot(NULL, NULL, xlim=xlimit, ylim=c(-0.05,1.05), 
       xlab = "", ylab="", axes = FALSE)
  # bin expression values of genes
  quantiles <- quantile(expressionValues, probs = seq(0.05, 0.95, 0.05), na.rm = T)
  for(i in 1:length(quantiles))
  {
    if(i == 1){
      tmp.id <- expressionValues < quantiles[i]
    }else if(i == length(quantiles)){
      tmp.id <- expressionValues > quantiles[i]
    }else{
      tmp.id <- expressionValues > quantiles[i] & expressionValues < quantiles[i + 1]
    }

    # plot quantiles
    means <- colMeans(codonCounts[tmp.id,], na.rm = T)
    std <- apply(codonCounts[tmp.id,], 2, sd, na.rm = T)
    for(k in 1:length(codons))
    {
      points(median(expressionValues[tmp.id]), means[k], 
             col=.codonColors[[ codons[k] ]] , pch=19, cex = 0.5)
      lines(rep(median(expressionValues[tmp.id]),2), c(means[k]-std[k], means[k]+std[k]), 
            col=.codonColors[[ codons[k] ]], lwd=0.8)
    }
  }
  
  # draw model fit
  #codonProbability <- do.call("rbind", codonProbability)
  for(i in 1:length(codons))
  {
    lines(phis, codonProbability[, i], col=.codonColors[[ codons[i] ]])
  }
  colors <- unlist(.codonColors[codons])
  
  # add indicator to optimal codon
  optim.codon.index <- which(min(c(selection, 0)) == c(selection, 0))
  codons[optim.codon.index] <- paste0(codons[optim.codon.index], "*") 
  legend("topleft", legend = codons, col=colors, bty = "n", lty=1, cex=0.75)
  
  return(xlimit)
}

