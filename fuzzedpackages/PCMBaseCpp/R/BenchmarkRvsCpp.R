#' Data for performing a benchmark
#'
#' A dataset containing three triplets trees, trait-values and models to
#' evaluate the likelihood calculation times for R and C++ implementations. 
#'
#' @format A data frame with 4 rows and 8 variables:
#' \describe{
#'   \item{tree}{phylogenetic tree (phylo) with set edge.regimes member}
#'   \item{model}{MGPM model used to simulate the data in X}
#'   \item{X}{trait values}
#'   \item{ll}{log-likelihood value}
#'   \item{modelBM}{a random BM model}
#'   \item{llBM}{log-likelihood value form modelBM}
#'   \item{modelOU}{a random OU model}
#'   \item{llOU}{log-likelihood value for modelOU}
#' }
"benchmarkData"

#' Results from running a performance benchmark on a personal computer 
#' including the time for parameter transformation
#' 
#' @format A data.table 
"benchmarkResults"


#' Results from running a performance benchmark on a personal computer 
#' excluding the time for parameter transformation
#' 
#' @format A data.table 
"benchmarkResultsNoTransform"

#' Evaluate the likelihood calculation times for example trees and data
#' @param data a `data.frame` with at least the following columns: 
#' \itemize{
#' \item{tree: }{a list column of phylo objects with an edge.part member set.}
#' \item{X: }{a list column of k x N numerical matrices.}
#' \item{model: }{a list column of PCM objects.}
#' }
#' Defaults: to `benchmarkData`, which is small data.table included
#' with the PCMBaseCpp package.
#' @param includeR logical (default TRUE) indicating if likelihood calculations
#' in R should be included in the benchmark (can be slow).
#' @param includeTransformationTime logical (default TRUE) indicating if the time for
#' \code{\link{PCMApplyTransformation}} should be included in the benchmark.
#' @param nRepsCpp : number of repetitions for the cpp likelihood calculation 
#' calls: a bigger value increases the precision of time estimation at the 
#' expense of longer running time for the benchmark. Defaults to 10.
#' @param listOptions options to set before measuring the calculation times. 
#' Defaults to `list(PCMBase.Lmr.mode = 11, PCMBase.Threshold.EV = 0, PCMBase.Threshold.SV = 0)`. 
#' `PCMBase.Lmr.mode` corresponds to the parallel traversal mode for the tree 
#' traversal algorithm (see 
#' \href{https://venelin.github.io/SPLITT/articles/SPLITTTraversalModes.html}{this page}
#' for possible values).
#' @param doProf logical indicating if profiling should be activated (see Rprof
#' from the utils R-package). Default: FALSE. Additional arguments to Rprof can 
#' be specified by assigning lists of arguments to the options 'PCMBaseCpp.ArgsRprofR'
#' and 'PCMBaseCpp.ArgsRprofCpp'. The default values for both options is
#' \code{list(append = TRUE, line.profiling = TRUE)}.
#' @param RprofR.out,RprofCpp.out character strings indicating Rprof.out files 
#' for the R and Cpp implementations; ignored if doProf is FALSE. Default values:
#' 'RprofR.out' and 'Rprofcpp.out'.
#' @return a data.frame.
#' @importFrom  PCMBase PCMInfo PCMLik PCMOptions MGPMDefaultModelTypes PCMTreeNumTips PCMTreeNumParts PCMApplyTransformation
#' @importFrom stats logLik
#' @examples
#' library(PCMBase)
#' library(PCMBaseCpp)
#' library(data.table)
#' 
#' testData <- PCMBaseCpp::benchmarkData[1]
#' # original MGPM model
#' MiniBenchmarkRvsCpp(data = testData)
#' 
#' # original MGPM model and parallel mode
#' MiniBenchmarkRvsCpp(
#' data = testData,
#' listOptions = list(PCMBase.Lmr.mode = 21, PCMBase.Threshold.EV = 1e-9, 
#' PCMBase.Threshold.SV = 1e-9))
#' 
#' # single-trait data, original MGPM model and single mode and enabled option 
#' # PCMBase.Use1DClasses
#' MiniBenchmarkRvsCpp(
#' data = PCMBaseCpp::benchmarkData[1, list(
#'  tree, 
#'  X = lapply(X, function(x) x[1,, drop=FALSE]), 
#'  model = lapply(model, function(m) PCMExtractDimensions(m, dims = 1)))],
#' listOptions = list(
#'   PCMBase.Lmr.mode = 11, 
#'   PCMBase.Threshold.EV = 1e-9, 
#'   PCMBase.Threshold.SV = 1e-9,
#'   PCMBase.Use1DClasses = FALSE))
#' 
#' @importFrom utils Rprof
#' @export
MiniBenchmarkRvsCpp <- function(
  data = PCMBaseCpp::benchmarkData, 
  includeR = TRUE,
  includeTransformationTime = TRUE,
  nRepsCpp = 10L, 
  listOptions = list(
    PCMBase.Lmr.mode = 11, PCMBase.Threshold.EV = 0, PCMBase.Threshold.SV = 0),
  doProf = FALSE, RprofR.out = "RprofR.out", RprofCpp.out = "RprofCpp.out") {
  
  listCurrentOptions <- options()
  
  do.call(options, listOptions)
  
  modelTypes <- MGPMDefaultModelTypes()
  
  res <- do.call(rbind, lapply(seq_len(nrow(data)), function(i) {
    tree <- data$tree[[i]]
    X <- data$X[[i]]
    model <- data$model[[i]]
    
    if(!includeTransformationTime) {
      # Apply the transformation first
      model <- PCMApplyTransformation(model) 
    }
    
    metaIR <- PCMInfo(X, tree, model)
    
    metaICpp <- PCMInfoCpp(X, tree, model, metaI = metaIR)
    
    valueR <- valueCpp <- NA_real_
  
    if(includeR) {
      if(doProf) {
        do.call(
          Rprof, 
          c(list(filename = RprofR.out),
            getOption("PCMBaseCpp.ArgsRprofR", 
                      list(append = TRUE, line.profiling = TRUE))))
      }
      timeR <- system.time({
        valueR <- PCMLik(X, tree, model, metaI = metaIR)
      })[3]
      if(doProf) {
        Rprof(NULL)
      }
    } else {
      timeR <- NA_real_
    }
    
    if(doProf) {
      do.call(
        Rprof, 
        c(list(filename = RprofCpp.out),
          getOption("PCMBaseCpp.ArgsRprofCpp", list(append = TRUE, line.profiling = TRUE))))
    }
    timeCpp <- system.time(
      valueCpp <-replicate(
        nRepsCpp, 
        PCMLik(X, tree, model, metaI = metaICpp))[1])[3] / nRepsCpp
    if(doProf) {
      Rprof(NULL)
    }
    
    data.frame(
      N = PCMTreeNumTips(tree), 
      R = PCMTreeNumParts(tree),
      mapping = I(list(names(attr(model, "modelTypes"))[attr(model, "mapping")])), 
      type = I(list(class(model))),
      PCMBase.Lmr.mode = PCMOptions()$PCMBase.Lmr.mode,
      logLik = valueR, 
      logLikCpp = valueCpp, 
      timeR = unname(timeR), 
      timeCpp = unname(timeCpp))
  }))
  resetOptionsStatus <- try(do.call(options, listCurrentOptions), silent = TRUE)
  if(inherits(resetOptionsStatus, "try-error")) {
    warning(
      paste0(
        "Error while resetting options at the end of a call to MiniBenchmarkRvsCpp: ",
        toString(resetOptionsStatus)))
  }
  res
}

#' A log-likelihood calculation time comparison for different numbers of traits 
#' and option-sets
#' @param ks a vector of positive integers, denoting different numbers of traits. 
#' Default: \code{c(1, 2, 4, 8)}.
#' @param optionSets a named list of lists of PCM-options. If NULL (the default) 
#' the option set is set to \code{DefaultBenchmarkOptions(k, includeParallelMode)}
#' for each \code{k} in \code{ks} (see the code in 
#' \code{PCMBaseCpp:::DefaultBenchmarkOptions}). 
#' @param includeParallelMode logical (default TRUE) indicating if the default 
#' optionSet should include parallel execution modes, i.e. setting the option 
#' PCMBase.Lmr.mode to 21 instead of 11. This argument is taken into account 
#' only with the argument \code{optionSets} set to NULL (the default). 
#' @param verbose logical indicating if log-messages should be printed to the console during the benchmark. Default FALSE.
#' @inheritParams MiniBenchmarkRvsCpp
#' @return a data.table for results similar to the data.table returned from \code{\link{MiniBenchmarkRvsCpp}} with 
#' additional columns for k, option-set and the type of model. 
#' @export
#' @importFrom data.table rbindlist
#' @importFrom PCMBase PCMExtractDimensions
BenchmarkRvsCpp <- function(
  ks = c(1, 2, 4, 8),
  includeR = TRUE,
  includeTransformationTime = TRUE, 
  optionSets = NULL,
  includeParallelMode = TRUE,
  doProf = FALSE, RprofR.out = "RprofR.out", RprofCpp.out = "RprofCpp.out",
  verbose = FALSE) {

  benchmarkData <- PCMBaseCpp::benchmarkData
  X <- tree <-  model <- modelBM <- modelOU <- modelType <- N <- R <- mapping <- 
    PCMBase.Lmr.mode <- logLik <- logLikCpp <- timeR <- timeCpp <- NULL
  
  resultList <- lapply(ks, function(k) {
    if(is.null(optionSets)) {
      optionSets <- DefaultBenchmarkOptions(k, includeParallelMode)
    } 
    
    
    resultList <- lapply(names(optionSets), function(oset) {
      
      if(verbose) {
        cat("Performing benchmark for k: ", k, "; optionSet: ", oset, "...\n")
      }
      fk <- findBiggestFactor(k, 8) 
      ds <- seq_len(fk)
      nRB <- k / fk
      
      testData <- rbindlist(list(
        benchmarkData[, list(
          k = k,
          modelType = "MGPM (A-F)",
          options = oset,
          tree, 
          X = lapply(X, function(x) x[rep(ds, nRB),, drop=FALSE]), 
          model = lapply(model, function(m) PCMExtractDimensions(m, dims = ds, nRepBlocks = nRB)))],
        benchmarkData[, list(
          k = k,
          modelType = "BM (B)",
          options = oset,
          tree, 
          X = lapply(X, function(x) x[rep(ds, nRB),, drop=FALSE]), 
          model = lapply(modelBM, function(m) PCMExtractDimensions(m, dims = ds, nRepBlocks = nRB)))],
        benchmarkData[, list(
          k = k,
          modelType = "OU (E)",
          options = oset,
          tree, 
          X = lapply(X, function(x) x[rep(ds, nRB),, drop=FALSE]), 
          model = lapply(modelOU, function(m) PCMExtractDimensions(m, dims = ds, nRepBlocks = nRB)))]))
      
      resultData <- MiniBenchmarkRvsCpp(
        data = testData,
        includeR = includeR,
        includeTransformationTime = includeTransformationTime,
        listOptions = optionSets[[oset]],
        doProf = doProf, RprofR.out = RprofR.out, RprofCpp.out = RprofCpp.out)
      
      resultData <- cbind(testData[, list(k, modelType, options)], resultData)
      if(verbose) {
        print(resultData[, list(
          k, modelType, N, R, mode = PCMBase.Lmr.mode, logLik, 
          logLikCpp, timeR, timeCpp)])
      }
      resultData
    })
    rbindlist(resultList, use.names = TRUE)
  })
  rbindlist(resultList, use.names = TRUE)
}

DefaultBenchmarkOptions <- function(k, includeParallelMode) {
  if(k == 1) {
    os <- list(
      `serial / 1D-multiv.` = list(PCMBase.Lmr.mode = 11, 
                                   PCMBase.Threshold.EV = 0, 
                                   PCMBase.Threshold.SV = 0, 
                                   PCMBase.Use1DClasses = FALSE),
      `serial / 1D-univar.` = list(PCMBase.Lmr.mode = 11, 
                                   PCMBase.Threshold.EV = 0, 
                                   PCMBase.Threshold.SV = 0, 
                                   PCMBase.Use1DClasses = TRUE))
    if(includeParallelMode) {
      os <- c(
        os, 
        list(
          `parallel / 1D-multiv.` = list(PCMBase.Lmr.mode = 21, 
                                         PCMBase.Threshold.EV = 0, 
                                         PCMBase.Threshold.SV = 0, 
                                         PCMBase.Use1DClasses = FALSE),
          `parallel / 1D-univar.` = list(PCMBase.Lmr.mode = 21, 
                                         PCMBase.Threshold.EV = 0, 
                                         PCMBase.Threshold.SV = 0, 
                                         PCMBase.Use1DClasses = TRUE)) )
    }
    os
  } else {
    os <- list(
      `serial / 1D-multiv.` = list(PCMBase.Lmr.mode = 11, 
                                   PCMBase.Threshold.EV = 0, 
                                   PCMBase.Threshold.SV = 0, 
                                   PCMBase.Use1DClasses = FALSE))
    if(includeParallelMode) {
      os <- c(
        os,
        list(`parallel / 1D-multiv.` = list(PCMBase.Lmr.mode = 21, 
                                            PCMBase.Threshold.EV = 0, 
                                            PCMBase.Threshold.SV = 0, 
                                            PCMBase.Use1DClasses = FALSE)))
    }
    os
  }
}



findBiggestFactor <- function(k, fMax) {
  f <- fMax
  while(f >= 1) {
    if(k %% f == 0) {
      break
    }
    f <- f - 1
  }
  f
}

