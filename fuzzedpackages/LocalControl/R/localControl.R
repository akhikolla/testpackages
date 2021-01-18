#' @importFrom Rcpp evalCpp sourceCpp
#' @importFrom utils data
#' @importFrom graphics box layout
#' @importFrom stats median quantile runif

#' @useDynLib LocalControl, .registration = TRUE

#' @name LocalControl
#'
#' @title Local Control
#'
#' @description Implements a non-parametric methodology for correcting biases when comparing the outcomes of two treatments in a
#'  cross-sectional or case control observational study. This implementation of Local Control uses nearest neighbors to each point within
#'  a given radius to compare treatment outcomes. Local Control matches along a continuum of similarity (radii), clustering the near neighbors
#'  to a given observation by variables thought to be sources of bias and confounding. This is analogous to combining a host of smaller studies
#'  that are each homogeneous within themselves, but represent the spectrum of variability of observations across diverse subpopulations. As the
#'  clusters get smaller, some of them can become noninformative, whereby all cluster members contain only one treatment, and there is no basis
#'  for comparison. Each observation has a unique set of near-neighbors, and the approach becomes more akin to a non-parametric
#'  density estimate using similar observations within a covariate hypersphere of a given radius. The global treatment difference is
#'  taken as the average of the treatment differences of the neighborhood around each observation.
#'
#'  While \code{\link{LocalControlClassic}} uses the number of clusters as a varying parameter to visualize treatment differences
#'  as a function of similarity of observations, this function instead uses a varying radius. The maximum radius enclosing all observations
#'  corresponds to the biased estimate which compares the outcome of all those with treatment A versus all those with treatment B.
#'  An easily interpretable graph can be created to illustrate the change in estimated outcome difference between two treatments, on average, across
#'  all clusters, as a function of using smaller and more homogenous clusters. The \code{\link{LocalControlNearestNeighborsConfidence}} procedure
#'  statistically resamples this Local Control process to generate confidence estimates.
#'  It is also helpful to plot a box-plot of the local treatment difference at a radius of zero, requiring that every observation has at
#'  least one perfect match on the other treatment. When perfect matches exist, one can estimate the treatment difference without making
#'  assumptions about the relative importance of the clustering variables. The \code{\link{plot.LocalControlCS}} function will plot both
#'  visualizations in a single graph.
#'
#'
#' @param data DataFrame containing all variables which will be used for the analysis.
#'
#' @param modelForm A formula containing the necessary variables for Local Control analysis. This can be used as an alternative
#'                  to the primary interface for cross-sectional studies. The formula should be in the following format:
#'                  "outcome ~ treatment | clusterVar1 ... clusterVarN".
#'
#' @param outcomeType Specifys the outcome type for the analysis. 
#'
#' @param treatmentColName A string containing the name of a column in data.
#'                         The column contains the treatment variable specifying the treatment groups.
#'
#' @param outcomeColName A string containing the name of a column in data.
#'                       The column contains the outcome variable to be compared between the treatment groups.
#'
#' @param cenCode A value specifying which of the outcome values corresponds to a censored observation.
#'
#' @param clusterVars A character vector containing column names in data.
#'                    Each column contains an X-variable, or covariate which will be used to form patient clusters.
#'
#' @param timeColName A string containing the name of a column in data.
#'                    The column contains the time to outcome for each of the observations in data.
#'
#' @param treatmentCode (optional) A string containing one of the factor levels from the treatment column.
#'                      If provided, the corresponding treatment will be considered "Treatment 1".
#'                      Otherwise, the first "level" of the column will be considered the primary treatment.
#'
#' @param labelColName (optional) A string containing the name of a column from data.
#'                     The column contains labels for each of the observations in data, defaults to the row indices.
#'
#' @param radStepType (optional) Used in the generation of correction radii.
#'                    The step type used to generate each correction radius after the maximum.
#'                    Currently accepts "unif" and "exp" (default).
#'                    "unif" for uniform decay ex: (radDecayRate = 0.1) (1, 0.9, 0.8, 0.7, ..., ~minRadFract, 0)
#'                    "exp" for exponential decay ex: (radDecayRate = 0.9) (1, 0.9, 0.81, 0.729, ..., ~minRadFract, 0)
#'
#' @param radDecayRate (optional) Used in the generation of correction radii.
#'                     The size of the "step" between each of the generated correction radii.
#'                     If radStepType == "exp", radDecayRate must be a value between (0,1).
#'                     This value defaults to 0.8.
#'
#' @param radMinFract (optional) Used in the generation of correction radii.
#'                    A floating point number representing the smallest fraction of the maximum radius to use as a correction radius.
#'
#' @param radiusLevels (optional) By default, Local Control builds a set of radii to fit data.
#'                     The radiusLevels parameter allows users to override the construction by explicitly providing a set of radii.
#'
#' @param normalize (optional) Logical value. Tells local control if it should or should not normalize the covariates. Default is TRUE.
#'
#' @param verbose (optional) Logical value. Display or suppress the console output during the call to Local Control. Default is FALSE.
#'
#' @param numThreads (optional) An integer value specifying the number of threads which will be assigned to the analysis.
#'                   The maximum number of threads varies depending on the system hardware. Defaults to 1 thread.
#'
#' @return A list containing the results from the call to LocalControl.
#' \itemize{
#'   \item{outcomes} {List containing two dataframes for the average T1 and T0 outcomes within each cluster at each radius.}
#'   \item{counts} {List containing two dataframes which hold the number of T1 and T0 patients within each cluster at each radius.}
#'   \item{ltds} {Dataframe containing the average LTD within each cluster at each radius.}
#'   \item{summary} {Dataframe containing summary statistics about the analysis for each radius.}
#'   \item{params} {List containing the parameters used to call LocalControl.}
#' }
#' @references
#' \itemize{
#'    \item Lauve NR, Nelson SJ, Young SS, Obenchain RL, Lambert CG. LocalControl: An R Package for Comparative Safety and Effectiveness Research. Journal of Statistical Software. 2020. p. 1–32. Available from: http://dx.doi.org/10.18637/jss.v096.i04
#'    \item Fischer K, Gartner B, Kutz M. Fast Smallest-Enclosing-Ball Computation in High Dimensions. In: Algorithms - ESA 2003. Springer, Berlin, Heidelberg; 2003:630-641.
#'    \item Martin Kutz, Kaspar Fischer, Bernd Gartner. miniball-1.0.3. \url{https://github.com/hbf/miniball}.
#' }
#' @examples
#'  # cross-sectional
#'
#'  data(lindner)
#'  linVars <- c("stent", "height", "female", "diabetic", "acutemi",
#'               "ejecfrac", "ves1proc")
#'  csresults = LocalControl(data = lindner,
#'                           clusterVars = linVars,
#'                           treatmentColName = "abcix",
#'                           outcomeColName = "cardbill",
#'                           treatmentCode = 1)
#'  plot(csresults)
#'
#'
#'  # survival / competing risks example
#'
#'  data(cardSim)
#'  crresults = LocalControl(data = cardSim, outcomeType = "survival",
#'                           outcomeColName = "status",
#'                           timeColName = "time",
#'                           treatmentColName = "drug",
#'                           treatmentCode = 1,
#'                           clusterVars = c("age", "bmi"))
#'  plot(crresults)
#' @export
LocalControl<- function(data,
                        modelForm = NULL,
                        outcomeType = "default",
                        treatmentColName,
                        outcomeColName,
                        cenCode = 0,
                        clusterVars,
                        timeColName="",
                        treatmentCode,
                        labelColName="",
                        radStepType = "exp",
                        radDecayRate = 0.8,
                        radMinFract = 0.01,
                        radiusLevels = numeric(),
                        normalize = TRUE,
                        verbose = FALSE,
                        numThreads=1)
{
  plist = list()

  if(!missing(modelForm)){
    if(!all(all.vars(modelForm) %in% names(data))){
      stop("All formula variables must be columns in data.")
    }
    else{
      mvars = all.vars(modelForm)
      outcomeColName =  mvars[1]
      treatmentColName =  mvars[2]
      clusterVars = mvars[3:length(mvars)]
    }
  }

  if(missing(data) || !inherits(data, "data.frame")){
    stop("Please provide a DataFrame with columns containing your Xs and Ys.")
  }
  plist$data = deparse(substitute(data))


  if(missing(outcomeColName) || !is.element(outcomeColName, dimnames(data)[[2]])){
    stop("outcomeColName must be a column in dataframe")
  }
  plist$outcomeColName = outcomeColName
  numRisks = length(unique(data[,outcomeColName]))
  if(cenCode %in% unique(data[,outcomeColName])){
    numRisks = numRisks - 1
  }

  if(missing(treatmentColName) || !is.element(treatmentColName, dimnames(data)[[2]])){
    stop("Please provide the name of the treatment column in your DataFrame.")
  }
  plist$treatmentColName = treatmentColName

  if(missing(treatmentCode)){
    treatmentCode = data[1, treatmentColName]
  }else if(!is.element(as.character(treatmentCode), levels(as.factor(data[,treatmentColName])))){
    stop("treatmentCode must be a level of treatmentCol")
  }
  plist$treatmentCode = treatmentCode
  treatments = as.numeric(data[,treatmentColName] == treatmentCode)

  if(radStepType == "exp" && (radDecayRate >= 1 || radDecayRate <= 0)){
    stop("Invalid decay rate, please choose a value in the range (0,1) (exclusive)")
  }
  plist$radDecayRate = radDecayRate

  if(radMinFract >= 1 || radMinFract <= 0){
    stop("Invalid minimum radius fraction. please choose a value in the range (0,1) (exclusive)")
  }
  plist$radMinFract = radMinFract

  if(outcomeType == "survival" && (missing(timeColName) ||
    !is.element(timeColName, dimnames(data)[[2]]))){
    stop("Please provide the name of time to reach outcome column in your DataFrame.")
  }
  plist$outcomeType = outcomeType

  if(missing(clusterVars) || !all(clusterVars %in% names(data))){
    stop("Please provide a list of clustering variables contained in data.")
  }
  plist$clusterVars = clusterVars

  if(verbose) message("Params checked...")
  clusterData = getClusterDF(data, clusterVars)
  cvData = clusterData$clusterDF
  svID = clusterData$svID
  svnames = clusterData$svnames

  if(normalize) cvData = normalizeData(cvData)

  if(verbose) message("Covars normalized...")

  #Calculate the vector of radii based on the covariate data
  if(length(radiusLevels) > 0){
    rads = radiusLevels
  }else{
    maxRad <- getMaxDist(cvData) #max dist between any two observations
    if(verbose) message("Got max rad...")
    rads <- getRadVect(maxRad,
                       radStepType = radStepType,
                       radDecayRate = radDecayRate,
                       radMinFract = radMinFract)
  }

  if(verbose) message("Got all rads...")

  patFrame = cbind(treatments, as.numeric(data[,outcomeColName]), cvData)

  if(outcomeType == "survival"){

    failTimes = c( as.numeric(data[,timeColName]), 0)

    returnFrame = newCRLC(patFrame, rads, failTimes, cenCode, numThreads)
    class(returnFrame) <- "LocalControlCR"

    if(numRisks > 1){
      maxTimes = unlist(by(data[,timeColName] , data[,treatmentColName], max))
      maxSharedTime = min(maxTimes)

      numT1 = nrow(data[which(data[,treatmentColName] == treatmentCode), ])
      numT0 = nrow(data[which(data[,treatmentColName] != treatmentCode), ])

      numEvents = length(returnFrame$Events$T0$rad_1[['Failcode_1']])

      pepe = lcHypothesis(returnFrame, maxSharedTime, numT1, numT0)
      extra = pepe[[2]]
      pepe = pepe[[1]]
      returnFrame$extra = extra
      returnFrame$pepe = pepe
    }
  }
  else{
    returnFrame <- newLC(patFrame, rads, numThreads)
    class(returnFrame) <- "LocalControlCS"
    returnFrame$ltds = returnFrame$outcomes$T1 - returnFrame$outcomes$T0
  }

  returnFrame$rads = rads
  if(normalize) returnFrame$normData = cvData

  returnFrame$params = plist

  returnFrame
}

#' @export
print.LocalControlCR = function(x, ...){
  lccr = x
  if(!inherits(lccr, "LocalControlCR")){
    stop('Please provide output from a call to LocalControl where outcomeType = "survival"')
  }
  message("LocalControl survival analysis results.\n")
  message(paste0("Data: ", lccr$params$data))
  message(paste0("Outcome variable: ", lccr$params$outcomeColName))
  message(paste0("Treatment variable: ", lccr$params$treatmentColName))
  message(paste0("Treatment 1: ", lccr$params$treatmentCode))
  message(paste0("Levels of correction: ", length(lccr$rads)))
  message(paste0("Clustering variables: ", paste0("\n\t- ", paste0(lccr$params$clusterVars, collapse = "\n\t- "))))
}

#' @export
print.LocalControlCS = function(x, ...){
  lccs = x
  if(!inherits(lccs, "LocalControlCS")){
    stop('Please provide output from a call to LocalControl where outcomeType = "default"')
  }
  message("LocalControl cross-sectional analysis results.\n")
  message(paste0("Data: ", lccs$params$data))
  message(paste0("Outcome variable: ", lccs$params$outcomeColName))
  message(paste0("Treatment variable: ", lccs$params$treatmentColName))
  message(paste0("Treatment 1: ", lccs$params$treatmentCode))
  message(paste0("Levels of correction: ", length(lccs$rads)))
  sumre = summary(lccs)
  message(paste0("Max-radius treatment difference: ", formatC(sumre[which.max(sumre$radius), "ltd"], digits = 2,format = "f")))
  hasPerfectMatches = sumre[which(sumre$radius == 0), "pct_data"] > 0
  message(paste0("Data contains perfect matches: ", hasPerfectMatches))
  if(hasPerfectMatches){
    message(paste0("Perfect match treatment difference: ", formatC(sumre[which(sumre$radius == 0), "ltd"], digits = 2,format = "f")))
  }
  message(paste0("Clustering variables: ", paste0("\n\t- ", paste0(lccs$params$clusterVars, collapse = "\n\t- "))))
}

#' @export
summary.LocalControlCR = function(object, ...){
  lccr = object
  sumre = data.frame("radius" = lccr$rads,
                     "pct_informative" = lccr$NumInf / lccr$NumInf[1],
                     "pct_radius" = lccr$rads / max(lccr$rads))
  if(!is.null(lccr$pepe)){
    sumre = cbind(sumre, lccr$pepe)
  }

  sumre
}

#' @export
summary.LocalControlCS = function(object, ...){
lccs = object
data.frame("radius" = lccs$rads,
           "ltd" = apply(X = lccs$ltds, MARGIN = 2, FUN = mean, na.rm = TRUE),
           "pct_radius" = lccs$rads / max(lccs$rads),
           "pct_data" = ((nrow(lccs$ltds) - apply(X = lccs$ltds, MARGIN = 2,
                          FUN = function(z) sum(is.na(z)))) / nrow(lccs$ltds)))
}

#' @name plot.LocalControlCR
#' @title Plot cumulative incidence functions (CIFs) from Local Control.
#' @description Given the results from LocalControl with outcomeType = "survival",
#'   plot a corrected and uncorrected cumulative incidence function (CIF) for both groups.
#'
#' @param x Return object from LocalControl with outcomeType = "survival".
#' @inheritDotParams graphics::plot -x -y
#' @param rad2plot The index or name ("rad_#") of the radius to plot. By default, the radius with pct_informative closest to
#'   0.8 will be selected.
#' @param xlim The x axis bounds. Defaults to c(0, max(lccrResults$Failtimes)).
#' @param ylim The y axis bounds. Defaults to c(0,1).
#' @param col1 The plot color for group 1.
#' @param col0 The plot color for group 0.
#' @param xlab The x axis label. Defaults to "Time".
#' @param ylab The y axis label. Defaults to "Cumulative incidence".
#' @param legendLocation The location to place the legend. Default "topleft".
#' @param main The main plot title. Default is empty.
#' @param group1 The name of the primary group (Treatment 1).
#' @param group0 The name of the secondary group (Treatment 0).
#' @references
#' \itemize{
#'    \item Lauve NR, Nelson SJ, Young SS, Obenchain RL, Lambert CG. LocalControl: An R Package for Comparative Safety and Effectiveness Research. Journal of Statistical Software. 2020. p. 1–32. Available from: http://dx.doi.org/10.18637/jss.v096.i04
#' }
#' @examples
#' data("cardSim")
#' results = LocalControl(data = cardSim,
#'                        outcomeType = "survival",
#'                        outcomeColName = "status",
#'                        timeColName = "time",
#'                        treatmentColName = "drug",
#'                        treatmentCode = 1,
#'                        clusterVars = c("age", "bmi"))
#' plot(results)
#'
#' @export
plot.LocalControlCR = function(x, ..., rad2plot, xlim, ylim = c(0,1),
                               col1 = "blue", col0 = "red",
                               xlab = "Time", ylab = "Cumulative incidence",
                               legendLocation = "topleft", main = "",
                               group1 = "Treatment 1", group0 = "Treatment 0"){
  lccrResults = x

  if(is.null(lccrResults$summary)){
    lccrResults$summary = summary(lccrResults)
  }
  if(missing(xlim)){
    xlim = c(0, max(lccrResults$Failtimes))
  }
  if(missing(rad2plot)){
    rad2plot = which.min(abs(lccrResults$summary[,"pct_informative"] - 0.8))
  }

  fcs = names(lccrResults$CIF)

  plot(NA, xlim = xlim, ylim = ylim, ylab = ylab, xlab = xlab, main = main, ...)
  grid()

  for(fc in fcs){
    lines(lccrResults$Failtimes, lccrResults$CIF[[fc]]$T1$rad_1, col = col1, type = "s", lwd = 2, lty = 3)
    lines(lccrResults$Failtimes, lccrResults$CIF[[fc]]$T0$rad_1, col = col0, type = "s", lwd = 2, lty = 3)
    lines(lccrResults$Failtimes, lccrResults$CIF[[fc]]$T1[,rad2plot], col = col1, type = "s", lwd = 2)
    lines(lccrResults$Failtimes, lccrResults$CIF[[fc]]$T0[,rad2plot], col = col0, type = "s", lwd = 2)
  }

  legend(legendLocation,
         col = c(col1, col1, col0, col0),
         lty = c(3,1,3,1), lwd = 2,
         legend = c(paste0(group1, " uncorrected"), paste0(group1, " corrected"),
                    paste0(group0, " uncorrected"), paste0(group0, " corrected")))

}


#' @name plot.LocalControlCS
#' @title Plots the local treatment difference as a function of radius for LocalControl.
#' @description Creates a plot where the y axis represents the local treatment difference,
#' while the x axis represents the percentage of the maximum radius. If the confidence summary (nnConfidence)
#' is provided, the 50\% and 95\% confidence estimates are also plotted.
#' @param x Return object from LocalControl with "default" outcomeType.
#' @inheritDotParams graphics::plot -x -y
#' @param nnConfidence Return object from LocalControlNearestNeighborsConfidence
#' @param ylim The y axis bounds. Defaults to c(0,1).
#' @param legendLocation The location to place the legend. Default "topleft".
#' @param ylab The y axis label. Defaults to "LTD".
#' @param xlab The x axis label. Defaults to "Fraction of maximum radius".
#' @param main The main plot title. Default is empty.
#'
#' @references
#' \itemize{
#'    \item Lauve NR, Nelson SJ, Young SS, Obenchain RL, Lambert CG. LocalControl: An R Package for Comparative Safety and Effectiveness Research. Journal of Statistical Software. 2020. p. 1–32. Available from: http://dx.doi.org/10.18637/jss.v096.i04
#' }
#' 
#' @examples
#' data(lindner)
#' # Specify clustering variables.
#' linVars <- c("stent", "height", "female", "diabetic",
#'              "acutemi", "ejecfrac", "ves1proc")
#'
#' # Call Local Control once.
#' linRes <- LocalControl(data = lindner,
#'                        clusterVars = linVars,
#'                        treatmentColName = "abcix",
#'                        outcomeColName = "cardbill",
#'                        treatmentCode = 1)
#'
#' # Plot the local treatment differences from Local Control without
#' # confidence intervals.
#' plot(linRes, ylim =  c(-6000, 3600))
#'
#' #If the confidence intervals are calculated:
#' #linConfidence = LocalControlNearestNeighborsConfidence(
#' #                                      data = lindner,
#' #                                      clusterVars = linVars,
#' #                                      treatmentColName = "abcix",
#' #                                      outcomeColName = "cardbill",
#' #                                      treatmentCode = 1, nBootstrap = 100)
#'
#' # Plot the local treatment difference with confidence intervals.
#' #plot(linRes, linConfidence)
#'
#' @export
plot.LocalControlCS = function(x, ..., nnConfidence, ylim,
                               legendLocation = "bottomleft", ylab = "LTD",
                               xlab = "Fraction of maximum radius", main =""){

  lcnnResults = x
  if(is.null(lcnnResults$summary)){
    lcnnResults$summary = summary(lcnnResults)
  }
  xaxi = lcnnResults$summary$pct_radius
  xl = c(1,min(lcnnResults$summary$pct_radius[which(lcnnResults$summary$pct_radius > 0)]))

  if(missing(ylim)){
    if(missing(nnConfidence)){
      ylim = c(min(lcnnResults$summary$ltd) - abs(min(lcnnResults$summary$ltd)/10),
               max(lcnnResults$summary$ltd) + max(lcnnResults$summary$ltd)/10)
    }
    else{
      ylim = c(min(nnConfidence) - abs(min(nnConfidence)/10), max(nnConfidence)+max(nnConfidence)/10)
    }
  }
  hasPMatches = (lcnnResults$summary$pct_data[nrow(lcnnResults$summary)] > 0) && (lcnnResults$summary$radius[nrow(lcnnResults$summary)] == 0)

  if(hasPMatches){
    layout(matrix(c(rep(1,9),2,rep(1,9),2), nrow = 2, ncol = 10, byrow = TRUE))
  }
  par(mar=c(5.1,4.1,4.1,0.5))
  plot(NA, xlim = xl, ylim = ylim,
       xlab = xlab,
       ylab = ylab, log = "x", panel.first=grid(equilogs=FALSE),
       main = main)

  box()
  if(!missing(nnConfidence)){
    lines(xaxi, nnConfidence$rs_median, lwd = 2, pch = 16, type = 'o')
    lines(xaxi, nnConfidence$rs_95L, lty = 1, lwd = 2, col = "green", pch = 16, type = 'o')
    lines(xaxi, nnConfidence$rs_95U, lty = 1, lwd = 2, col = "green", pch = 16, type = 'o')
    legend(legendLocation,
           lty = rep(1,2),
           lwd = rep(2,2),
           pch = rep(16,2),
           col = c("black", "green"),
           legend = c("Resampling median LTD",
                      "Resampling 95% confidence"))
  }
  else{
    lines(xaxi, lcnnResults$summary$ltd, lwd = 2, pch = 16, type = 'o')
    legend(legendLocation,
           lty = 1,
           lwd = 2,
           pch = 16,
           col = "black",
           legend = "Average LTD")
  }



  if(hasPMatches){
    if(!missing(nnConfidence)){
      pmatches = as.numeric(nnConfidence[nrow(nnConfidence),])
      par(mar=c(5.1,0.5,4.1,2.1))
      plot(NA, ylim = ylim, panel.first = grid(nx = 2, ny = NULL), xaxt = "n", yaxt = "n", xlab = NA)
      par(new = T)
      boxplot(pmatches, ylim = ylim, show.names = TRUE, names = "0", yaxt = "n",notch = F,add = T,lwd = 2)
      par(new = F)
    }
    else{
      pmatches = as.numeric(lcnnResults$summary[nrow(lcnnResults$summary),"ltd"])
      par(mar=c(5.1,0.5,4.1,2.1))
      plot( x = 0, y = pmatches, ylim = ylim, panel.first = grid(nx = 2, ny = NULL), yaxt = "n", xlab = NA, lwd = 2, pch = 16, xaxt = "n")
      lines(x = c(-.3, .3), y = c(pmatches, pmatches), lwd = 2)
      axis(side = 1, at = 0)
    }
  }
  mtext(side = 1, text = "Perfect\n matches", cex = 0.6, line = 3)
  par(mfrow=c(1,1))
  par(mar = c(5.1,4.1,4.1,2.1))
}



#' @name LocalControlCompetingRisksConfidence
#' @title Calculate confidence intervals around the cumulative incidence functions (CIFs) generated by LocalControl when outcomeType = "survival".
#' @description Given the output of \code{\link{LocalControl}}, this function produces pointwise standard error estimates
#' for the cumulative incidence functions (CIFs) using a modified version of Choudhury's approach (2002). This function currently supports
#' the creation of 90\%, 95\%, 98\%, and 99\% confidence intervals with linear, log(-log), and arcsine transformations of the estimates.
#' @param LCCompRisk Output from a successful call to LocalControl with outcomeType = "survival".
#' @param confLevel Level of confidence with which the confidence intervals will be formed. Choices are: "90\%", "95\%", "98\%", "99\%".
#' @param confTransform Transformation of the confidence intervals, defaults to arcsin ("asin").
#'   "log" and "linear" are also implemented.
#'
#' @references
#' \itemize{
#'    \item Lauve NR, Nelson SJ, Young SS, Obenchain RL, Lambert CG. LocalControl: An R Package for Comparative Safety and Effectiveness Research. Journal of Statistical Software. 2020. p. 1–32. Available from: http://dx.doi.org/10.18637/jss.v096.i04
#'    \item Choudhury JB (2002) Non-parametric confidence interval estimation for competing risks analysis: application to contraceptive data. Stat Med 21:1129-1144. doi: 10.1002/sim.1070
#' }
#'
#' @examples
#'  data(cardSim)
#'  results = LocalControl(data = cardSim,
#'                         outcomeType = "survival",
#'                         outcomeColName = "status",
#'                         timeColName = "time",
#'                         treatmentColName = "drug",
#'                         treatmentCode = 1,
#'                         clusterVars = c("age", "bmi"))
#'
#'  conf = LocalControlCompetingRisksConfidence(results)
#'
#' @export
LocalControlCompetingRisksConfidence <- function(LCCompRisk, confLevel = "95%", confTransform = "asin"){

  if(class(LCCompRisk) != "LocalControlCR"){
    stop("LCCompRisk param must be type 'LocalControlCR' returned from the function 'LocalControl()' where outcomeType = 'survival'")
  }

  numTimes = length(LCCompRisk$Failtimes)

  rads = names(LCCompRisk$Events$T0)
  numRads = length(rads)

  failCodes = names(LCCompRisk$CIF)
  numFailcodes = length(failCodes)

  numCIs = numRads * numFailcodes * 2

  message(paste0("Building ", numCIs, " confidence intervals.",
                " (NumberOfRads * NumberOfFailCodes * NumberOfTreatments = ", numRads, ' * ', numFailcodes, ' * 2)'))

  LCConfidence = list()

  for(treatment in c("T1", "T0")){
    fcList = list()

    for(failCode in failCodes){

      lowerConf = data.frame(matrix(NA, nrow = numTimes, numRads, dimnames = list(NULL, rads)))
      upperConf = data.frame(matrix(NA, nrow = numTimes, numRads, dimnames = list(NULL, rads)))

      for(rad in rads){

        ciDF = getConfidence(cifVect = LCCompRisk[["CIF"]][[failCode]][[treatment]][,rad],
                             sErr    = LCCompRisk[["SDR"]][[failCode]][[treatment]][,rad],
                             confLevel = confLevel,
                             type = confTransform)

        lowerConf[,rad] = ciDF$Lower
        upperConf[,rad] = ciDF$Upper
      }

      confInts = list(LOWER_CONF = lowerConf, UPPER_CONF = upperConf)
      fcList[[failCode]] = confInts

    }

    names(fcList) = failCodes
    LCConfidence[[treatment]] = fcList

  }

  class(LCConfidence) = "LocalControlConf"
  LCConfidence
}


#' @name LocalControlNearestNeighborsConfidence
#' @title Provides a bootstrapped confidence interval estimate for LocalControl LTDs.
#' @description Given a number of bootstrap iterations and the params used to call
#'   \code{\link{LocalControl}} with outcomeType = "default", this function calls LocalControl nBootstrap times.
#'   The 50\% and 95\% quantiles are drawn from the distribution of results to produce the LTD confidence intervals.
#'
#' @param data DataFrame containing all variables which will be used for the analysis.
#'
#' @param treatmentColName A string containing the name of a column in data.
#'                         The column contains the treatment variable specifying the treatment groups.
#'
#' @param outcomeColName A string containing the name of a column in data.
#'                       The column contains the outcome variable to be compared between the treatment groups.
#'                       If outcomeType = "survival", the outcome column holds the failure/censor assignments.
#'
#' @param clusterVars A character vector containing column names in data.
#'                    Each column contains an X-variable, or covariate which will be used to form patient clusters.
#'
#' @param treatmentCode (optional) A string containing one of the factor levels from the treatment column.
#'                      If provided, the corresponding treatment will be considered "Treatment 1".
#'                      Otherwise, the first "level" of the column will be considered the primary treatment.
#'
#' @param labelColName (optional) A string containing the name of a column from data.
#'                     The column contains labels for each of the observations in data, defaults to the row indices.
#'
#' @param radStepType (optional) Used in the generation of correction radii.
#'                    The step type used to generate each correction radius after the maximum.
#'                    Currently accepts "unif" and "exp" (default).
#'                    "unif" for uniform decay ex: (radDecayRate = 0.1) (1, 0.9, 0.8, 0.7, ..., ~minRadFract, 0)
#'                    "exp" for exponential decay ex: (radDecayRate = 0.9) (1, 0.9, 0.81, 0.729, ..., ~minRadFract, 0)
#'
#' @param radDecayRate (optional) Used in the generation of correction radii.
#'                     The size of the "step" between each of the generated correction radii.
#'                     If radStepType == "exp", radDecayRate must be a value between (0,1).
#'                     This value defaults to 0.8.
#'
#' @param radMinFract (optional) Used in the generation of correction radii.
#'                    A floating point number representing the smallest fraction of the maximum radius to use as a correction radius.
#'
#' @param radiusLevels (optional) By default, Local Control builds a set of radii to fit data.
#'                     The radiusLevels parameter allows users to override the construction by explicitly providing a set of radii.
#'
#' @param normalize (optional) Logical value. Tells local control if it should or should not normalize the covariates. Default is TRUE.
#'
#' @param verbose (optional) Logical value. Display or suppress the console output during the call to Local Control. Default is FALSE.
#'
#' @param numThreads (optional) An integer value specifying the number of threads which will be assigned to the analysis.
#'                   The maximum number of threads varies depending on the system hardware. Defaults to 1 thread.
#'
#' @param nBootstrap The number of times to resample and run LocalControl for the confidence intervals.
#' @param randSeed The seed used to set random number generator state prior to resampling. No default value, provide one for reproducible results.
#' @references
#' \itemize{
#'    \item Lauve NR, Nelson SJ, Young SS, Obenchain RL, Lambert CG. LocalControl: An R Package for Comparative Safety and Effectiveness Research. Journal of Statistical Software. 2020. p. 1–32. Available from: http://dx.doi.org/10.18637/jss.v096.i04
#'    \item Kereiakes DJ, Obenchain RL, Barber BL, Smith A, McDonald M, Broderick TM, Runyon JP, Shimshak TM, Schneider JF, Hattemer CR, Roth EM, Whang DD, Cocks D, Abbottsmith CW. Abciximab provides cost-effective survival advantage in high-volume interventional practice. Am Heart J. 2000 Oct;140(4):603–610. PMID: 11011333
#' }
#' 
#' @examples
#' \dontrun{
#' #input the abciximab study data of Kereiakes et al. (2000).
#' data(lindner)
#'
#' linVars <- c("stent", "height", "female", "diabetic", "acutemi",
#'              "ejecfrac", "ves1proc")
#' results <- LocalControl(data = lindner,
#'                         clusterVars = linVars,
#'                         treatmentColName = "abcix",
#'                         outcomeColName = "cardbill",
#'                         treatmentCode = 1)
#'
#' #Calculate the confidence intervals via resampling.
#' confResults = LocalControlNearestNeighborsConfidence(
#'                                         data = lindner,
#'                                         clusterVars = linVars,
#'                                         treatmentColName = "abcix",
#'                                         outcomeColName = "cardbill",
#'                                         treatmentCode = 1, nBootstrap = 20)
#'
#' # Plot the local treatment difference with confidence intervals.
#' plot(results, confResults)
#' }
#'
#' @export
LocalControlNearestNeighborsConfidence = function(data, nBootstrap, randSeed,
                                                  treatmentColName,
                                                  treatmentCode="",
                                                  outcomeColName,
                                                  clusterVars,
                                                  labelColName="",
                                                  numThreads=1,
                                                  radiusLevels = numeric(),
                                                  radStepType = "exp",
                                                  radDecayRate = 0.8,
                                                  radMinFract = 0.01,
                                                  normalize = TRUE, verbose = FALSE){

  ltds = list()
  if(!missing(randSeed)){
    set.seed(randSeed)
  }
  for(i in 1:nBootstrap){

    newGuys = as.integer(runif(n = nrow(data), min = 1, max = nrow(data)))

    xSults = LocalControl(data = data[newGuys, ],
                          treatmentColName = treatmentColName,
                          treatmentCode = treatmentCode,
                          outcomeColName = outcomeColName,
                          clusterVars = clusterVars,
                          labelColName = labelColName,
                          numThreads = numThreads,
                          radiusLevels = radiusLevels,
                          radStepType = radStepType,
                          radDecayRate = radDecayRate,
                          radMinFract = radMinFract,
                          normalize = normalize, verbose = verbose)

    ltds[[i]] = summary(xSults)$ltd

  }

  ltdf = data.frame(ltds)
  vsumm = varSumm(ltdf)

  vsumm
}


varSumm = function(vdf){
  dfRows = nrow(vdf)
  dfCols = 5
  sdf = data.frame(matrix(nrow = dfRows, ncol = dfCols))
  names(sdf) = c("rs_median", "rs_95U", "rs_95L")
  for(i in 1:dfRows){
    sdf[i, "rs_median"] = median(as.numeric(vdf[i,]), na.rm = T)
    sdf[i, 2:dfCols] = quantile(as.numeric(vdf[i,]), probs = c(0.025, 0.975), na.rm = T)
  }
  sdf
}

# The Pepe and Mori hypothesis testing code was derived from the code given in :
# Pintilie, M. (2006) Competing Risks - Definitions, in Competing Risks: A Practical Perspective,
# John Wiley & Sons, Ltd, Chichester, UK. doi: 10.1002/9780470870709.ch3
lcHypothesis = function(retFrame, maxTime, nT1, nT0){
  fTimes = retFrame$Failtimes[2:length(retFrame$Failtimes)]
  sTimes = c(fTimes[2:length(fTimes)], NA)

  numTimes = table(sTimes <= maxTime)["TRUE"]

  timeDiffs = fTimes - sTimes
  timeDiffs = timeDiffs[1:numTimes]
  fTimes = fTimes[1:numTimes]

  numLimits = length(retFrame$Events$T1)
  pepeFrame = data.frame(matrix(nrow = numLimits, ncol = 2))
  names(pepeFrame) = c("chisq", "p_value")

  for( i in 1:numLimits){
    radIdx = names(retFrame$Events$T1)[i]
    numInf = retFrame$NumInf[i]

    t1 = data.frame(Failcode_1 = retFrame$Events$T1[[radIdx]]$Failcode_1[2:(numTimes+1)],
                    Failcode_2 = retFrame$Events$T1[[radIdx]]$Failcode_2[2:(numTimes+1)],
                    Censored   = retFrame$Events$T1[[radIdx]]$Censored[2:(numTimes+1)])

    t0 = data.frame(Failcode_1 = retFrame$Events$T0[[radIdx]]$Failcode_1[2:(numTimes+1)],
                    Failcode_2 = retFrame$Events$T0[[radIdx]]$Failcode_2[2:(numTimes+1)],
                    Censored   = retFrame$Events$T0[[radIdx]]$Censored[2:(numTimes+1)])

    t1 = round(t1 * nT1 / numInf)
    t0 = round(t0 * nT0 / numInf)

    t1 = lcSurvival(t1, nT1, numTimes)
    t0 = lcSurvival(t0, nT0, numTimes)

    timeWeights = (t1$censKM * t0$censKM * (nT1 + nT0)) / (nT1 * t1$censKM + nT0 * t0$censKM)

    si = timeWeights * (t1$mainCIF - t0$mainCIF) * timeDiffs
    s = sqrt(nT1*nT0/(nT1 + nT0))*sum(si, na.rm = T)

    sig1 = pmSigma(t1, timeWeights, timeDiffs)
    sig0 = pmSigma(t0, timeWeights, timeDiffs)

    sigma = nT0 * nT1 * (sig1+sig0)/(nT0 + nT1)
    z = s^2 / sigma
    pvalue = pchisq(z, 1, lower.tail = F)
    if(i == 1){
      extraFrame = t1
    }
    pepeFrame[i, ] = c(z, pvalue)

  }

  r = list()
  r[[1]] = pepeFrame
  r[[2]] = extraFrame

  r
}

lcSurvival <- function(df, nStart, nTimes){
  df$allEvents = apply(X = df, MARGIN = 1, FUN = function(x) sum(x))
  df$atRisk = nStart - cumsum(c(0, df$allEvents[1:(nTimes-1)]))
  df$mainEvents = df$Failcode_1
  df$compEvents = df$Failcode_2
  df$censEvents = df$Failcode_2 + df$Censored
  df$allFailures = df$Failcode_1 + df$Failcode_2
  #km
  df$kapHazard = (df$atRisk - df$allFailures) / df$atRisk
  df$km = cumprod(df$kapHazard)
  df$km = c(1, df$km[1:(length(df$km)-1)])
  #cifs
  df$mainHazard = df$Failcode_1 / df$atRisk * df$km
  df$mainCIF = cumsum(df$mainHazard)
  df$compHazard = df$Failcode_2 / df$atRisk * df$km
  df$compCIF = cumsum(df$compHazard)
  #Kaplan-Meier of the censoring distribution
  df$censHazard = (df$atRisk - df$censEvents) / df$atRisk
  df$censKM = cumprod(df$censHazard)
  df$censKM = c(1, df$censKM[1:(length(df$censKM)-1)])

  df
}

pmSigma = function(df, wt, dt){
  temp = wt * (1-df$mainCIF) * dt
  temp[is.na(temp)] = 0
  v1 = rev(cumsum(rev(temp)))

  temp = wt * dt
  temp[is.na(temp)] = 0
  v2 = df$compCIF * rev(cumsum(rev(temp)))

  temp = wt * df$mainCIF * dt
  temp[is.na(temp)] = 0
  v3 = rev(cumsum(rev(temp)))

  brackets = v1 - v2
  sigt = (df$mainEvents * brackets^2 + (df$allFailures - df$mainEvents) * v3^2)/(df$atRisk * (df$atRisk - 1))
  sig = sum(sigt, na.rm=T)

  sig
}

getConfidence = function(cifVect, sErr, confLevel = "95%", type = "asin"){

  zVal = switch(confLevel, "90%" = 1.645, "95%" = 1.96, "98%" = 2.326, "99%" = 2.576)

  n = length(cifVect)

  if(type == "linear"){
    up = cifVect + zVal*sErr
    low = cifVect - zVal*sErr

  }
  else if(type == "log"){
    up = cifVect ^ exp(zVal*sErr/(cifVect*log(cifVect)))
    low = cifVect ^ exp(-zVal*sErr/(cifVect*log(cifVect)))

    low[is.na(low)] = 0
    up[is.na(up)] = 0
  }
  else if(type == "asin"){
    cifVect = as.double(cifVect)
    lowSin = asin(sqrt(cifVect)) - .5*zVal*sErr*sqrt(1/(cifVect*(1-cifVect)))
    upSin = asin(sqrt(cifVect)) + .5*zVal*sErr*sqrt(1/(cifVect*(1-cifVect)))

    lowSin[is.na(lowSin)] = 0
    upSin[is.na(upSin)] = 0
    up = rep(0,n)
    low = rep(0,n)

    for(d in 1:n){
     up[d] = (sin(min(pi/2, upSin[d])))^2
     low[d] = (sin(max(0, lowSin[d])))^2

    }
  }
  data.frame("Upper" = up, "Lower" = low)
}

getRadVect <- function(maxRad, radStepType = "exp", radDecayRate = .8, radMinFract = 0.01){
  rads <- c(maxRad)
  x <- maxRad
  repeat {
    if(radStepType == "exp"){
      x <- x * radDecayRate
    }else if(radStepType == "unif"){
      x <- x - radDecayRate
    }else{
      stop("Invalid step type, please choose between \"exp\", and \"unif\".")
    }
    rads <- c(rads, x)
    if(x < (radMinFract * maxRad)) {
      break
    }
  }
  rads <- c(rads, 0)
  rads
}

normalizeData <- function(normData){
  for(j in 1:ncol(normData)){
    tmpSD <- sd(normData[, j])
    tmpM <- mean(normData[, j])
    normData[, j] <- (normData[, j] - tmpM) / tmpSD
  }
  normData
}

getClusterDF <- function(data, clusterVars){
  cvData = numeric()
  for(i in clusterVars){
    if(is.factor(data[, i])){
      if(length(levels(data[, i])) == 2){
        cvData = cbind(cvData, as.integer(data[, i]) - 1)
      }
      else{
        cvData = cbind(cvData, model.matrix(~data[, i] - 1))
      }
    }
    else{
      cvData = cbind(cvData, data[, i])
    }
  }
  retList <- list()
  retList$clusterDF <- cvData
  retList
}

