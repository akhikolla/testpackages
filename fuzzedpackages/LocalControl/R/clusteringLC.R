#' @importFrom cluster diana agnes
#' @importFrom gss ssden dssden qssden
#' @importFrom grDevices dev.interactive
#' @importFrom lattice bwplot
#' @importFrom graphics abline barplot hist lines par plot points symbols title axis boxplot grid legend mtext
#' @importFrom stats aov as.formula binomial chisq.test complete.cases rnorm
#' @importFrom stats cutree density dist fitted.values glm hclust kmeans smooth.spline
#' @importFrom stats lm loess model.matrix na.omit pchisq prcomp predict sd var

#' @name LocalControlClassic
#' @title Local Control Classic
#'
#' @description LocalControlClassic was originally contained in the deprecated
#' CRAN package USPS, this function is a combination
#' of three of the original USPS functions, UPShclus, UPSaccum, and UPSnnltd.
#' This replicates the original implementation of the Local Control functionality in
#' Robert Obenchain's USPS package. Some of the features have been
#' removed due to deprecation of R packages distributed through CRAN.
#' For a given number of patient clusters in baseline X-covariate space, LocalControlClassic()
#' characterizes the distribution of Nearest Neighbor "Local Treatement Differences" (LTDs) on
#' a specified Y-outcome variable.
#'
#' @param data The data frame containing all baseline X covariates.
#' @param clusterVars List of names of X variable(s).
#' @param treatmentColName Name of treatment factor variable.
#' @param outcomeColName Name of outcome Y variable.
#' @param faclev Maximum number of different numerical values an outcome variable can assume
#' without automatically being converted into a "factor" variable; faclev=1 causes a binary
#' indicator to be treated as a continuous variable determining an average or proportion.
#' @param scedas Scedasticity assumption: "homo" or "hete".
#' @param clusterMethod Type of clustering method, defaults to "complete". Currently implemented methods:  "ward", "single", "complete" or "average".
#' @param clusterDist Distance type to use, defaults to "euclidean". Currently implemented: "euclidiean", "manhattan", "maximum", or "minkowski".
#' @param clusterCounts A vector containing different number of clusters in baseline X-covariate space which Local Control will iterate over.
#'
#' @return Returns a list containing several elements.
#'   \item{hiclus}{Name of clustering object created by UPShclus().}
#'   \item{dframe}{Name of data.frame containing X, t & Y variables.}
#'   \item{trtm}{Name of treatment factor variable.}
#'   \item{yvar}{Name of outcome Y variable.}
#'   \item{numclust}{Number of clusters requested.}
#'   \item{actclust}{Number of clusters actually produced.}
#'   \item{scedas}{Scedasticity assumption: "homo" or "hete"}
#'   \item{PStdif}{Character string describing the treatment difference.}
#'   \item{nnhbindf}{Vector containing cluster number for each patient. }
#'   \item{rawmean}{Unadjusted outcome mean by treatment group.}
#'   \item{rawvars}{Unadjusted outcome variance by treatment group.}
#'   \item{rawfreq}{Number of patients by treatment group.}
#'   \item{ratdif}{Unadjusted mean outcome difference between treatments.}
#'   \item{ratsde}{Standard error of unadjusted mean treatment difference.}
#'   \item{binmean}{Unadjusted mean outcome by cluster and treatment.}
#'   \item{binvars}{Unadjusted variance by cluster and treatment.}
#'   \item{binfreq}{Number of patients by bin and treatment.}
#'   \item{awbdif}{Across cluster average difference with cluster size weights.}
#'   \item{awbsde}{Standard error of awbdif.}
#'   \item{wwbdif}{Across cluster average difference, inverse variance weights.}
#'   \item{wwbsde}{Standard error of wwbdif.}
#'   \item{faclev}{Maximum number of different numerical values an outcome variable can assume without
#'    automatically being converted into a "factor" variable; faclev=1 causes a binary indicator to be
#'    treated as a continuous variable determining an average or proportion.}
#'   \item{youtype}{"continuous" => only next eight outputs; "factor" => only last three outputs.}
#'   \item{aovdiff}{ANOVA summary for treatment main effect only.}
#'   \item{form2}{Formula for outcome differences due to bins and to treatment nested within bins.}
#'   \item{bindiff}{ANOVA summary for treatment nested within cluster.}
#'   \item{sig2}{Estimate of error mean square in nested model.}
#'   \item{pbindif}{Unadjusted treatment difference by cluster.}
#'   \item{pbinsde}{Standard error of the unadjusted difference by cluster.}
#'   \item{pbinsiz}{Cluster radii measure: square root of total number of patients.}
#'   \item{symsiz}{Symbol size of largest possible Snowball in a UPSnnltd() plot with 1 cluster.}
#'   \item{factab}{Marginal table of counts by Y-factor level and treatment.}
#'   \item{cumchi}{Cumulative Chi-Square statistic for interaction in the three-way, nested table.}
#'   \item{cumdf}{Degrees of-Freedom for the Cumulative Chi-Squared.}
#'
#' @references
#' \itemize{
#'  \item Obenchain, RL. \emph{USPS package: Unsupervised and Supervised Propensity Scoring in R.} \url{https://cran.r-project.org/src/contrib/Archive/USPS/} 2005.
#'  \item Obenchain, RL. \emph{The ''Local Control'' Approach to Adjustment for Treatment Selection Bias and Confounding (illustrated with JMP Scripts)}. Observational Studies. Cary, NC: SAS Press. 2009.
#'  \item Obenchain RL. The local control approach using JMP. In: Faries D, Leon AC, Haro JM, Obenchain RL, eds. Analysis of Observational Health Care Data Using SAS. Cary, NC: SAS Institute; 2010:151-194.
#'  \item Obenchain RL, Young SS. Advancing statistical thinking in observational health care research. J Stat Theory Pract. 2013;7(2):456-506.
#'  \item Faries DE, Chen Y, Lipkovich I, Zagar A, Liu X, Obenchain RL. Local control for identifying subgroups of interest in observational research: persistence of treatment for major depressive disorder. Int J Methods Psychiatr Res. 2013;22(3):185-194.
#'  \item Lopiano KK, Obenchain RL, Young SS. Fair treatment comparisons in observational research. Stat Anal Data Min. 2014;7(5):376-384.
#'  \item Young SS, Obenchain RL, Lambert CG (2016) A problem of bias and response heterogeneity. In: Alan Moghissi A, Ross G (eds) Standing with giants: A collection of public health essays in memoriam to Dr. Elizabeth M. Whelan. American Council on Science and Health, New York, NY, pp 153-169.
#'  }
#'
#' @examples
#'  data(lindner)
#'
#'  cvars <- c("stent","height","female","diabetic","acutemi",
#'             "ejecfrac","ves1proc")
#'  numClusters <- c(1, 2, 10, 15, 20, 25, 30, 35, 40, 45, 50)
#'  results <- LocalControlClassic( data = lindner,
#'                                 clusterVars = cvars,
#'                                 treatmentColName = "abcix",
#'                                 outcomeColName = "cardbill",
#'                                 clusterCounts = numClusters)
#'  UPSLTDdist(results,ylim=c(-15000,15000))
#'
#' @export
LocalControlClassic <- function(data,
                                clusterVars,
                                treatmentColName,
                                outcomeColName,
                                faclev = 3,
                                scedas = "homo",
                                clusterMethod = "ward",
                                clusterDist = "euclidean",
                                clusterCounts = c(50,100,200))
{
  dfName = deparse(substitute(data))
  # We deparse here to get the name of the dataframe in the global environment.
  # deparse(substitute()) only works one level deep, so calling it in newLCCenv
  # means we would be naming it after the parameter in this function, 'data'.
  LCenv = newLCCenv(data, dfName)

  UPShclus(dframe = dfName, xvars = clusterVars, method = clusterMethod, metric = clusterDist, envir = LCenv)

  UPSaccum(dframe = dfName, trtm = treatmentColName, yvar = outcomeColName, faclev = faclev, scedas=scedas, envir = LCenv)

  for(cc in clusterCounts){
    UPSnnltd(LCenv, cc)
  }

  LCenv
}

newLCCenv <-function(data, dfName = ""){
  if(dfName == ""){
    dfName = deparse(substitute(data))
  }

  LCenv = new.env(parent = emptyenv())
  class(LCenv) = "LocalControlC"
  LCenv[[dfName]] = data
  LCenv
}

"UPlinint" <- function(q, xmin, n, x, w)
{
  if( q < 0 || q > 1 )
      stop("Desired quantiles must be between 0 and 1, inclusive.")
  cdf <- xmin
  while( TRUE ) {
      if( q == 0 ) break
      if( q == 1 ) {
          cdf <- x[n]
          break
          }
      j <- 0
      qlo <- 0
      qhi <- w[1] + 0.0000001
      while( q > qhi ) {
           j <- j+1
           if( j > (n-1) ) {
               j <- j-1
               break
               }
           qlo <- qhi
           qhi <- qhi + w[(j+1)]
           }
      if( j== 0 )
          cdf <- xmin + (q-qlo)*(x[1]-xmin)/(qhi-qlo)
      else
          cdf <- x[j] + (q-qlo)*(x[(j+1)]-x[j])/(qhi-qlo)
      break
      }
  cdf
}

#'
#'
#' @name UPShclus
#' @title Hierarchical Clustering of Patients on X-covariates for Unsupervised Propensiy Scoring
#' @description Derive a full, hierarchical clustering tree (dendrogram) for all patients (regardless
#'  of treatment received) using Mahalonobis between-patient distances computed from specified
#'  baseline X-covariate characteristics.
#'
#' @param dframe {Name of data.frame containing baseline X covariates.}
#' @param xvars {List of names of X variable(s).}
#' @param method {Hierarchical Clustering Method: "diana", "agnes" or "hclus".}
#' @param envir name of the working local control classic environment.
#' @param metric A valid distance metric for clustering.
#'
#' @details The first step in an Unsupervised Propensity Scoring alalysis is always
#'  to hierarchically cluster patients in baseline X-covariate space.  UPShclus uses
#'  a Mahalabobis metric and clustering methods from the R "cluster" library for this
#'  key initial step.
#'
#' @return An output list object of class UPShclus:
#' \itemize{
#'  \item{dframe}{Name of data.frame containing baseline X covariates.}
#'  \item{xvars}{List of names of X variable(s).}
#'  \item{method}{Hierarchical Clustering Method: "diana", "agnes" or "hclus".}
#'  \item{upshcl}{Hierarchical clustering object created by choice between three possible methods.}
#' }
#'
#' @references {
#'  Kaufman L, Rousseeuw PJ.  (1990) \bold{Finding Groups in Data.  An Introduction to
#'   Cluster Analysis}.  New York: John Wiley and Sons.
#'
#'  Kereiakes DJ, Obenchain RL, Barber BL, et al. (2000) Abciximab provides
#'   cost effective survival advantage in high volume interventional practice.
#'   \emph{Am Heart J} \bold{140}: 603-610.
#'
#'  Obenchain RL. (2004) Unsupervised Propensity Scoring: NN and IV Plots.
#'   \emph{Proceedings of the American Statistical Association (on CD)} 8 pages.
#'
#'  Obenchain RL. (2011) \bold{USPSinR.pdf}  USPS R-package vignette, 40 pages.
#'
#'  Rubin DB. (1980) Bias reduction using Mahalanobis metric matching.
#'  \emph{Biometrics} \bold{36}: 293-298.
#' }
#'
#' @seealso \code{\link{UPSaccum}}, \code{\link{UPSnnltd}} and \code{\link{UPSgraph}}.
#' @keywords cluster design
#' @author Bob Obenchain <wizbob@att.net>
#' @export
"UPShclus" <- function (envir, dframe, xvars, method, metric)
{
    if (missing(envir))
        stop("First argument to UPShclus must be an environment to work in.\nUse ?UPShclus for more information.")
    if (missing(dframe) || !exists(dframe, envir = envir))
        stop("Second argument to UPShclus must be the name of a Data Frame in the working environment.")
    if (missing(xvars))
        stop("Third argument to UPShclus must be a list of X variables which are columns of dframe.")
    if (missing(method))
        stop("Fourth argument to UPShclus must be a valid clustering method.")
    if (missing(metric))
        stop("Last argument to UPShclus must be a valid clustering metric.")

    data = envir[[dframe]]
    xpc <- prcomp(data[, xvars], retx = TRUE, center = TRUE, scale. = TRUE)

    if (method == "diana") {
        upshcl <- diana(dist(xpc$x/xpc$sdev), metric = "euclidean",stand = TRUE, keep.diss = FALSE, keep.data = FALSE)
    }
    else if (method == "agnes") {
        upshcl <- agnes(dist(xpc$x/xpc$sdev), diss = FALSE, metric = "euclidean",method = "complete", stand = TRUE, keep.diss = FALSE, keep.data = FALSE)
    }
    else if(method == "ward"){
        upshcl <- hclust(dist(xpc$x/xpc$sdev), method = "ward.D2")
    }
    else{
        upshcl <- hclust(dist(xpc$x/xpc$sdev), method = method)
    }

    HCLolist <- list(dframe = dframe, xvars = xvars, method = method, upshcl = upshcl)

    class(HCLolist) <- "UPShclus"

    assign("HCLolist", HCLolist, envir = envir)

    HCLolist
}


#' @name UPSaccum
#' @title Prepare for Accumulation of (Outcome,Treatment) Results in Unsupervised Propensity Scoring
#' @description Specify key result accumulation parameters: Treatment t-Factor, Outcome
#'  Y-variable, faclev setting, scedasticity assumption, and name of the UPSgraph() data
#'  accumulation object.
#'
#' @param envir name of the working local control classic environment.
#' @param dframe {Name of data.frame containing the X, t & Y variables.}
#' @param trtm {Name of treatment factor variable.}
#' @param yvar {Name of outcome Y variable.}
#' @param faclev {Maximum number of different numerical values an outcome variable can assume
#'   without automatically being converted into a "factor" variable; faclev=1 causes a binary
#'   indicator to be treated as a continuous variable determining an average or proportion.}
#' @param scedas {Scedasticity assumption: "homo" or "hete"}
#'
#' @details The second phase in an Unsupervised Propensity Scoring analysis is to prepare to
#'  accumulate results over a wide range of values for "Number of Clusters."  As the number of
#'  such clusters increases, individual clusters will tend to become smaller and smaller and,
#'  thus, more and more compact in covariate X-space.
#'
#' @return
#' \itemize{
#'  \item{hiclus}{Name of a diana, agnes or hclust object created by UPShclus().}
#'  \item{dframe}{Name of data.frame containing the X, t & Y variables.}
#'  \item{trtm}{Name of treatment factor variable.}
#'  \item{yvar}{Name of outcome Y variable.}
#'  \item{faclev}{Maximum number of different numerical values an outcome variable can assume
#'   without automatically being converted into a "factor" variable; faclev=1 causes a binary
#'   indicator to be treated as a continuous variable determining a proportion.}
#'  \item{scedas}{Scedasticity assumption: "homo" or "hete"}
#'  \item{accobj}{Name of the object for accumulation of I-plots to be ultimately displayed
#'     using UPSgraph().}
#'  \item{nnymax}{Maximum NN LTD Standard Error observed; Upper NN plot limit;
#'     initialized to zero.}
#'  \item{nnxmin}{Minimum NN LTD observed; Left NN plot limit; initialized to zero.}
#'  \item{nnxmax}{Maximum NN LTD observed; Right NN plot limit; initialized to zero.}
#'}
#'
#' @references {
#'  Obenchain RL. (2004) Unsupervised Propensity Scoring: NN and IV Plots.
#'  \emph{Proceedings of the American Statistical Association (on CD)} 8 pages.
#'
#'  Obenchain RL. (2011) \bold{USPSinR.pdf}  USPS R-package vignette, 40 pages.
#' }
#' @seealso \code{\link{UPSnnltd}}, \code{\link{UPSivadj}} and \code{\link{UPShclus}}.
#'
#' @keywords univar design
#' @author Bob Obenchain <wizbob@att.net>
#' @export
"UPSaccum" <- function (envir, dframe, trtm, yvar, faclev = 3, scedas = "homo")
{
  if (missing(envir))
        stop("First argument to UPSaccum must be an environment to work in.\nUse ?UPSaccum for more information.")
  if (missing(dframe) || !exists(dframe, envir = envir))
        stop("Second argument to UPSaccum must be the name of a Data Frame in the working environment.")
  if (missing(trtm))
    stop("Third argument to UPSaccum must name the Treatment factor variable in dframe.")
  if (missing(yvar))
    stop("Fourth argument to UPSaccum must name the Target variable in dframe.")
  if(!exists("HCLolist", envir = envir))
    stop("The environment provided has no HCLolist, please call UPShclus first.")

  data = envir[[dframe]]

  if (!is.element(yvar, dimnames(data)[[2]]))
    stop("Target Outcome or Covariate must be an existing Data Frame variable.")
  if (!is.element(trtm, dimnames(data)[[2]]))
    stop("Treatment factor must be an existing Data Frame variable.")
  if (length(table(data[, trtm])) != 2)
    stop("Treatment factor must assume exactly two different levels.")


  if (scedas != "homo")
    scedas <- "hete"

  nnymax <- 0
  nnxmin <- 0
  nnxmax <- 0

  pars <- cbind("HCLolist", dframe, trtm, yvar, faclev, scedas, "LCacc", nnymax, nnxmin, nnxmax)
  assign("UPSaccum.pars", pars, envir = envir)

  accdf <- as.data.frame(cbind("NONE", "HCLolist", dframe, trtm, yvar, 1, 0, 0))
  names(accdf) <- c("NNIV", "hicl", "dfrm", "trtm", "yvar", "bins", "tdif", "tdse")

  assign("LCacc", accdf, envir = envir)
}


#'
#' @name UPSnnltd
#' @title Nearest Neighbor Distribution of LTDs in Unsupervised Propensiy Scoring
#' @description For a given number of patient clusters in baseline X-covariate space, UPSnnltd()
#'  characterizes the distribution of Nearest Neighbor "Local Treatemnt Differences" (LTDs) on
#'  a specified Y-outcome variable.
#'
#' @param envir name of the working local control classic environment.
#' @param numclust {Number of clusters in baseline X-covariate space.}
#' @details Multiple calls to UPSnnltd(n) for varying numbers of clusters, n, are typically made
#'  after first invoking UPShclus() to hierarchically cluster patients in X-space and then
#'  invoking UPSaccum() to specify a Y outcome variable and a two-level treatment factor t.
#'  UPSnnltd(n) then determines the LTD Distribution corresponding to n clusters and,
#'  optionally, displays this distribution in a "Snowball" plot.
#'
#' @return An output list object of class UPSnnltd:
#' \itemize{
#'  \item{hiclus}{Name of clustering object created by UPShclus().}
#'  \item{dframe}{Name of data.frame containing X, t & Y variables.}
#'  \item{trtm}{Name of treatment factor variable.}
#'  \item{yvar}{Name of outcome Y variable.}
#'  \item{numclust}{Number of clusters requested.}
#'  \item{actclust}{Number of clusters actually produced.}
#'  \item{scedas}{Scedasticity assumption: "homo" or "hete"}
#'  \item{PStdif}{Character string describing the treatment difference.}
#'  \item{nnhbindf}{Vector containing cluster number for each patient. }
#'  \item{rawmean}{Unadjusted outcome mean by treatment group.}
#'  \item{rawvars}{Unadjusted outcome variance by treatment group.}
#'  \item{rawfreq}{Number of patients by treatment group.}
#'  \item{ratdif}{Unadjusted mean outcome difference between treatments.}
#'  \item{ratsde}{Standard error of unadjusted mean treatment difference.}
#'  \item{binmean}{Unadjusted mean outcome by cluster and treatment.}
#'  \item{binvars}{Unadjusted variance by cluster and treatment.}
#'  \item{binfreq}{Number of patients by bin and treatment.}
#'  \item{awbdif}{Across cluster average difference with cluster size weights.}
#'  \item{awbsde}{Standard error of awbdif.}
#'  \item{wwbdif}{Across cluster average difference, inverse variance weights.}
#'  \item{wwbsde}{Standard error of wwbdif.}
#'  \item{faclev}{Maximum number of different numerical values an outcome variable can assume without
#'     automatically being converted into a "factor" variable; faclev=1 causes a binary indicator to be
#'     treated as a continuous variable determining an average or proportion.}
#'  \item{youtype}{"contin"uous => only next eight outputs; "factor" => only last three outputs.}
#'  \item{aovdiff}{ANOVA summary for treatment main effect only.}
#'  \item{form2}{Formula for outcome differences due to bins and to treatment nested within bins.}
#'  \item{bindiff}{ANOVA summary for treatment nested within cluster.}
#'  \item{sig2}{Estimate of error mean square in nested model.}
#'  \item{pbindif}{Unadjusted treatment difference by cluster.}
#'  \item{pbinsde}{Standard error of the unadjusted difference by cluster.}
#'  \item{pbinsiz}{Cluster radii measure: square root of total number of patients.}
#'  \item{symsiz}{Symbol size of largest possible Snowball in a UPSnnltd() plot with 1 cluster.}
#'  \item{factab}{Marginal table of counts by Y-factor level and treatment.}
#'  \item{cumchi}{Cumulative Chi-Square statistic for interaction in the three-way, nested table.}
#'  \item{cumdf}{Degrees of-Freedom for the Cumulative Chi-Squared.}
#' }
#'
#' @references {
#'   Obenchain RL. (2004) Unsupervised Propensity Scoring: NN and IV Plots.
#'  \emph{Proceedings of the American Statistical Association (on CD)} 8 pages.
#'
#'  Obenchain RL. (2011) \bold{USPSinR.pdf}  USPS R-package vignette, 40 pages.
#'
#'  Rosenbaum PR, Rubin RB. (1983) The Central Role of the Propensity Score
#'  in Observational Studies for Causal Effects. \emph{Biometrika} \bold{70}:
#'  41--55.
#'
#'  Rubin DB. (1980) Bias reduction using Mahalanobis metric matching.
#'  \emph{Biometrics} \bold{36}: 293-298.
#' }
#' @seealso \code{\link{UPSivadj}}, \code{\link{UPSaccum}} and \code{\link{UPSgraph}}.
#' @keywords nonparametric
#' @author Bob Obenchain <wizbob@att.net>
#' @export
"UPSnnltd" <- function (envir, numclust)
{
  if (missing(envir))
    stop("First argument to UPSnnltd must be an environment to work in.\nUse ?UPSnnltd for more information.")
  if (missing(numclust))
    stop("The argument to UPSnnltd must specify a Number of Clusters.")

  if (numclust < 1)
    numclust <- 1

  UPSpars <-get("UPSaccum.pars", envir = envir)
  hiclus <- get( UPSpars[1],  envir = envir)
  dframe <- get( UPSpars[2],  envir = envir)
  trtm <- UPSpars[3]
  yvar <- UPSpars[4]
  faclev <- as.numeric(UPSpars[5])
  scedas <- UPSpars[6]
  UPSdf  <- get(UPSpars[7], envir = envir)
  hclbin <- "HclusBin"

  hbins <- as.data.frame(cutree(hiclus$upshcl, k = numclust))

  names(hbins) <- hclbin
  dfnnltd <- merge(dframe, hbins, by.x = "row.names", by.y = "row.names",all.x = TRUE)
  dfnnltd[, hclbin] <- as.factor(dfnnltd[, hclbin])
  bins <- length(table(dfnnltd[, hclbin]))

  NNolist <- list(hiclus = UPSpars[1],
                  dframe = UPSpars[2],
                  trtm = trtm,
                  yvar = yvar,
                  numclust = numclust,
                  actclust = bins,
                  scedas = scedas)

  PSmean <- as.matrix(tapply(dfnnltd[, yvar], dfnnltd[, trtm], na.rm = TRUE, mean))
  PStrtm <- paste(trtm, "=", dimnames(PSmean)[[1]])
  PStdif <- paste(PStrtm[2], "minus", PStrtm[1])
  PSvars <- as.matrix(tapply(dfnnltd[, yvar], dfnnltd[, trtm], na.rm = TRUE, var))
  PSfreq <- as.matrix(table(dfnnltd[, trtm]))
  PSvars <- PSvars/PSfreq


  NNolist <- c(NNolist, list(PStdif = PStdif,
                             nnhbindf = hbins,
                             rawmean = PSmean,
                             rawvars = PSvars,
                             rawfreq = PSfreq))


  RATdif <- sum(PSmean[2, ] - PSmean[1, ])
  RATsde <- sum(sqrt(PSvars[2, ] + PSvars[1, ]))


  NNolist <- c(NNolist, list(ratdif = RATdif,
                             ratsde = RATsde))


  PSmean <- tapply(dfnnltd[, yvar], list(dfnnltd[, hclbin],dfnnltd[, trtm]), na.rm = TRUE, mean)
  PSmean <- cbind(matrix(1:bins, bins, 1), PSmean)
  dimnames(PSmean) <- list(1:bins, c("BIN", PStrtm[1], PStrtm[2]))
  PSmean <- as.data.frame(PSmean)
  PSmean[, "BIN"] <- as.factor(PSmean[, "BIN"])
  PSvars <- as.matrix(tapply(dfnnltd[, yvar], list(dfnnltd[,hclbin], dfnnltd[, trtm]), na.rm = TRUE, var))
  PSfreq <- as.matrix(table(dfnnltd[, hclbin], dfnnltd[, trtm]))
  PSvars <- PSvars/PSfreq
  PSvars[PSvars == Inf] <- NA
  PSfreq <- cbind(matrix(1:bins, bins, 1), PSfreq)
  dimnames(PSfreq) <- list(1:bins, c("BIN", PStrtm[1], PStrtm[2]))
  PSfreq <- as.data.frame(PSfreq)
  PSfreq[, "BIN"] <- as.factor(PSfreq[, "BIN"])


  NNolist <- c(NNolist, list(binmean = PSmean,
                             binvars = PSvars,
                             binfreq = PSfreq))

  sig2 <- 0
  if (length(table(dfnnltd[, yvar])) > faclev) {

    youtype <- "contin"
    form <- as.formula(paste(yvar, "~", trtm))
    aovdiff <- aov(form, dfnnltd, na.action = na.omit)


    NNolist <- c(NNolist, list(youtype = youtype,
                               faclev = faclev,
                               aovdiff = invisible(summary(aovdiff))))


    if (bins > 1) {
      form <- as.formula(paste(yvar, "~", hclbin, "+", trtm, "%in%", hclbin))
      aovdiff <- aov(form, dfnnltd, na.action = na.omit)
      raov <- as.matrix(aovdiff$residuals)
      sig2 <- apply(raov, 2, na.rm = TRUE, var) * (length(raov) -1)/aovdiff$df.residual


      NNolist <- c(NNolist, list(form2 = form,
                                 bindiff = invisible(summary(aovdiff)),
                                 sig2 = sig2))


    }
    pbindif <- PSmean[, 3] - PSmean[, 2]
    nnxmin <- min(pbindif, na.rm = TRUE)
    nnxmax <- max(pbindif, na.rm = TRUE)
    if (scedas == "homo" && sig2 > 0) {
      PSvars <- sig2/as.matrix(table(dfnnltd[, hclbin], dfnnltd[, trtm]))
      PSvars[PSvars == Inf] <- NA
    }
    pbinsde <- sqrt(PSvars[, 2] + PSvars[, 1])
    nnymax <- max(pbinsde, na.rm = TRUE)
    pbinsiz <- sqrt(PSfreq[, 3] + PSfreq[, 2])
    symsiz <- length(dframe[, yvar])


    NNolist <- c(NNolist, list(pbindif = pbindif,
                               pbinsde = pbinsde,
                               pbinsiz = pbinsiz,
                               symsiz = symsiz))


  }
  else {
    youtype <- "factor"
    df3 <- as.data.frame(cbind(dfnnltd[, yvar], dfnnltd[,trtm]))
    df3[, 3] <- dfnnltd[, hclbin]
    df3 <- na.omit(df3)
    names(df3) <- c(yvar, trtm, hclbin)
    df3[, 1] <- as.factor(df3[, 1])
    tab <- table(df3[, 1], df3[, 2])


    NNolist <- c(NNolist, list(youtype = youtype,
                               faclev = faclev,
                               factab = tab))


    tab <- table(df3[, 1], df3[, 2], df3[, 3])
    cumchi <- 0
    cumdf <- 0
    for (i in 1:bins) {
      ht <- chisq.test(tab[, , i])
      if (!is.na(ht$statistic) && is.finite(ht$statistic)) {
        cumchi <- cumchi + ht$statistic
        cumdf <- cumdf + ht$parameters
      }
    }
    nnymax <- 0
    nnxmin <- 0
    nnxmax <- 0


    NNolist <- c(NNolist, list(cumchi = cumchi,cumdf = cumdf))


  }

  #Grab only the clusters which are informative
  infRows = which(PSfreq[,2] > 0 & PSfreq[,3] > 0)
  finfRows = which(PSfreq[,2] > 1 & PSfreq[,3] > 1)

  #use only info from informative clusters.
  INfreq = PSfreq[infRows,]
  INmean = PSmean[infRows,]
  INmean[,"LTD"] = INmean[,3] - INmean[,2]

  if(bins > 1){
    FINvars = data.frame(varT0 = PSvars[,1], varT1 = PSvars[,2])
    FINvars[,"BIN"] = seq(1,nrow(FINvars))
    FINvars[,"poolvar"] = FINvars$varT0 + FINvars$varT1
    FINvars = FINvars[finfRows,]

    FINmean = PSmean[finfRows,]
    FINfreq = PSfreq[finfRows,]

    awbdif <- sum(na.omit(((INmean[, 3] - INmean[, 2]) * (INfreq[,3] + INfreq[, 2])) / sum(PSfreq[, 3] + PSfreq[, 2])))
    #awbsde <- sqrt(sum(na.omit((FINvars[, 2] + FINvars[, 1]) * ((FINfreq[, 3] + FINfreq[, 2])^2))/(sum(FINfreq[, 3] + FINfreq[,2]))^2 ))

    awbsde = sqrt(sig2 * sum(((FINfreq[,3] + FINfreq[,2])^3 / (FINfreq[,2] * FINfreq[,3])) / sum(FINfreq[,3] + FINfreq[,2])^2))

  }else{
    awbdif <- sum(na.omit(((PSmean[, 3] - PSmean[, 2]) * (PSfreq[,3] + PSfreq[, 2])) / sum(PSfreq[, 3] + PSfreq[, 2])))
    awbsde <- sqrt(sum(na.omit((PSvars[, 2] + PSvars[, 1]) * (PSfreq[, 3] + PSfreq[, 2])^2))/ (sum(PSfreq[, 3] + PSfreq[,2]))^2)
  }

  wwbdif <- sum(na.omit( (PSmean[, 3] - PSmean[, 2])/(PSvars[, 2] + PSvars[, 1])/sum(1/na.omit(PSvars[, 2] + PSvars[,1]))))
  wwbsde <- sqrt(1/sum(1/na.omit(PSvars[, 2] + PSvars[, 1])))

  form <- as.formula(paste(yvar, "~", trtm))


  NNolist <- c(NNolist, list(awbdif = awbdif,
                             awbsde = awbsde,
                             wwbdif = wwbdif,
                             wwbsde = wwbsde,
                             form = form))


  accnew <- as.data.frame(cbind("NN", UPSpars[1], UPSpars[2],trtm, yvar, bins, awbdif, awbsde))
  names(accnew) <- c("NNIV", "hicl", "dfrm", "trtm", "yvar","bins", "tdif", "tdse")
  UPSdf <- as.data.frame(rbind(UPSdf, accnew))
  accnew <- as.data.frame(cbind("NW", UPSpars[1], UPSpars[2], trtm, yvar, bins, wwbdif, wwbsde))
  names(accnew) <- c("NNIV", "hicl", "dfrm", "trtm", "yvar","bins", "tdif", "tdse")
  UPSdf <- as.data.frame(rbind(UPSdf, accnew))
  assign(UPSpars[7], UPSdf, envir = envir)

  if (as.numeric(UPSpars[8]) < nnymax)
    envir$UPSaccum.pars[8] <- nnymax
  if (as.numeric(UPSpars[9]) > nnxmin)
    envir$UPSaccum.pars[9] <- nnxmin
  if (as.numeric(UPSpars[10]) < nnxmax)
    envir$UPSaccum.pars[10] <- nnxmax

  class(NNolist) <- "UPSnnltd"
  UPSSummary(NNolist, envir)

  envir[[paste0("UPSnnltd", numclust)]] = NNolist
}


#' @name UPSgraph
#' @title Display Sensitivity Analysis Graphic in Unsupervised Propensiy Scoring
#' @description Plot summary of results from multiple calls to UPSnnltd() and/or UPSivadj() after an
#'  initial setup call to UPSaccum().  The UPSgraph() plot displays any sensitivity of the LTD and
#'  LOA Distributions to choice of Number of Clusters in X-space.
#'
#' @param nncol {optional; string specifying color for display of the Mean of the LTD
#'   distribution when weighted by cluster size from any calls to UPSnnltd().}
#' @param nwcol {optional; string specifying color for display of the Mean of the LTD
#'   distribution when weighted inversely proportional to variance from any calls to UPSnnltd().}
#' @param ivcol {optional; string specifying color for display of the Difference in LOA
#'   predictions, at PS = 100\% minus that at PS = 0\%, from any calls to UPSivadj().}
#' @param envir name of the working local control classic environment.
#' @param ... Additional arguments to pass to the plotting function.
#'
#' @details The third phase of Unsupervised Propensity Scoring is a graphical Sensitivity
#'  Analysis that depicts how the Overall Means of the LTD and LOA distributions change with
#'  the number of clusters.
#'
#' @references   {
#' Kaufman L, Rousseeuw PJ. (1990) \bold{Finding Groups in Data.  An Introduction to
#'  Cluster Analysis}.  \emph{New York: John Wiley and Sons}.
#'
#'  Obenchain RL. (2004) Unsupervised Propensity Scoring: NN and IV Plots.
#'  \emph{Proceedings of the American Statistical Association (on CD)} 8 pages.
#'
#'  Obenchain RL. (2011) \bold{USPSinR.pdf}  USPS R-package vignette, 40 pages.
#'
#'  Rubin DB. (1980) Bias reduction using Mahalanobis metric matching.
#'  \emph{Biometrics} \bold{36}: 293-298.
#' }
#'
#' @seealso \code{\link{UPSnnltd}}, \code{\link{UPSivadj}} and \code{\link{UPSaccum}}.
#'
#' @author Bob Obenchain <wizbob@att.net>
#'
#' @keywords hplot
#' @export
"UPSgraph" <- function (envir, nncol = "red", nwcol = "green3", ivcol = "blue", ...)
{

    if (missing(envir))
        stop("First argument to UPSgraph must be an environment to work in.\nUse ?UPSgraph for more information.")

    UPSpars <- get("UPSaccum.pars", envir = envir)
    hiclus <- UPSpars[1]
    dframe <- UPSpars[2]
    trtm <- UPSpars[3]
    yvar <- UPSpars[4]
    faclev <- as.numeric(UPSpars[5])
    scedas <- UPSpars[6]
    UPSdf <- get(UPSpars[7], envir = envir)

    if (UPSdf$NNIV[1] == "NONE")
        UPSdf <- UPSdf[2:length(UPSdf$NNIV), ]

    UPSdf$bins <- as.numeric(as.character(UPSdf$bins))
    UPSdf$tdif <- as.numeric(as.character(UPSdf$tdif))
    UPSdf$tdse <- as.numeric(as.character(UPSdf$tdse))
    m <- order(UPSdf$hicl, UPSdf$dfrm, UPSdf$trtm, UPSdf$yvar, UPSdf$NNIV, UPSdf$bins)

    UPSdf <- UPSdf[m, ]
    assign(UPSpars[7], UPSdf, envir = envir)
    plot(UPSdf$bins, UPSdf$tdif, ylim = c(min(0, min(UPSdf$tdif -
        2 * UPSdf$tdse)), max(0, max(UPSdf$tdif + 2 * UPSdf$tdse))),
        log = "x", ann = FALSE, type = "n")

    symbols(UPSdf$bins, UPSdf$tdif, boxplots = cbind(0, 0, 2 *
        UPSdf$tdse, 2 * UPSdf$tdse, 0), inches = FALSE, add = TRUE)
    x <- UPSdf[UPSdf$NNIV == "NW", ]
    points(x$bins, x$tdif, pch = 21, col = nwcol)
    x <- UPSdf[UPSdf$NNIV == "IV", ]
    points(x$bins, x$tdif, pch = 21, col = ivcol)
    x <- UPSdf[UPSdf$NNIV == "NN", ]
    points(x$bins, x$tdif, pch = 21, col = nncol)
    abline(h = UPSdf$tdif[1], lty = 2, col = nncol)
    title(
          main = "Unsupervised LTD Distributiun Sensitivity",
          xlab = "Number of Clusters", ylab = "Mean LTD +/-2 Sigma LTD")
}

#' @author Bob Obenchain <wizbob@att.net>
#' @export
"plot.UPSnnltd" <- function (x, envir, pballs = TRUE, nnplot = "snob", nnalpha = 1.4, ...)
{
  if (missing(envir)){
    stop("You must provide a valid local control environment to plot this object.")
  }

  if (x$youtype == "contin") {
    UPSpars <- get("UPSaccum.pars", envir = envir)
    if (nnplot != "snob" && nnplot != "dens" && nnplot != "cdf" && nnplot != "seq")
      nnplot <- "all"
    opar <- par(no.readonly = TRUE, ask = dev.interactive(orNone = TRUE))
    on.exit(par(opar))
    if (nnplot == "all")
      par(mfrow=c(2,2))
    else
      par(mfrow=c(1,1))
    if (nnplot == "all" || nnplot == "dens" || nnplot == "cdf" || nnplot == "seq") {
      nm <- complete.cases(x$pbinsde)
      nn <- as.vector(x$pbinsiz[nm]^2)
      xx <- as.vector(x$pbindif[nm])
      ww <- as.vector(x$pbinsde[nm]^-2)
      ww <- round( sum(nn)* ww / sum(ww) )
      ltdfit <- ssden(~xx, weights = ww, alpha = nnalpha, maxiter = 30)
      xp <- seq(as.numeric(UPSpars[9]),as.numeric(UPSpars[10]), len=101)
    }
    if (nnplot == "all" || nnplot == "seq" || nnplot == "snob") {
      inchsz <- 2 * max(x$pbinsiz)^2/x$symsiz
      plot(x$pbindif, x$pbinsde, ann = FALSE, type = "n",
           ylim = c(0, as.numeric(UPSpars[8])),
           xlim = c(as.numeric(UPSpars[9]), as.numeric(UPSpars[10])))
      symbols(x$pbindif, x$pbinsde, circles = x$pbinsiz, inches = inchsz,
              add = TRUE)
      par(lty = 2, col = "red")
      abline(v = x$awbdif)
      par(lty = 3, col = "green3")
      abline(v = x$wwbdif)
      par(lty = 1, col = "black")
      abline(v = 0)
      if (pballs) {
        symbols(x$awbdif, x$awbsde, circles = x$awbsde, inches = 0.25,
                bg = "red", add = TRUE)
        symbols(x$wwbdif, x$wwbsde, circles = x$wwbsde, inches = 0.25,
                bg = "green3", add = TRUE)
      }
      title(main = paste("Unsupervised NN/LTD Snow Balls for",
                         length(na.omit(x$pbinsde)), "of", x$numclust, "Clusters"),
            xlab = "Local Treatment Difference (LTD)", ylab = "LTD Standard Error")
      if (nnplot == "seq") {
        cat("\nPress the Enter key to view the DENS nnplot...")
        scan()
      }
    }
    if (nnplot == "all" || nnplot == "seq" || nnplot == "dens") {
      plot(xp, dssden(ltdfit, xp), ann=FALSE, type="l", lty=1, col="blue")
      par(lty=2, col="red")
      abline(v=x$awbdif)
      par(lty=3, col="green3")
      abline(v=x$wwbdif)
      par(lty=1, col="black")
      abline(v=0)
      title(main=paste("NN/LTD Weighted Density for", length(na.omit(x$pbinsde)),
                       "Informative Clusters"), xlab = "Local Treatment Difference (LTD)",
            ylab = "gss Probability Density")
      if (nnplot == "seq") {
        cat("\nPress the Enter key to view the CDF nnplot...")
        scan()
      }
    }
    if (nnplot == "all" || nnplot == "seq" || nnplot == "cdf") {
      qq <- seq(0,1,len=51)
      cdf <- qssden(ltdfit, qq)
      plot(cdf, qq, ann=FALSE, type="l", lty=1, col="blue")
      par(lty=2, col="red")
      abline(v=x$awbdif)
      par(lty=3, col="green3")
      abline(v=x$wwbdif)
      par(lty=1, col="black")
      abline(v=0)
      title(main=paste("NN/LTD Weighted CDF for", length(na.omit(x$pbinsde)),
                       "Informative Clusters"), xlab = "Local Treatment Difference (LTD)",
            ylab = "gss Probability Less Than")
    }
  }
}

#' @author Bob Obenchain <wizbob@att.net>
#' @export
"print.UPSnnltd" <- function (x, ...)
{
  cat("\nUPSnnltd Object: Nearest Neighbor Local Treatment Differences\n")
  cat("Hierarchical Clustering object:", x$hiclus, "\n")
  cat("Data Frame input:", x$dframe, "\n")
  cat("Outcome variable:", x$yvar, "\n")
  cat("Treatment difference:", x$PStdif, "\n")
  cat("Scedasticity assumption:", x$scedas, "\n")

  if (length(x$sig2) == 1)
    cat("Homoscedastic Sigma:", sqrt(x$sig2), "\n")

  cat("Maximum Number of Clusters:", x$numclust, "\n")

  if (x$actclust != x$numclust)
    cat("Actual Number of Clusters =", x$actclust, "\n")

  cat("Number of Informative Clusters =", length(na.omit(x$pbinsde)),"\n")
  cat("\n    Overall Raw Average Treatment Difference =", x$ratdif)
  cat("\n    Standard Deviation of this Raw Difference =", x$ratsde)
  cat("\n    Average Within Cluster Treatment Difference =", x$awbdif)
  cat("\n    Standard Deviation of Average Within Cluster Difference =", x$awbsde)
  cat("\n    Inverse Variance Weighted Difference =", x$wwbdif)
  cat("\n    Standard Deviation of this Weighted Difference =", x$wwbsde)
  cat("\n\nTest for Raw / Unadjusted Treatment Difference:\n")
  print(x$form)
  if (x$youtype == "contin") {
    print(x$aovdiff)
    if (x$actclust > 1) {
      cat("\nTest for Treatment Difference within Clusters:\n")
      print(x$form2)
      print(x$bindiff)
    }
  }
  else {
    print(x$tab)
    print(chisq.test(x$tab))
    cat("\nTest for Treatment Difference within Clusters")
    cat("\n    Cumulative ChiSquare = ", round(x$cumchi, digits = 2))
    cat("\n    Cumulative df = ", x$cumdf)
    cat("\n    Cumulative p.value = ", round(1 - pchisq(x$cumchi, x$cumdf), digits = 4), "\n\n")
  }
}

"SPSsmoot" <-
function (envir, dframe, trtm, pscr, yvar, faclev = 3, df = 5, spar = NULL, cv = FALSE, penalty = 1)
{
    if (missing(envir))
      stop("First argument to SPSsmoot must be the name of an existing Data Frame.")
    if (missing(dframe) || !exists(dframe, envir = envir))
      stop("Second argument to SPSsmoot must be the name of a Data Frame in the working environment.")
    if (missing(trtm))
      stop("Third argument to SPSsmoot must name the Treatment factor in dframe.")
    if (missing(yvar))
      stop("Fifth argument to SPSsmoot must name the Target variable in dframe.")

    data = envir[[dframe]]
    if (!is.element(yvar, dimnames(data)[[2]]))
      stop("Target Outcome or Covariate must be an existing Data Frame variable.")
    if (!is.element(trtm, dimnames(data)[[2]]))
      stop("Treatment factor must be an existing Data Frame variable.")
    if (length(table(data[, trtm])) != 2)
      stop("Treatment factor must assume exactly two different levels.")

    if (missing(pscr))
        stop("Fourth argument to SPSsmoot must be fitted Propensity Scores.")

    if (!is.element(pscr, dimnames(data)[[2]]))
        stop("Propensity Score must be an existing Data Frame variable.")

    dframe = data

    SPSolist <- list(dframe = deparse(substitute(dframe)), trtm = trtm,
        pscr = pscr, yvar = yvar)
    PSmean <- as.matrix(tapply(dframe[, yvar], dframe[, trtm],
        na.rm = TRUE, mean))
    PSvars <- as.matrix(tapply(dframe[, yvar], dframe[, trtm],
        na.rm = TRUE, var))
    PSfreq <- as.matrix(table(dframe[, trtm]))
    PSvars <- PSvars/PSfreq
    form <- as.formula(paste(yvar, "~", trtm))
    aovdiff <- aov(form, dframe, na.action = na.omit)
    SPSolist <- c(SPSolist, list(rawmean = PSmean, rawvars = PSvars,
        rawfreq = PSfreq, form = form, aovdiff = aovdiff))
    ssfit <- as.data.frame(cbind(dframe[, pscr], dframe[, yvar],
        as.numeric(as.character(dframe[, trtm])), 0, 0))
    ssfit <- na.omit(ssfit)
    names(ssfit) <- c("PS", "YVAR", "TRTM", "FIT", "SEP")
    ssgrid <- as.data.frame(cbind(seq(0.005, 0.995, length = 100),
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
    names(ssgrid) <- c("PS", "F0", "S0", "C0", "F1", "S1", "C1",
        "DIF", "SED", "HST", "DEN")
    spsub0 <- subset.data.frame(ssfit, ssfit$TRTM == 0)
    spsub1 <- subset.data.frame(ssfit, ssfit$TRTM == 1)
    ssobj0 <- smooth.spline(spsub0$PS, spsub0$YVAR, df = df,
        spar = spar, cv = cv, penalty = penalty)
    res <- (ssobj0$yin - ssobj0$y)/(1 - ssobj0$lev)
    sigma <- sqrt(var(res))
    spsub0 <- as.data.frame(cbind(ssobj0$x, ssobj0$yin, 0, ssobj0$y,
        sigma * sqrt(ssobj0$lev)))
    names(spsub0) <- c("PS", "YAVG", "TRTM", "FIT", "SEP")
    ssgrid$F0 <- predict(ssobj0, ssgrid$PS)$y
    ssobj0 <- smooth.spline(spsub0$PS, spsub0$SEP, df = df, spar = spar,
        cv = cv, penalty = penalty)
    ssgrid$S0 <- predict(ssobj0, ssgrid$PS)$y
    ssobj1 <- smooth.spline(spsub1$PS, spsub1$YVAR, df = df,
        spar = spar, cv = cv, penalty = penalty)
    res <- (ssobj1$yin - ssobj1$y)/(1 - ssobj1$lev)
    sigma <- sqrt(var(res))
    spsub1 <- as.data.frame(cbind(ssobj1$x, ssobj1$yin, 1, ssobj1$y,
        sigma * sqrt(ssobj1$lev)))
    names(spsub1) <- c("PS", "YAVG", "TRTM", "FIT", "SEP")
    ssgrid$F1 <- predict(ssobj1, ssgrid$PS)$y
    ssobj1 <- smooth.spline(spsub1$PS, spsub1$SEP, df = df, spar = spar,
        cv = cv, penalty = penalty)
    ssgrid$S1 <- predict(ssobj1, ssgrid$PS)$y
    ssgrid$DIF <- ssgrid$F1 - ssgrid$F0
    ssgrid$SED <- sqrt(ssgrid$S0^2 + ssgrid$S1^2)
    ssgrid$C0 <- hist(spsub0$PS, breaks = seq(0, 1, length = 101),
        plot = FALSE, probability = FALSE)$counts
    ssgrid$C1 <- hist(spsub1$PS, breaks = seq(0, 1, length = 101),
        plot = FALSE, probability = FALSE)$counts
    ssgrid$HST <- hist(ssfit$PS, breaks = seq(0, 1, length = 101),
        plot = FALSE, probability = TRUE)$counts
    ssgrid$DEN <- density(ssfit$PS, bw = "nrd0", kernel = "gaussian",
        n = 100, from = 0.005, to = 0.995)$y
    grid <- na.omit(ssgrid)
    if (sum(grid$DEN) > 0)
        grid$DEN <- grid$DEN/sum(grid$DEN)
    PSmean <- sum(grid$DIF * grid$DEN)
    PSster <- sum(grid$SED * grid$DEN)
    SPSolist <- c(SPSolist, list(ssgrid = ssgrid, spsub0 = spsub0,
        spsub1 = spsub1, df = df, sptdif = PSmean, sptsde = PSster))
    class(SPSolist) <- "SPSsmoot"
    SPSolist
}

#'
#' @name UPSivadj
#' @title Instrumental Variable LATE Linear Fitting in Unsupervised Propensiy Scoring
#' @description For a given number of patient clusters in baseline X-covariate space and a
#'  specified Y-outcome variable, linearly smooth the distribution of Local Average
#'  Treatment Effects (LATEs) plotted versus Within-Cluster Treatment Selection (PS)
#'  Percentages.
#'
#' @param envir name of the working local control classic environment.
#' @param numclust {Number of clusters in baseline X-covariate space.}
#'
#' @details Multiple calls to UPSivadj(n) for varying numbers of clusters n are made after first
#'  invoking UPShclus() to hierarchically cluster patients in X-space and then invoking UPSaccum()
#'  to specify a Y outcome variable and a two-level treatment factor t.  UPSivadj(n) linearly
#'  smoothes the LATE distribution when plotted versus within cluster propensity score percentages.
#'
#' @return An output list object of class UPSivadj:
#'  \itemize{
#'   \item{hiclus}{Name of clustering object created by UPShclus().}
#'   \item{dframe}{Name of data.frame containing X, t & Y variables.}
#'   \item{trtm}{Name of treatment factor variable.}
#'   \item{yvar}{Name of outcome Y variable.}
#'   \item{numclust}{Number of clusters requested.}
#'   \item{actclust}{Number of clusters actually produced.}
#'   \item{scedas}{Scedasticity assumption: "homo" or "hete"}
#'   \item{PStdif}{Character string describing the treatment difference.}
#'   \item{ivhbindf}{Vector containing cluster number for each patient. }
#'   \item{rawmean}{Unadjusted outcome mean by treatment group.}
#'   \item{rawvars}{Unadjusted outcome variance by treatment group.}
#'   \item{rawfreq}{Number of patients by treatment group.}
#'   \item{ratdif}{Unadjusted mean outcome difference between treatments.}
#'   \item{ratsde}{Standard error of unadjusted mean treatment difference.}
#'   \item{binmean}{Unadjusted mean outcome by cluster and treatment.}
#'   \item{binfreq}{Number of patients by bin and treatment.}
#'   \item{faclev}{Maximum number of different numerical values an outcome variable can assume without
#'    automatically being converted into a "factor" variable; faclev=1 causes a binary indicator to be
#'    treated as a continuous variable determining an average or proportion.}
#'   \item{youtype}{"contin"uous => next eleven outputs; "factor" => no additional output items.}
#'   \item{pbinout}{LATE regardless of treatment by cluster.}
#'   \item{pbinpsp}{Within-Cluster Treatment Percentage = non-parametric Propensity Score.}
#'   \item{pbinsiz}{Cluster radii measure: square root of total number of patients.}
#'   \item{symsiz}{Symbol size of largest possible Snowball in a UPSivadj() plot with 1 cluster.}
#'   \item{ivfit}{lm() output for linear smooth across clusters.}
#'   \item{ivtzero}{Predicted outcome at PS percentage zero.}
#'   \item{ivtxsde}{Standard deviation of outcome prediction at PS percentage zero.}
#'   \item{ivtdiff}{Predicted outcome difference for PS percentage 100 minus that at zero.}
#'   \item{ivtdsde}{Standard deviation of outcome difference.}
#'   \item{ivt100p}{Predicted outcome at PS percentage 100.}
#'   \item{ivt1pse}{Standard deviation of outcome prediction at PS percentage 100.}
#' }
#' @references {
#'  Imbens GW, Angrist JD. (1994) Identification and Estimation  of
#'   Local Average Treatment Effects (LATEs). \emph{Econometrica} \bold{62}: 467-475.
#'
#'  Obenchain RL. (2004) Unsupervised Propensity Scoring: NN and IV Plots.
#'   \emph{Proceedings of the American Statistical Association (on CD)} 8 pages.
#'
#'  Obenchain RL. (2011) \bold{USPSinR.pdf}  USPS R-package vignette, 40 pages.-
#'
#'  McClellan M, McNeil BJ, Newhouse JP. (1994) Does More Intensive Treatment of
#'   Myocardial Infarction in the Elderly Reduce Mortality?: Analysis Using Instrumental
#'   Variables. \emph{JAMA} \bold{272}: 859-866.
#'
#'  Rosenbaum PR, Rubin RB. (1983) The Central Role of the Propensity Score
#'   in Observational Studies for Causal Effects. \emph{Biometrika} \bold{70}: 41-55.
#' }
#'
#' @seealso \code{\link{UPSnnltd}}, \code{\link{UPSaccum}} and \code{\link{UPSgraph}}.
#' @keywords hplot nonparametric
#' @author Bob Obenchain <wizbob@att.net>
#' @export
"UPSivadj" <-function (envir, numclust)
{

    if(missing(envir)){
      stop("You must provide a valid environment to UPSivadj")
    }

    if (missing(numclust))
        stop("You must specify a number of clusters to create.")

    if (numclust < 3)
        numclust <- 3

    UPSpars <- get("UPSaccum.pars", envir = envir)
    hiclus <- get(UPSpars[1], envir = envir)
    dframe <- get(UPSpars[2], envir = envir)
    trtm <- UPSpars[3]
    yvar <- UPSpars[4]
    faclev <- as.numeric(UPSpars[5])
    scedas <- UPSpars[6]
    UPSdf <- get(UPSpars[7], envir = envir)
    hclbin <- "HclusBin"
    hbins <- as.data.frame(cutree(hiclus$upshcl, k = numclust))
    names(hbins) <- hclbin
    dfivadj <- merge(dframe, hbins, by.x = "row.names", by.y = "row.names", all.x = TRUE)
    dfivadj[, hclbin] <- as.factor(dfivadj[, hclbin])
    bins <- length(table(dfivadj[, hclbin]))
    IVolist <- list(hiclus = UPSpars[1], dframe = UPSpars[2],
                    trtm = trtm, yvar = yvar, numclust = numclust, actclust = bins,
                    scedas = scedas)
    PSmean <- as.matrix(tapply(dfivadj[, yvar], dfivadj[, trtm], na.rm = TRUE, mean))
    PStrtm <- paste(trtm, "=", dimnames(PSmean)[[1]])
    PSfreq <- as.matrix(table(dfivadj[, trtm]))
    IVolist <- c(IVolist, list(ivhbindf = hbins, rawmean = PSmean, rawfreq = PSfreq))
    PSmean <- tapply(dfivadj[, yvar], dfivadj[, hclbin], na.rm = TRUE,mean)
    PSmean <- cbind(matrix(1:bins, bins, 1), PSmean)
    dimnames(PSmean) <- list(1:bins, c("BIN", "OUT"))
    PSmean <- as.data.frame(PSmean)
    PSmean[, "BIN"] <- as.factor(PSmean[, "BIN"])
    PSfreq <- as.matrix(table(dfivadj[, hclbin], dfivadj[, trtm]))
    PSfreq <- cbind(matrix(1:bins, bins, 1), PSfreq)
    dimnames(PSfreq) <- list(1:bins, c("BIN", PStrtm[1], PStrtm[2]))
    PSfreq <- as.data.frame(PSfreq)
    PSfreq[, "BIN"] <- as.factor(PSfreq[, "BIN"])
    IVolist <- c(IVolist, list(binmean = PSmean, binfreq = PSfreq))
    if (length(table(dfivadj[, yvar])) > faclev) {
        youtype <- "contin"
        IVolist <- c(IVolist, list(youtype = youtype, faclev = faclev))
        pbinout <- PSmean[, 2]
        pbinpsp <- PSfreq[, 3] + PSfreq[, 2]
        pbinpsp <- 100 * PSfreq[, 3]/pbinpsp
        pbinsiz <- sqrt(PSfreq[, 3] + PSfreq[, 2])
        symsiz <- length(dframe[, yvar])
        IVolist <- c(IVolist, list(pbinout = pbinout, pbinpsp = pbinpsp,
                     pbinsiz = pbinsiz, symsiz = symsiz))
        if (bins > 1) {
            ivfit <- lm(pbinout ~ pbinpsp, weights = pbinsiz^2)
            ivsum <- summary(ivfit)
            ivtzero <- ivsum$coefficients[[1]]
            ivtzsde <- ivsum$coefficients[[3]]
            ivtdiff <- 100 * ivsum$coefficients[[2]]
            ivtdsde <- 100 * ivsum$coefficients[[4]]
            ivsum <- predict(ivfit, data.frame(pbinpsp = 100), se.fit = TRUE)
            ivt100p <- ivsum$fit
            ivt1pse <- ivsum$se.fit
            IVolist <- c(IVolist, list(ivfit = ivfit, ivtzero = ivtzero,
                         ivtzsde = ivtzsde, ivtdiff = ivtdiff, ivtdsde = ivtdsde,
                         ivt100p = ivt100p, ivt1pse = ivt1pse))
        }
    }
    else {
        youtype <- "factor"
        IVolist <- c(IVolist, list(youtype = youtype, faclev = faclev))
    }
    accnew <- as.data.frame(cbind("IV", UPSpars[1], UPSpars[2],
                                  trtm, yvar, bins, ivtdiff, ivtdsde))
    names(accnew) <- c("NNIV", "hicl", "dfrm", "trtm", "yvar",
                       "bins", "tdif", "tdse")
    UPSdf <- as.data.frame(rbind(UPSdf, accnew))
    assign(UPSpars[7], UPSdf, envir = envir)
    class(IVolist) <- "UPSivadj"
    IVolist
}



#'
#' @name UPSaltdd
#' @title Artificial Distribution of LTDs from Random Clusters
#' @description For a given number of clusters, UPSaltdd() characterizes the potentially biased
#'  distribution of "Local Treatment Differences" (LTDs) in a continuous outcome y-variable between
#'  two treatment groups due to Random Clusterings.  When the NNobj argument is not NA and specifies
#'  an existing UPSnnltd() object, UPSaltdd() also computes a smoothed CDF for the NN/LTD
#'  distribution for direct comparison with the Artificial LTD distribution.
#'
#' @param dframe {Name of data.frame containing a treatment-factor and the outcome y-variable.}
#' @param trtm {Name of treatment factor variable with two levels.}
#' @param yvar {Name of continuous outcome variable.}
#' @param faclev {Maximum number of different numerical values an outcome variable can assume
#'   without automatically being converted into a "factor" variable; faclev=1 causes a binary
#'   indicator to be treated as a continuous variable determining an average or proportion.}
#' @param scedas {Scedasticity assumption: "homo" or "hete"}
#' @param NNobj {Name of an existing UPSnnltd object or NA.}
#' @param clus {Number of Random Clusters requested per Replication; ignored when NNobj is not NA.}
#' @param reps {Number of overall Replications, each with the same number of requested clusters.}
#' @param seed {Seed for Monte Carlo random number generator.}
#' @param envir name of the working local control classic environment.
#'
#' @details Multiple calls to UPSaltdd() for different UPSnnltd objects or different numbers of
#'  clusters are typically made after first invoking UPSgraph().
#'
#' @return
#' \itemize{
#'   \item{dframe}{Name of data.frame containing X, t & Y variables.}
#'   \item{trtm}{Name of treatment factor variable.}
#'   \item{yvar}{Name of outcome Y variable.}
#'   \item{faclev}{Maximum number of different numerical values an outcome variable can assume without
#'    automatically being converted into a "factor" variable; faclev=1 causes a binary indicator to be
#'    treated as a continuous variable determining an average or proportion.}
#'   \item{scedas}{Scedasticity assumption: "homo" or "hete"}
#'   \item{NNobj}{Name of an existing UPSnnltd object or NA.}
#'   \item{clus}{Number of Random Clusters requested per Replication.}
#'   \item{reps}{Number of overall Replications, each with the same number of requested clusters.}
#'   \item{pats}{Number of patients with no NAs in their yvar outcome and trtm factor.}
#'   \item{seed}{Seed for Monte Carlo random number generator.}
#'   \item{altdd}{Matrix of LTDs and relative weights from artificial clusters.}
#'   \item{alxmin}{Minimum artificial LTD value.}
#'   \item{alxmax}{Maximum artificial LTD value.}
#'   \item{alymax}{Maximum weight among artificial LTDs.}
#'   \item{altdcdf}{Vector of artificial LTD x-coordinates for smoothed CDF.}
#'   \item{qq}{Vector of equally spaced CDF values from 0.0 to 1.0.}
#'   \item{nnltdd}{Optional matrix of relevant NN/LTDs and relative weights.}
#'   \item{nnlxmin}{Optional minimum NN/LTD value.}
#'   \item{nnlxmax}{Optional maximum NN/LTD value.}
#'   \item{nnlymax}{Optional maximum weight among NN/LTDs.}
#'   \item{nnltdcdf}{Optional vector of NN/LTD x-coordinates for smoothed CDF.}
#'   \item{nq}{Optional vector of equally spaced CDF values from 0.0 to 1.0.}
#'
#' }
#'
#' @seealso \code{\link{UPSnnltd}}, \code{\link{UPSaccum}} and \code{\link{UPSgraph}}.
#'
#' @references {
#'   Obenchain RL. (2004) Unsupervised Propensity Scoring: NN and IV Plots.
#'  \emph{Proceedings of the American Statistical Association (on CD)} 8 pages.
#'
#'  Obenchain RL. (2011) \bold{USPSinR.pdf}  USPS R-package vignette, 40 pages.
#'
#'  Rosenbaum PR, Rubin RB. (1983) The Central Role of the Propensity Score
#'  in Observational Studies for Causal Effects. \emph{Biometrika} \bold{70}:
#'  41-55.
#'
#'  Rubin DB. (1980) Bias reduction using Mahalanobis metric matching.
#'  \emph{Biometrics} \bold{36}: 293-298.
#'
#' }
#'
#' @author Bob Obenchain <wizbob@att.net>
#' @export
"UPSaltdd" <- function(envir, dframe, trtm, yvar, faclev=3, scedas="homo", NNobj=NA, clus=50, reps=10, seed=12345)
{
    # Compute the Artificial LTD Distribution for random patient clusterings  ...i.e. include
    # the least or less relevant comparisons as well as those neurtal, more or most relevant.
    if (missing(envir))
      stop("First argument to UPSaltdd must be the name of an existing Data Frame.")
    if (missing(dframe) || !exists(dframe, envir = envir))
      stop("Second argument to UPSaltdd must be the name of a Data Frame in the working environment.")
    if (missing(trtm))
      stop("Third argument to UPSaltdd must name the Treatment factor.")
    if (missing(yvar))
      stop("Fourth argument to UPSaltdd must name the Target variable.")

    data = envir[[dframe]]
    if (!is.element(yvar, dimnames(data)[[2]]))
      stop("Target Outcome or Covariate must be an existing Data Frame variable.")
    if (!is.element(trtm, dimnames(data)[[2]]))
      stop("Treatment factor must be an existing Data Frame variable.")
    if (length(table(data[, trtm])) != 2)
      stop("Treatment factor must assume exactly two different levels.")


    if(scedas!="hete") scedas <- "homo" # variances either homoscedastic or heteroscedastic

    NNobjnam <- deparse(substitute(NNobj))
    if(NNobjnam!="NA"&&!inherits(NNobj,"UPSnnltd"))
        stop("The NNobj argument to UPSaltdd must either be NA or an existing UPSnnltd object.")

    if(NNobjnam!="NA")
       clus <- NNobj$numclust
    if(clus < 2)
       clus <- 2

    ytvars <- c(yvar,trtm)
    ytdata <- data[,ytvars]
    if(length(table(ytdata[,1])) <= faclev)
        stop("To be considered continuous, Outcome Variable must assume more than faclev values.")

    pats <- length(ytdata[,1])
    altdg <- "Agroups"
    ALmean <- as.matrix(tapply(ytdata[,1], ytdata[,2], na.rm = TRUE, mean))
    ALtrtm <- paste(trtm, "=", dimnames(ALmean)[[1]])
    # Start forming UPSaltdd output list...

    ALolist <- list(dframe=dframe, trtm=trtm, yvar=yvar, faclev=faclev, scedas=scedas,
                    NNobj=NNobjnam, clus=clus, reps=reps, pats=pats, seed=seed)
    set.seed(seed)   # Set seed for Monte Carlo pseudo random sequence...

    for(i in 1:reps) {

       xrand <- as.data.frame(cbind(rnorm(pats),rnorm(pats)))
       crand <- kmeans(xrand, clus)
       crand <- as.data.frame(crand$cluster)

       ALmean <- as.matrix(tapply(ytdata[,1], list(crand[,1], ytdata[,2]), na.rm = TRUE, mean))
       ALmean <- cbind(matrix(1:clus, clus, 1), ALmean)
       dimnames(ALmean) <- list(1:clus, c("ALC", ALtrtm[1], ALtrtm[2]))
       ALmean <- as.data.frame(ALmean)
       ALmean[,"ALC"] <- as.factor(ALmean[,"ALC"])
       ALvars <- as.matrix(tapply(ytdata[,1], list(crand[,1], ytdata[,2]), na.rm = TRUE, var))
       ALfreq <- as.matrix(table(crand[,1], ytdata[,2]))
       ALvars <- ALvars / ALfreq  # redefined below when scedas=="homo"
       ALvars[ALvars==Inf] <- NA

       altdif <- na.omit( ALmean[,3] - ALmean[,2] )
       if( scedas=="homo" ) {
           ALvars <- 1 / as.matrix(table(crand[,1], ytdata[,2]))
           ALvars[ALvars==Inf] <- NA
           }
       altdwt <- ALvars[,2] + ALvars[,1]
       altdwt <- na.omit( 1 / altdwt )
       altdwt <- altdwt / sum(altdwt)
       if( i == 1 ) {
           altdacc <- as.matrix(na.omit(cbind(altdif,altdwt)))
           }
       else {
           altdacc <- rbind( altdacc, as.matrix(na.omit(cbind(altdif,altdwt))))
           }
       }
    alxmin <- min(altdacc[,1])
    alxmax <- max(altdacc[,1])
    alymax <- max(altdacc[,2])
    ALolist <- c(ALolist, list(altdd=altdacc, alxmin=alxmin, alxmax=alxmax, alymax=alymax))

    xx <- as.vector(altdacc[,1])
    yy <- as.vector(altdacc[,2])
    oo <- order(xx)
    xx <- xx[oo]
    yy <- yy[oo]
    nx <- length(xx)
    # Locate and eliminate duplicate xx values...
    j <- 0
    if( nx > 1 ) {
        for(i in seq(nx, 2, by=-1) ) {
            if( xx[(i-1)]==xx[i] ) {
                yy[(i-1)] <- yy[(i-1)]+yy[i]
                xx[i] <- alxmax + 9999
                j <- j+1
                }
            }
        }
    if( j > 0 ) {
        oo <- order(xx)
        xx <- xx[oo]
        yy <- yy[oo]
        nx <- nx - j
        }
    yy <- yy/sum(yy[1:nx])
    # Locate mid-points between jumps...
    if( nx > 2 ) {
        for(i in 2:(nx-1)) {
            xx[i] <- (xx[i]+xx[(i+1)])/2
            }
        }
    qq <- seq(0,1,len=min(1001,1+round(nx*4)))
    cdf <- qq
    for(i in 1:length(qq)) {
        cdf[i] <- UPlinint(qq[i], alxmin, nx, xx, yy)
        }

    ALolist <- c(ALolist, list(altdcdf=cdf, qq=qq))

    if(NNobjnam!="NA") {
       nnltdif <- NNobj$pbindif
       nnltdwt <- NNobj$pbinsde
       nnltdsm <- min(mean(nnltdwt, na.rm=TRUE), min(nnltdwt[nnltdwt!=0.0]))
       nnltdwt[nnltdwt==0.0] <- nnltdsm
       nnltdwt <- 1 / nnltdwt^2
       nnltdacc <- as.matrix(na.omit(cbind(nnltdif,nnltdwt)))
       nnltdacc[,2] <- nnltdacc[,2] / sum(nnltdacc[,2])
       nnlxmin <- min(nnltdacc[,1])
       nnlxmax <- max(nnltdacc[,1])
       nnlymax <- max(nnltdacc[,2])
       ALolist <- c(ALolist, list(nnltdd=nnltdacc, nnlxmin=nnlxmin, nnlxmax=nnlxmax, nnlymax=nnlymax))

       xx <- as.vector(nnltdacc[,1])
       yy <- as.vector(nnltdacc[,2])
       oo <- order(xx)
       xx <- xx[oo]
       yy <- yy[oo]
       nx <- length(xx)
       # Locate and eliminate duplicate xx values...
       j <- 0
       if( nx > 1 ) {
          for(i in seq(nx, 2, by=-1) ) {
               if( xx[(i-1)]==xx[i] ) {
                   yy[(i-1)] <- yy[(i-1)]+yy[i]
                   xx[i] <- alxmax + 9999
                   j <- j+1
                   }
               }
           }
       if( j > 0 ) {
           oo <- order(xx)
           xx <- xx[oo]
           yy <- yy[oo]
           nx <- nx - j
           }
       yy <- yy/sum(yy[1:nx])
       # Locate mid-points between jumps...
       if( nx > 2 ) {
           for(i in 2:(nx-1)) {
               xx[i] <- (xx[i]+xx[(i+1)])/2
               }
           }
       qq <- seq(0,1,len=min(501,1+round(nx*4)))
       cdf <- qq
       for(i in 1:length(qq)) {
           cdf[i] <- UPlinint(qq[i], alxmin, nx, xx, yy)
           }

       ALolist <- c(ALolist, list(nnltdcdf=cdf, nq=qq))
       }

    class(ALolist) <- "UPSaltdd"
    ALolist
}

#' @author Bob Obenchain <wizbob@att.net>
#' @export
"plot.SPSbalan" <- function (x, ...)
{
    if (x$youtype == "contin" && x$bins > 1) {
        form <- as.formula(paste(x$trtm, "~", x$xvar, "|", x$qbin))
        bwplot(form, data = x$df3)
    }
}

#' @author Bob Obenchain <wizbob@att.net>
#' @export
"plot.SPSloess" <- function (x, tcol = "blue", ucol = "red", dcol = "green3", ...)
{
    # First plot...
    plot(x$losub0[, 1], x$losub0[, 2], type = "p", ann = FALSE,
        col = tcol, pch = 21)
    lines(x$losub0[, 1], x$losub0[, 4], type = "l", col = tcol,
        lty = "solid")
    points(x$losub1[, 1], x$losub1[, 2], type = "p", col = ucol,
        pch = 24)
    lines(x$losub1[, 1], x$losub1[, 4], type = "l", col = ucol,
        lty = "dashed")
    title(main = paste("Outcomes and PS Loess Fit, Span =", x$span),
        ylab = paste("Observed and Smoothed", x$yvar), xlab = "Estimated Propensity Score")
    opar <- par(ask = dev.interactive(orNone = TRUE))
    # Second plot...
    plot(x$losub0[, 1], x$losub0[, 4], type = "l", ann = FALSE,
        col = tcol, lty = "solid")
    lines(x$losub1[, 1], x$losub1[, 4], type = "l", col = ucol,
        lty = "dashed")
    title(main = paste("PS Loess Fit, Span =", x$span), ylab = paste("Smoothed",
        x$yvar), xlab = "Estimated Propensity Score")
    # Third plot...
    plot(x$logrid[, 1], 0.1 * x$logrid[, 10], type = "p", ann = FALSE,
        pch = 21)
    lines(x$logrid[, 1], x$logrid[, 11], type = "l", col = dcol,
        lty = "solid")
    title(main = "PS Probability Density", xlab = "Estimated Propensity Score")
    par(opar)
    NULL
}

#' @author Bob Obenchain <wizbob@att.net>
#' @export
"plot.SPSoutco" <- function (x, ...)
{
    # First plot...
    barplot(t(as.matrix(x$binfreq[, 2:3])), beside = TRUE, names.arg = as.vector(x$binfreq[,
        1]), main = "Frequencies by PS Bin")
    if (x$youtype == "contin") {
        opar <- par(ask = dev.interactive(orNone = TRUE))
        # Second plot...
        barplot(t(as.matrix(x$binmean[, 2:3])), beside = TRUE,
            names.arg = as.vector(x$binmean[, 1]), main = "Outcome Means by PS Bin")
        # Third plot...
        plot(x$pbindif, x$pbinsde, ann = FALSE, type = "n", ylim = c(0,
            max(x$pbinsde, na.rm = TRUE)))
        symbols(x$pbindif, x$pbinsde, circles = x$pbinsiz, inches = 0.25,
            add = TRUE)
        par(lty = 1)
        abline(v = 0)
        par(lty = 2)
        abline(v = x$awbdif)
        par(lty = 3)
        abline(v = x$wwbdif)
        title(main = paste("Supervised Propensity Scoring: BINS =",
            x$bins), xlab = "Within Bin Treatment Difference",
            ylab = "Difference Standard Deviation", sub = paste("Number of Informative PS Bins =",
                length(na.omit(x$pbinsde))))
        par(opar)
    }
}

#' @author Bob Obenchain <wizbob@att.net>
#' @export
"plot.SPSsmoot" <- function (x, tcol = "blue", ucol = "red", dcol = "green3", ...)
{
    # First plot...
    plot(x$spsub0[, 1], x$spsub0[, 2], type = "p", ann = FALSE,
        col = tcol, pch = 21)
    lines(x$spsub0[, 1], x$spsub0[, 4], type = "l", col = tcol,
        lty = "solid")
    points(x$spsub1[, 1], x$spsub1[, 2], type = "p", col = ucol,
        pch = 24)
    lines(x$spsub1[, 1], x$spsub1[, 4], type = "l", col = ucol,
        lty = "dashed")
    title(main = paste("Outcomes and PS Smoothing Splines, df =",
        x$df), ylab = paste("Observed and Smoothed", x$yvar),
        xlab = "Estimated Propensity Score")
    opar <- par(ask = dev.interactive(orNone = TRUE))
    # Second plot...
    plot(x$spsub0[, 1], x$spsub0[, 4], type = "l", ann = FALSE,
        col = tcol, lty = "solid")
    lines(x$spsub1[, 1], x$spsub1[, 4], type = "l", col = ucol,
        lty = "dashed")
    title(main = paste("PS Smoothing Splines, df =", x$df), ylab = paste("Smoothed",
        x$yvar), xlab = "Estimated Propensity Score")
    # Third plot...
    plot(x$ssgrid[, 1], 0.1 * x$ssgrid[, 10], type = "p", ann = FALSE,
        pch = 21)
    lines(x$ssgrid[, 1], x$ssgrid[, 11], type = "l", col = dcol,
        lty = "solid")
    title(main = "PS Probability Density", xlab = "Estimated Propensity Score")
    par(opar)
    NULL
}

#' @author Bob Obenchain <wizbob@att.net>
#' @export
"plot.UPSaltdd" <- function(x, breaks="Sturges", ...)
{
  opar <- par(no.readonly = TRUE, ask = dev.interactive(orNone = TRUE))
  on.exit(par(opar))
  if( x$NNobj=="NA" ) par(mfrow = c(2, 1)) else par(mfrow = c(2, 2))
  nx <- length(x$altdd[,1])
  if( x$NNobj!="NA") {
      xmin <- min(x$alxmin, x$nnlxmin)
      xmax <- max(x$alxmax, x$nnlxmax)
      }
  else {
      xmin <- x$alxmin
      xmax <- x$alxmax
      }
  plot(x$altdcdf, x$qq, ann = FALSE, type = "l", xlim = c(xmin,xmax))
  par(lty=1)
  title(main=paste("Artificial LTD CDF"),
      ylab = "Probability Less Than",
      xlab =paste("From ", nx, "Random, Informative Clusters within", x$reps, "reps of", x$clus))
  hist(x$altdcdf, breaks=breaks, xlim=c(xmin,xmax),
      main=paste("Artificial LTD Histogram"), ylab = "Probability Density",
      xlab = paste("Artificial Effects of", x$trtm, "on", x$yvar, "using Random Clusters.")
      )
  if( x$NNobj!="NA" ) {
      nx <- length(x$nnltdd[,1])
      plot(x$nnltdcdf, x$nq, ann = FALSE, type = "l", xlim = c(xmin,xmax))
      par(lty=1)
      title(main=paste("Nearest Neighbor LTD CDF"),
          ylab = "Probability Less Than",
          xlab =paste("From ", nx, "Informative Clusters within", x$clus))
      hist(x$nnltdcdf, breaks=breaks, xlim=c(xmin,xmax),
          main=paste("Nearest Neighbor LTD Histogram"), ylab = "Probability Density",
          xlab = paste("Local Effects of", x$trtm, "on", x$yvar, "using Relevant Clusters.")
          )
    }
}

#' @author Bob Obenchain <wizbob@att.net>
#' @export
"plot.UPShclus" <- function (x, ...)
{
    if (x$method == "diana") {
        plot(x$upshcl, main = "Unsupervised Divisive Hierarchy",
            sub = paste("Divisive Coefficient = ", round(x$upshcl$dc,
                digits = 2)))
    }
    else if (x$method == "agnes") {
        plot(x$upshcl, main = "Unsupervised Agglomerative Hierarchy",
            sub = paste("Agglomerative Coefficient = ", round(x$upshcl$ac,
                digits = 2)))
    }
    else {
        plot(x$upshcl, main = "Unsupervised Clustering Hierarchy")
    }
}

#' @author Bob Obenchain <wizbob@att.net>
#' @export
"plot.UPSivadj" <- function (x, ...)
{
    if (x$youtype == "contin") {
        inchsz <- 2 * max(x$pbinsiz)^2/x$symsiz
        plot(x$pbinpsp, x$pbinout, ann = FALSE, type = "n")
        symbols(x$pbinpsp, x$pbinout, circles = x$pbinsiz, inches = inchsz,
            add = TRUE)
        par(lty = 1)
        if (x$actclust > 1)
            abline(x$ivfit)
        title(main = paste("Unsupervised IV Clusters =", x$actclust),
            xlab = "Within-Cluster Treatment Percentage (PS)",
            ylab = "Observed LATE", sub = "Symbol Area proportional to Cluster Size")
    }
}

#' @author Bob Obenchain <wizbob@att.net>
#' @export
"print.SPSbalan" <- function (x, ...)
{
    cat("\nSPSbalan: Checking for Within-Bin BALANCE on X Covariates\n")
    cat("Data Frame input:", x$dframe, "\n")
    cat("Treatment factor:", x$trtm, "\n")
    cat("Number of PS Bins:", x$bins, "\n")
    cat("Covariate X variable:", x$xvar, "\n")
    cat("\nTest for Overall Covariate Differences between Treatment Groups:\n")
    print(x$form)
    if (x$youtype == "contin") {
        print(summary(x$aovdiff))
        if (x$bins > 1) {
            cat("\nTest for Covariate Differences between Treatments within Bins:\n")
            print(x$form2)
            print(summary(x$bindiff))
        }
    }
    else {
        print(x$tab)
        print(chisq.test(x$tab))
        cat("\nTest for Covariate Differences between Treatments within Bins")
        cat("\n    Cumulative ChiSquare = ", round(x$cumchi,
            digits = 2))
        cat("\n    Cumulative df = ", x$cumdf)
        cat("\n    Cumulative p.value = ", round(1 - pchisq(x$cumchi,
            x$cumdf), digits = 4), "\n\n")
    }
}

#' @author Bob Obenchain <wizbob@att.net>
#' @export
"print.SPSloess" <- function (x, ...)
{
    cat("\n\nSPSloess Object: Supervised PS LOESS Smoothing Analysis\n")
    cat("Data Frame input:", x$dframe, "\n")
    cat("Treatment factor:", x$trtm, "\n")
    cat("Propensity Score variable:", x$pscr, "\n")
    cat("Outcome variable:", x$yvar, "\n")
    cat("\nOverall Average Outcome Means by Treatment")
    print(x$rawmean)
    cat("\nOverall Treatment Frequencies")
    print(x$rawfreq)
    cat("\nRaw Average Outcome Difference =", x$rawmean[2, ] -
        x$rawmean[1, ])
    cat("\nStandard Deviation of Difference =", sqrt(x$rawvars[2,
        ] + x$rawvars[1, ]))
    cat("\nTest of Raw Treatment Difference:\n")
    print(summary(x$aovdiff))
    cat("\n\nLOESS Smoothing Span =", x$span)
    cat("\nMean Spline Weighted Difference =", x$lotdif)
    cat("\nStandard Error of Spline Difference =", x$lotsde,
        "\n\n")
}

#' @author Bob Obenchain <wizbob@att.net>
#' @export
"print.SPSlogit" <- function (x, ...)
{
    cat("\nSPSlogit Object: Supervised Estimation of Propensity Scores\n")
    cat("Data Frame  input:", x$dfname, "\n")
    cat("Data Frame output:", x$dfoutnam, "\n")
    cat("Treatment Factor:", x$trtm, "\n")
    cat("Logistic Regression Formula:\n")
    print(x$form)
    cat("\n    Variable containing fitted PScores =", x$pfit)
    cat("\n    Variable containing PS Ranks =", x$prnk)
    cat("\n    PScore Bin Number Factor =", x$qbin)
    cat("\n    Number of PS Bins =", x$bins, "\n")
    print(summary(x$glmobj))
}

#' @author Bob Obenchain <wizbob@att.net>
#' @export
"print.SPSoutco" <- function (x, ...)
{
    cat("\nSPSoutco Object: Within Bin Outcome LTD\n")
    cat("Data Frame input:", x$dframe, "\n")
    cat("Outcome Y variable:", x$yvar, "\n")
    cat("Treatment difference:", x$PStdif, "\n")
    cat("Number of PS Bins:", x$bins, "\n")
    cat("Number of Informative PS Bins =", length(na.omit(x$pbinsde)),
        "\n")
    cat("\n    Overall Raw Average Treatment Difference =", x$ratdif)
    cat("\n    Standard Deviation of this Raw Difference =",
        x$ratsde)
    cat("\n    Average Within Bin Treatment Difference =", x$awbdif)
    cat("\n    Standard Deviation of Average Within Bin Difference =",
        x$awbsde)
    cat("\n    Inverse Variance Weighted Difference =", x$wwbdif)
    cat("\n    Standard Deviation of this Weighted Difference =",
        x$wwbsde)
    cat("\n\nTest for Raw / Unadjusted Treatment Difference:\n")
    print(x$form)
    if (x$youtype == "contin") {
        print(summary(x$aovdiff))
        if (x$bins > 1) {
            cat("\nTest for Treatment Difference within Bins:\n")
            print(x$form2)
            print(summary(x$bindiff))
        }
    }
    else {
        print(x$tab)
        print(chisq.test(x$tab))
        cat("\nTest for Treatment Difference within Bins")
        cat("\n    Cumulative ChiSquare = ", round(x$cumchi,
            digits = 2))
        cat("\n    Cumulative df = ", x$cumdf)
        cat("\n    Cumulative p.value = ", round(1 - pchisq(x$cumchi,
            x$cumdf), digits = 4), "\n\n")
    }
}

#' @author Bob Obenchain <wizbob@att.net>
#' @export
"print.SPSsmoot" <- function (x, ...)
{
    cat("\n\nSPSsmoot Object: Supervised PS Smoothing Spline Analysis\n")
    cat("Data Frame input:", x$dframe, "\n")
    cat("Treatment factor:", x$trtm, "\n")
    cat("Propensity Score variable:", x$pscr, "\n")
    cat("Outcome variable:", x$yvar, "\n")
    cat("\nOverall Average Outcome Means by Treatment")
    print(x$rawmean)
    cat("\nOverall Treatment Frequencies")
    print(x$rawfreq)
    cat("\nRaw Average Outcome Difference =", x$rawmean[2, ] -
        x$rawmean[1, ])
    cat("\nStandard Deviation of Difference =", sqrt(x$rawvars[2,
        ] + x$rawvars[1, ]))
    cat("\nTest of Raw Treatment Difference:\n")
    print(summary(x$aovdiff))
    cat("\n\nSmoothing Spline Degrees-of-Freedom =", x$df)
    cat("\nMean Spline Weighted Difference =", x$sptdif)
    cat("\nStandard Error of Spline Difference =", x$sptsde,
        "\n\n")
}

#' @author Bob Obenchain <wizbob@att.net>
#' @export
"print.UPSaltdd" <- function(x, ...)
{
    cat("\nUPSaltdd Object: Artificial Distribution of LTDs for random clusters...\n")
    cat("Data Frame:", x$dframe, "\n")
    cat("Outcome Variable:", x$yvar, "\n")
    cat("Treatment Factor:", x$trtm, "\n")
    cat("Scedasticity assumption:", x$scedas, "\n")
    cat("Number of Replications:", x$reps, "\n")
    cat("Number of Clusters per Replication:", x$clus, "\n")
    cat("Total Number of Informative Clusters =", length(x$altdd[,1]), "\n" )
    cat("\n    Mean Artificial Treatment Difference =", mean(x$altdd[,1]))
    cat("\n    Number of Smoothed Sample Quantiles  =", length(x$altdcdf))
    cat("\n    Mean of Smoothed Sample Quantiles    =", mean(x$altdcdf))
    cat("\n    Std. Deviation of Sample Quantiles   =", var(x$altdcdf)^0.5, "\n")
    if( x$NNobj!="NA" ) {
        cat("\nUPSnnltd Object:", x$NNobj, "\n")
        cat("Number of Informative Clusters =", length(x$nnltdd[,1]), "\n" )
        cat("\n    Mean of Observed LTD Distribution    =", mean(x$nnltdd[,1]))
        cat("\n    Number of Smoothed Sample Quantiles  =", length(x$nnltdcdf))
        cat("\n    Mean of Smoothed Sample Quantiles    =", mean(x$nnltdcdf))
        cat("\n    Std. Deviation of Sample Quantiles   =", var(x$nnltdcdf)^0.5, "\n")
        }
    }

#' @author Bob Obenchain <wizbob@att.net>
#' @export
"print.UPShclus" <- function (x, ...)
{
    cat("\n\nUPShclus object: Unsupervised Hierarchical Clustering\n")
    cat("\nData Frame input:", x$dframe)
    cat("\nClustering algorithm used:", x$method)
    cat("\nCovariate X variables:")
    print(x$xvars, quote = FALSE)
}

#' @author Bob Obenchain <wizbob@att.net>
#' @export
"print.UPSivadj" <- function (x, ...)
{
    cat("\nUPSivadj Object: Clustering Instrumental Variable (IV) Adjustment\n")
    cat("Hierarchical Clustering object:", x$hiclus, "\n")
    cat("Data Frame input:", x$dframe, "\n")
    cat("Treatment variable:", x$trtm, "\n")
    cat("Outcome variable:", x$yvar, "\n")
    cat("Maximum Number of Clusters:", x$numclust, "\n")
    if (x$actclust != x$numclust)
        cat("Actual Number of Clusters =", x$actclust, "\n")
    cat("\n    Overall Raw Average Outcome by Treatment\n")
    print(x$rawmean)
    cat("\n    Overall Treatment Frequency\n")
    print(x$rawfreq)
    cat("\n    Predicted Outcome at Treament Percentage Zero =",
        x$ivtzero)
    cat("\n    Standard Error at Treament Percentage Zero =",
        x$ivtzsde)
    cat("\n    Predicted Outcome at Treament Percentage 100 = ",
        x$ivt100p)
    cat("\n    Standard Error at Treament Percentage 100 =",
        x$ivt1pse)
    cat("\n    Predicted Overall Treatment Outcome Difference =",
        x$ivtdiff)
    cat("\n    Standard Error of Overall Outcome Difference =",
        x$ivtdsde)
    cat("\n\nWithin Cluster Average Outcomes\n")
    print(x$binmean)
    cat("\nWithin Cluster Treatment Frequencies\n")
    print(x$binfreq)
    if (x$youtype == "contin" && x$actclust > 1) {
        cat("\nWithin Cluster Treatment Percentages\n")
        print(x$pbinpsp)
        cat("\nCluster Symbol Radii (Root Total Sizes)\n")
        print(x$pbinsiz)
    }
}

#' Test for Within-Bin X-covariate Balance in Supervised Propensiy Scoring
#'
#' @description Test for Conditional Independence of X-covariate Distributions from Treatment
#'   Selection within Given, Adjacent PS Bins. The second step in Supervised Propensity Scoring analyses is to verify that baseline
#'   X-covariates have the same distribution, regardless of treatment, within each fitted PS bin.
#'
#' @param dframe Name of augmented data.frame written to the appn="" argument of SPSlogit().
#' @param trtm Name of the two-level treatment factor variable.
#' @param qbin Name of variable containing bin numbers.
#' @param envir The local control environment
#' @param yvar The outcome variable.
#' @param xvar Name of one baseline covariate X variable used in the SPSlogit() PS model.
#' @param faclev Maximum number of different numerical values an X-covariate can assume without
#'    automatically being converted into a "factor" variable; faclev=1 causes a binary indicator
#'    to be treated as a continuous variable determining a proportion.
#' @return An output list object of class SPSbalan. The first four are returned with a continuous x-variable. The next 4 are used if it is a factor variable.
#' \itemize{
#'   \item {aovdiff}{ANOVA output for marginal test.}
#'   \item {form2}{Formula for differences in X due to bins and to treatment nested within bins.}
#'   \item {bindiff}{ANOVA output for the nested within bin model.}
#'   \item {df3}{Output data.frame containing 3 variables: X-covariate, treatment and bin.}
#'
#'   \item{factab}{Marginal table of counts by X-factor level and treatment.}
#'   \item{tab}{Three-way table of counts by X-factor level, treatment and bin.}
#'   \item{cumchi}{Cumulative Chi-Square statistic for interaction in the three-way, nested table.}
#'   \item{cumdf}{Degrees of-Freedom for the Cumulative Chi-Squared.}
#' }
#'
#' @references
#' \itemize{
#'   \item Cochran WG. (1968) The effectiveness of adjustment by subclassification
#'   in removing bias in observational studies. \emph{Biometrics} \bold{24}:
#'   205-213.
#'
#'   \item Obenchain RL. (2011) \bold{USPSinR.pdf}  USPS R-package vignette, 40 pages.
#'
#'   \item Rosenbaum PR, Rubin RB. (1983) The Central Role of the Propensity Score
#'   in Observational Studies for Causal Effects. \emph{Biometrika} \bold{70}:
#'   41-55.
#'
#'   \item Rosenbaum PR, Rubin DB. (1984) Reducing Bias in Observational Studies
#'   Using Subclassification on a Propensity Score. \emph{J Amer Stat Assoc}
#'   \bold{79}: 516-524.
#' }
#'
#' @author Bob Obenchain <wizbob@att.net>
#' @export
"SPSbalan" <- function (envir, dframe, trtm, yvar, qbin, xvar, faclev = 3)
{
    if (missing(envir))
      stop("First argument to SPSbalan must be an environment to work in.\nUse ?SPSbalan for more information.")
    if (missing(dframe) || !exists(dframe, envir = envir))
      stop("Second argument to SPSbalan must be the name of a Data Frame in the working environment.")
    if (missing(trtm))
      stop("Third argument to SPSbalan must name the Treatment factor.")
    if (missing(yvar))
      stop("Fourth argument to SPSbalan must name the Target variable.")
    if(!exists("HCLolist", envir = envir))
      stop("The environment provided has no HCLolist, please call UPShclus first.")
    if (missing(qbin))
        stop("Third argument to SPSbalan must be the PS Bin Number variable.")

    data = envir[[dframe]]
    if (!is.element(yvar, dimnames(data)[[2]]))
      stop("Target Outcome or Covariate must be an existing Data Frame variable.")
    if (!is.element(trtm, dimnames(data)[[2]]))
      stop("Treatment factor must be an existing Data Frame variable.")
    if (length(table(data[, trtm])) != 2)
      stop("Treatment factor must assume exactly two different levels.")
    if (!is.element(qbin, dimnames(data)[[2]]))
        stop("PS Bin number factor must be an existing Data Frame variable.")
    if (missing(xvar))
        stop("Sixth argument to SPSbalan must name a X predictor variable.")
    if (!is.element(xvar, dimnames(data)[[2]]))
        stop("Predictor X variable must be an existing Data Frame variable.")

    bins <- length(table(data[, qbin]))
    SPSolist <- list(dframe = dframe, trtm = trtm,
        qbin = qbin, bins = bins, xvar = xvar)
    form <- as.formula(paste(xvar, "~", trtm))
    SPSolist <- c(SPSolist, list(form = form))
    if (length(table(data[, xvar])) > faclev) {
        youtype <- "contin"
        aovdiff <- aov(form, data, na.action = na.omit)
        SPSolist <- c(SPSolist, list(youtype = youtype, faclev = faclev,
            aovdiff = aovdiff))
        if (bins > 1) {
            form <- as.formula(paste(xvar, "~", qbin, "+", trtm,
                "%in%", qbin))
            aovdiff <- aov(form, data, na.action = na.omit)
            df3 <- as.data.frame(cbind(data[, xvar], as.numeric(as.character(data[,
                trtm])), data[, qbin]))
            df3[, 2] <- as.factor(df3[, 2])
            df3[, 3] <- as.factor(df3[, 3])
            names(df3) <- c(xvar, trtm, qbin)
            SPSolist <- c(SPSolist, list(form2 = form, bindiff = aovdiff,
                df3 = df3))
        }
    }
    else {
        youtype <- "factor"
        df3 <- na.omit(as.data.frame(cbind(as.factor(data[,
            xvar]), data[, trtm], data[, qbin])))
        names(df3) <- c(xvar, trtm, qbin)
        tab <- table(df3[, 1], df3[, 2])
        SPSolist <- c(SPSolist, list(youtype = youtype, faclev = faclev,
            factab = tab))
        tab <- table(df3[, 1], df3[, 2], df3[, 3])
        cumchi <- 0
        cumdf <- 0
        for (i in 1:bins) {
            ht <- chisq.test(tab[, , i])
            if (is.finite(ht$statistic)) {
                cumchi <- cumchi + ht$statistic
                cumdf <- cumdf + ht$parameters
            }
        }
        SPSolist <- c(SPSolist, list(tab = tab, cumchi = cumchi,
            cumdf = cumdf))
    }
    class(SPSolist) <- "SPSbalan"
    SPSolist
}

#
#' SPSloess
#'
#' @title LOESS Smoothing of Outcome by Treatment in Supervised Propensiy Scoring
#' @description Express Expected Outcome by Treatment as LOESS Smooths of Fitted Propensity Scores.
#' @param dframe data.frame of the form returned by SPSlogit().
#' @param trtm the two-level factor on the left-hand-side in the formula argument to SPSlogit().
#' @param pscr fitted propensity scores of the form returned by SPSlogit().
#' @param yvar continuous outcome measure or result unknown at the time patient was assigned
#'             (possibly non-randomly) to treatment; "NA"s are allowed in yvar.
#' @param faclev optional; maximum number of distinct numerical values a variable can assume
#'               and yet still be converted into a factor variable; faclev=1 causes a binary indicator to be
#'               treated as a continuous variable determining a proportion.
#' @param deg optional; degree (1=linear or 2=quadratic) of the local fit.
#' @param span optional; span (0 to 2) argument for the loess() function.
#' @param fam optional; "gaussian" or "symmetric".
#' @param envir Local control classic environment.
#'
#' @details {Once one has fitted a somewhat smooth curve through scatters of observed outcomes, Y,
#'  versus the fitted propensity scores, X, for the patients in each of the two treatment groups,
#'  one can consider the question: "Over the range where both smooth curves are defined (i.e. their
#'  common support), what is the (weighted) average signed difference between these two curves?"
#'
#'  If the distribution of patients (either treated or untreated) were UNIFORM over this range, the
#'  (unweighted) average signed difference (treated minus untreated) would be an appropriate
#'  estimate of the overall difference in outcome due to choice of treatment.
#'
#'  Histogram patient counts within 100 cells of width 0.01 provide a naive "non-parametric density
#'  estimate" for the distribution of total patients (treated or untreated) along the propensity
#'  score axis.  The weighted average difference (and standard error) displayed by SPSsmoot() are
#'  based on an R density() smooth of these counts.
#'
#'  In situations where the propensity scoring distribution for all patients in a therapeutic class
#'  is known to differ from that of the patients within the current study, that population weighted
#'  average would also be of interest.  Thus the SPSloess() output object contains two data frames,
#'  logrid and lofit, useful in further computations.
#' }
#' \itemize{
#'   \item{logrid}{loess grid data.frame containing 11 variables and 100 observations. The PS
#'    variable contains propensity score "cell means" of 0.005 to 0.995 in steps of 0.010.
#'    Variables F0, S0 and C0 for treatment 0 and variables F1, S1 and C1 for treatment 1 contain
#'    fitted smooth spline values, standard error estimates and patient counts, respectively.  The
#'    DIF variable is simply (F1\-F0), the SED variable is sqrt(S1\^2+S0\^2), the HST variable is
#'    proportional to (C0+C1), and the DEN variable is the estimated probability density of patients
#'    along the PS axis.  Observations with "NA" for variables F0, S0, F1 or S1 represent "extremes"
#'    where the lowess fits could not be extrapolated because no observed outcomes were available.}
#'   \item{losub0, losub1}{loess fit data.frame contains 4 variables for each distinct PS value in lofit.
#'    These 4 variables are named PS, YAVG, TRT==0 and 1, respectively, and FIT = spline prediction for
#'    the specified degrees-of-freedom (default df=1.)}
#'   \item{span}{loess span setting.}
#'   \item{lotdif}{outcome treatment difference mean.}
#'   \item{lotsde}{outcome treatment difference standard deviation.}
#' }
#' @references{
#'   Cleveland WS, Devlin SJ. (1988) Locally-weighted regression: an
#'   approach to regression analysis by local fitting. \emph{J Amer Stat Assoc}
#'   \bold{83}: 596-610.
#'
#'   Cleveland WS, Grosse E, Shyu WM. (1992) Local regression models. Chapter 8 of
#'   \bold{Statistical Models in S} eds Chambers JM and Hastie TJ. \emph{Wadsworth & Brooks/Cole}.
#'
#'   Obenchain RL. (2011) \bold{USPSinR.pdf}  USPS R-package vignette, 40 pages.
#'
#'   Ripley BD, loess() based on the 'cloess' package of Cleveland, Grosse and Shyu.
#' }
#'
#' @author Bob Obenchain <wizbob@att.net>
#' @export
"SPSloess" <- function (envir, dframe, trtm, pscr, yvar, faclev = 3, deg = 2, span = 0.75,fam = "symmetric")
{
    if (missing(envir))
      stop("First argument to SPSloess must be the name of an existing Data Frame.")
    if (missing(dframe) || !exists(dframe, envir = envir))
      stop("Second argument to SPSloess must be the name of a Data Frame in the working environment.")
    if (missing(trtm))
      stop("Third argument to SPSloess must name the Treatment factor.")
    if (missing(yvar))
      stop("Fourth argument to SPSloess must name the Target variable.")

    data = envir[[dframe]]
    if (!is.element(yvar, dimnames(data)[[2]]))
      stop("Target Outcome or Covariate must be an existing Data Frame variable.")
    if (!is.element(trtm, dimnames(data)[[2]]))
      stop("Treatment factor must be an existing Data Frame variable.")
    if (length(table(data[, trtm])) != 2)
      stop("Treatment factor must assume exactly two different levels.")

    if (missing(pscr))
        stop("Third argument to SPSlowes must name the PScore variable.")
    if (!is.element(pscr, dimnames(data)[[2]]))
        stop("Propensity Score must be an existing Data Frame variable.")


    SPSolist <- list(dframe = dframe, trtm = trtm,
        pscr = pscr, yvar = yvar)

    dframe = data
    PSmean <- as.matrix(tapply(dframe[, yvar], dframe[, trtm],
        na.rm = TRUE, mean))
    PStrtm <- paste(trtm, "=", dimnames(PSmean)[[1]])
    PSvars <- as.matrix(tapply(dframe[, yvar], dframe[, trtm],
        na.rm = TRUE, var))
    PSfreq <- as.matrix(table(dframe[, trtm]))
    PSvars <- PSvars/PSfreq
    form <- as.formula(paste(yvar, "~", trtm))
    aovdiff <- aov(form, dframe, na.action = na.omit)
    SPSolist <- c(SPSolist, list(rawmean = PSmean, rawvars = PSvars,
        rawfreq = PSfreq, form = form, aovdiff = aovdiff))
    m <- order(dframe[, pscr])
    lofit <- as.data.frame(cbind(dframe[m, pscr], dframe[m, yvar],
        as.numeric(as.character(dframe[m, trtm])), 0))
    lofit <- na.omit(lofit)
    names(lofit) <- c("PS", "YVAR", "TRTM", "FIT")
    logrid <- as.data.frame(cbind(seq(0.005, 0.995, length = 100),
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
    names(logrid) <- c("PS", "F0", "S0", "C0", "F1", "S1", "C1",
        "DIF", "SED", "HST", "DEN")
    losub0 <- subset.data.frame(lofit, lofit$TRTM == 0)
    losub1 <- subset.data.frame(lofit, lofit$TRTM == 1)
    loobj0 <- loess(YVAR ~ PS, losub0, family = fam, degree = deg,
        span = span)
    losub0$FIT <- fitted.values(loobj0)
    fit0 <- predict(loobj0, logrid$PS, se = TRUE)
    logrid$F0 <- fit0$fit
    logrid$S0 <- fit0$se.fit
    logrid$C0 <- hist(losub0$PS, breaks = seq(0, 1, length = 101),
        plot = FALSE)$counts
    loobj1 <- loess(YVAR ~ PS, losub1, family = fam, degree = deg,
        span = span)
    losub1$FIT <- fitted.values(loobj1)
    fit1 <- predict(loobj1, logrid$PS, se = TRUE)
    logrid$F1 <- fit1$fit
    logrid$S1 <- fit1$se.fit
    logrid$C1 <- hist(losub1$PS, breaks = seq(0, 1, length = 101),
        plot = FALSE)$counts
    logrid$DIF <- logrid$F1 - logrid$F0
    logrid$SED <- sqrt(logrid$S0^2 + logrid$S1^2)
    logrid$HST <- hist(lofit$PS, breaks = seq(0, 1, length = 101),
        plot = FALSE)$counts
    logrid$DEN <- density(lofit$PS, bw = "nrd0", kernel = "gaussian",
        n = 100, from = 0.005, to = 0.995)$y
    grid <- na.omit(logrid)
    if (sum(grid$DEN) > 0)
        grid$DEN <- grid$DEN/sum(grid$DEN)
    PSmean <- sum(grid$DIF * grid$DEN)
    PSster <- sum(grid$SED * grid$DEN)
    SPSolist <- c(SPSolist, list(logrid = logrid, losub0 = losub0,
        losub1 = losub1, span = span, lotdif = PSmean, lotsde = PSster))
    class(SPSolist) <- "SPSloess"
    SPSolist
}


#' @name SPSlogit
#' @title Propensity Score prediction of Treatment Selection from Patient Baseline X-covariates
#' @description Use a logistic regression model to predict Treatment Selection from Patient
#'  Baseline X-covariates in Supervised Propensity Scoring.
#'
#' @param envir name of the working local control classic environment.
#' @param dframe data.frame containing X, t and Y variables.
#' @param form Valid formula for glm()with family = binomial(), with the two-level treatment
#'  factor variable as the left-hand-side of the formula.
#' @param pfit Name of variable to store PS predictions.
#' @param prnk Name of variable to store tied-ranks of PS predictions.
#' @param qbin Name of variable to store the assigned bin number for each patient.
#' @param bins optional; number of adjacent PS bins desired; default to 5.
#' @param appn optional; append the pfit, prank and qbin variables to the input dfname when
#'   appn=="", else save augmented data.frame to name specified within a non-blank appn string.
#'
#'
#' @details {The first phase of Supervised Propensity Scoring is to develop a logit (or probit) model
#'  predicting treatment choice from patient baseline X characteristics.  SPSlogit uses a call to
#'  glm()with family = binomial() to fit a logistic regression.}
#' @return An output list object of class SPSlogit:
#' \itemize{
#'  \item{dframe}{Name of input data.frame containing X, t & Y variables.}
#'  \item{dfoutnam}{Name of output data.frame augmented by pfit, prank and qbin variables.}
#'  \item{trtm}{Name of two-level treatment factor variable.}
#'  \item{form}{glm() formula for logistic regression.}
#'  \item{pfit}{Name of predicted PS variable.}
#'  \item{prank}{Name of variable containing PS tied-ranks.}
#'  \item{qbin}{Name of variable containing assigned PS bin number for each patient.}
#'  \item{bins}{Number of adjacent PS bins desired.}
#'  \item{glmobj}{Output object from invocation of glm() with family = binomial().}
#' }
#' @references{
#' Cochran WG. (1968) The effectiveness of adjustment by subclassification
#'  in removing bias in observational studies. \emph{Biometrics} \bold{24}:
#'  205-213.
#'
#'  Kereiakes DJ, Obenchain RL, Barber BL, et al. (2000) Abciximab provides
#'  cost effective survival advantage in high volume interventional practice.
#'  \emph{Am Heart J} \bold{140}: 603-610.
#'
#'  Obenchain RL. (2011) \bold{USPSinR.pdf}  USPS R-package vignette, 40 pages.
#'
#'  Rosenbaum PR, Rubin RB. (1983) The Central Role of the Propensity Score
#'  in Observational Studies for Causal Effects. \emph{Biometrika} \bold{70}:
#'  41-55.
#'
#'  Rosenbaum PR, Rubin DB. (1984) Reducing Bias in Observational Studies
#'  Using Subclassification on a Propensity Score. \emph{J Amer Stat Assoc}
#'  \bold{79}: 516-524.
#' }
#' @seealso \code{\link{SPSbalan}}, \code{\link{SPSnbins}} and \code{\link{SPSoutco}}.
#' @keywords models
#' @author Bob Obenchain <wizbob@att.net>
#' @export
"SPSlogit" <- function (envir, dframe, form, pfit, prnk, qbin, bins = 5, appn = "")
{
    if (missing(envir))
      stop("First argument to SPSlogit must be the name of an existing Data Frame.")
    if (missing(dframe) || !exists(dframe, envir = envir))
      stop("Second argument to SPSlogit must be the name of a Data Frame in the working environment.")
    if (missing(trtm))
      stop("Third argument to SPSlogit must name the Treatment factor.")
    if (missing(form) || class(form) != "formula")
        stop("Second argument to SPSlogit must be a valid formula.")
    trtm <- deparse(form[[2]])
    dfname = deparse(substitute(dframe))
    if (missing(pfit))
        stop("Third argument to SPSlogit must name the PScore variable.")
    if (missing(prnk))
        stop("Fourth argument to SPSlogit must name the PS Rank variable.")
    if (missing(qbin))
        stop("Fifth argument to SPSlogit must name the Bin Number variable.")

    glmobj <- glm(form, family = binomial(), data = dframe)
    df3 <- as.data.frame(fitted.values(glmobj))
    pfit <- deparse(substitute(pfit))
    names(df3) <- pfit
    prnk <- deparse(substitute(prnk))
    df3[, prnk] <- rank(df3[, pfit], na.last = TRUE)
    qbin <- deparse(substitute(qbin))
    df3[, qbin] <- factor(1 + floor((bins * df3[, prnk])/(1 +
        length(df3[, prnk]))))
    dframe <- merge(dframe, df3, by.x = "row.names", by.y = "row.names",
        all.x = TRUE)
    if (appn == "")
        dfoutnam <- dfname
    else dfoutnam <- appn
    assign(dfoutnam, dframe, envir = envir)
    SPSolist <- list(dfname = dfname, dfoutnam = dfoutnam, trtm = trtm,
        form = form, pfit = pfit, prnk = prnk, qbin = qbin, bins = bins,
        glmobj = glmobj)
    class(SPSolist) <- "SPSlogit"
    SPSolist
}



#' @name SPSnbins
#' @title Change the Number of Bins in Supervised Propensiy Scoring
#' @description {Change the Number of Bins in Supervised Propensiy Scoring}
#'
#' @param dframe {Name of data.frame of the form output by SPSlogit().}
#' @param envir name of the working local control classic environment.
#' @param prnk {Name of PS tied-rank variable from previous call to SPSlogit().}
#' @param qbin {Name of variable to contain the re-assigned bin number for each patient.}
#' @param bins {Number of PS bins desired.}
#'
#' @details Part or all of the first phase of Supervised Propensity Scoring will need to be redone
#'   if SPSbalan() detects dependence of within-bin X-covariate distributions upon treatment choice.
#'   Use SPSnbins() to change (increase) the number of adjacent PS bins.  If this does not achieve
#'   balance, invoke SPSlogit() again to modify the form of your PS logistic model, typically by
#'   adding interaction and/or curvature terms in continuous X-covariates.
#' @return An output data.frame with new variables inserted:
#' \itemize{
#'   \item{dframe2}{Modified version of the data.frame specified as the first argument to SPSnbins().}
#' }
#' @references{
#'   Cochran WG. (1968) The effectiveness of adjustment by subclassification
#'   in removing bias in observational studies. \emph{Biometrics} \bold{24}:
#'   205-213.
#'
#'   Obenchain RL. (2011) \bold{USPSinR.pdf}  USPS R-package vignette, 40 pages.
#'
#'   Rosenbaum PR, Rubin DB. (1984) Reducing Bias in Observational Studies
#'   Using Subclassification on a Propensity Score. \emph{J Amer Stat Assoc}
#'   \bold{79}: 516-524.
#' }
#' @seealso \code{\link{SPSlogit}}, \code{\link{SPSbalan}} and \code{\link{SPSoutco}}.
#' @keywords design
#' @author Bob Obenchain <wizbob@att.net>
#' @export
"SPSnbins" <- function (envir, dframe, prnk, qbin, bins = 8)
{
    if (missing(envir))
      stop("First argument to SPSnbins must be the name of an existing Data Frame.")
    if (missing(dframe) || !exists(dframe, envir = envir))
      stop("Second argument to SPSnbins must be the name of a Data Frame in the working environment.")
    dframe <- envir[[dframe]]
    if (missing(prnk))
        stop("Second argument to SPSnbins must name the patient rank variable.")
    if (!is.element(prnk, dimnames(dframe)[[2]]))
        stop("SPSnbins patient ranks must be stored in an existing variable.")
    if (missing(qbin))
        stop("Third argument to SPSnbins must name the bin number variable.")

    subjects <- length(dframe[, prnk])
    if (bins < 2)
        bins <- 2
    if (bins > floor(subjects/2))
        bins <- floor(subjects/2)
    dframe[, qbin] <- factor(1 + floor((bins * dframe[, prnk])/(1 +
        subjects)))
    dframe
}



#' @name SPSoutco
#' @title {Examine Treatment Differences on an Outcome Measure in Supervised Propensiy Scoring}
#' @description {Examine Within-Bin Treatment Differences on an Outcome Measure and Average these
#'  Differences across Bins.}
#'
#' @param envir name of the working local control classic environment.
#' @param dframe {Name of augmented data.frame written to the appn="" argument of SPSlogit().}
#' @param trtm {Name of treatment factor variable.}
#' @param yvar {Name of an outcome Y variable.}
#' @param qbin {Name of variable containing the PS bin number for each patient.}
#' @param faclev {Maximum number of different numerical values an X-covariate can assume without
#'   automatically being converted into a "factor" variable; faclev=1 causes a binary indicator to
#'   be treated as a continuous variable determining an average or proportion.}
#'
#' @details Once the second phase of Supervised Propensity Scoring confirms, using SPSbalan(), that
#'  X-covariate Distributions have been Balanced Within-Bins, the third phase can start: Examining
#'  Within-Bin Outcome Difference due to Treatment and Averaging these Differences across Bins.
#'  Graphical displays of SPSoutco() results feature R barplot() invocations.
#' @return An output list object of class SPSoutco:
#' \itemize{
#'  \item{dframe}{Name of augmented data.frame written to the appn="" argument of SPSlogit().}
#'  \item{trtm}{Name of the two-level treatment factor variable.}
#'  \item{yvar}{Name of an outcome Y variable.}
#'  \item{bins}{Number of variable containing bin numbers.}
#'  \item{PStdif}{Character string describing the treatment difference.}
#'  \item{rawmean}{Unadjusted outcome mean by treatment group.}
#'  \item{rawvars}{Unadjusted outcome variance by treatment group.}
#'  \item{rawfreq}{Number of patients by treatment group.}
#'  \item{ratdif}{Unadjusted mean outcome difference between treatments.}
#'  \item{ratsde}{Standard error of unadjusted mean treatment difference.}
#'  \item{binmean}{Unadjusted mean outcome by cluster and treatment.}
#'  \item{binvars}{Unadjusted variance by cluster and treatment.}
#'  \item{binfreq}{Number of patients by bin and treatment.}
#'  \item{awbdif}{Across cluster average difference with cluster size weights.}
#'  \item{awbsde}{Standard error of awbdif.}
#'  \item{wwbdif}{Across cluster average difference, inverse variance weights.}
#'  \item{wwbsde}{Standard error of wwbdif.}
#'  \item{form}{Formula for overall, marginal treatment difference on X-covariate.}
#'  \item{faclev}{Maximum number of different numerical values an X-covariate can assume without
#'   automatically being converted into a "factor" variable; faclev=1 causes a binary indicator to be
#'   treated as a continuous variable determining an average or proportion.}
#'  \item{youtype}{"contin"uous => only next six outputs; "factor" => only last four outputs.}
#'  \item{aovdiff}{ANOVA output for marginal test.}
#'  \item{form2}{Formula for differences in X due to bins and to treatment nested within bins.}
#'  \item{bindiff}{ANOVA summary for treatment nested within bin.}
#'  \item{pbindif}{Unadjusted treatment difference by cluster.}
#'  \item{pbinsde}{Standard error of the unadjusted difference by cluster.}
#'  \item{pbinsiz}{Cluster radii measure: square root of total number of patients.}
#'  \item{factab}{Marginal table of counts by Y-factor level and treatment.}
#'  \item{tab}{Three-way table of counts by Y-factor level, treatment and bin.}
#'  \item{cumchi}{Cumulative Chi-Square statistic for interaction in the three-way, nested table.}
#'  \item{cumdf}{Degrees of-Freedom for the Cumulative Chi-Squared.}
#' }
#' @references{
#'  Cochran WG. (1968) The effectiveness of adjustment by subclassification
#'  in removing bias in observational studies. \emph{Biometrics} \bold{24}:
#'  205-213.
#'
#'  Obenchain RL. (2011) \bold{USPSinR.pdf}  USPS R-package vignette, 40 pages.
#'
#'  Rosenbaum PR, Rubin RB. (1983) The Central Role of the Propensity Score
#'  in Observational Studies for Causal Effects. \emph{Biometrika} \bold{70}:
#'  41-55.
#'
#'  Rosenbaum PR, Rubin DB. (1984) Reducing Bias in Observational Studies
#'  Using Subclassification on a Propensity Score. \emph{J Amer Stat Assoc}
#'  \bold{79}: 516-524.
#' }
#' @seealso {\code{\link{SPSlogit}}, \code{\link{SPSbalan}} and \code{\link{SPSnbins}}.}
#' @keywords nonparametric hplot
#' @author Bob Obenchain <wizbob@att.net>
#' @export
"SPSoutco" <- function (envir, dframe, trtm, qbin, yvar, faclev = 3)
{

    if (missing(envir))
      stop("First argument to SPSoutco must be the name of an existing Data Frame.")
    if (missing(dframe) || !exists(dframe, envir = envir))
      stop("Second argument to SPSoutco must be the name of a Data Frame in the working environment.")
    if (missing(trtm))
      stop("Third argument to SPSoutco must name the Treatment factor.")
    if (missing(yvar))
      stop("Fourth argument to SPSoutco must name the Target variable.")

    data = envir[[dframe]]
    if (!is.element(yvar, dimnames(data)[[2]]))
      stop("Target Outcome or Covariate must be an existing Data Frame variable.")
    if (!is.element(trtm, dimnames(data)[[2]]))
      stop("Treatment factor must be an existing Data Frame variable.")
    if (length(table(data[, trtm])) != 2)
      stop("Treatment factor must assume exactly two different levels.")
    if (missing(qbin))
        stop("Third argument to SPSoutco must be the PS Bin Number variable.")
    if (!is.element(qbin, dimnames(data)[[2]]))
        stop("PS Bin number must be an existing Data Frame variable.")


    dframe = data
    bins <- length(table(dframe[, qbin]))
    SPSolist <- list(dframe = deparse(substitute(dframe)), trtm = trtm,
        yvar = yvar, bins = bins)
    PSmean <- as.matrix(tapply(dframe[, yvar], dframe[, trtm],
        na.rm = TRUE, mean))
    PStrtm <- paste(trtm, "=", dimnames(PSmean)[[1]])
    PStdif <- paste(PStrtm[2], "minus", PStrtm[1])
    PSvars <- as.matrix(tapply(dframe[, yvar], dframe[, trtm],
        na.rm = TRUE, var))
    PSfreq <- as.matrix(table(dframe[, trtm]))
    PSvars <- PSvars/PSfreq
    SPSolist <- c(SPSolist, list(PStdif = PStdif, rawmean = PSmean,
        rawvars = PSvars, rawfreq = PSfreq))
    RATdif <- sum(PSmean[2, ] - PSmean[1, ])
    RATsde <- sum(sqrt(PSvars[2, ] + PSvars[1, ]))
    SPSolist <- c(SPSolist, list(ratdif = RATdif, ratsde = RATsde))
    PSmean <- tapply(dframe[, yvar], list(dframe[, qbin], dframe[,
        trtm]), na.rm = TRUE, mean)
    PSmean <- cbind(matrix(1:bins, bins, 1), PSmean)
    dimnames(PSmean) <- list(1:bins, c("BIN", PStrtm[1], PStrtm[2]))
    PSmean <- as.data.frame(PSmean)
    PSmean[, "BIN"] <- as.factor(PSmean[, "BIN"])
    PSvars <- as.matrix(tapply(dframe[, yvar], list(dframe[,
        qbin], dframe[, trtm]), na.rm = TRUE, var))
    PSfreq <- as.matrix(table(dframe[, qbin], dframe[, trtm]))
    PSvars <- PSvars/PSfreq
    PSfreq <- cbind(matrix(1:bins, bins, 1), PSfreq)
    dimnames(PSfreq) <- list(1:bins, c("BIN", PStrtm[1], PStrtm[2]))
    PSfreq <- as.data.frame(PSfreq)
    PSfreq[, "BIN"] <- as.factor(PSfreq[, "BIN"])
    SPSolist <- c(SPSolist, list(binmean = PSmean, binvars = PSvars,
        binfreq = PSfreq))
    awbdif <- sum(na.omit(((PSmean[, 3] - PSmean[, 2]) * (PSfreq[,
        3] + PSfreq[, 2]))/sum(PSfreq[, 3] + PSfreq[, 2])))
    awbsde <- sqrt(sum(na.omit((PSvars[, 2] + PSvars[, 1]) *
        (PSfreq[, 3] + PSfreq[, 2])^2))/(sum(PSfreq[, 3] + PSfreq[,
        2]))^2)
    wwbdif <- sum(na.omit((PSmean[, 3] - PSmean[, 2])/(PSvars[,
        2] + PSvars[, 1])/sum(1/na.omit(PSvars[, 2] + PSvars[,
        1]))))
    wwbsde <- sqrt(1/sum(1/na.omit(PSvars[, 2] + PSvars[, 1])))
    form <- as.formula(paste(yvar, "~", trtm))
    SPSolist <- c(SPSolist, list(awbdif = awbdif, awbsde = awbsde,
        wwbdif = wwbdif, wwbsde = wwbsde, form = form))
    if (length(table(dframe[, yvar])) > faclev) {
        youtype <- "contin"
        aovdiff <- aov(form, dframe, na.action = na.omit)
        SPSolist <- c(SPSolist, list(youtype = youtype, faclev = faclev,
            aovdiff = aovdiff))
        if (bins > 1) {
            form <- as.formula(paste(yvar, "~", qbin, "+", trtm,
                "%in%", qbin))
            aovdiff <- aov(form, dframe, na.action = na.omit)
            SPSolist <- c(SPSolist, list(form2 = form, bindiff = aovdiff))
        }
        pbindif <- PSmean[, 3] - PSmean[, 2]
        pbinsde <- sqrt(PSvars[, 2] + PSvars[, 1])
        pbinsiz <- sqrt(PSfreq[, 3] + PSfreq[, 2])
        SPSolist <- c(SPSolist, list(pbindif = pbindif, pbinsde = pbinsde,
            pbinsiz = pbinsiz))
    }
    else {
        youtype <- "factor"
        df3 <- as.data.frame(cbind(dframe[, yvar], dframe[, trtm]))
        df3[, 3] <- dframe[, qbin]
        df3 <- na.omit(df3)
        names(df3) <- c(yvar, trtm, qbin)
        df3[, 1] <- as.factor(df3[, 1])
        tab <- table(df3[, 1], df3[, 2])
        SPSolist <- c(SPSolist, list(youtype = youtype, faclev = faclev,
            factab = tab))
        tab <- table(df3[, 1], df3[, 2], df3[, 3])
        cumchi <- 0
        cumdf <- 0
        for (i in 1:bins) {
            ht <- chisq.test(tab[, , i])
            if (!is.na(ht$statistic) && is.finite(ht$statistic)) {
                cumchi <- cumchi + ht$statistic
                cumdf <- cumdf + ht$parameters
            }
        }
        SPSolist <- c(SPSolist, list(tab = tab, cumchi = cumchi,
            cumdf = cumdf))
    }
    class(SPSolist) <- "SPSoutco"
    SPSolist
}

#Creates a summary frame for this Local Control given an output object. Makes for some easier plotting.
"UPSSummary" <- function(LCobj, envir)
{
  lcsumm = numeric(9)
  names(lcsumm) = c("Requested_clusters", "Clusters_created", "Informative_clusters",
                    "Informative_patients", "Informative_fraction", "Across_cluster_LTD",
                    "LTD_standard_error", "LTD_lower_confidence", "LTD_upper_confidence")

  lcsumm["Requested_clusters"] = LCobj$numclust
  lcsumm["Clusters_created"] = LCobj$actclust
  lcsumm["Across_cluster_LTD"] = LCobj$awbdif
  #lcsumm["pbindif"] = LCobj$pbindif
  #lcsumm["ratdif"] = LCobj$ratdif
  lcsumm["wwbdif"] = LCobj$wwbdif
  lcsumm["LTD_standard_error"] = LCobj$awbsde
  lcsumm["Informative_clusters"] = length(which(!is.na(LCobj$pbindif)))
  infPatFrame = LCobj$binfreq
  infPatFrame = infPatFrame[which(infPatFrame[,2] > 0 & infPatFrame[,3] > 0),] #Count the patients in informative clusters.
  lcsumm["Informative_patients"] = sum(infPatFrame[2] + infPatFrame[3])
  lcsumm["Informative_fraction"] = lcsumm["Informative_patients"] / sum(LCobj$rawfreq)
  lcsumm["LTD_lower_confidence"] = LCobj$awbdif - 2*sqrt(LCobj$awbsde)
  lcsumm["LTD_upper_confidence"] = LCobj$awbdif + 2*sqrt(LCobj$awbsde)

  if(exists("summary", envir = envir)){
    envir[["summary"]] = rbind(envir[["summary"]], lcsumm)
  }
  else{
    envir[["summary"]] = data.frame(t(lcsumm))
  }

  envir[["summary"]]
}


#' @name UPSLTDdist
#' @title Plot the LTD distribution as a function of the number of clusters.
#' @description This function creates a plot displaying the distribution of
#'  Local Treatment Differences (LTDs) as a function of the number of clusters
#'  created for all UPSnnltd objects in the provided environment. The hinges and
#'  whiskers are generated using \code{\link[grDevices]{boxplot.stats}}.
#' @param envir A LocalControlClassic environment containing UPSnnltd objects.
#' @param legloc Where to place the legend in the returned plot. Defaults to "bottomleft".
#' @inheritDotParams graphics::plot -x -y
#'
#' @return Returns the LTD distribution plot.
#' @return Adds the "ltdds" object to envir.
#'
#' @examples
#'
#'  data(lindner)
#'  cvars <- c("stent","height","female","diabetic","acutemi",
#'             "ejecfrac","ves1proc")
#'  numClusters <- c(1, 2, 10, 15, 20, 25, 30, 35, 40, 45, 50)
#'  results <- LocalControlClassic(data = lindner,
#'                                 clusterVars = cvars,
#'                                 treatmentColName = "abcix",
#'                                 outcomeColName = "cardbill",
#'                                 clusterCounts = numClusters)
#'  UPSLTDdist(results,ylim=c(-15000,15000))
#'
#' @export
"UPSLTDdist" <- function(envir, legloc = "bottomleft", ...){

  clusts = envir[["summary"]]$Requested_clusters

  ltdds = list()
  for(cc in clusts){
      objName = paste0("UPSnnltd", cc)
      nnObj = envir[[objName]]
      ltdds[[objName]] = nnObj$pbindif[nnObj$nnhbindf[,1]]
  }

  envir[["ltdds"]] = ltdds

  infPatFrac = envir[["summary"]]$Informative_fraction

  actClusterCounts = envir[["summary"]]$Clusters_created
  ltdm = sapply(ltdds, function(x) mean(x,na.rm=T))
  boxes = lapply(ltdds, FUN = function(x) boxplot(x, plot=FALSE))
  bStats = sapply(boxes, FUN = function(x) x$stats)

  #add the boxplot stats for each of calls to upsnnltd
  envir[["boxstats"]] = bStats

  yname = envir[["UPSaccum.pars"]][,"yvar"]

  globMean = ltdm[1] ## Appropriately Weighted by Cluster Size = Number of Patients

  ## New LTD Distribution Graphic...
  par(mar = c(5,5,2,5))
  plot(..., x = actClusterCounts, y = ltdm,
  	lwd = 3, type = "l", col = "black",
  	xlab = "Number of Clusters",
  	ylab = paste0(yname," LTD Distribution"),
  	main = paste0(yname," Local Treatment Differences vs. Number of Clusters"))
  grid()

  lines(actClusterCounts, bStats[1,], col = "blue", lty = 3, lwd = 2) #lower whisker
  lines(actClusterCounts, bStats[2,], col = "black", lty = 3, lwd = 2)#lower hinge
  lines(actClusterCounts, bStats[3,], col = "red", lty = 3, lwd = 2)  #median
  lines(actClusterCounts, bStats[4,], col = "black", lty = 3, lwd = 2)#upper hinge
  lines(actClusterCounts, bStats[5,], col = "blue", lty = 3, lwd = 2) #upper whisker

  par(new = T)
  plot( x = actClusterCounts, y = infPatFrac,
      	type = "l", lwd = 3, lty = 3,
      	axes=F, xlab=NA, ylab=NA,
      	cex=1.2, col = "green", ylim = c(0,1))

  axis(side = 4)
  mtext(side = 4, line = 3, 'Fraction of Patients Informative')
  legend( legloc,
        	col = c("black", "green", "red", "black", "blue"),
        	lty = c(1, 3, 3, 3, 3), lwd = c(3, 2, 2, 2, 2),
        	legend = c( "LTD Main-Effect",
                  		"Fraction of Patients Informative",
                  		"LTD Median",
                  		"LTD Hinges",
                  		"LTD Whiskers"))
}

#' @name UPSboxplot
#'
#' @title Returns a series of boxplots comparing LTD distributions given different numbers of clusters.
#' @description Given the output of \code{\link{LocalControlClassic}}, this function uses all or some of the
#' UPSnnltd objects contained to create a series of boxplots of the local treatment difference at each of the
#' different numbers of requested clusters.
#' @param envir A LocalControlClassic environment containing UPSnnltd objects.
#' @param clusterSubset (optional) A vector containing requested cluster counts.
#'  If provided, the boxplot is created using only the UPSnnltd objects corresponding to the requested cluster counts.
#'
#' @return Returns the call to boxplot with the formula: "ltd ~ numclst".
#' @return Adds the "ltdds" object to the Local Control environment.
#'
#' @examples
#'
#' data(lindner)
#' cvars <- c("stent","height","female","diabetic","acutemi",
#'            "ejecfrac","ves1proc")
#' numClusters <- c(1, 5, 10, 20, 40, 50)
#'
#' results <- LocalControlClassic(data = lindner,
#'                                clusterVars = cvars,
#'                                treatmentColName = "abcix",
#'                                outcomeColName = "cardbill",
#'                                clusterCounts = numClusters)
#'
#' bxp <- UPSboxplot(results)
#'
#' @export
"UPSboxplot" <- function(envir, clusterSubset = c()){

  if(length(clusterSubset > 0)) {
    clusts = clusterSubset
  }
  else {
    clusts = envir[["summary"]]$Requested_clusters
  }

  ltdds = list()
  for(cc in clusts){
      objName = paste0("UPSnnltd", cc)
      nnObj = envir[[objName]]
      ltdds[[objName]] = nnObj$pbindif[nnObj$nnhbindf[,1]]
  }

  envir[["ltdds"]] = ltdds

  bxps = do.call("rbind",
                 mapply(function(X,Y) {data.frame(ltd = X, numclst = Y)},
                        X=ltdds, Y=clusts, SIMPLIFY = F, USE.NAMES = F))

  bxpstats <- boxplot(ltd ~ numclst,
                      data = bxps, xaxt = 'n',
                      col = "lightgray",
                      main="LTD Distribution depends upon Number of Clusters",
                      ylab="LTD Distribution",
                      xlab="Number of Clusters")

  axis(side = 1, labels = clusts, at = 1:length(clusts))

  bxpstats
}
