#' @title Identify shifts in the rate of trait diversification
#' @description Summarises phenotypic rate variation on phylogenies.
#' @param object Output of a traitMedusa analysis in \code{\link{transformPhylo.ML}}.
#' @param ... Additional arguments to be passed to \code{\link{summary}}.
#' @param cutoff Cutoff value for differences in AIC scores when comparing models. More complex models with an AIC score more than this number of units lower than simpler models are retained (as per runMedusa in the \pkg{geiger}).
#' @param AICc If true, AICc is used instead of AIC.
#' @param lowerBound Minimum value for parameter estimates.
#' @param upperBound Maximum value for parameter estimates.
#' @param print.warnings Logical. If TRUE, warnings are issued if confidence intervals fall outside upper or lower bounds.
#' @details This functions summarises the output of a "medusa" model in transformPhylo.ML (see below). The best overall model is chosen based on AIC (or AICc if AICc=TRUE). The cut-off point for improvement in AIC score between successively more complex models can be defined using cutoff. The default cutoff is 4 but this is somewhat arbitrary and a "good" cut-off may well vary between data sets so it may well be worth exploring different cutoffs.
#' @return ModelFit Summary of the best optimal rate shift model.
#' @return Rates Summary of the rate parameters from the best rate shift model.
#' @return optimalTree A phylo object with branch lengths scaled relative to rate.
#' @seealso \code{\link{transformPhylo.ML}}
#' @references Alfaro ME, Santini F, Brock CD, Alamillo H, Dornburg A, Carnevale G, Rabosky D & Harmon LJ. 2009. Nine exceptional radiations plus high turnover explain species diversity in jawed vertebrates. PNAS 106, 13410-13414.
#' @references O'Meara BC, Ane C, Sanderson MJ & Wainwright PC. 2006. Testing for different rates of continuous trait evolution using likelihood. Evolution 60, 922-933
#' @references Thomas GH, Meiri S, & Phillimore AB. 2009. Body size diversification in Anolis: novel environments and island effects. Evolution 63, 2017-2030.
#' @author Gavin Thomas
#' @examples 
#' data(anolis.tree)
#' data(anolis.data)
#' attach(anolis.data)
#' male.length <- matrix(Male_SVL, dimnames=list(rownames(anolis.data)))
#' sortedData <- sortTraitData(anolis.tree, male.length)
#' phy <- sortedData$phy
#' male.length <- sortedData$trait
#' phy.clade <- extract.clade(phy, 182)
#' male.length.clade <- as.matrix(male.length[match(phy.clade$tip.label, rownames(male.length)),])
#' tm1 <- transformPhylo.ML(male.length.clade, phy=phy.clade, model="tm1", minCladeSize=10, nSplits=1)
#' tm1_out <- summary(tm1, cutoff=1)
#' @export summary.traitMedusa
#' @export

summary.traitMedusa <- function (object, ..., cutoff=4, AICc=TRUE, lowerBound=1e-8, upperBound=200, print.warnings=FALSE) {
	
  if (class(object) != "traitMedusa")
    stop("Please sapply output from transformPhylo.ML with model = 'tm1' or model = 'tm2'.")
  
  traitMedusaObject <- object
  y <- traitMedusaObject[[3]]
  original.phy <- phy <- traitMedusaObject[[4]]
  breaks <- numeric()
  output.silent <- summary(y, ...)
  
  if (AICc == TRUE) {
    datCol = 6 
  } else {
  	datCol = 5
  }
  
  bestModel <- traitMedusaObject[[1]][1, ]
  rateType <- traitMedusaObject[[1]][2:nrow(traitMedusaObject[[1]]), 2]
  
  for (i in 2:dim(traitMedusaObject[[1]])[1]) {
  	if ((traitMedusaObject[[1]][i - 1, datCol] - traitMedusaObject[[1]][i, datCol]) < cutoff)
  	  break
  	bestModel <- traitMedusaObject[[1]][i,]
  	where.model <- i
  	}
  
  out <- vector(mode="list", length = 3)
  names(out) <- c("ModelFit", "Rates", "optimalTree")
  
  if (bestModel$node == 0) {
  	out$optimalTree <- phy 
  	out$ModelFit <- bestModel
  	out$Rates <- "Single rate"
  } else {
  	foo <- function(param) {
  	  ll <- transformPhylo.ll(y, phyClade, model = "clade", nodeIDs = SingleNode, cladeRates = param, rateType = whichRateType)$logLikelihood
  	  return(as.numeric(ll - as.numeric(bestModel["lnL"]) + 1.92))
  	  }
    nodeIDs <- traitMedusaObject[[1]][2:where.model, "node"]
    cladeRates <- traitMedusaObject[[1]][where.model, -c(1:6)]
    if(any(is.na(cladeRates)))
      cladeRates <- cladeRates[-which(is.na(cladeRates))]
    cladeRates <- as.numeric(cladeRates)
    rateType <- traitMedusaObject[[1]][2:where.model, "shiftPos"]
    optimalTree <- transformPhylo(phy, model="clade", nodeIDs=sort(nodeIDs), cladeRates=as.numeric(cladeRates), rateType=rateType)
    out$Rates <- matrix(NA, length(nodeIDs), 5, byrow = TRUE)
    out$Rates <- as.data.frame(out$Rates)
    colnames(out$Rates) <- c("node", "shiftPos", "MLRate", "LowerCI", "UpperCI")
    out$ModelFit <- bestModel[ ,3:6]
    out$optimalTree <- optimalTree
    for (i in 1:length(nodeIDs)) {
      SingleNode <- sort(nodeIDs[i])
      whichRateType <- rateType[i]
      phyClade <- transformPhylo(optimalTree, model = "clade", nodeIDs = SingleNode, cladeRates = 1 / cladeRates[i], rateType = whichRateType)
      LCI <- NULL
      UCI <- NULL
      out.lower <- try(uniroot(foo, interval = c(lowerBound, cladeRates[i]))$root, silent = TRUE)
      if (is.numeric(out.lower)) {
      	LCI <- uniroot(foo, interval = c(lowerBound, cladeRates[i]))$root 				
      } else {
      	LCI <- NA 
      }
      out.upper <- try(uniroot(foo, interval = c(cladeRates[i], upperBound))$root, silent = TRUE)
      if(is.numeric(out.upper)) {
      	UCI <- uniroot(foo, interval = c(cladeRates[i], upperBound))$root
      } else {
      	UCI <- NA
      }
      out$Rates[i,] <- aaa <- c(nodeIDs[i], rateType[i], cladeRates[i], LCI, UCI)
    }
		
	if (any(is.na(out$Rates[,3:4])) && print.warnings)
	  warning("Confidence limits fall outside the current parameter bounds - consider changing lowerBound and/or upperBound")
  }
  out$original.phy <- original.phy
  class(out) <- c("traitMedusa.model")
  return(out)
}