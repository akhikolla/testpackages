
#' Fits a LOLOG model via Monte Carlo Generalized Method of Moments
#'
#'
#' @description
#' \code{lolog} is used to fit Latent Order Logistic Graph (LOLOG) models. LOLOG models are
#' motivated by the idea of network growth where the network begins empty, and edge variables
#' are sequentially 'added' to the network with an either unobserved, or partially observed
#' order \eqn{s}. Conditional upon the inclusion order, the probability of an edge has a
#' logistic relationship with the change in network statistics.
#'
#'
#' @param formula A lolog formula for the sufficient statistics (see details).
#' @param auxFormula A lolog formula of statistics to use for moment matching.
#' @param theta Initial parameters values. Estimated via \code{\link{lologVariational}} if NULL.
#' @param nsamp The number of sample networks to draw at each iteration.
#' @param includeOrderIndependent If TRUE, all order independent terms in formula are used for 
#' moment matching.
#' @param targetStats A vector of network statistics to use as the target for the moment equations.
#' If \code{NULL}, the observed statistics for the network are used.
#' @param weights The type of weights to use in the GMM objective. Either 'full' for the inverse 
#' of the full covariance matrix or 'diagonal' for the inverse of the diagonal of the covariance matrix.
#' @param tol The Hotelling's T^2 p-value tolerance for convergence for the transformed moment conditions.
#' @param nHalfSteps The maximum number of half steps to take when the objective is not improved 
#' in an iteration.
#' @param maxIter The maximum number of iterations.
#' @param minIter The minimum number of iterations.
#' @param startingStepSize The starting dampening of the parameter update.
#' @param maxStepSize The largest allowed value for dampening.
#' @param cluster A parallel cluster to use for graph simulation.
#' @param verbose Level of verbosity 0-3.
#'
#'
#' @details
#' LOLOG represents the probability of a tie, given the network grown up to a time point as
#' \deqn{
#'   \textrm{logit}\big(p(y_{s_t}=1 | \eta, y^{t-1}, s_{ \leq t})\big) = \theta \cdot c(y_{s_t}=1 | y^{t-1}, s_{ \leq t})
#' }
#' where \eqn{s_{\leq t}} is the growth order of the network up to time \eqn{t}, \eqn{y^{t-1}} is the
#' state of the graph at time \eqn{t-1}. \eqn{c(y_{s_t} | y^{t-1}, s_{ \leq t})} is a vector
#' representing the change in graph statistics from time \eqn{t-1} to \eqn{t} if an edge is present, and
#' \eqn{\theta} is a vector of parameters.
#'
#' The motivating growth order proceeds 'by vertex.' The network begins 'empty' and then vertices are 'added'
#' to the network sequentially. The order of vertex inclusion may be random or fixed. When a vertex 'enters' the
#' network, each of the edge variables connecting it and vertices already in the network are considered for
#' edge creation in a completely random order.
#'
#' LOLOG formulas contain a network, DirectedNet or UndirectedNet object on the left hand side.
#' the right hand side contains the model terms used. for example,
#'
#' \code{net ~ edges}
#'
#' represents and Erdos-Renyi model and
#'
#' \code{net ~ edges + preferentialAttachment()}
#'
#' represents a Barabasi-Albert model. See \code{\link{lolog-terms}} for a list of allowed model statistics
#'
#' Conditioning on (partial) vertex order can be done by
#' placing an ordering variable on the right hand side of the '|' operator, as in
#'
#' \code{net ~ edges + preferentialAttachment() | order}
#'
#' 'order' should be a numeric vector with as many elements as there are vertices in the network.
#' Ties are allowed. Vertices with higher order values will always be included later. Those with the same
#' values will be included in a random order in each simulated network.
#'
#' offsets and constraints are specified by wrapping them with either \code{offset()} or \code{constraint()},
#' for example, the following specifies an Erdos-Renyi model with the constraint that degrees must be less
#' that 10
#'
#' \code{net ~ edges + constraint(boundedDegree(0L, 10L))}
#'
#' If the model contains any order dependent statistics, additional moment constraints
#' must be specified in \code{auxFormula}. Ideally these should be chosen to capture
#' the features modeled by the order dependent statistics. For example, \code{preferentialAttachment}
#' models the degree structure, so we might choose two-stars as a moment constraint.
#'
#'  \code{lolog(net ~ edges + preferentialAttachment(), net ~ star(2))}
#'
#' will fit a Barabasi-Albert model with the number of edges and number of two-stars as moment constraints.
#'
#'
#' @return An object of class 'lolog'. If the model is dyad independent, the returned object will
#' also be of class "lologVariational" (see \code{\link{lologVariational}}, otherwise it will
#' also be a "lologGmm" object.
#'
#' lologGmm objects contain:
#'
#' \item{method}{"Method of Moments" for order independent models, otherwise "Generalized Method of Moments"}
#' \item{formula}{The model formula}
#' \item{auxFormula}{The formula containing additional moment conditions}
#' \item{theta}{The parameter estimates}
#' \item{stats}{The statistics for each network in the last iteration}
#' \item{estats}{The expected stats (G(y,s)) for each network in the last iteration}
#' \item{obsStats}{The observed h(y) network statistics}
#' \item{targetStats}{The target network statistics}
#' \item{obsModelStats}{The observed g(y,s) network statistics}
#' \item{net}{A network simulated from the fit model}
#' \item{grad}{The gradient of the moment conditions (D)}
#' \item{vcov}{The asymptotic covariance matrix of the parameter estimates}
#' \item{likelihoodModel}{An object of class *LatentOrderLikelihood at the fit parameters}
#'
#'
#' @examples
#' library(network)
#' set.seed(1)
#' data(flo)
#' flomarriage <- network(flo,directed=FALSE)
#' flomarriage %v% "wealth" <- c(10,36,27,146,55,44,20,8,42,103,48,49,10,48,32,3)
#'
#' # A dyad independent model
#' fit <- lolog(flomarriage ~ edges + nodeCov("wealth"))
#' summary(fit)
#'
#' # A dyad dependent model with 2-stars and triangles
#' fit2 <- lolog(flomarriage ~ edges + nodeCov("wealth") + star(2) + triangles, verbose=FALSE)
#' summary(fit2)
#'
#' \dontrun{
#'
#' # An order dependent model
#' fit3 <- lolog(flomarriage ~ edges + nodeCov("wealth") + preferentialAttachment(),
#'               flomarriage ~ star(2:3), verbose=FALSE)
#' summary(fit3)
#' 
#' # Try something a bit more real
#' data(ukFaculty)
#' 
#' # Delete vertices missing group
#' delete.vertices(ukFaculty, which(is.na(ukFaculty %v% "Group")))
#' 
#' fituk <- lolog(ukFaculty ~ edges() + nodeMatch("GroupC") + nodeCov("GroupC") + triangles + star(2))
#' summary(fituk)
#' plot(fituk$net, vertex.col= ukFaculty %v% "Group" + 2)
#'
#' }
#'
lolog <- function(formula,
                  auxFormula = NULL,
                  theta = NULL,
                  nsamp = 1000,
                  includeOrderIndependent = TRUE,
                  targetStats = NULL,
                  weights = "full",
                  tol = .1,
                  nHalfSteps = 10,
                  maxIter = 100,
                  minIter = 2,
                  startingStepSize = .1,
                  maxStepSize = .5,
                  cluster = NULL,
                  verbose = TRUE) {
  vcat <- function(..., vl=1){ if(verbose >= vl) cat(...) }
  vprint <- function(..., vl=1){ if(verbose >= vl) print(...) }
  
  #initialize theta via variational inference
  if (is.null(theta)) {
    vcat("Initializing Using Variational Fit\n")
    varFit <- lologVariational(formula, dyadInclusionRate = 1)
    if (varFit$allDyadIndependent) {
      vcat("Model is dyad independent. Returning maximum likelihood estimate.\n")
      return(varFit)
    }
    theta <- varFit$theta
    vcat("Initial Theta:\n", drop(theta),"\n")
  }
  
  lolik <- createLatentOrderLikelihood(formula, theta = theta)
  obsModelStats <- lolik$getModel()$statistics()
  statNames <- names(obsModelStats)
  
  orderIndependent <- lolik$getModel()$isIndependent(FALSE, TRUE)
  dyadIndependent <- lolik$getModel()$isIndependent(TRUE, TRUE)
  dyadIndependentOffsets <-
    lolik$getModel()$isIndependent(TRUE, FALSE)
  if (all(dyadIndependent) && all(dyadIndependentOffsets) && is.null(targetStats)) {
    vcat("Model is dyad independent. Returning maximum likelihood estimate.\n")
    varFit <- lologVariational(formula, dyadInclusionRate = 1)
    return(varFit)
  }
  obsModelStats[!orderIndependent] <- NA
  
  terms <- .prepModelTerms(formula)
  auxTerms <- .prepModelTerms(auxFormula)
  samp <- NULL
  obsStats <- NULL
  if (!is.null(auxFormula)) {
    auxModel <- createCppModel(auxFormula)
    auxModel$setNetwork(lolik$getModel()$getNetwork())
    auxModel$calculate()
    obsStats <- auxModel$statistics()
  }
  if (includeOrderIndependent) {
    obsStats <-
      c(lolik$getModel()$statistics()[orderIndependent], obsStats)
  }
  if(!is.null(targetStats)){
    if(any(is.na(targetStats)))
      stop("targetStats may not have missing values")
    if(length(obsStats)!=length(targetStats))
      stop("Incorrect length of the targetStats vector: should be ", 
           length(obsStats), " but is ",length(targetStats),".")
  }else{
    targetStats <- obsStats
  }
  if(length(targetStats) < length(theta))
    stop("Too few moment conditions. Specify more in auxFormula.")
  
  stepSize <- startingStepSize
  lastTheta <- NULL
  lastObjective <- Inf
  hsCount <- 0
  iter <- 0
  if (!is.null(cluster)) {
    tmpNet <- lolik$getModel()$getNetwork()$clone()
    tmpNet$emptyGraph()
    network <- as.network(tmpNet)
    clusterExport(cluster, "terms", envir = environment())
    clusterExport(cluster, "auxTerms", envir = environment())
    clusterExport(cluster, "network", envir = environment())
    clusterExport(cluster, "orderIndependent", envir = environment())
    clusterExport(cluster, "includeOrderIndependent", envir = environment())
    clusterEvalQ(cluster, {
      # Load lolog on each node
      library(lolog)
      
      # Pull from package while making R CMD check happy
      .createLatentOrderLikelihoodFromTerms <- eval(parse(text="lolog:::.createLatentOrderLikelihoodFromTerms"))
      .makeCppModelFromTerms <- eval(parse(text="lolog:::.makeCppModelFromTerms"))
      
      # Assign variables with external pointers
      enet <- as.BinaryNet(network)
      lolik2 <-
        .createLatentOrderLikelihoodFromTerms(terms, enet)
      if (!is.null(auxTerms))
        auxModel2 <- .makeCppModelFromTerms(auxTerms, enet)
      NULL
    })
  }
  while (iter < maxIter) {
    iter <- iter + 1
    vcat("\n\nIteration", iter,"\n") 
    
    #generate networks
    lolik$setThetas(theta)
    stats <- matrix(0, ncol = length(theta), nrow = nsamp)
    estats <- matrix(0, ncol = length(theta), nrow = nsamp)
    if (includeOrderIndependent)
      auxStats <-
      matrix(0,
             ncol = length(targetStats) - sum(orderIndependent),
             nrow = nsamp)
    else
      auxStats <- matrix(0, ncol = length(targetStats), nrow = nsamp)
    if (is.null(cluster)) {
      vcat("Drawing", nsamp, "Monte Carlo Samples:\n")
      if(verbose)
        pb <- utils::txtProgressBar(min = 0, max = nsamp, style = ifelse(interactive(),3,1))
      for (i in 1:nsamp) {
        if (verbose)
          utils::setTxtProgressBar(pb, i)
        samp <- lolik$generateNetwork()
        if (!is.null(auxFormula)) {
          auxModel$setNetwork(samp$network)
          auxModel$calculate()
          auxStats[i, ] <- auxModel$statistics()
        }
        stats[i, ] <- samp$stats + samp$emptyNetworkStats
        estats[i, ] <- samp$expectedStats + samp$emptyNetworkStats
      }
      if(verbose)
        close(pb)
      vcat("\n")
      if (includeOrderIndependent)
        auxStats <- cbind(stats[, orderIndependent], auxStats)
    } else{
      worker <- function(i, theta) {
        lolik2$setThetas(theta)
        samp <- lolik2$generateNetwork()
        if (!is.null(auxTerms)) {
          auxModel2$setNetwork(samp$network)
          auxModel2$calculate()
          as <- auxModel2$statistics()
        } else{
          as <- numeric()
        }
        list(
          stats = samp$stats + samp$emptyNetworkStats,
          estats = samp$expectedStats + samp$emptyNetworkStats,
          auxStats = as
        )
      }
      results <-
        parallel::parLapply(cluster, 1:nsamp, worker, theta = theta)
      stats <- t(sapply(results, function(x)
        x$stats))
      estats <- t(sapply(results, function(x)
        x$estats))
      colnames(stats) <- colnames(estats) <- statNames
      if (!is.null(auxFormula))
        auxStats <- t(sapply(results, function(x)
          x$auxStats))
      else
        auxStats <- NULL
      if (includeOrderIndependent)
        auxStats <- cbind(stats[, orderIndependent], drop(auxStats))
      colnames(auxStats) <- names(obsStats)
    }
    
    # Calculate gradient of moment conditions
    grad <- matrix(0, ncol = length(theta), nrow = length(targetStats))
    for (i in 1:length(targetStats)) {
      for (j in 1:length(theta)) {
        grad[i, j] <-
          -(cov(auxStats[, i], stats[, j]) - cov(auxStats[, i], estats[, j]))
      }
    }
    
    if (weights == "diagonal")
      W <- diag(1 / (diag(var(auxStats))))
    else{
      tr <- try(W <- solve(var(auxStats)))
      if(inherits(tr, "try-error")){
        warning("Singular statistic covariance matrix. Using diagnoal.")
        W <- diag(1 / (diag(var(auxStats))))
      }
    }
    
    # Calculate moment conditions and stat/observed stat differences transformed by W.
    mh <- colMeans(auxStats)
    diffs <- -sweep(auxStats, 2, targetStats)
    transformedDiffs <- t(t(grad) %*% W %*% t(diffs))
    momentCondition <- colMeans(transformedDiffs)
    
    objective <- colMeans(diffs) %*% W %*% colMeans(diffs)
    vcat("Objective: ", drop(objective), "\n")
    
    objCrit <-
      max(-1000000, objective - lastObjective) / (lastObjective + 1)
    
    # Calculate inverse
    invFailed <-
      inherits(try(gradInv <-
                     solve(t(grad) %*% W %*% grad), silent = TRUE)
               ,"try-error")
    
    if (verbose >= 3){
      ns <- nrow(stats)
      os <- obsModelStats
      pairs(rbind(stats, os), pch='.', cex=c(rep(1, ns), 10),
	    col=c(rep("black", ns), "red"), diag.panel = .panelHist,
            main=paste("Iteration",iter),cex.main=0.9)
    }
    vcat("Statistic GM Skewness Coef :", apply(stats,2, .gmSkewness),"\n", vl=2)
    
    # If inverse failed, or the objective has increased significantly, initiate half stepping
    if (hsCount < nHalfSteps &&
        !is.null(lastTheta) && (invFailed || objCrit > .3)) {
      vcat("Half Step Back\n")
      theta <- (lastTheta + theta) / 2
      hsCount <- hsCount + 1
      stepSize <- stepSize / 2
      vcat("Theta:", drop(theta), "\n", vl=2)
      next
    } else{
      stepSize <- min(maxStepSize, stepSize * 1.25)
      hsCount <- 0
    }
    
    vcat("Step size:", stepSize, "\n", vl=2)
    
    # Update theta
    lastTheta <- theta
    theta <- theta - stepSize * gradInv %*% momentCondition
    lastObjective <- objective
    
    # Print diagnostic Information
    algoState <- data.frame(lastTheta, theta,momentCondition)
    colnames(algoState) <- c("Theta","Next Theta","Moment Conditions")
    rownames(algoState) <- statNames
    if(!is.null(auxFormula)){
      momState <- data.frame(colMeans(diffs)/sqrt(diag(var(diffs))))
      colnames(momState) <- "(h(y) - E(h(Y))) / sd(h(Y))"
      rownames(momState) <- names(obsStats)
      vprint(algoState, vl=2)
      vprint(momState, vl=2)
    }else{
      algoState[["(h(y) - E(h(Y))) / sd(h(Y))"]] <- colMeans(diffs)/sqrt(diag(var(diffs)))
      vprint(algoState, vl=2)
    }
    
    
    #Hotelling's T^2 test
    hotT <-
      momentCondition %*% solve(var(transformedDiffs) / nrow(transformedDiffs)) %*% momentCondition
    pvalue <- pchisq(hotT, df = length(theta), lower.tail = FALSE)
    
    vcat("Hotelling's T2 p-value: ", format.pval(pvalue,digits=5,eps=1e-5), "\n")
    
    if (pvalue > tol && iter >= minIter) {
      break
    } else if (iter < maxIter) {
      
    }
  }
  
  if (is.null(samp)) {
    samp <- lolik$generateNetwork()
  }
  
  # Calculate parameter covariances
  omega <- var(auxStats)
  vcov <- solve(t(grad) %*% W %*% grad) %*%
    t(grad) %*% W %*% omega %*% t(W) %*% grad %*%
    solve(t(grad) %*% t(W) %*% grad)
  
  # Some formatting of return items
  lastTheta <- drop(lastTheta)
  rownames(grad) <- colnames(auxStats) <- names(obsStats)
  rownames(vcov) <- colnames(vcov) <- colnames(grad) <- 
    colnames(stats) <- colnames(estats) <- names(lastTheta) <- statNames
  
  method <-
    if (is.null(auxFormula))
      "Method of Moments"
  else
    "Generalized Method Of Moments"
  
  result <- list(
    method = method,
    formula = formula,
    auxFromula = auxFormula,
    theta = lastTheta,
    stats = stats,
    estats = estats,
    auxStats = auxStats,
    obsStats = obsStats,
    targetStats = targetStats,
    obsModelStats = obsModelStats,
    net = samp$network,
    grad = grad,
    vcov = vcov,
    likelihoodModel = lolik
  )
  class(result) <- c("lologGmm", "lolog", "list")
  result
}

#' Print a `lolog` object
#' @param x the object
#' @param ... additional parameters (unused)
#' @method print lolog
print.lolog <- function(x, ...) {
  cat(x$method, "Coefficients:\n")
  print(x$theta)
}


#' Summary of a `lolog` object
#' @param object the object
#' @param ... additional parameters (unused)
#' @method summary lolog
#' @examples
#' data(lazega)
#' fit <- lologVariational(lazega ~ edges() + nodeMatch("office") + triangles, 
#'                         nReplicates=50L, dyadInclusionRate=1)
#' summary(fit)
#' @method summary lolog
summary.lolog <- function(object, ...) {
  x <- object
  theta <- x$theta
  se <- sqrt(diag(x$vcov))
  pvalue <- 2 * pnorm(abs(theta / se), lower.tail = FALSE)
  stats <- x$likelihoodModel$getModel()$statistics()
  orderInd <- x$likelihoodModel$getModel()$isIndependent(FALSE, TRUE)
  stats[!orderInd] <- NA
  result <-
    data.frame(
      observed_statistics = stats,
      theta = theta,
      se = se,
      pvalue = round(pvalue, 4)
    )
  rownames(result) <- names(stats)
  result
}


#' Generates BinaryNetworks from a fit lolog object
#'
#'
#' @param object A `lolog` object.
#' @param nsim The number of simulated networks
#' @param seed Either NULL or an integer that will be used in a call to set.seed before simulating
#' @param convert convert to a network object#'
#' @param ... unused
#' 
#' @return A list of BinaryNet (or network if convert=TRUE) objects. Networks contain an additional
#' vertex covariate "__order__" that indicates the sequence order in which the vertex was 'added' 
#' into the network.
#'
#'
#' @examples
#' library(network)
#' data(flo)
#' flomarriage <- network(flo,directed=FALSE)
#' flomarriage %v% "wealth" <- c(10,36,27,146,55,44,20,8,42,103,48,49,10,48,32,3)
#' fit <- lolog(flomarriage ~ edges + nodeCov("wealth"))
#' net <- simulate(fit)[[1]]
#' plot(net)
#'
#' @method simulate lolog
simulate.lolog <- function(object, nsim = 1, seed = NULL, convert = FALSE, ...) {
  if (!is.null(seed))
    set.seed(seed)
  l <- list()
  for (i in 1:nsim) {
    l[[i]] <- object$likelihoodModel$generateNetwork()$network
    if (convert)
      l[[i]] <- as.network(l[[i]])
  }
  l
}


#' Extracts estimated model coefficients.
#' 
#' @param object A `lolog` object.
#' @param ... unused
#' @examples
#' # Extract parameter estimates as a numeric vector:
#' data(ukFaculty)
#' fit <- lolog(ukFaculty ~ edges)
#' coef(fit)
#' @method coef lolog
coef.lolog <- function(object, ...){
  object$theta
}



#' Conduct Monte Carlo diagnostics on a lolog model fit
#' 
#' This function creates simple diagnostic
#' plots for MC sampled statistics produced from a lolog fit.
#' 
#' Plots are produced that represent the distributions of the 
#' output sampled statistic values or the target statistics values.
#' The values of the observed target statistics for the networks are
#' also represented for comparison with the sampled statistics.
#'
#' @param x A model fit object to be diagnosed.
#' @param type The type of diagnostic plot. "histograms", the default, produces histograms of the sampled
#' output statistic values with the observed statistics represented by vertical lines. "target" produces a pairs plot of the 
#' target output statistic values with the pairs of observed target statistics represented by red squares.
#' output statistic values with the observed statistics represented by vertical lines. "model" produces a pairs plot of the 
#' sampled output statistic values with the pairs of observed statistics represented by red squares.
#' @param ... Additional parameters. Passed to \link{geom_histogram} if type="histogram" 
#' and \link{pairs} otherwise.
#' 
#' @examples 
#' library(network)
#' set.seed(1)
#' data(flo)
#' flomarriage <- network(flo,directed=FALSE)
#' flomarriage %v% "wealth" <- c(10,36,27,146,55,44,20,8,42,103,48,49,10,48,32,3)
#'
#'
#' # An order dependent model
#' fit3 <- lolog(flomarriage ~ edges + nodeCov("wealth") + preferentialAttachment(),
#'               flomarriage ~ star(2:3), verbose=FALSE)
#' plot(fit3)
#' plot(fit3, "target")
#' plot(fit3, "model")
plot.lologGmm <- function(x, type=c("histograms", "target","model"), ...) {
  type <- match.arg(type, c("histograms", "target","model"))
  
  if(type == "target"){
    stats <- x$auxStats
    ns <- nrow(stats)
    stats <- rbind(stats, x$targetStats)
    pch <- '.'
    cex <- c(rep(1, ns), 10)
    col <- c(rep("black", ns), "red")
    pairs(stats, pch=pch, cex=cex, col=col, ...)
  }else if(type == "model"){
    stats <- x$stats
    ns <- nrow(stats)
    stats <- rbind(stats, x$obsModelStats)
    pch <- '.'
    cex <- c(rep(1, ns), 10)
    col <- c(rep("black", ns), "red")
    pairs(stats, pch=pch, cex=cex, col=col, diag.panel = .panelHist, ...)    
  }else if(type == "histograms"){
    stats <- x$auxStats
    obs <- x$targetStats
    nin <- !(colnames(x$stats) %in% colnames(stats))
    stats <- cbind(x$stats[,nin,drop=FALSE], stats)
    obs <- c(x$obsModelStats[nin], obs)
    
    Var2 <- value <- NULL # for R CMD check
    ss <- reshape2::melt(stats)
    oo <- na.omit(data.frame(Var2=names(obs),value=obs))
    pp <- ggplot2::ggplot(data=ss, ggplot2::aes(x=value)) + 
      ggplot2::geom_histogram(...) + 
      ggplot2::geom_vline(ggplot2::aes(xintercept=value), data=oo, color=I("red"), size=I(1)) + 
      ggplot2::facet_wrap(~Var2, scales="free") + ggplot2::xlab("") + ggplot2::theme_bw()
    invisible(print(pp))
  }else
    stop("unknown plot type")
  invisible()
}
