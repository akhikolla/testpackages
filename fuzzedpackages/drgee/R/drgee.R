drgee <-
    function(outcome,
             exposure,
             oformula,
             eformula,
             iaformula = formula(~1),
             olink = c("identity", "log", "logit"),
             elink = c("identity", "log", "logit"),
             data,
             subset = NULL, 
             estimation.method = c("dr", "o", "e"),
             cond = FALSE,
             clusterid,
             clusterid.vcov,
             rootFinder = findRoots,
             intercept = TRUE, 
             ...
             ) {

        call <- match.call()
        olink <- match.arg(olink)
        elink <- match.arg(elink)
        estimation.method <- match.arg(estimation.method)

        if (estimation.method != "o" & olink == "logit" & elink != "logit") {
            warning("\nFor dr- and e-estimation, olink=\"logit\" can only be combined with elink=\"logit\"\nelink changed to \"logit\"\n")
            elink <- "logit"
        }

        if (cond & missing(clusterid)) {
            stop("\nFor conditional estimation, a clusterid is required")
        }

        if (missing(eformula) & missing(exposure)) {
            stop("An exposure needs to be specified\n\n")
        }

        if (missing(oformula) & missing(outcome)) {
            stop("An outcome needs to be specified\n\n")
        }


        m <- match(c("oformula", "olink", "outcome",
                     "iaformula",
                     "eformula", "elink", "exposure",
                     "cond",
                     "data",
		     "subset",
                     "estimation.method", 
                     "clusterid",
                     "clusterid.vcov"),
                   names(call), 0L)

        dD <- call[c(1L, m)]

        dD[[1L]] <- quote(drgeeData)

        drgee.data <- eval(dD, parent.frame())

        n.obs <- drgee.data$n.obs

        if (drgee.data$n.obs == 0) {
            stop("\nNo data\n")
        }
        
        if (length(drgee.data$used.rows) == 0) {
            stop("\nNo data can be used for estimation\n")
        }
        
        if (estimation.method == "o") {
            fit <- oFit(drgee.data)
        } else if (estimation.method == "e") {
            fit <- eFit(drgee.data, rootFinder, ...)
        } else if (estimation.method == "dr") {
            fit <- drFit(drgee.data, rootFinder, intercept, ...)
        }

        ## Calculate variance based on the estimating equations
        vcov <- as.matrix( robustVcov(fit$U, fit$d.U.sum, drgee.data$id.vcov) )

        dimnames(vcov) <- list( fit$coef.names, fit$coef.names)
        
        fit$call <- call

        fit$coefficients.all <- fit$coefficients
        fit$coefficients <- fit$coefficients[1:ncol(drgee.data$ax)]

        fit$vcov.all <- vcov
        fit$vcov <- vcov[1:ncol(drgee.data$ax), 1:ncol(drgee.data$ax), drop = F]

        fit$drgee.data <- drgee.data

        fit$estimation.method <- estimation.method

        if (!missing(eformula)) {
            fit$eformula <- eformula
        }

        if (!missing(oformula)) {
            fit$oformula <- oformula
        }

        if (!missing(iaformula)) {
            fit$iaformula <- iaformula
        }

        if (!missing(data) ) {
            fit$data <- data
        }
        
        class(fit) <- c("drgee")

        return(fit)
    }

print.drgee <-
    function(x, digits = max(3L, getOption("digits") - 3L), ...){
        if(length(x$coefficients)) {
            cat("\nCoefficients for main effect:\n")
            print.default(format(coef(x), digits = digits),
                          print.gap = 2, quote = FALSE)
            cat("\n")
        } else {
            cat("No coefficients\n\n")
        }
    }

summary.drgee <-
    function(object, digits = max(3L, getOption("digits") - 3L), ...){

	s.err <- sqrt(diag(as.matrix(vcov(object))))
	zvalue <- coef(object) / s.err
	pvalue <- 2 * pnorm(-abs(zvalue))

	coef.table <- as.matrix(cbind(coef(object), s.err, zvalue, pvalue))

	dimnames(coef.table) <- list(names(coef(object)), c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))

	model.info <- summary(object$drgee.data)

	ans <- list(call = object$call, coefficients = coef.table,
                    vcov = vcov(object), estimation.method = object$estimation.method,
                    model.info = model.info, n.obs = model.info$n.obs,
                    n.clust = model.info$n.clust)

        class(ans) <- "summary.drgee"
        return(ans)
    }

print.summary.drgee <-
    function(x, digits = max(3L, getOption("digits") - 3L),
             signif.stars = getOption("show.signif.stars"), ...){
        cat("\nCall:  ",
            paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n", sep = "")
        cat("\nOutcome: ", x$model.info$outcome, "\n")
        cat("\nExposure: ", x$model.info$exposure, "\n")
        cat("\nCovariates: ", paste( x$model.info$covariates, collapse  = ","), "\n")

        cat("\nMain model: ", x$model.info$main.model,"\n")

        if (x$estimation.method != "e") {
            cat("\nOutcome nuisance model: ",
                x$model.info$outcome.nuisance.model, "\n")
        }

        cat("\nOutcome link function: ", x$model.info$olink, "\n")

        if (x$estimation.method != "o") {
            cat("\nExposure nuisance model: ",
                x$model.info$exposure.nuisance.model, "\n")
            cat("\nExposure link function: ", x$model.info$elink, "\n")
        }

        if (length(x$coefficients)) {

            cat("\n")

            printCoefmat(x$coefficients, digits = digits, signif.stars = signif.stars,
                         na.print = "NA", ...)

            cat("\n(Note: The estimated parameters quantify the conditional\nexposure-outcome association, given the covariates\nincluded in the nuisance models)\n")

            cat("\n", x$n.obs, " complete observations used\n")

            if (x$n.clust < x$n.obs) {
                cat("\nCluster-robust Std. errors\nusing ", x$n.clust,
                    " clusters defined by levels of ",  x$model.info$clustname, "\n")
            }
        } else {
            cat("No coefficients estimated\n\n")
        }

    }

coef.drgee <- function(object, ...) {
    return(object$coefficients)
}

vcov.drgee <- function(object, ...) {
    
    return(object$vcov)

}

## naiveVcov.drgee <- function(object) {
    
##     d.U <- object$d.U.sum / object$drgee.data$n.obs
##     return( -solve( d.U ) )
    
## }

clusterRobustVcov.drgee <- function(object, clusterid = NULL){

    n.params <- length( object$coefficients )
    
    if( is.null(clusterid) ){
        
        return( object$vcov.all[1:n.params, 1:n.params])
            
    } else {
            
        if ( is.character(clusterid) ) {
                
            clusterid <- get(clusterid, envir = parent.frame())
                
        }

        clusterid.reordered <- clusterid[object$drgee.data$used.rows]
        
        return( robustVcov(object$U,
                           object$d.U.sum,
                           id = clusterid.reordered)[1:n.params, 1:n.params] )
    }
}


