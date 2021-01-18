gee <-
    function (formula,
              link = c("identity","log","logit"),
              data,
              subset, 
              cond = FALSE,
              clusterid,
              clusterid.vcov,
              rootFinder = findRoots,
              ...
              ) {

        call <- match.call()

        link <- match.arg(link)

        m <- match(c("data", "subset", "cond", "clusterid", "clusterid.vcov"), names(call), 0L)

        dD <- call[c(1L, m)]
        dD$oformula <- formula
        dD$olink <- link
        dD$estimation.method <- "o"

        dD[[1L]] <- quote(drgeeData)

        gee.data <- eval(dD, parent.frame())

        n.obs <- gee.data$n.obs

        if (gee.data$n.obs == 0) {
            stop("\nNo data\n")
        }
        
        if (length(gee.data$used.rows) == 0) {
            stop("\nNo data can be used for estimation\n")
        }
        
        n.vars <- ncol(gee.data$v)

        if (n.vars == 0) {
            stop("\nNo parameters to estimate\n")
        }
        
        if (cond) {
            
            fit <- geeFitCond(gee.data$y,
                              gee.data$v,
                              gee.data$y.names,
                              gee.data$v.names, 
                              link = link,
                              gee.data$id,
                              rootFinder,
                              ...)
            
        } else {
        
            fit <- geeFit(gee.data$y,
                          gee.data$v,
                          link = link)

        }

        coefficients = fit$coefficients
        names(coefficients) <- gee.data$v.names

        x <- gee.data$v[gee.data$orig.order,, drop = F]

        obs.names <- rownames(x)
        
        y <- as.vector(gee.data$y[gee.data$orig.order])
        names(y) <- obs.names

        res <- as.vector(fit$res[gee.data$orig.order])
        names(res) <- obs.names

        d.res <- fit$d.res[gee.data$orig.order,, drop = F]
        
        ## Calculate asymptotic variance        
        U <- fit$eq.x * fit$res
        d.U.sum <- crossprod( fit$eq.x , fit$d.res )
        vcov <- as.matrix( robustVcov(U, d.U.sum, gee.data$id.vcov) )
        dimnames(vcov) <- list(gee.data$v.names, gee.data$v.names)

        result <- list(coefficients = coefficients,
                       vcov = vcov,
                       call = call,
                       cond = cond,
                       y = y,
                       x = x,
                       gee.data = gee.data, 
                       optim.object = fit$optim.object,
                       U = U,
                       d.U.sum = d.U.sum,
                       res = res, 
                       d.res = d.res, 
                       formula = formula)

        if (!missing(data) ) {
            result$data <- data
        }
        
        class(result) <- "gee"

        return(result)
    }

print.gee <-
    function(x, digits = max(3L, getOption("digits") - 3L), ...) {
        if (length(x$coefficients)) {
            cat("\nCoefficients:\n")
            print.default(format(coef(x), digits = digits),
                          print.gap = 2, quote = FALSE)
            cat("\n")
        } else {
            cat("No coefficients\n\n")
        }

    }

summary.gee <-
    function(object, digits = max(3L, getOption("digits") - 3L), ...) {

	s.err <- sqrt(diag(as.matrix(vcov(object))))
	zvalue <- coef(object) / s.err
	pvalue <- 2 * pnorm(-abs(zvalue))

	coef.table <- as.matrix(cbind(coef(object), s.err, zvalue, pvalue))

	dimnames(coef.table) <- list(names(coef(object)), c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))

        summ <- summary(object$gee.data)

        if ( ( object$gee.data$cond & ncol(object$gee.data$v) == 0 ) | ( !object$gee.data$cond & ncol(object$gee.data$v) == 1 ) ) {
            model.formula <- paste( object$gee.data$y.names, " ~ 1", sep = "")
        } else {
            model.formula = summ$outcome.nuisance.model
        }
        
	ans <- list(call = object$call,
                    coefficients = coef.table,
                    vcov = vcov(object),
                    model.formula = model.formula,
                    link = summ$olink,
                    n.obs = summ$n.obs,
                    n.clust = summ$n.clust)

        class(ans) <- "summary.gee"
        return(ans)
    }

print.summary.gee <-
    function(x, digits = max(3L, getOption("digits") - 3L),
             signif.stars = getOption("show.signif.stars"), ...){
        cat("\nCall:  ",
            paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n", sep = "")

        

        cat("\nModel: ",x$model.formula,"\n")

        cat("\nLink function: ", x$link,"\n")

        if (length(x$coefficients)) {
            cat("\n")
            printCoefmat(x$coefficients, digits = digits, signif.stars = signif.stars,
                         na.print = "NA", ...)
            cat("\n", x$n.obs, " complete observations used\n")

            if (x$n.clust < x$n.obs){
                cat("\nCluster-robust Std. errors\n", x$n.clust, " clusters\n")
            }
        } else {
            cat("No coefficients estimated\n\n")
        }

    }

coef.gee <- function(object, ...) {
    return(object$coefficients)
}

vcov.gee <- function(object, ...) {
    
    return(object$vcov)
}

naiveVcov.gee <- function(object) {
    
    return( -solve( object$d.U.sum ) )
    
}

clusterRobustVcov.gee <- function(object, clusterid = NULL){

    if(is.null(clusterid)){
        
        return( object$vcov )
            
    }else{
            
        if ( is.character(clusterid) ) {
                
            clusterid <- get(clusterid, envir = parent.frame())
                
        }

        clusterid.reordered <- clusterid[object$gee.data$used.rows]
        
        return( robustVcov(object$U,
                           object$d.U.sum,
                           id = clusterid.reordered) )
    }
}

