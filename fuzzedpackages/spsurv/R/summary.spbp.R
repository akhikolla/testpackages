#' Bernstein Polynomial Based Regression Object Summary
#'
#' @export
#' @param object an object of class spbp
#' @param interval interval coverage (confidence or credibility)
#' @param ... further arguments passed to or from other methods
#' @method summary spbp
#' @return An object of class analogous to for e.g. 'summary.bppo.bayes'.


summary.spbp <- function(object, interval = 0.95, ...){
    ## mle approach
    if(object$call$approach == "mle"){
        beta <- object$coefficients[1:object$q]
        var <- object$var[1:object$q, 1:object$q]

        ### Error handling ###
        # Null model
        if (is.null(object$coefficients)) return(object)

        beta2 <- beta[!(is.na(beta))] #non-missing coefs
        if(is.null(beta) | is.null(var )) stop("Input is not valid")

        se <- suppressWarnings(sqrt(diag(object$var)[1:object$q]))

        output <- list(call = object$call,
                       return_code = object$return_code,
                       n = object$n,
                       loglik = object$loglik)

        if (!is.null(object$nevent))
        output$nevent <- object$nevent

        output$coefficients  <- cbind(beta, exp(beta), se, beta/se,
                         pchisq((beta/ se)^2, 1, lower.tail=FALSE))
        dimnames(output$coefficients) <- list(names(beta), c("coef", "exp(coef)",
                                                 "se(coef)", "z", "Pr(>|z|)"))
        if (interval) {
            z <- qnorm((1 + interval)/2, 0, 1)
            output$interval <- cbind(exp(beta), exp(-beta), exp(beta - z * se),
                         exp(beta + z * se))
            dimnames(output$interval) <- list(names(beta), c("exp(coef)", "exp(-coef)",
                                                 paste("lower .", round(100 * interval, 2), sep = ""),
                                                 paste("upper .", round(100 * interval, 2), sep = "")))
        }
        df <- length(beta2)
        logtest <- -2 * (object$loglik[1] - object$loglik[2])
        output$logtest <- c(test=logtest,
                            df=df,
                            pvalue= pchisq(logtest, df, lower.tail=FALSE))
        output$rsq<-c(rsq=1-exp(-logtest/object$n),
                      maxrsq=1-exp(2*object$loglik[1]/object$n))
        output$waldtest<-c(test=as.vector(round(object$wald.test, 2)),
                           df=df,
                           pvalue= pchisq(as.vector(object$wald.test), df,
                                          lower.tail=FALSE))


        class(output) <- switch (object$call$model, "po"  = "summary.bppo.mle",
                                             "ph"  = "summary.bpph.mle",
                                             "aft" = "summary.bpaft.mle")
    }
    else{
        status <- eval(object$call$data)$status
        output <- list(call = object$call,
                       n = length(status),
                       loglik = rstan::extract(object$stanfit, "log_lik")$log_lik,
                       nevent = sum(status)
                       )

        output$coef_names <- colnames(model.matrix(object))
        aux <- rstan::summary(object$stanfit, probs = .5, pars = "beta")$summary
        exp_samp <- coda::mcmc(exp(rstan::extract(object$stanfit, "beta")$beta))

        output$summary_chain <- cbind(object$pmode[1:length( output$coef_names)],
                                      aux,
                                coda::HPDinterval(log(exp_samp), prob = interval))

        colnames(output$summary_chain) <- c("mode", colnames(aux), "lowerHPD", "upperHPD")
        rownames(output$summary_chain) <- output$coef_names
        #####

        output$summary_exp <- cbind(apply(exp_samp, 2, mean),
                                 apply(exp_samp, 2, median),
                                 apply(exp_samp, 2, sd),
                                 coda::HPDinterval(exp_samp, prob = interval)
                                 )
        rownames(output$summary_exp) <- rownames(output$summary_chain)
        colnames(output$summary_exp) <-  c("mean_exp", "median_exp", "sd_exp",
                                        "lowerHPD_exp", "upperHPD_exp")

        output$waic <- object$waic
        output$loo <- object$loo

        class(output) <- switch (object$call$model, "po"  = "summary.bppo.bayes",
                                 "ph"  = "summary.bpph.bayes",
                                 "aft" = "summary.bpaft.bayes")

    }
    return(output)
}
