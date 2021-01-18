geeFit <-
    function(y,
             x,
             link = c("identity", "log", "logit")) {

	link <- match.arg(link)

        eq.x <- x

	if (link == "identity") {
            family <- gaussian()
	} else if (link == "log") {
            family = quasipoisson()
	} else if(link == "logit") {
            family = quasibinomial()
	}
        
        fit <- try( glm.fit(x = x,
                            y = y,
                            family = family) )

        if (inherits(fit, "try-error")) { 

            d.res = matrix(rep(NA,nrow(x) * ncol(x)), nrow = nrow(x))

            return( list(coefficents = rep(NA, ncol(x)),
                         res = rep(NA,nrow(y)),
                         d.res = d.res,
                         eq.x = eq.x,
                         optim.object = NULL))

        } else {

            d.res <- -x * family$mu.eta( eta = fit$linear.predictors )
            
            colnames(d.res) <- colnames(x)

            return(list(coefficients = fit$coefficients,
                        res = as.vector(y - fit$fitted.values),
                        d.res = d.res,
                        eq.x = x,
                        optim.object = NULL))
        }

    }
