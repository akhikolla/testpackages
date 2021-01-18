oFit <-
    function(object, inv = FALSE, ...) {

	if (!inherits(object, "drgeeData")) {
            stop("An object of class \"drgeeData\" is expected")
        }

        ## For retrospective logistic regression
        if(inv){

            ## For e-estimation with logit link, 
            ## let y and a switch place and replace v with z
            ## and run retrospective logistic regression

            if (object$olink != "logit" | object$elink != "logit")
                
                stop("\nReverse regression only possible in the logit-logit case\n")
            
            if (object$cond) {
                
                fit <- conditFit(y = object$a,
                                 x = cbind(object$yx, object$z),
                                 y.names = object$a.names,
                                 x.names = c(object$yx.names, object$z.names), 
                                 id = object$id)
                
            } else {
                
                fit <- geeFit(y = object$a,
                              x = cbind(object$yx, object$z),
                              link = object$olink)
                
            }

            coef.names <- c(object$yx.names, object$z.names)
            
        } else {
            
            if (object$cond) {
                
                fit <- geeFitCond(y = object$y, x = cbind(object$ax, object$v),
                                  link = object$olink, id = object$id, ...)
                
            } else {
                
                fit <- geeFit(y = object$y, x = cbind(object$ax, object$v), link = object$olink)
                
            }

            coef.names <- c(object$ax.names, object$v.names)

        }

        U <- fit$eq.x * fit$res
        d.U.sum <- crossprod( fit$eq.x , fit$d.res ) 

        coefficients <- fit$coefficients
        names(coefficients) <- coef.names

        result <- list(coefficients = coefficients,
                       coef.names = coef.names, 
                       U = U,
                       d.U.sum = d.U.sum,
                       optim.object = fit$optim.object,
                       optim.object.o = fit$optim.object,
                       optim.object.e = NULL,
                       id = object$id,
                       id.vcov = object$id.vcov)

        return(result)
    }
