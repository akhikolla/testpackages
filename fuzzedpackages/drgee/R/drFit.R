drFit <-
    function(object, rootFinder = findRoots, intercept, ...){

	if (!inherits(object, "drgeeData")) {
            stop("An object of class \"drgeeData\" is expected")
	}
        
        if (object$cond & object$olink == "logit") {
            
            return( drConditFit(object, rootFinder, intercept) )
            
        } else if (object$cond) {
            
            return( dreFitCond(object,
                               omodel = TRUE,
                               rootFinder = rootFinder, ...) )
            
        } else {
            
            return( dreFit(object, omodel = TRUE,
                           rootFinder = rootFinder, ...) )
            
        }
    }
