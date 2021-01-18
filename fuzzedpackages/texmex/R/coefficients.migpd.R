#' @export
coef.migpd <-
function(object, ...){
    co <- sapply(object$models, coef)
    up <- sapply(object$models, endPoint, verbose=FALSE)
    co <- rbind(object$mth, object$mqu, co, up)
    dimnames(co) <- list(c("Threshold", "P(X < threshold)", "sigma", "xi", "Upper end point"),
                           names(object$models))
    co[3, ] <- exp(co[3, ])
    co
}


