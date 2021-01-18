logLik.JMcmprsk <-
function (object, ...) {
    if (!inherits(object, "JMcmprsk"))
        stop("Use only with 'JMcmprsk' objects.\n")
    out <- object$loglike
    #attr(out, "df") <- nrow(object$Hessian)
    attr(out, "nobs") <- object$k
    class(out) <- "logLik"
    out
}
