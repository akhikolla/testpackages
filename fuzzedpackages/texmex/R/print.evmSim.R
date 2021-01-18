#' @export
print.evmSim <- function(x , print.seed=FALSE, ...){
    print(x$call)
    print(x$map$family, verbose=FALSE)

    if (print.seed){
        cat("Random seed:", x$seed, "\n")
    }
    acc <- sapply(x$chains, function(X) signif(attr(X, "acceptance"), 3))
    acc <- paste(acc, collapse = ", ")
    cat("Acceptance rate: ", acc)

    cat("\n\nMAP estimates:\n")
    co <- coef(x)
    map <- x$map$coefficients
    names(map) <- names(co)

    print(map , ...)
    cat("\nPosterior means:\n")
    m <- coef(x)
    print(m , ...)

    invisible(x)
}
