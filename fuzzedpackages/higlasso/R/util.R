generate_design_matrices <- function(X, degree)
{
    generate.Xm <- function(i)
    {
        m <- splines::bs(X[, i], degree = degree)
        apply(m, 2, function(x) (x - mean(x)) / stats::sd(x))
    }

    Xm <- purrr::map(1:ncol(X), generate.Xm)
    Xi <- generate_Xi(Xm)
    decomp <- function(Xmi)
    {
        if (ncol(Xmi) > 0)
            apply(qr.Q(qr(Xmi)), 2, function(x) x / stats::sd(x))
        else
            Xmi
    }

    Xm <- purrr::map(Xm, decomp)
    Xi <- purrr::map(Xi, decomp)

    if (is.null(colnames(X)))
        names(Xm) <- paste0("V", seq(length(Xm)))
    else
        names(Xm) <- colnames(X)

    j <- purrr::map_lgl(Xi, ~ ncol(.x) > 0)

    ngroups <- length(Xm)
    groups <- purrr::flatten_dbl(c(
        purrr::imap(unname(Xm),    function(Xm.i, i) rep(i, ncol(Xm.i))),
        purrr::imap(Xi[j], function(Xi.i, i) rep(ngroups + i, ncol(Xi.i)))
    ))

    # construct "inverse" of groups
    igroups <- vector("list", max(groups))
    for (i in seq(length(groups)))
        igroups[[groups[i]]] <- c(igroups[[groups[i]]], i)

    X.xp <- do.call("cbind", c(unname(Xm), Xi[j]))


    list(Xm = Xm, Xi = Xi, X.xp = X.xp, groups = groups, igroups = igroups)
}

# Type checking functions to save space.
check.Y <- function(Y)
{
    name <- deparse(substitute(Y))
    if (!is.vector(Y) || !is.numeric(Y))
        stop("'", name ,"' must be a numeric vector.")
        if (any(is.na(Y)))
        stop("'", name, "' cannot contain missing values.")
}

check.XZ <- function(XZ, Y)
{
    name.XZ <- deparse(substitute(XZ))
    name.Y <- deparse(substitute(Y))
    if (!is.matrix(XZ) || !is.numeric(XZ))
        stop("'", name.XZ, "' must be a numeric matrix.")
    if (nrow(XZ) != length(Y))
        stop("The number of rows of '", name.XZ, "' does not match the length",
             " of '", name.Y, "'.")
    if (any(is.na(XZ)))
        stop("'", name.XZ, "' cannot contain missing values.")

}
