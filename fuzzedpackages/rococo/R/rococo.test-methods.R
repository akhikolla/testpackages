rococo.test.numeric <- function(x, y, similarity=c("linear", "exp", "gauss",
                                                   "epstol", "classical"),
                                tnorm="min", r=0, numtests=1000,
                                storeValues=FALSE, exact=FALSE,
                                alternative=c("two.sided", "less", "greater"),
                                noVarReturnZero=TRUE)
{
    if (!is.numeric(x) || !is.numeric(y) || length(x) != length(y))
        stop("'x' and 'y' need to be numeric vectors of the same length")

    if (!is.numeric(r) || any(r < 0))
        stop("'r' must be a (vector of) non-negative real number(s)")

    if (missing(similarity))
        similarity <- "linear"

    if (!is.logical(exact) || length(exact) > 1)
        stop("'exact' must be single logical")

    if (exact)
    {
        if (length(x) > 10)
        {
            warning("exact test intractable for more than 10 samples;",
                    " setting 'exact=FALSE'")
            exact <- FALSE
        }
        else
            numtests <- factorial(length(x))
    }

    similarity <- match.arg(similarity, several.ok=TRUE)
    alternative <- match.arg(alternative)
    altId <- switch(alternative, two.sided=0, less=1, greater=2)

    if (length(r) == 1)
        r[2] <- r[1]

    if (!is.numeric(numtests) || length(numtests) > 1 || numtests < 1 ||
        numtests != floor(numtests))
        stop("'numtests' must by a single positive integer")
    else if (numtests < 100 && !exact)
        warning("for the sake of significance, it is not recommend to use",
                "numtests < 100")

    if (!is.logical(storeValues) || length(storeValues) > 1)
        stop("'storeValues' must be single logical")

    tnormSel <- -1

    if (is.character(tnorm))
    {
        tnorm <- match.arg(tnorm, c("min", "prod", "lukasiewicz"))
        tnormSel <- switch(tnorm, min=1, prod=2, lukasiewicz=3)
        tnlist <- list(name=tnorm)
    }
    else if (is.function(tnorm))
    {
        if (length(formals(tnorm)) != 2)
            stop("'tnorm' should be a function of two arguments, ",
                 "e.g. 'tnorm=function(a, b) a * b' for the product t-norm")

        if (abs(tnorm(1, 0.5) - 0.5) > .Machine$double.eps ||
            abs(tnorm(0, 0.25)) > .Machine$double.eps)
            stop("supplied function does not appear to be a valid t-norm")

        if (!is.null(attr(tnorm, "name")))
            tnlist <- list(name=attr(tnorm, "name"), def=tnorm)
        else
            tnlist <- list(name="user-defined t-norm", def=tnorm)

        if (requireNamespace("compiler", quietly=TRUE))
            tnorm <- compiler::cmpfun(tnorm)
    }
    else
        stop("'tnorm' should be valid string or a function of two arguments, ",
             "e.g. 'tnorm=function(a, b) a * b' for the product t-norm")

    if (length(similarity) > 1 && similarity[1] != similarity[2])
    {
        if (similarity[1] != "classical" && r[1] == 0)
        {
            r[1] <- 0.1 * IQR(x)

            if (r[1] == 0)
            {
                warning("IQR(x) = 0 => ",
                        "could not determine tolerance; ",
                        "reverting to classical crisp similarity")
                similarity[1] <- "classical"
            }
        }

        if (similarity[2] != "classical" && r[2] == 0)
        {
            r[2] <- 0.1 * IQR(y)

            if (r[2] == 0)
            {
                warning("IQR(y) = 0 => ",
                        "could not determine tolerance; ",
                        "reverting to classical crisp similarity")
                similarity[2] <- "classical"
            }
        }

        Rx <- switch(similarity[1],
                     linear=.Call("rcor_matrix_linear",
                                  vx=as.double(x), as.double(r[1])),
                     exp=.Call("rcor_matrix_exp",
                               vx=as.double(x), as.double(r[1])),
                     gauss=.Call("rcor_matrix_gauss",
                                 vx=as.double(x), as.double(r[1])),
                     epstol=.Call("rcor_matrix_epstol",
                                  vx=as.double(x), as.double(r[1])),
                     classical=.Call("rcor_matrix_classical",
                                     vx=as.double(x), as.double(r[1])))
        Ry <- switch(similarity[2],
                     linear=.Call("rcor_matrix_linear",
                                  vx=as.double(y), as.double(r[2])),
                     exp=.Call("rcor_matrix_exp",
                               vx=as.double(y), as.double(r[2])),
                     gauss=.Call("rcor_matrix_gauss",
                                 vx=as.double(y), as.double(r[2])),
                     epstol=.Call("rcor_matrix_epstol",
                                  vx=as.double(y), as.double(r[2])),
                     classical=.Call("rcor_matrix_classical",
                                     vx=as.double(y), as.double(r[2])))
    }
    else
    {
        if (similarity[1] != "classical")
        {
            if (r[1] == 0)
            {
                r[1] <- 0.1 * IQR(x)

                if (r[1] == 0)
                {
                    warning("IQR(x) = 0 => ",
                            "could not determine tolerance; ",
                            "reverting to classical crisp similarity")
                    similarity[1] <- "classical"
                }
            }

            if (r[2] == 0)
            {
                r[2] <- 0.1 * IQR(y)

                if (r[2] == 0)
                {
                    warning("IQR(y) = 0 => ",
                            "could not determine tolerance; ",
                            "reverting to classical crisp similarity")
                    similarity[2] <- "classical"
                }
            }
        }

        if (length(similarity) > 1 && similarity[1] != similarity[2])
        {
            Rx <- switch(similarity[1],
                         linear=.Call("rcor_matrix_linear",
                                      vx=as.double(x), as.double(r[1])),
                         exp=.Call("rcor_matrix_exp",
                                   vx=as.double(x), as.double(r[1])),
                         gauss=.Call("rcor_matrix_gauss",
                                     vx=as.double(x), as.double(r[1])),
                         epstol=.Call("rcor_matrix_epstol",
                                      vx=as.double(x), as.double(r[1])),
                         classical=.Call("rcor_matrix_classical",
                                         vx=as.double(x), as.double(r[1])))
            Ry <- switch(similarity[2],
                         linear=.Call("rcor_matrix_linear",
                                      vx=as.double(y), as.double(r[2])),
                         exp=.Call("rcor_matrix_exp",
                                   vx=as.double(y), as.double(r[2])),
                         gauss=.Call("rcor_matrix_gauss",
                                     vx=as.double(y), as.double(r[2])),
                         epstol=.Call("rcor_matrix_epstol",
                                      vx=as.double(y), as.double(r[2])),
                         classical=.Call("rcor_matrix_classical",
                                         vx=as.double(y), as.double(r[2])))
        }
        else
        {
            matrices <- switch(similarity[1],
                               linear=.Call("rcor_matrices_linear",
                                            vx=as.double(x),
                                            vy=as.double(y),
                                            as.double(r[1]),
                                            as.double(r[2])),
                               exp=.Call("rcor_matrices_exp",
                                         vx=as.double(x),
                                         vy=as.double(y),
                                         as.double(r[1]),
                                         as.double(r[2])),
                               gauss=.Call("rcor_matrices_gauss",
                                           vx=as.double(x),
                                           vy=as.double(y),
                                           as.double(r[1]),
                                           as.double(r[2])),
                               epstol=.Call("rcor_matrices_epstol",
                                            vx=as.double(x),
                                            vy=as.double(y),
                                            as.double(r[1]),
                                            as.double(r[2])),
                               classical=.Call("rcor_matrices_classical",
                                               vx=as.double(x),
                                               vy=as.double(y),
                                               as.double(r[1]),
                                               as.double(r[2])))

            Rx <- matrices$Rx
            Ry <- matrices$Ry
        }
    }

    if (tnormSel > 0)
    {
        res <- .Call("rcor", Rx, Ry, as.integer(tnormSel))
        c <- res$c
        d <- res$d
    }
    else
    {
        c <- sum(mapply(tnorm, Rx, Ry))
        d <- sum(mapply(tnorm, Rx, t(Ry)))
    }

    if (!is.na(c) && !is.na(d) && c + d > 0)
        oldgamma <- (c - d) / (c + d)
    else if (identical(noVarReturnZero, TRUE))
        oldgamma <- 0
    else
    {
        oldgamma <- NA
        warning("no variation in at least one of the two observables")
    }

    if (is.na(oldgamma))
    {
        out <- new("RococoTestResults",
                   count=as.integer(0),
                   tnorm=tnlist,
                   input=paste(deparse(substitute(x, env=parent.frame())), "and",
                               deparse(substitute(y, env=parent.frame()))),
                   length=length(x),
                   p.value=1,
                   p.value.approx=1,
                   r.values=r[1:2],
                   numtests=as.integer(numtests),
                   exact=as.logical(exact),
                   similarity=similarity,
                   sample.gamma=as.numeric(NA),
                   H0gamma.mu=as.numeric(NA),
                   H0gamma.sd=as.numeric(NA),
                   perm.gamma=as.numeric(NA),
                   alternative=alternative)

        return(out)
    }

    ## Run tests
    cnt <- 0

    if (tnormSel < -1)
    {
        samples <- vector(mode="numeric", length=numtests)

        if (exact)
        {
            perm <- 1:length(x)
            sigt <- rep(as.integer(-1), length(x))

            i <- 1
            repeat
            {
                c <- sum(mapply(tnorm, Rx, Ry[perm, perm]))
                d <- sum(mapply(tnorm, Rx, t(Ry[perm, perm])))

                if (!is.na(c) && !is.na(d) && c + d > 0)
                    newgamma <- (c - d) / (c + d)
                else
                    newgamma <- 0

                samples[i] <- newgamma

                if (identical(altId, 0) && abs(newgamma) >= abs(oldgamma))
                    cnt <- cnt + 1
                else if (identical(altId, 1) && newgamma <= oldgamma)
                    cnt <- cnt + 1
                else if (identical(altId, 2) && newgamma >= oldgamma)
                    cnt <- cnt + 1

                ret <- .Call("permNextWrapper", perm, sigt)

                if (is.null(ret))
                    break
                else
                {
                    perm <- ret$perm
                    sigt <- ret$sign
                }

                i <- i + 1
            }
        }
        else
        {
            i <- 1
            while (i <= numtests)
            {
                perm <- sample.int(length(x))
                c <- sum(mapply(tnorm, Rx, Ry[perm, perm]))
                d <- sum(mapply(tnorm, Rx, t(Ry[perm, perm])))

                if (!is.na(c) && !is.na(d) && c + d > 0)
                    newgamma <- (c - d) / (c + d)
                else
                    newgamma <- 0

                samples[i] <- newgamma

                if (identical(altId, 0) && abs(newgamma) >= abs(oldgamma))
                    cnt <- cnt + 1
                else if (identical(altId, 1) && newgamma <= oldgamma)
                    cnt <- cnt + 1
                else if (identical(altId, 2) && newgamma >= oldgamma)
                    cnt <- cnt + 1

                i <- i + 1
            }
        }

        sampleMU <- mean(samples)
        sampleSD <- sd(samples)
    }
    else
    {
        if (exact)
            res <- .Call("rcor_exacttest", Rx, Ry, as.integer(tnormSel),
                         as.integer(numtests), as.double(oldgamma),
                         as.integer(altId), storeValues)
        else
            res <- .Call("rcor_permtest", Rx, Ry, as.integer(tnormSel),
                         as.integer(numtests), as.double(oldgamma),
                         as.integer(altId), storeValues)

        cnt <- res$cnt
        sampleMU <- res$H0mu
        sampleSD <- res$H0sd
        if (storeValues)
            samples <- res$values
    }

    pval <- cnt / numtests

    if (alternative == "greater")
        pval2 <- pnorm(oldgamma, mean=sampleMU, sd=sampleSD, lower.tail=FALSE)
    else if (alternative == "less")
        pval2 <- pnorm(oldgamma, mean=sampleMU, sd=sampleSD, lower.tail=TRUE)
    else
        pval2 <- 2 * pnorm(abs(oldgamma), mean=sampleMU, sd=sampleSD,
                           lower.tail=FALSE)

    new("RococoTestResults",
        count=as.integer(cnt),
        tnorm=tnlist,
        input=paste(deparse(substitute(x, env=parent.frame())), "and",
                    deparse(substitute(y, env=parent.frame()))),
        length=length(x),
        p.value=pval,
        p.value.approx=pval2,
        r.values=r[1:2],
        numtests=as.integer(numtests),
        exact=as.logical(exact),
        similarity=similarity,
        sample.gamma=oldgamma,
        H0gamma.mu=sampleMU,
        H0gamma.sd=sampleSD,
        perm.gamma=if (storeValues) samples
                   else vector(mode="numeric", length=0),
        alternative=alternative)
}

setMethod("rococo.test", signature(x="numeric", y="numeric"),
          rococo.test.numeric)


rococo.test.formula <- function(x, y, na.action, ...)
{
    if (length(x) != 2L)
        stop("formula invalid")
    m <- match.call(expand.dots = FALSE)
    if (is.matrix(eval(m$y, parent.frame())))
        m$y <- as.data.frame(y)
    m[[1L]] <- as.name("model.frame")
    m$... <- NULL
    m$formula <- m$x
    m$data <- m$y
    m$x <- NULL
    m$y <- NULL
    mf <- eval(m, parent.frame())
    if (length(mf) != 2L)
        stop("formula invalid")
    DNAME <- paste(names(mf), collapse = " and ")

    ret <- rococo.test(mf[[1]], mf[[2]], ...)
    ret@input <- DNAME
    ret
}

setMethod("rococo.test", signature(x="formula", y="data.frame"),
          rococo.test.formula)
