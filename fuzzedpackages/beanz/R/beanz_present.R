
##-----------------------------------------------------------------------------------
##                     comparison
##-----------------------------------------------------------------------------------
#' Comparison of posterior treatment effects
#'
#' Present the difference in the posterior treatment effects
#' between subgroups
#'
#'
#' @name bzComp
#'
#' @return \code{bzSummaryComp} generates a data frame with summary statistics
#'     of the difference of treatment effects between the selected subgroups.
#'     \code{bzPlotComp} generates the density plot of the difference in the
#'     posterior treatment effects between subgroups. \code{bzForestComp}
#'     generates the forest plot of the difference in the posterior treatment
#'     effects between subgroups.
#'
#'
#' @examples
#' \dontrun{
#' var.cov    <- c("sodium", "lvef", "any.vasodilator.use");
#' var.resp   <- "y";
#' var.trt    <- "trt";
#' var.censor <- "censor";
#' resptype   <- "survival";
#' var.estvar <- c("Estimate", "Variance");
#'
#' subgrp.effect <- bzGetSubgrpRaw(solvd.sub,
#'                              var.resp   = var.resp,
#'                              var.trt    = var.trt,
#'                              var.cov    = var.cov,
#'                              var.censor = var.censor,
#'                              resptype   = resptype);
#'
#' rst.sr     <- bzCallStan("sr", dat.sub=subgrp.effect,
#'                          var.estvar=var.estvar, var.cov = var.cov,
#'                          par.pri=c(B=1000, C=1000),
#'                          chains=4, iter=500,
#'                          warmup=100, thin=2, seed=1000);
#'
#' sel.grps <- c(1,4,5);
#' tbl.sub <- bzSummaryComp(rst.sr, sel.grps=sel.grps);
#' bzPlot(rst.sr, sel.grps = sel.grps);
#' bzForest(rst.sr, sel.grps = sel.grps);}
#'
#' @seealso \code{\link{bzCallStan}}
#'
#'
NULL


#' @rdname bzComp
#'
#' @inheritParams bzSummary
#'
#' @param seed random seed
#'
#' @export
#'
bzSummaryComp <- function(stan.rst, sel.grps=NULL, cut=0, digits=3, seed = NULL) {

    if (!is.null(seed))
        set.seed(seed);

    stopifnot(is(stan.rst, "beanz.stan"));

    mus      <- stan.rst$get.mus();
    sel.grps <- get.sel.subgrp(mus, sel.grps);
    if (length(sel.grps) < 2)
        return(NULL);

    mus  <- mus[,sel.grps, drop=FALSE];
    fsum <- function(x) {
        crst <- c(mean(x),
                  sd(x),
                  quantile(x, c(0.025, 0.25, 0.5, 0.75, 0.975)),
                  mean(x<cut));
        crst <- round(crst, digits);
    }

    rst <- NULL;
    for (i in 1:(ncol(mus)-1)) {
        for (j in (i+1):ncol(mus)) {
            cur.j <- sample(mus[,j], nrow(mus), TRUE);
            cur.i <- sample(mus[,i], nrow(mus), TRUE);
            rst   <- rbind(rst, fsum(cur.j-cur.i));
            rownames(rst)[nrow(rst)] <- paste("Subgroup ", sel.grps[j], "-", sel.grps[i], sep="");
        }
    }


    rst <- cbind(rownames(rst), rst);
    colnames(rst) <- c("Comparison", "Mean", "SD", "Q025", "Q25",
                       "Median", "Q75", "Q975", paste("ProbLT", cut, sep=""));
    rownames(rst) <- NULL;
    data.frame(rst);
}



#' @rdname bzComp
#'
#'
#' @inheritParams bzPlot
#'
#'
#' @export
#'
bzPlotComp <- function(stan.rst, sel.grps=NULL, ..., seed = NULL) {

    if (!is.null(seed))
        set.seed(seed);

    stopifnot(is(stan.rst, "beanz.stan"));
    mus      <- stan.rst$get.mus();
    sel.grps <- get.sel.subgrp(mus, sel.grps);

    if (length(sel.grps) > 5 | length(sel.grps) < 2)
        return(NULL);

    mus     <- mus[,sel.grps, drop=FALSE];
    lst.den <- NULL;
    cmp.mus <- NULL;
    for (i in 1:(ncol(mus)-1)) {
        for (j in (i+1):ncol(mus)) {
            cur.j <- sample(mus[,j], nrow(mus), TRUE);
            cur.i <- sample(mus[,i], nrow(mus), TRUE);
            cmp.mus <- cbind(cmp.mus, cur.j-cur.i);
            lst.den[[length(lst.den)+1]]    <- density(cur.j - cur.i);
            names(lst.den)[length(lst.den)] <- paste("Subgroup ", sel.grps[j], "-", sel.grps[i], sep="");
        }
    }

    plot.densities(list(lst.den=lst.den, mus=cmp.mus),
                   main="Comparison of Subgroup Effects",
                   xlab="Difference of Treatment Effects",
                   ...);
}

#' @rdname bzComp
#'
#' @inheritParams bzForest
#'
#' @export
#'
bzForestComp <- function(stan.rst, sel.grps=NULL, ..., quants=c(0.025,0.975), seed = NULL) {
    stopifnot(is(stan.rst, "beanz.stan"));

    if (!is.null(seed))
        set.seed(seed);

    mus      <- stan.rst$get.mus();
    sel.grps <- get.sel.subgrp(mus, sel.grps);

    if (length(sel.grps) > 5 | length(sel.grps) < 2)
        return(NULL);

    mus      <- mus[,sel.grps];
    mu.qmean <- NULL;
    for (i in 1:(ncol(mus)-1)) {
        for (j in (i+1):ncol(mus)) {
            cur.j <- sample(mus[,j], nrow(mus), TRUE);
            cur.i <- sample(mus[,i], nrow(mus), TRUE);
            cur.m <- cur.j-cur.i;
            cq    <- quantile(cur.m, quants);

            mu.qmean <- rbind(mu.qmean, c(cq[1], mean(cq), cq[2]));
            rownames(mu.qmean)[nrow(mu.qmean)] <- paste("Subgroup ", sel.grps[j], "-", sel.grps[i], sep="");
        }
    }

    plot.forest(mu.qmean, main="Difference of Subgroup Effects Forest Plot", ...);
}


##-----------------------------------------------------------------------------------
##                     posterior summary
##-----------------------------------------------------------------------------------

#' Posterior subgroup treatment effects
#'
#' Present the posterior subgroup treatment effects
#'
#' @name bzSummary
#'
#' @return \code{bzSummary} generates a dataframe with summary statistics
#'     of the posterior treatment effect for the selected subgroups.
#'     \code{bzPlot} generates the density plot of the posterior treatment
#'     effects for the selected subgroups. \code{bzForest}
#'     generates the forest plot of the posterior treatment
#'     effects.
#'
#'@examples
#' \dontrun{
#' sel.grps <- c(1,4,5);
#' tbl.sub <- bzSummary(rst.sr, ref.stan.rst=rst.nse, ref.sel.grps=1);
#' bzPlot(rst.sr, sel.grps = sel.grps, ref.stan.rst=rst.nse, ref.sel.grps=1);
#' bzForest(rst.sr, sel.grps = sel.grps, ref.stan.rst=rst.nse, ref.sel.grps=1);}
#'
#' @seealso \code{\link{bzCallStan}}
#'
#'
NULL


#' @rdname bzSummary
#'
#' @param cut cut point to compute the probabiliby that the posterior subgroup
#'     treatment effects is below
#'
#' @param digits number of digits in the summary result table
#'
#' @inheritParams bzPlot
#'
#' @export
#'
bzSummary <- function(stan.rst, sel.grps=NULL, ref.stan.rst=NULL, ref.sel.grps=1, cut=0, digits=3) {

    stopifnot(is(stan.rst, "beanz.stan"));

    s.mus <- get.all.mus(stan.rst, sel.grps, ref.stan.rst, ref.sel.grps);

    rst <- apply(s.mus, 2,
                 function(x) {
        crst <- c(mean(x),
                  sd(x),
                  quantile(x, c(0.025, 0.25, 0.5, 0.75, 0.975)),
                  mean(x<cut)
                  );
        crst <- round(crst, digits);
    });

    rst <- t(rst);
    rst <- cbind(colnames(s.mus), rst);

    colnames(rst) <- c("Subgroup", "Mean", "SD", "Q025", "Q25", "Median", "Q75",
                       "Q975", paste("ProbLT", cut, sep=""));
    rownames(rst) <- NULL;
    data.frame(rst)
}

#' @rdname bzSummary
#'
#' @param stan.rst a class \code{beanz.stan} object generated by
#'     \code{\link{bzCallStan}}
#'
#' @param sel.grps an array of subgroup numbers to be included in the summary results
#'
#' @param ref.stan.rst a class \code{beanz.stan} object from \code{\link{bzCallStan}} that
#'     is used as the reference
#'
#' @param ref.sel.grps subgroups from the reference model to be included in the
#'     summary table
#'
#' @param ... options for \code{plot} function
#'
#' @export
#'
bzPlot <- function(stan.rst, sel.grps=NULL,
                        ref.stan.rst=NULL, ref.sel.grps=1,
                        ... ) {

    stopifnot(is(stan.rst, "beanz.stan"));

    s.mus          <- get.all.mus(stan.rst, sel.grps, ref.stan.rst, ref.sel.grps);
    lst.den        <- apply(s.mus, 2, density);
    names(lst.den) <- colnames(s.mus);

    plot.densities(list(mus=s.mus,lst.den=lst.den),
                   main="Posterior Distribution of Subgroup Effects",
                   xlab="Treatment Effect",
                   ...);
}

#' @rdname bzSummary
#' @inheritParams bzPlot
#' @param quants lower and upper quantiles of the credible intervals in the
#'     forest plot
#'
#' @export
#'
bzForest <- function(stan.rst, sel.grps=NULL,
                     ref.stan.rst=NULL, ref.sel.grps=1,
                     ..., quants=c(0.025,0.975)) {

    stopifnot(is(stan.rst, "beanz.stan"));

    s.mus    <- get.all.mus(stan.rst, sel.grps, ref.stan.rst, ref.sel.grps);
    mu.q     <- apply(s.mus, 2, quantile, quants);
    mu.mean  <- apply(s.mus, 2, mean);
    mu.qmean <- cbind(mu.q[1,], mu.mean, mu.q[2,]);
    rownames(mu.qmean) <- colnames(s.mus);

    plot.forest(mu.qmean, main="Subgroup Effects Forest Plot", ...);
}

#' Predictive Distribution
#'
#' Get the predictive distribution of the subgroup treatment effects
#'
#' @inheritParams bzCallStan
#' @inheritParams bzSummary
#'
#' @return A dataframe of predicted subgroup treament effects. That is, the
#'     distribution of \deqn{ \theta_g | \widehat{\theta}_1, \widehat{\sigma}^2_1, \ldots,
#'     \widehat{\theta}_G, \widehat{\sigma}^2_G.}
#'
#' @examples
#' \dontrun{
#' var.cov    <- c("sodium", "lvef", "any.vasodilator.use");
#' var.resp   <- "y";
#' var.trt    <- "trt";
#' var.censor <- "censor";
#' resptype   <- "survival";
#' var.estvar <- c("Estimate", "Variance");
#'
#' subgrp.effect <- bzGetSubgrp(solvd.sub,
#'                                   var.resp   = var.resp,
#'                                   var.trt    = var.trt,
#'                                   var.cov    = var.cov,
#'                                   var.censor = var.censor,
#'                                   resptype   = resptype);
#'
#' rst.nse    <- bzCallStan("nse", dat.sub=subgrp.effect,
#'                          var.estvar = var.estvar, var.cov = var.cov,
#'                          par.pri = c(B=1000),
#'                          chains=4, iter=4000,
#'                          warmup=2000, thin=2, seed=1000);
#'
#' pred.effect <- bzPredSubgrp(rst.nes,
#'                             dat.sub = solvd.sub,
#'                             var.estvar = var.estvar);}
#' @export
#'
bzPredSubgrp <- function(stan.rst, dat.sub, var.estvar) {

    stopifnot(class(stan.rst) == "beanz.stan");

    mus       <- stan.rst$get.mus();
    prior.sig <- stan.rst$prior.sig;
    delta     <- stan.rst$delta;
    sigma2    <- dat.sub[, var.estvar[2]];
    nsub      <- nrow(dat.sub);

    f.sig <- function() {
        if (0 == prior.sig) {
            eps <- runif(nsub, -delta, delta);
        } else {
            eps <- rnorm(nsub, 0, sqrt(delta))
        }
        lsig <- log(sqrt(sigma2)) + eps;
        exp(lsig);
    }

    rst <- apply(mus, 1, function(x) {
        cur.sig <- f.sig();
        cur.rst <- rnorm(nsub, x, cur.sig);
        cur.rst
    })

    t(rst)
}

##-----------------------------------------------------------------------------------
##                     SUMMARY REPORT
##-----------------------------------------------------------------------------------

#' Summary table of treatment effects
#'
#' Compare the DIC from different models and report the summary of treatment effects
#' based on the model with the smallest DIC value
#'
#' @param lst.stan.rst list of class \code{beanz.stan} results from
#'     \code{\link{bzCallStan}} for different models
#'
#' @inheritParams bzCallStan
#' @inheritParams bzSummary
#'
#' @return A dataframe with summary statistics of the model selected by DIC
#'
#' @export
#'
bzRptTbl <- function(lst.stan.rst, dat.sub, var.cov, cut=0, digits=3) {

    if (is.null(dat.sub) | is.null(lst.stan.rst))
        return(NULL);

    looic <- NULL;
    for (i in 1:length(lst.stan.rst)) {
        looic <- c(looic, lst.stan.rst[[i]]$looic$looic);
    }

    min.mdl <- which.min(looic);
    mus     <- lst.stan.rst[[min.mdl]]$get.mus();
    rst     <- apply(mus, 2,
                     function(x) {
        crst <- c(mean(x),
                  sd(x),
                  mean(x < cut)
                  );
        crst <- round(crst, digits);
    });

    rst           <- t(rst);
    colnames(rst) <- c("Mean", "SD", paste("Prob < ", cut, sep=""));
    rst           <- cbind(Model = lst.stan.rst[[min.mdl]]$mdl,
                           dat.sub[, c("Subgroup", var.cov)],
                           rst);

    rst
}

##-----------------------------------------------------------------------------------
##                             Frequentist GailSimon Test
##-----------------------------------------------------------------------------------
#' Gail-Simon Test
#'
#' Gail-Simon qualitative interaction test.
#'
#' @param effects subgroup treatment effects
#'
#' @param sderr standard deviation of the estimated treatment effects
#'
#' @param d clinically meaningful difference
#'
#'
#' @examples
#' \dontrun{
#' var.cov    <- c("sodium", "lvef", "any.vasodilator.use");
#' var.resp   <- "y";
#' var.trt    <- "trt";
#' var.censor <- "censor";
#' resptype   <- "survival";
#' subgrp.effect <- bzGetSubgrp(solvd.sub,
#'                                   var.resp   = var.resp,
#'                                   var.trt    = var.trt,
#'                                   var.cov    = var.cov,
#'                                   var.censor = var.censor,
#'                                   resptype   = resptype);
#'
#' gs.pval <- bzGailSimon(subgrp.effect$Estimate,
#'                        subgrp.effect$Variance); }
#'
#'
#' @export
#'
bzGailSimon <- function(effects, sderr, d=0) {
    d  <- abs(d);
    I  <- length(effects);
    Qm <- sum((effects > d) * ((effects-d)/sderr)^2);
    Qp <- sum((effects < -d) * ((effects+d)/sderr)^2);
    test.stats <- min(Qm, Qp);
    pval       <- sum(dbinom(1:(I-1), I-1, .5) * (1 - pchisq(test.stats, df=1:(I-1))));
    pval
}



##------------------------------------------------------------------
##------------tool kit----------------------------------------------
##------------------------------------------------------------------

plot.densities <- function(lst.mus.den,
                           main="", xlab="",
                           ylims=NULL,
                           xlims=NULL,
                           cols=c(rep(c("black", "red", "green", "brown", "gray", "yellow",
                                        "cyan", "blue", "ivory", "wheat"), 4),
                                  colors()),
                           ltys=c(rep(1:4, each=10),
                                  rep(1:6, each=100)),
                           quants=c(0.025,0.975)) {

    mus     <- lst.mus.den$mus;
    lst.den <- lst.mus.den$lst.den;

    ##get min max of x and y in dens
    f.tmp <- function(cur.den, cur.mu) {
        max.y <- max(cur.den$y);
        min.x <- quantile(cur.mu, quants[1]);
        max.x <- quantile(cur.mu, quants[2]);
        rst   <- c(max.y, min.x, max.x);
    }

    ##boundary
    lims <- NULL;
    for (i in 1:length(lst.den)) {
        cur.lim <- f.tmp(lst.den[[i]], mus[,i]);
        lims    <- cbind(lims, cur.lim);
    }

    if (is.null(ylims)) {
        max.y <- max(lims[1,]);
        ylims <- c(0,1.1*max.y);
    }

    if (is.null(xlims)) {
        min.x <- min(lims[2,]);
        max.x <- max(lims[3,]);
        xlims <- c(min.x,max.x)
    }

    plot(NULL, ylim=ylims, xlim=xlims,
         main=main,
         xlab=xlab,
         ylab="Density");

    ns <- length(lst.den);
    for (i in 1:ns) {
        lines(lst.den[[i]]$x,
              lst.den[[i]]$y,
              col=cols[i],  lty=ltys[i], lwd=2);
    }

    legend("topleft",
           legend=names(lst.den),
           lty=ltys[1:ns],
           col=cols[1:ns], bty="n", cex=1.2);
}




##plot forest plot
plot.forest <- function(quants,
                        cut=0,
                        y.bottom=0.1, y.top=0.95,
                        x.labs=0.2,
                        main="", xlab="") {

    ##where to put labels
    ng    <- nrow(quants);
    ys    <- seq(y.top, y.bottom, length.out=ng);
    lab.g <- rownames(quants);

    ##convert
    xrs <- range(c(quants, cut));
    xrs <- range(xrs, 0.1*(xrs[2]-xrs[1]), -0.1*(xrs[2]-xrs[1]));

    xs <- c(cut, seq(xrs[1], xrs[2], length.out=4));
    xs <- round(xs, digits=1);

    f.con <- function(x) {
        (x - xrs[1])/(xrs[2] - xrs[1])*(1-x.labs) + x.labs;
    }

    ##plot
    plot(NULL, ylim=c(0,1), xlim=c(0,1), main=main, xlab=xlab, ylab="", axes=FALSE);
    lines(f.con(c(cut,cut)), c(0,1), lty=2, lwd=2, col="gray");
    axis(1, at=f.con(xs), labels=xs);

    for (i in 1:nrow(quants)) {
        cur.q <- quants[i,];
        text(x.labs, ys[i], lab.g[i], pos=2);
        lines(f.con(cur.q[c(1,3)]), c(ys[i], ys[i]), lwd=2, lty=1, col="gray");
        points(f.con(cur.q[2]), ys[i], pch=23);
    }
}

get.sel.subgrp <- function(mus, sel.grps=NULL) {

    smus <- as.numeric(1:ncol(mus));

    if (is.null(sel.grps) |
        !is.numeric(sel.grps)) {
        sel.grps <- smus;
    } else {
        sel.grps <- as.numeric(sel.grps);
        sel.grps <- sel.grps[sel.grps %in% smus];
        if (0 == length(sel.grps))
            sel.grps <- smus;
    }

    sel.grps;
}


##comibine all mus for presentation
get.all.mus <- function(stan.rst, sel.grps=NULL, ref.stan.rst=NULL, ref.sel.grps=1) {

    stopifnot(class(stan.rst) == "beanz.stan");
    stopifnot(is.null(ref.stan.rst) | class(ref.stan.rst) == "beanz.stan");

    mus      <- stan.rst$get.mus();
    sel.grps <- get.sel.subgrp(mus, sel.grps);
    s.mus    <- mus[, sel.grps, drop=FALSE];

    if (!is.null(ref.stan.rst)) {

        ref.mus <- ref.stan.rst$get.mus();
        if (is.null(ref.sel.grps))
            ref.sel.grps <- sel.grps;

        r.mus <- ref.mus[, ref.sel.grps, drop=FALSE];
        colnames(r.mus) <- paste(ref.stan.rst$mdl, "(", ref.sel.grps, ")", sep="");
        s.mus <- cbind(s.mus, r.mus);
    }

    s.mus
}
