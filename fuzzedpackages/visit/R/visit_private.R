##constants in the package
get.const <- function() {
    rst <- list(REGIONS   = c("Toxic", "Ineffective", "Safe,Effective", "Effective,Safety concern"),
                THETA     = c("No DLT, No Response", "No DLT, Response",
                              "DLT, No Response", "DLT, Response"),
                CLSPRIOR  = "VTPRIOR",
                CLSPOST   = "VTPOST",
                CLSTRUEPS = "VTTRUEPS",
                CLSSIMU   = "VTSIMU",
                CLSDEC    = "VTDEC");

    rst$DENLEGEND <- c("Toxicity Rate", "Response Rate", rst$THETA);
    rst
}

##set options
set.option <- function(x, opt.x) {
    if (is.null(x)) {
        x <- opt.x;
    } else {
        x <- x[x %in% opt.x];
        if (0 == length(x))  {
            x <- opt.x;
        }
    }

    x
}

## p     = P(tox);
## q     = P(response);
## rho   = odds ratio p00p11/p10/p01
get.p11 <- function(pqr) {
    p   <- pqr[1];
    q   <- pqr[2];
    rho <- pqr[3];

    if (1 == rho) {
        p11 = p*q;
        p01 = (1-p)*q;
        p10 = p*(1-q);
        p00 = (1-p)*(1-q);
    } else {
        p11 = -(sqrt((p+q-p*rho-q*rho-1)^2-4*(rho-1)*p*q*rho)+(p+q-p*rho-q*rho-1))/2/(rho-1);
        p01 = q-p11;
        p10 = p-p11;
        p00 = p01*p10*rho/p11;
    }

    ##check
    if (0) {
        cp   <- p11 + p10;
        cq   <- p11 + p01;
        crho <- p11*p00/p01/p10;
        print(c(cp, cq, crho));
    }

    rst <- c(p00,p01,p10,p11);
    rst <- rst/sum(rst);

    rst
}

rdirichlet <- function (n, alpha) {
    l <- length(alpha);
    x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE);
    sm <- x %*% rep(1, l);

    return(x/as.vector(sm))
}

##get samples using the priors
get.prior.smp <- function(p.prior, model = 0, n=10000) {

    stopifnot(get.const()$CLSPRIOR %in% class(p.prior));

    ndose      <- nrow(p.prior);
    last.alpha <- 0;
    rst        <- list(NULL);
    for (dose.level in 1:ndose) {
        beta <- p.prior[dose.level,];

        ## toxcity rate
        tau   <- beta["TAU"];
        if (0 == model) {
            alpha      <- rnorm(n, last.alpha, beta["SDALPHA"]);
            last.alpha <- alpha;
        } else {
            alpha <- rnorm(n, beta["MEANALPHA"], beta["SDALPHA"]);
        }
        smp.tox <- tau^exp(alpha);

        ## response
        smp.rep <- rbeta(n, beta["A"], beta["B"]);

        ## odds ratio
        smp.r <- exp(rnorm(n, beta["C"], beta["D"]));

        ## theta ij
        smp.theta <- apply(cbind(smp.tox, smp.rep, smp.r), 1, get.p11);

        rst[[dose.level]] <- list(smp.tox=smp.tox,
                                  smp.rep=smp.rep,
                                  smp.r=smp.r,
                                  smp.theta=t(smp.theta));
    }
    rst
}

##combine simu results
lst.combine <- function(lst, var) {
    simplify2array(lapply(lst, function(x) {x[[var]]}));
}

##---------------------------------------------------
## project specic functions
##---------------------------------------------------

plot.densities <- function(lst.dens, draw.curves = 1:6, draw.levels=NULL,
                           legends=NULL, mfrow=NULL,
                           ltys = c(1,1,2,2,2,2),
                           cols = c("red", "blue", "brown", "black", "gray", "green"),
                           max.y.adjust = 1.1
                           ) {

    if (is.null(legends)) {
        legends <- get.const()$DENLEGEND;
    }

    draw.levels <- set.option(draw.levels, 1:length(lst.dens));
    draw.curves <- set.option(draw.curves, 1:6);

    if (is.null(mfrow)) {
        mfrow <- c(ceiling(length(draw.levels)/2), min(2, length(draw.levels)));
    }

    ##get all densities
    max.y <- 0;
    for (l in draw.levels) {
        cur.den  <- lst.dens[[l]];
        for (k in draw.curves) {
            max.y <- max(max.y, cur.den[[k]]$y);
        }
    }

    ## plot
    par(mfrow = mfrow);
    for (l in draw.levels) {
        plot(NULL, xlim=c(0,1), ylim=c(0, max.y.adjust*max(max.y)),
             xlab="Rate", ylab="Density", main=paste("Level", l, sep=" "),
             bty="n");
        for (k in draw.curves) {
            lines(lst.dens[[l]][[k]], col = cols[k], lty = ltys[k]);
        }

        if (draw.levels[1] == l) {
            legend("topright", bty="n", legend = legends[draw.curves],
                   lty = ltys[draw.curves], col = cols[draw.curves]);
        }
    }
}


plot.VTPRIOR <- function(x, model = 0, n = 10000, draw.levels = NULL, adjust = 1.2,
                         ...) {
    ##get samples
    priors <- x;
    draw.levels <- set.option(draw.levels, 1:nrow(priors));
    prior.smps  <- get.prior.smp(priors, model = model, n = n);

    ##get all densities
    lst.dens <- list(NULL);
    for (l in draw.levels) {
        cur.den      <- list(NULL);
        cur.den[[1]] <- density(prior.smps[[l]]$smp.tox, adjust = adjust, na.rm = TRUE);
        cur.den[[2]] <- density(prior.smps[[l]]$smp.rep, adjust = adjust, na.rm = TRUE);
        for (k in 1:4) {
            cur.den[[k+2]] <- density(prior.smps[[l]]$smp.theta[,k], adjust=adjust, na.rm = TRUE);
        }
        lst.dens[[l]] <- cur.den;
    }

    ## plot
    plot.densities(lst.dens, draw.levels = draw.levels, ...);
}

summary.VTPRIOR <- function(object, qs = c(0.025,0.5, 0.975), model = 0, n = 10000, ...) {
    f.sum <- function(l, lbl, smp) {
        c("Level" = l,
          "Prob"  = lbl,
          "Mean"  = mean(smp, na.rm = TRUE),
          quantile(smp, qs, na.rm = TRUE));
    }

    legends <- get.const()$DENLEGEND;

    ##get samples
    priors     <- object;
    prior.smps <- get.prior.smp(priors, model = model, n = n);

    ##get all densities
    rst <- NULL;
    for (l in 1:nrow(priors)) {
        rst <- rbind(rst, f.sum(l, legends[1], prior.smps[[l]]$smp.tox));
        rst <- rbind(rst, f.sum(l, legends[2], prior.smps[[l]]$smp.rep));
        for (k in 1:4) {
            rst <- rbind(rst,
                         f.sum(l, legends[k+2], prior.smps[[l]]$smp.theta[,k]));
        }
    }

    ## plot
    data.frame(rst);
}
