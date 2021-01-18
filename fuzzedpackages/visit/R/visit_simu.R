#' S3 Summary function
#'
#' @param x object
#' @param ... reserved parameters
#' 
#' @export
summary2 <- function (x, ...) {
    UseMethod("summary2", x)
}

#' Set simulation scenario
#'
#' Simulation function. Get true \eqn{\theta}'s using marginal probabilities and
#' odds ratio \eqn{\rho} for all dose levels.
#'
#' @rdname vtScenario
#'
#' @param tox Vector of marginal DLT risk rates for all levels
#' @param res Vector of marginal immune response rates for all levels
#' @param rho Vector of odds ratio for all levels. If length of \code{rho} is
#'     shorter than the length of \code{tox} or \code{res}, vector \code{rho} is
#'     repeated to have the same length as \code{tox} and \code{res}.
#'
#' @details
#'
#' The calculation is as following. If \eqn{\rho = 1}, then \eqn{\theta_{11} =
#' pq}, \eqn{\theta_{01} = (1-p)q}, \eqn{\theta_{10} = p(1-q)}, and
#' \eqn{\theta_{00} = (1-p)(1-q)}. Otherwise, \eqn{ \theta_{11} = -(\sqrt{A+B}},
#' \eqn{\theta_{01} = q-\theta_{11}}, \eqn{\theta_{10} = p-\theta_{11}}, and
#' \eqn{\theta_{00} = \theta_{01}\theta_{10}\rho/\theta_{11}}, where
#' \eqn{A=(p+q-p \rho-q\rho-1)^2-4(\rho-1)pq\rho)} and
#' \eqn{B=(p+q-p\rho-q\rho-1))/2/(\rho-1)}.
#'
#'
#' @return a \code{VTTRUEPS} object containing all \eqn{\theta}'s in a matrix
#'     with its number of rows equaling the number of dose levels and its number
#'     of columns being 4.
#'
#' @examples
#'  rst.sce <- vtScenario(tox=c(0.05, 0.05, 0.08), res=c(0.2, 0.3, 0.5), rho=1)
#'
#' @export
#'
vtScenario <- function(tox = c(0.05, 0.05, 0.08),
                       res = c(0.2, 0.3, 0.5),
                       rho = 1) {

    if (!is.vector(rho)) {
        rho <- rep(rho, length(tox));
    }

    true.theta <- apply(cbind(tox, res, rho),
                        1,
                        function(x) {
        tt      <- get.p11(x);
        tt[4]   <- 1 - sum(tt[1:3]);
        tt
    });


    rst           <- t(true.theta);
    colnames(rst) <- get.const()$THETA;
    class(rst)    <- get.const()$CLSTRUEPS;

    invisible(rst)
}


#' Get prior distribution parameters
#'
#' Get prior distribution parameters for partially parametric or partially parametric+ models
#'
#' @rdname vtPriorPar
#'
#' @param prior.y Historical data for generating prior parameters. It has the
#'     same structure as \code{obs.y} in \code{\link{vtPost}}.
#' @param tau Vector of \eqn{\tau} values. See \code{\link{visit}} for details.
#'     Can not be \code{NULL} if \code{prior.y} is \code{NULL}.
#' @param sdalpha \eqn{\sigma_\alpha}. See \code{\link{visit}} for details.
#' @param sdrho \eqn{\sigma_\rho}.
#' @param vtheta Additional variance term for eliciting prior parameters from
#'     \code{prior.y}
#'
#' @details
#'
#' The priors are specified as \eqn{q^{(l)} \sim Beta(a_q^{(l)}, b_q^{(l)})},
#' and \eqn{\log\rho^{(l)} \sim N(0, \sigma_\rho^2)}.
#'
#' @return A \code{VTPRIOR} list with
#'
#' \itemize{
#'
#' \item{TAU:}{vector of \eqn{\tau}'s for each level}
#'
#' \item{ABCD:}{A matrix of 4 columns: \eqn{a_q}, \eqn{b_q}, \eqn{a_\rho},
#'     \eqn{{\sigma_\rho}}. Each row represents a dose level.}}
#'
#' @examples
#'
#' par.prior <- vtPriorPar(tau = c(0.2, 0.4, 0.6), sdalpha = 10);
#'
#' @export
#'
vtPriorPar <- function(prior.y = NULL, tau = NULL, sdalpha=10, sdrho=10, vtheta=NULL) {
    if (!is.null(prior.y)) {
        beta <- NULL;
        for (i in 1:nrow(prior.y)) {
            pt <- prior.y[i,];
            if (!is.null(vtheta))
                pt <- vtheta * pt/sum(pt);

            ## toxcity rate
            p.tox <- sum(pt[3:4])/sum(pt);

            ## response
            a.rep <- sum(pt[c(2,4)]);
            b.rep <- sum(pt[c(1,3)]);

            ## odds ratio
            a.lr  <- log(pt[1]*pt[4]/pt[2]/pt[3]);
            b.lr  <- sqrt(sum(1/pt));
            beta  <- rbind(beta, c(p.tox, a.rep, b.rep, a.lr, b.lr));
        }
        colnames(beta) <- c("TAU","A", "B","C","D");
        rst            <- list(TAU       = beta[,1],
                               ABCD      = beta[,2:5]);
    } else {
        if (is.null(tau))
            stop("Please provide prior.y or tau");

        abcd       <- array(0, dim=c(length(tau), 4));
        abcd[,1:2] <- 0.5;
        abcd[,4]   <- sdrho;
        rst        <- list(TAU  = tau,
                           ABCD = abcd);
    }

    rst$SDALPHA <- sdalpha;
    ##return
    class(rst) <- get.const()$CLSPRIOR;
    rst
}
#' Simulate a single trial
#'
#' Simulation function for simulating a single trial
#'
#' @rdname vtSingleTrial
#'
#' @inheritParams parameters
#'
#' @param trueps True \eqn{\theta}'s. A \code{VTTRUEPS} object made from
#'     \code{\link{vtScenario}}
#' @param size.cohort Size of each cohort
#' @param size.level Maximum number of patients for each dose level
#' @param ... Optional arguments for \code{vtPost}
#'
#'
#' @return
#' \itemize{
#' \item{\code{dose}: }{Optimal dose level}
#' \item{\code{n.patients}: }{Number of patients for each dose level and each cohort}
#' \item{\code{ptox}: }{Posterior mean of DLT risk rate after each interim analysis}
#' \item{\code{pres}: }{Posterior mean of immune response rate after each interim analysis}
#' \item{\code{region}: }{Identified region in the decision map after each interim analysis}
#' \item{\code{prob}: }{Posterior mean of \eqn{\theta}'s after each interim analysis}
#' \item{\code{smps}: }{Observed data after each cohort}
#' }
#'
#' @examples
#' rst.sce <- vtScenario(tox = c(0.05, 0.05, 0.08),
#'                       res = c(0.2, 0.3, 0.5),
#'                       rho = 1)
#' rst.simu  <- vtSingleTrial(trueps = rst.sce, size.cohort=3, size.level=12,
#'                            prob.mdl="NONPARA");
#'
#' @export
#'
vtSingleTrial <- function(trueps,
                          size.cohort = 3,
                          size.level  = NULL,
                          etas        = c(0.1,0.3),
                          dec.cut     = 0.65,
                          prob.mdl    = c("NONPARA", "NONPARA+", "PARA", "PARA+"),
                          priors      = NULL,
                          ...) {

    ## parameters
    ndose <- nrow(trueps);
    if (is.null(size.level))
        size.level <- rep(size.cohort * 2, ndose);

    if (length(size.level) < ndose) {
        size.level <- c(size.level,
                        rep(size.level[length(size.level)], ndose - length(size.level)));
    }

    prob.mdl <- match.arg(prob.mdl);
    nstages  <- ceiling(max(size.level)/size.cohort);

    ## prepare returns
    rst.np     <- array(0,  dim=c(ndose, nstages));    ##n patients treated on each level
    rst.smp    <- array(0,  dim=c(ndose, nstages, 4)); ##simulated data
    rst.ptox   <- array(NA, dim=c(ndose, nstages));    ##n dlt on each level
    rst.pres   <- array(NA, dim=c(ndose, nstages));    ##n responses on each level
    rst.prob   <- array(NA, dim=c(ndose, nstages, 4)); ##region prob
    rst.region <- array(NA, dim=c(ndose, nstages));    ##region selected

    ## simulation
    smps      <- NULL;
    action    <- 1; ##1 escalate; 0 stay;
    while (action >= 0) {
        if (1 == action) {
            smps      <- rbind(smps, rep(0,4));
            cur.dose  <- nrow(smps);
            cur.stage <- 1;

            if (1 == cur.dose) {
                prev.res <- 0;
            } else {
                prev.res  <- apply(post.smp[,c(2,4)],1,sum);
            }
        } else {
            cur.stage <- cur.stage + 1;
        }

        ##sample next cohort
        cur.smp <- rmultinom(1,
                             min(size.cohort, size.level[cur.dose] - sum(smps[cur.dose,])),
                             trueps[cur.dose,]);
        smps[cur.dose,] <- smps[cur.dose,] + cur.smp;

        ##bayesian
        post.smp <- vtPost(smps, prob.mdl, priors, ...);

        ##decision
        cur.dec  <- vtDecMap(post.smp, etas, prev.res=prev.res, dec.cut=dec.cut);

        ##keep records
        rst.np[cur.dose,     cur.stage]  <- sum(cur.smp);
        rst.smp[cur.dose,    cur.stage,] <- cur.smp;
        rst.ptox[cur.dose,   cur.stage]  <- cur.dec$ptox;
        rst.pres[cur.dose,   cur.stage]  <- cur.dec$pres;
        rst.prob[cur.dose,   cur.stage,] <- cur.dec$prob;
        rst.region[cur.dose, cur.stage]  <- cur.region <- cur.dec$region;

        ##decision
        action <- c(-1,-1,1,0)[cur.region];

        ##action
        if (0 == action &
            size.level[cur.dose] == sum(smps[cur.dose,])) {
            action <- 1;
        }

        if (ndose == cur.dose & 1 == action) {
            rst.dose <- ndose;
            action   <- -2;
        } else if (-1 == action) {
            rst.dose <- cur.dose - 1;
        }
    }

    ##return
    list(dose=rst.dose,
         n.patients=rst.np,
         ptox=rst.ptox,
         pres=rst.pres,
         region=rst.region,
         prob=rst.prob,
         smps=rst.smp);
}

#' Conduct simulation study
#'
#' Simulate clinical trials with given settings for multiple times to evaluate
#' the study operating characteristics.
#'
#' @rdname vtSimu
#'
#' @param n.rep Number of repetitions
#' @param seed Seed
#' @param ... Optional parameters for \code{\link{vtSingleTrial}}
#' @param n.cores Number of cores for parallel computations
#' @param update.progress Reserved parameter for Shiny GUI
#'
#' @return
#'
#' A class \code{VTSIMU} list with length \code{n.rep} of results. Each item is
#'     a list return from \code{\link{vtSingleTrial}}.
#'
#' @examples
#'
#' rst.sce <- vtScenario(tox = c(0.05, 0.05, 0.08),
#'                       res = c(0.2, 0.3, 0.5),
#'                       rho = 1)
#'
#' rst.simu  <- vtSimu(n.rep = 100, n.cors = 2, trueps = rst.sce,
#'                     size.cohort=3, size.level=12, prob.mdl="NONPARA");
#'
#'
#' @export
#'
vtSimu <- function(n.rep=100, seed=NULL, ..., n.cores=1, update.progress=NULL) {

    if (!is.null(seed)) {
        old_seed <- .Random.seed;
        set.seed(seed)
    }
    
    if ("PROGRESS" %in% toupper(class(update.progress)))
        update.progress$set(value=1, detail=paste(""));

    rst <- parallel::mclapply(1:n.rep,
                              function(x) {
                         if ("PROGRESS" %in% toupper(class(update.progress))) {
                             update.progress$set(value=x/n.rep,
                                                 detail=paste("Replication", x, sep=" "));
                         } else {
                             print(x);
                         }

                         vtSingleTrial(...);
                     }, mc.cores=n.cores);

    class(rst) <- get.const()$CLSSIMU;

    if (!is.null(seed)) {
        set.seed(old_seed);
    }
    
    rst
}


#' Summarize simulation results
#'
#' Summarize simulation results to get the frequency of a dose level is identified
#' as the optimal dose level and the number of DLT's and responses
#'
#' @rdname summary2.VTSIMU
#'
#' @param x A class \code{VTSIMU} list generated by \code{\link{vtSimu}}
#' @param ... Reserved parameters
#'
#'
#' @return A numeric array that shows 1: number of times each level is selected,
#'     2. total number of times any level is selected, 3. frequency each
#'     level is selected, 4. frequency any level is selected, 5. average number
#'     of DLT's and responders for each level, 6. average total number of DLT's
#'     and responders
#'
#' @examples
#'
#' rst.sce <- vtScenario(tox = c(0.05, 0.05, 0.08),
#'                       res = c(0.2, 0.3, 0.5),
#'                       rho = 1)
#' rst.simu  <- vtSimu(n.rep = 20, n.cors = 2, trueps = rst.sce,
#'                     size.cohort=3, size.level=12,
#'                     prob.mdl="NONPARA");
#' sum.simu <- summary2(rst.simu)
#'
#'
#' @method summary2 VTSIMU
#'
#' @export
#'
summary2.VTSIMU <- function(x, ...) {
    cur.rst  <- summary(x);
    rst      <- NULL;

    ##dose selection percentage
    rst <- c(rst, cur.rst$dose[2,]);
    rst <- c(rst, sum(cur.rst$dose[2,]));
    rst <- c(rst, cur.rst$npat[,"total"]);
    rst <- c(rst, sum(cur.rst$npat[,"total"]));

    ss <- 0;
    for (j in 1:ncol(cur.rst$dose)) {
        cur.ss   <- cur.rst$samples[paste("level ",j," total", sep = ""), c("T", "R")];
        ss       <- ss + cur.ss;
        rst      <- c(rst, cur.ss);
    }
    rst <- c(rst, ss);
    invisible(rst);
}


#' Plot true parameters
#'
#' Plot true DLT risk rates and response rates.
#'
#' @rdname plot.VTTRUEPS
#'
#' @param x A class \code{VTTRUEPS} matrix generated by \code{\link{vtScenario}}
#' @param draw.levels Select dose levels to draw. Default \code{NULL} draws all
#'     levels.
#' @param draw.curves Indicate which curves to plot. The options are
#'
#' \itemize{
#' \item{1:}{p: DLT risk rate}
#' \item{2:}{q: Response rate}
#' \item{3:}{\eqn{\theta_{00}}}
#' \item{4:}{\eqn{\theta_{01}}}
#' \item{5:}{\eqn{\theta_{10}}}
#' \item{6:}{\eqn{\theta_{11}}}
#' }
#'
#' See \code{\link{visit}} for details.
#'
#' @param legends Line legends
#' @param ltys Line types
#' @param pch Line PCH
#' @param ylim Y limits
#' @param cols Line colors
#' @param add.legend Include legends (TRUE) or not (FALSE)
#' @param ... optional arguments for plot
#'
#'
#'
#' @examples
#' rst.sce <- vtScenario(tox = c(0.05, 0.05, 0.08),
#'                       res = c(0.2, 0.3, 0.5),
#'                       rho = 1)
#' plot(rst.sce, draw.levels = 1:2, draw.curves=1:6)
#'
#' @method plot VTTRUEPS
#'
#' @export
#'
plot.VTTRUEPS <- function(x, draw.levels = NULL, draw.curves = 1:6,
                          legends = NULL, ltys = c(1,1,2,2,2,2), pch=19:24, ylim = c(0,1),
                          cols = c("red", "blue", "brown", "black", "gray", "green"),
                          add.legend = TRUE, ...) {
    ## parameter checking
    stopifnot(get.const()$CLSTRUEPS %in% class(x));
    draw.levels <- set.option(draw.levels, 1:nrow(x));
    draw.curves <- set.option(draw.curves, 1:6);

    if (is.null(legends)) {
        legends <- get.const()$DENLEGEND;
    }

    ## append trueps
    x <- cbind(x[,3] + x[,4],
               x[,2] + x[,4],
               x);

    ## plot
    plot(NULL, xlim=range(draw.levels), xlab="Levels", ylab="Rate", ylim = ylim, ...);
    for (k in draw.curves) {
        lines(draw.levels, x[draw.levels,k], type = "b",
              col = cols[k], lty = ltys[k], pch = pch[k]) ;
    }

    if (add.legend) {
        legend("topleft", bty="n",
               legend = legends[draw.curves],
               pch = pch[draw.curves],
               lty = ltys[draw.curves], col = cols[draw.curves]);
    }
}

#' Print true probabilities
#'
#' Print the true probabilities, with probabilities of toxicity and resistance,
#'      and \eqn{\rho}.
#'
#' @rdname summary.VTTRUEPS
#'
#' @param object A class \code{VTTRUEPS} matrix generated by \code{\link{vtScenario}}
#' @inheritParams parameters
#'
#'
#' @return
#'
#' A table showing the summary of the \code{VTTRUEPS} object. The first
#'     four columns are individual probability, fifth and sixth are probability
#'     for toxicity and resistance, and seventh is rho, the correlation.
#'
#' @examples
#' rst.sce <- vtScenario(tox = c(0.05, 0.05, 0.08),
#'                       res = c(0.2, 0.3, 0.5),
#'                       rho = 1)
#' summary(rst.sce)
#'
#' @method summary VTTRUEPS
#'
#' @export
#'
summary.VTTRUEPS <- function(object, digits = 2, ...) {

    stopifnot(get.const()$CLSTRUEPS %in% class(object));
    rst        <- round(object, digits = digits);

    rst2 <- cbind(object[,3]+object[,4],
                  object[,2]+object[,4]);

    colnames(rst2) <- c("Toxicity Rate", "Response Rate");

    rst3 <- object[,1]*object[,4]/object[,2]/object[,3];
    rst3 <- cbind(rst3);
    colnames(rst3) <- "Rho";

    cbind(rst, rst2, rst3);
}

#' Print true probabilities in latex format
#'
#' Print the true probabilities, with probabilities of toxicity and resistance,
#'      and \eqn{\rho}, in latex format
#'
#' @rdname summary2.VTTRUEPS
#'
#' @param x A class \code{VTTRUEPS} matrix generated by \code{\link{vtScenario}}
#' @inheritParams parameters
#' @param rp2d Columns to be in bold font
#'
#' @return
#'
#' A summary of the true probabilities in latex format.
#'
#' @examples
#'
#' rst.sce <- vtScenario(tox = c(0.05, 0.05, 0.08),
#'                       res = c(0.2, 0.3, 0.5),
#'                       rho = 1)
#' ltx.ps <- summary2(rst.sce)
#'
#'
#' @method summary2 VTTRUEPS
#'
#' @export
#'
summary2.VTTRUEPS <- function(x, rp2d = -1, digits = 2, ...) {
    fill.pq <- NULL;
    cur.pq  <- summary(x)[,5:6];
    for (j in 1:nrow(cur.pq)) {
        cur.pqnum <- round(cur.pq[j,], digits = digits);
        cur.pqtxt <- paste("(", cur.pqnum[1], ",", cur.pqnum[2], ")", sep="");
        if (j %in% rp2d) {
            cur.pqtxt <- paste("\\\\textbf{", cur.pqtxt,"}");
        }
        fill.pq <- c(fill.pq, cur.pqtxt);
    }

    invisible(fill.pq);
}

#' Summarize simulation results
#'
#' Summarize the simulation results with numerous statistical measures
#'
#' @rdname summary.VTSIMU
#'
#' @param object A class \code{VTSIMU} list generated by \code{\link{vtSimu}}
#' @inheritParams parameters
#'
#' @return
#'
#' A list containing
#' \itemize{
#' \item{dose}{: Frequency for each dose level being selected as the optimal dose level}
#' \item{npat}{: Average number of patients for each cohort and each dose level}
#' \item{samples}{: Average number of DLT risks and responses for each cohort on each dose level}
#' \item{decision}{: Frequency each region in the decision map is selected for each cohort on each dose level}
#'
#' \item{prob}{: Average conditional probabilities corresponding to each region in
#' the decision map for each cohort on each dose level}
#'
#' \item{ptox}{: Mean and credible interval of DLT risk rates for each cohort on each dose level}
#'
#' \item{pres}{: Mean and credible interval of immune response rates for each cohort on each dose level}
#' }
#'
#' @examples
#' rst.sce <- vtScenario(tox = c(0.05, 0.05, 0.08),
#'                       res = c(0.2, 0.3, 0.5),
#'                       rho = 1)
#' rst.simu  <- vtSimu(n.rep = 50, n.cors = 2, trueps = rst.sce,
#'                     size.cohort=3, size.level=12,
#'                     prob.mdl="NONPARA");
#' sum.simu <- summary(rst.simu)
#'
#' @method summary VTSIMU
#'
#' @export
#'
summary.VTSIMU <- function(object, ...) {
    f.tp <- function(var) {
        rst.p <- lst.combine(object, var);
        rst.p <- apply(rst.p, 1:2, function(x) {
            if (all(is.na(x))) {
                rst <- rep(NA, 3);
            } else {
                rst <- c(mean(x, na.rm = TRUE),
                         quantile(x, c(0.025,0.975), na.rm = TRUE))
            };
            rst
        });

        rst  <- NULL;
        for (i in 1:n.dose) {
            for (j in 1:n.stage) {
                rst  <- rbind(rst, rst.p[,i,j]);
            }
        }
        colnames(rst) <- c("Mean", "LB", "UB");
        rst
    }

    chk.pts <- object[[1]]$n.patients;
    n.dose  <- nrow(chk.pts);
    n.stage <- ncol(chk.pts);
    n.reps  <- length(object);

    ##selected dose
    rst.dose <- lst.combine(object, "dose");
    rst.dose <- table(factor(rst.dose, levels=1:n.dose));
    rst.dose <- rbind(rst.dose, rst.dose/n.reps*100);
    rownames(rst.dose) <- c("Frequency", "%");

    ##enrolled patients
    rst.np   <- lst.combine(object, "n.patients");
    rst.np   <- apply(rst.np, 1:2, mean);
    rst.np   <- cbind(rst.np, apply(rst.np, 1, sum));
    rownames(rst.np) <- paste("level", 1:n.dose,  sep = " ");
    colnames(rst.np) <- c(paste("stage", 1:n.stage, sep = " "), "total");

    ##patients details
    rst.smps     <- lst.combine(object, "smps");
    rst.smps     <- apply(rst.smps, 1:3, mean);
    s.smps       <- NULL;
    s.smps.names <- NULL;
    for (i in 1:n.dose) {
        for (j in 1:n.stage) {
            s.smps.names <- c(s.smps.names,
                              paste("level", i, "stage", j, sep = " "));

            cur.smps <- rst.smps[i,j,];
            s.smps   <- rbind(s.smps,
                              c(cur.smps,
                                cur.smps[3] + cur.smps[4],
                                cur.smps[2] + cur.smps[4]));
        }
        s.smps.names <- c(s.smps.names,
                          paste("level", i, "total", sep = " "));
        cur.level    <- apply(rst.smps[i,,], 2, sum);
        s.smps       <- rbind(s.smps,
                              c(cur.level,
                                cur.level[3] + cur.level[4],
                                cur.level[2] + cur.level[4]));
    }

    rownames(s.smps) <- s.smps.names;
    colnames(s.smps) <- c("No T, No R", "No T, R", "T, No R", "T,R", "T", "R");

    ##probabilities
    rst.prob <- lst.combine(object, "prob");
    rst.prob <- apply(rst.prob, 1:3, mean, na.rm=TRUE);
    s.probs       <- NULL;
    s.probs.names <- NULL;
    for (i in 1:n.dose) {
        for (j in 1:n.stage) {
            s.probs.names <- c(s.probs.names,
                              paste("level", i, "stage", j, sep = " "));
            s.probs   <- rbind(s.probs, rst.prob[i,j,]);
        }
    }
    rownames(s.probs) <- s.probs.names;
    colnames(s.probs) <- get.const()$REGIONS;

    ##decision
    rst.decision <- lst.combine(object, "region");
    rst.decision <- apply(rst.decision, 1:2, function(x) {
        if (all(is.na(x))) {
            rst <- rep(NA, 4);
        } else {
            rst <- table(factor(x, levels=1:4))/sum(!is.na(x)) * 100;
        };
        c(sum(!is.na(x)), rst);
    });

    s.dec       <- NULL;
    s.dec.names <- NULL;
    for (i in 1:n.dose) {
        for (j in 1:n.stage) {
            s.dec.names <- c(s.dec.names,
                              paste("level", i, "stage", j, sep = " "));

            s.dec  <- rbind(s.dec, rst.decision[,i,j]);
        }
    }

    rownames(s.dec) <- s.dec.names;
    colnames(s.dec) <- c("N", get.const()$REGIONS);

    ##p.tox
    p.tox <- f.tp("ptox");
    p.res <- f.tp("pres");
    rownames(p.tox) <- s.dec.names;
    rownames(p.res) <- s.dec.names;

    ##return
    rst <- list(dose    = rst.dose,
                npat    = rst.np,
                samples = s.smps,
                decison = s.dec,
                prob    = s.probs,
                ptox    = p.tox,
                pres    = p.res);
}
