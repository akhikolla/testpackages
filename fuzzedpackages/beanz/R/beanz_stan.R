##----------------------------------------------------------------
##                  STAN MODEL
##----------------------------------------------------------------
#' Call STAN models
#'
#' Call STAN to draw posterior samples for Bayesian HTE models.
#'
#' @param mdls name of the Bayesian HTE model. The options are:
#'
#' \describe{
#'   \item{nse}{No subgroup effect model}
#'   \item{fs}{Full stratification model}
#'   \item{sr}{Simple regression model}
#'   \item{bs}{Basic shrinkage model}
#'   \item{srs}{Simple regression with shrinkage model}
#'   \item{ds}{Dixon-Simon model}
#'   \item{eds}{Extended Dixon-Simon model}
#' }
#'
#' @param dat.sub dataset with subgroup treatment effect summary data
#'
#' @param var.estvar column names in dat.sub that corresponds to treatment effect
#'     estimation and the estimated variance
#'
#' @param var.cov array of column names in dat.sub that corresponds to binary or
#'     ordinal baseline covariates
#'
#' @param var.nom array of column names in dat.sub that corresponds to nominal
#'     baseline covariates
#'
#' @param delta parameter for specifying the informative priors of \eqn{\sigma_g}
#'
#' @param prior.sig option for the informative prior on \eqn{\sigma_g}. 0: uniform prior and
#'      1: log-normal prior
#'
#' @param par.pri vector of prior parameters for each model. See
#'     \code{\link{beanz-package}} for the details of model specification.
#'
#' \describe{
#'   \item{nse, fs}{\code{B}}
#'   \item{sr}{\code{B}, \code{C}}
#'   \item{bs, ds, eds}{\code{B}, \code{D}}
#'   \item{srs}{\code{B}, \code{C}, \code{D}}
#' }
#'
#' @param chains STAN options. Number of chains.
#'
#' @param ... options to call STAN sampling. These options include
#'     \code{iter}, \code{warmup}, \code{thin}, \code{algorithm}.
#'     See \code{rstan::sampling} for details.
#'
#' @return A class \code{beanz.stan} list containing
#'  \describe{
#'    \item{mdl}{name of the Bayesian HTE model}
#'    \item{stan.rst}{raw \code{rstan} \code{sampling} results}
#'    \item{smps}{matrix of the posterior samples}
#'    \item{get.mus}{method to return the posterior sample of the subgroup treatment effects}
#'    \item{DIC}{DIC value}
#'    \item{looic}{leave-one-out cross-validation information criterion}
#'    \item{rhat}{Gelman and Rubin potential scale reduction statistic}
#'    \item{prior.sig}{option for the informative prior on \eqn{\sigma_g}}
#'    \item{delta}{parameter for specifying the informative priors of \eqn{\sigma_g}}
#' }
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
#'                                   var.resp   = var.resp,
#'                                   var.trt    = var.trt,
#'                                   var.cov    = var.cov,
#'                                   var.censor = var.censor,
#'                                   resptype   = resptype);
#'
#' rst.nse    <- bzCallStan("nse", dat.sub=subgrp.effect,
#'                          var.estvar = var.estvar, var.cov = var.cov,
#'                          par.pri = c(B=1000),
#'                          chains=4, iter=600,
#'                          warmup=200, thin=2, seed=1000);
#'
#' rst.sr     <- bzCallStan("sr", dat.sub=subgrp.effect,
#'                         var.estvar=var.estvar, var.cov = var.cov,
#'                         par.pri=c(B=1000, C=1000),
#'                         chains=4, iter=600,
#'                         warmup=200, thin=2, seed=1000);}
#' @export
#'
#'
#'
bzCallStan <- function(mdls = c("nse", "fs", "sr", "bs", "srs", "ds", "eds"),
                      dat.sub,
                      var.estvar,
                      var.cov,
                      par.pri=c(B=1000.0,C=1000.0,D=1.0),
                      var.nom=NULL,
                      delta = 0.0,
                      prior.sig=1,
                      chains=4,
                      ...) {

    mdls <- match.arg(mdls);

    ##check data
    if (!all(c(var.estvar, var.cov) %in% colnames(dat.sub)))
        stop("Variables specified are not in the dataset.");

    ##check par.pri
    if (!all(names(par.pri) %in% c("B","C","D")))
        stop("Prior parameters are not recognized.");

    ##check number of chains
    if (chains < 2)
        stop("At least 2 Markov Chains are required in order to check convergence.")

    stopifnot(all(par.pri >= 0));

    ##check delta
    stopifnot(delta >= 0);

    ##check priorsig
    stopifnot(prior.sig %in% c(0,1));

    ##basic data
    lst.basic <- list(SIZE     = nrow(dat.sub),
                      Y        = dat.sub[,var.estvar[1]],
                      SIGY     = sqrt(dat.sub[,var.estvar[2]]),
                      DELTA    = delta,
                      PRIORSIG = prior.sig);

    ##model specific
    vlist <- NULL;
    eval(parse(text=paste("vlist <- stan.model.", mdls,
               "(dat.sub, var.estvar, var.cov, var.nom)", sep="")));

    ## prior data
    names(par.pri) <- toupper(names(par.pri));

    ##call stan
    stan.rst <- rstan::sampling(stanmodels[[mdls]],
                                data=c(lst.basic, vlist, par.pri),
                                pars=c("uvs", "nvs", "nomega", "nphi"),
                                include=FALSE,
                                chains=chains,
                                ...);

    smps <- rstan::extract(stan.rst, permuted=FALSE);

    ##diagnositics-DIC
    dic   <- get.dic(dat.sub[,var.estvar[1]], smps);
    ##diagnositics-LOOIC
    looic <- loo::loo(loo::extract_log_lik(stan.rst), cores = 1);
    ##rhat
    rhat  <- rstan::summary(stan.rst)$summary[,"Rhat"];

    ##return
    rst <- list(mdl       = get.mdl.name(mdls),
                stan.rst  = stan.rst,
                smps      = smps,
                dic       = dic,
                looic     = looic,
                prior.sig = prior.sig,
                delta     = delta,
                rhat      = rhat,
                get.mus=function() {
        n.mu <- stan.rst@par_dims$mu;
        mj    <- paste("mu[", 1:n.mu, "]", sep="");
        s.mus <- smps[,,mj];

        mus <- NULL;
        for (k in 1:dim(s.mus)[2]) {
            mus <- rbind(mus, s.mus[,k,]);
        }
        colnames(mus) <- paste("Subgroup", 1:n.mu);
        mus
    })

    class(rst) <- "beanz.stan";
    rst
}

##no subgroup effect model
stan.model.nse <- function(dat.sub, var.estvar, var.cov, var.nom) {
    NULL;
}
##full stratification
stan.model.fs <- stan.model.nse;
##basic shrinkage
stan.model.bs <- stan.model.fs;

##simple regression
stan.model.sr <- function(dat.sub, var.estvar, var.cov, var.nom) {
    ##design matrix
    dx     <- dat.sub[,var.cov, drop=FALSE];
    dx     <- df.convert(dx, var.nom);
    fml    <- paste("~", paste(var.cov, collapse="+"), sep="");
    des.x  <- model.matrix(formula(fml), dx);
    X      <- des.x[,-1]; #remove intercept
    NX     <- ncol(X);

    vname  <- make.list(environment());
    vname
}

##simple regression + shrinkage
stan.model.srs <- stan.model.sr;
##Dixon and Simon
stan.model.ds <- stan.model.sr;

##extended dixon and simon
stan.model.eds <- function(dat.sub, var.estvar, var.cov, var.nom) {
    ##design matrix
    dx     <- dat.sub[,var.cov, drop=FALSE];
    dx     <- df.convert(dx, var.nom);
    NTAU   <- ncol(dx);

    fml    <- paste("~", paste(var.cov, collapse="+"), "^", NTAU, sep="");
    des.x  <- model.matrix(formula(fml), dx);
    X      <- des.x[,-1]; #remove intercept
    NX     <- ncol(X);

    ##use the number of : for the order of interaction
    d.name <- colnames(X);
    TAUINX <- 1 + nchar(d.name) - nchar(gsub(":", "", d.name));

    ##return
    vname <- make.list(environment());
    vname
}

##convert dataframe into numeric values
##except the var.nom (nominal covariates)
df.convert <- function(dat, var.nom=NULL) {
    rst <- as.data.frame(lapply(dat, as.numeric));
    if (!is.null(var.nom)) {
        for (i in 1:length(var.nom)) {
            rst[,var.nom[i]] <- as.character(rst[,var.nom[i]]);
        }
    }
    rst
}

make.list <- function(par.env) {
    objlst  <- objects(par.env);
    rst     <- list();
    for (o in 1:length(objlst)) {
        curo <- objlst[o];
        if (toupper(curo) == curo) {
            rst[[curo]] <- get(curo, envir=par.env);
        }
    }
    rst
}

##-----------------------------------------------------------------------------------
##                             DIC
##-----------------------------------------------------------------------------------
get.log.ll <- function(y, mu, vs) {
    rst <- dnorm(y, mu, vs, log=TRUE);
    sum(rst);
}

get.dic <- function(y, stan.smps) {
    ny    <- length(y);
    l.mus <- paste("mu[", 1:ny, "]", sep="");
    l.vs  <- paste("vs[", 1:ny, "]", sep="");
    l.ll  <- paste("log_lik[", 1:ny, "]", sep="");

    mus     <- stan.smps[, , l.mus];
    vs      <- stan.smps[, , l.vs];
    ll      <- stan.smps[, , l.ll];
    mean.mu <- apply(mus, 3, mean);
    mean.vs <- 1/(apply(1/vs, 3, mean));
    d.theta.bar <- -2*get.log.ll(y, mean.mu, mean.vs);
    all.d       <- -2*apply(ll, c(1,2), sum);
    bar.d       <- mean(all.d);

    ##dic
    DIC <- 2*bar.d - d.theta.bar;
    DIC
}

get.mdl.name <- function(mdl) {

    ALL.MODELS  <- c("No subgroup effect",
                     "Full stratification",
                     "Simple regression",
                     "Basic shrinkage",
                     "Regression and shrinkage",
                     "Dixon and Simon",
                     "Extended Dixon and Simon"
                     );
    STAN.NAME   <- c("nse", "fs", "sr", "bs", "srs", "ds", "eds");

    rst <- ALL.MODELS[which(mdl == STAN.NAME)];
}

