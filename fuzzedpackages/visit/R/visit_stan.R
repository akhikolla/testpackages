#' Call STAN models for MCMC sampling
#'
#' Call STAN to draw posterior samples of the joint distribution of
#' immunogenicity rate and toxicity risk for parametric and parametric+ model
#'
#' @rdname vtStan
#'
#' @inheritParams parameters
#'
#' @param model option of the probability models:
#' \describe{
#'   \item{0:}{parametric model}
#'   \item{1:}{parametric+ model}
#' }
#' See \code{\link{visit}} for details.
#'
#' @param iter STAN option: number of iterations
#' @param chains STAN option: number of chains
#' @param warmup STAN option: number of warmup
#' @param ... additional parameters for package rstan's sampling method. These
#'     options include \code{iter}, \code{warmup}, \code{thin},
#'     \code{algorithm}. See \code{rstan::sampling} for details.
#'
#'
#' @return A \code{rstan} object that contains the posterior sampling results
#'
#'
vtStan <- function(obs.y,
                   priors,
                   model  = 0,
                   iter   = 4000,
                   chains = 4,
                   warmup = 2000,
                   ...) {

    ##model
    if (0 == model) {
        ##PARA
        FIXRHOAT1   <- 0;
        SINGLERHO   <- 1;
        SINGLEALPHA <- 1;
    } else if (1 == model) {
        ##PARA+
        FIXRHOAT1   <- 1;
        SINGLERHO   <- 1;
        SINGLEALPHA <- 1;
    } else if (2 == model) {
        FIXRHOAT1   <- 0;
        SINGLERHO   <- 1;
        SINGLEALPHA <- 0;
    } else if (3 == model) {
        FIXRHOAT1   <- 1;
        SINGLERHO   <- 1;
        SINGLEALPHA <- 0;
    }

    ##data
    NDOSE     <- length(priors$TAU);
    SDALPHA   <- priors$SDALPHA;
    TAU       <- priors$TAU;
    PAR       <- priors$ABCD;

    OBSY <- array(0, dim=c(1,4)); ##pseudo record for stan
    if (is.null(obs.y)) {
        NOBSLEVEL <- 0;
    } else {
        NOBSLEVEL <- nrow(obs.y);
        OBSY      <- rbind(obs.y, OBSY);
    }

    ##call stan
    stan.rst <- rstan::sampling(stanmodels[["visit"]],
                                data=list(NDOSE=NDOSE,
                                          TAU=as.array(TAU),
                                          PAR=PAR,
                                          SDALPHA=SDALPHA,
                                          NOBSLEVEL=NOBSLEVEL,
                                          OBSY=OBSY,
                                          SINGLEALPHA=SINGLEALPHA,
                                          SINGLERHO=SINGLERHO,
                                          FIXRHOAT1=FIXRHOAT1),
                                iter=iter,
                                chains=chains,
                                warmup=warmup,
                                ...);
    ##return
    stan.rst
}


#' Postetrior sampling for given observed samples
#'
#' Call STAN to draw posterior samples of the joint distribution of
#' immunogenicity rate and toxicity risk
#'
#' @rdname vtPost
#'
#' @inheritParams parameters
#'
#'
#' @param nsmp number of iterations
#'
#' @param ... additional parameters for package rstan's sampling method. These
#'     options include \code{warmup}, \code{thin}, \code{algorithm}. See
#'     \code{rstan::sampling} for details.
#'
#' @param prior.const Specify \eqn{\alpha} for a Beta(\eqn{\alpha},
#'     \eqn{\alpha}) prior. The Beta prior is used for \code{NONPARA} and
#'     \code{NONPARA+} models. Default value \eqn{0.5}.
#'
#' @return A class \code{VTPOST} matrix of posterior samples with \code{nsmp}
#'     rows and 4 columns. Columns 1-4 correspond to\eqn{\theta^{(l)}_{00},
#'     \theta^{(l)}_{01}, \theta^{(l)}_{10}, \theta^{(l)}_{11}}. See
#'     \code{\link{visit}} for details about \eqn{\theta}'s.
#'
#'
#'
#' @examples
#' obs.y    <- rbind(c(5, 2, 0, 0), c(3, 4, 0, 0), c(1, 6, 0, 0))
#' prior <- vtPriorPar(prior.y = NULL, tau = c(0.1, 0.3, 0.6), sdalpha=10, sdrho=10, vtheta=NULL)
#' rst.post <- vtPost(obs.y, priors = prior, warmup = 100, prob.mdl = "PARA", nsmp = 200)
#'
#' @export
#'
vtPost <- function(obs.y,
                   prob.mdl = c("NONPARA", "NONPARA+", "PARA", "PARA+"),
                   priors = NULL,
                   ...,
                   nsmp = 4000,
                   prior.const = 0.5) {

    stopifnot(4 == ncol(obs.y));

    ## parameter
    prob.mdl <- match.arg(prob.mdl);
    if ("PARA"  == prob.mdl |
        "PARA+" == prob.mdl) {
        if (is.null(priors))
            stop("Priors need to be specified for parametric models.")
    }
    cur.level <- nrow(obs.y);

    ##posterior sampling
    if ("NONPARA" == prob.mdl) {
        post.smp <- rdirichlet(nsmp, obs.y[cur.level,] + prior.const);
    } else if ("NONPARA+" == prob.mdl) {
        post.dlt <- rbeta(nsmp,
                          sum(obs.y[cur.level,c(3,4)])+prior.const,
                          sum(obs.y[cur.level,c(1,2)])+prior.const);
        post.rep <- rbeta(nsmp,
                          sum(obs.y[cur.level,c(2,4)])+prior.const,
                          sum(obs.y[cur.level,c(1,3)])+prior.const);

        post.smp <- cbind((1-post.dlt)*(1-post.rep), (1-post.dlt)*post.rep,
                          post.dlt*(1-post.rep), post.dlt*post.rep);

    } else {
        if ("PARA" == prob.mdl) {
            post.rst <- vtStan(obs.y, priors, model=0, iter = nsmp, ...);
        } else if ("PARA+" == prob.mdl){
            post.rst <- vtStan(obs.y, priors, model=1, iter = nsmp, ...);
        }
        post.theta <- extract(post.rst)$theta;
        post.smp   <- post.theta[,cur.level,];
    }

    ##return
    class(post.smp) <- get.const()$CLSPOST;

    post.smp
}
