
#' @useDynLib BayesLN, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom Matrix t
#' @importFrom methods as
NULL


#'Numerical evaluation of the log-normal conditioned means posterior moments
#'
#'Function that evaluates the existence conditions for moments of useful quantities in the original data scale
#'when a log-normal linear mixed model is estimated.
#'
#'
#'@param X Design matrix for fixed effects.
#'@param Z Design matrix for random effects.
#'@param Xtilde Covariate patterns used for the leverage computation.
#'@param order_moment Order of the posterior moments required to be finite.
#'@param s Number of variances of the random effects.
#'@param m Vector of size \code{s} (if s>1) that indicates the dimensions of the random effect vectors.
#'
#'
#'@details
#'
#' This function computes the existence conditions for the moments up to order fixed by \code{order_moment} of the log-normal
#' linear mixed model specified by the design matrices \code{X} and \code{Z}. It considers the prediction based on multiple
#' covariate patterns stored in the rows of the \code{Xtilde} matrix.
#'
#'@return
#' Both the values of the factors determining the existence condition and the values of the gamma parameters for the different
#' variance components are provided.
#'
#'
#'
#' @export


LN_hier_existence <-
  function(X,
           Z,
           Xtilde,
           order_moment = 2,
           s = 1,
           m = NULL) {
    lev <- numeric(dim(Xtilde)[1])
    invXtX<-solve(t(X) %*% X)
    invZtZ<-MASS::ginv(t(Z) %*% Z)
    for (j in 1:dim(Xtilde)[1]) {
      lev[j] <- Xtilde[j, ] %*% invXtX %*% t(Xtilde)[, j]
    }
    pred_cond_cond <- max(lev)
    if (s == 1) {
      pred_cond <- 0
      mat_inv <- solve(
        1e6 * t(X) %*% (diag(rep(1, dim(
          Z
        )[1])) - Z %*% invZtZ %*% t(Z)) %*% X +
          (t(X) %*% (
            Z %*% invZtZ %*% invZtZ %*% t(Z)
          ) %*% X)
      )
      for (j in 1:dim(Xtilde)[1]) {
        pred_cond <-
          max(c(
            pred_cond,
            Xtilde[j, ] %*% mat_inv %*% t(Xtilde)[,j]
          ))
      }
    } else{
      pred_cond <- numeric(s)
      m_cum <- c(0, cumsum(m))
      for (i in 1:s) {
        diag_vec <- rep(0, m_cum[s + 1])
        diag_vec[(m_cum[i] + 1):(m_cum[i + 1])] <- 1
        D <- diag(diag_vec)
        mat_inv <- solve(
          1e6 * t(X) %*% (diag(rep(1, dim(
            Z
          )[1])) - Z %*% invZtZ %*% t(Z)) %*% X +
            (
              t(X) %*% (
                Z %*% invZtZ %*% D %*% invZtZ %*% t(Z)
              ) %*% X
            )
        )
        for (j in 1:dim(Xtilde)[1]) {
          pred_cond[i] <-
            max(c(
              pred_cond[i],
              Xtilde[j, ] %*% mat_inv %*% t(Xtilde)[, j]
            ))
        }
      }
    }

    gamma_cond_re <-
      sqrt(order_moment + 1 + (order_moment + 1) ^ 2 * pred_cond_cond)
    gamma_red <-
      sqrt((order_moment + 1) ^ 2 + (order_moment + 1) ^ 2 * pred_cond_cond)
    gamma_cond_mean <- sqrt(order_moment + 1 + (order_moment + 1) ^ 2 * pred_cond)

    gamma <- c(gamma_cond_re, gamma_red, gamma_cond_mean)
    factor_cond <- c(pred_cond_cond, pred_cond_cond, pred_cond)
    out <- cbind(factor_cond, gamma)
    colnames(out) <- c("Factor condition", "Gamma parameter")
    rownames(out) <- c("Sigma2", "Sigma2 - predictive", rep("Tau2", s))
    return(out)

  }



#' Bayesian estimation of a log - normal hierarchical model
#'
#'Function that estimates a log-normal linear mixed model with GIG priors on the variance components,
#'in order to assure the existence of the posterior moments of key functionals in the original data scale like conditioned means
#'or the posterior predictive distribution.
#'
#'@param formula_lme A two-sided linear formula object describing
#'both the fixed-effects and random-effects part of the model is required. For details see \code{\link{lmer}}.
#'@param data_lme Optional data frame containing the variables named in \code{formula_lme}.
#'@param y_transf Logical. If \code{TRUE}, the response variable is assumed already as log-transformed.
#'@param functional Functionals of interest: \code{"Subject"} for subject-specific conditional mean,
#'     \code{"Marginal"} for the overall expectation and \code{"PostPredictive"} for the posterior predictive distribution.
#'@param data_pred Data frame with the covariate patterns of interest for prediction. All the covariates present in the \code{data_lme} object must be included. If \code{NULL} the design matrix of the model is used.
#'@param order_moment Order of the posterior moments that are required to be finite.
#'@param nsamp Number of Monte Carlo iterations.
#'@param par_tau List of vectors defining the triplets of hyperparaemters for each random effect variance (as many vectors as the number of specified random effects variances).
#'@param par_sigma Vector containing the tiplet of hyperparameters for the prior of the data variance.
#'@param var_pri_beta Prior variance for the model coefficients.
#'@param inits List of object for initializing the chains. Objects with compatible dimensions must be named with \code{beta}, \code{sigma2} and \code{tau2}.
#'@param verbose Logical. If \code{FALSE}, the messages from the Gibbs sampler are not shown.
#'@param burnin Number of iterations to consider as burn-in.
#'@param n_thin Number of thinning observations.
#'
#'@details
#'The function allows to estimate a log-normal linear mixed model through a Gibbs sampler. The model equation is specified as in \code{\link{lmer}} model and the target functionals to estimate need to be declared.
#'A weakly informative prior setting is automatically assumed, always keeping the finiteness of the posterior moments of the target functionals.
#'
#'
#'
#'@return
#'The output list provided is composed of three parts. The object \code{$par_prior} contains the parameters fixed for the variance components priors. The object \code{$samples} contains the posterior samples for all the paramters.
#'They are returned as a \code{\link{mcmc}} object and they can be analysed trough the functions contained in the
#'\code{coda} package in order to check for the convergence of the algorithm. Finally, in \code{$summaries} an overview of the posteriors of the model parameters and of the target functionals is provided.
#'
#'
#'@examples
#' \donttest{
#' library(BayesLN)
#' # Load the dataset included in the package
#' data("laminators")
#' data_pred_new <- data.frame(Worker = unique(laminators$Worker))
#' Mod_est<-LN_hierarchical(formula_lme = log_Y~(1|Worker),
#'                          data_lme = laminators,
#'                          data_pred = data_pred_new,
#'                          functional = c("Subject","Marginal"),
#'                          order_moment = 2, nsamp = 50000, burnin = 10000)
#'}
#'
#' @export



LN_hierarchical <- function(formula_lme,
                            data_lme,
                            y_transf = TRUE,
                            functional = c("Subject", "Marginal", "PostPredictive"),
                            data_pred = NULL,
                            order_moment = 2,
                            nsamp = 10000,
                            par_tau = NULL,
                            par_sigma = NULL,
                            var_pri_beta = 1e3,
                            inits = list(NULL),
                            verbose = TRUE,
                            burnin = 0.1 * nsamp,
                            n_thin = 1) {
  if(!is.data.frame(data_lme)){
    stop("'data_lme' must be a data.frame")
  }
  if(!is.data.frame(data_lme)){
    stop("'data_lme' must be a data.frame")
  }
  if(!is.null(data_pred)){
    if(!is.data.frame(data_pred)){
      stop("'data_pred' must be a data.frame")
    }
  n<-nrow(data_lme)
  if(!as.character(formula_lme[[2]])%in%colnames(data_pred)){
    data_pred<-cbind(1,data_pred)
    colnames(data_pred)[1]<-as.character(formula_lme[[2]])
  }

  if(ncol(data_pred)!=ncol(data_lme)){
    stop("'data_pred' must have the same number of columns of 'data_lme'")
  }
  if(!all(colnames(data_pred)%in%colnames(data_lme))){
    stop("Covariates in 'data_pred' must have the same names than 'data_lme'")
  }
  data_lme<-data.table::rbindlist(list(data_lme,data_pred), use.names = T)
  }
  lmer_fitted <-
    suppressMessages(suppressWarnings(lme4::lmer(formula_lme, data = data_lme)))
  X <- as.matrix(lme4::getME(lmer_fitted, "X")[, ])
  X <- array(X, dim = dim(X),dimnames = list(NULL, colnames(X)))
  y <- as.numeric(lme4::getME(lmer_fitted, "y"))
  Zlist <- lme4::getME(lmer_fitted, "Ztlist")
  s <- length(Zlist)
  Z <- NULL
  m <- numeric(s)
  for (i in 1:s) {
    Zlist[[i]]<-Matrix::t(Zlist[[i]])
    Zprovv<- as.matrix(Zlist[[i]])
    m[i] <- dim(Zprovv)[2]
    colnames(Zprovv)<-paste(names(Zlist)[i],colnames(Zprovv), sep="_")
    Z <- cbind(Z, Zprovv)
  }
  Zprovv <- NULL
  if(is.null(data_pred)){
    Xtilde <- X
    Ztilde <- Z
    Ztilde_list <- Zlist
  }else{
    Xtilde <- matrix(X[-c(1:n),], ncol = ncol(X), dimnames = list(NULL, colnames(X)))
    X <- matrix(X[1:n,], ncol = ncol(X), dimnames = list(NULL, colnames(X)))
    Ztilde <- matrix(Z[-c(1:n),], ncol = ncol(Z), dimnames = list(NULL, colnames(Z)))
    Z <- matrix(Z[1:n,], ncol = ncol(Z), dimnames = list(NULL, colnames(Z)))
    y <- y[1:n]
    Ztilde_list<-list(NULL)
    for (i in 1:s) {
      Zlist[[i]]<-Zlist[[i]][1:n,]
      Ztilde_list[[i]]<-Zlist[[i]][-c(1:n),]
    }
  }

  if (y_transf == FALSE) {
    y_log <- log(y)
  } else {
    y_log <- y
  }

  ###parameter
  if (is.null(par_tau) && is.null(par_sigma)) {
    cond <- LN_hier_existence(
      X = X,
      Z = Z,
      Xtilde = Xtilde,
      order_moment = order_moment,
      s = s,
      m = m
    )
    gammas <-
      list(
        "Subject" = cond[1, 2],
        "Marginal" = max(cond[-2, 2]),
        "PostPredictive" = cond[2, 2]
      )
    gamma_pri <-
      max(unlist(gammas[names(gammas)[names(gammas) %in% functional]]))
    par_tau <- list(NULL)
    for (i in 1:s) {
      par_tau[[i]] <- c(1, 0.01, gamma_pri)
    }
    par_sigma <- c(1, 0.01, gamma_pri)
  } else{
    for (i in 1:s) {
      if (length(par_tau[[i]]) != 3) {
        stop("Each element of par_tau must be a vector with dimension 3")
      }
    }
    if (length(par_sigma) != 3) {
      stop("par_tau must be a vector with dimension 3")
    }
  }
    l_s <- par_sigma[1]
    d_s <- par_sigma[2]
    g_s <- par_sigma[3]
    par_tau_mat <- matrix(unlist(par_tau), ncol = 3, byrow = T)
    l_t <- par_tau_mat[, 1]
    d_t <- par_tau_mat[, 2]
    g_t <- par_tau_mat[, 3]
    m_s_cum <- cumsum(m)
    S_beta_pri<-var_pri_beta * diag(1, ncol(X))

    Precmat_Eff<-list(NULL)

    for(i in 1:s){
      Precmat_Eff[[i]] <- as(diag(1,nrow = m[i]), "dgCMatrix")
    }

    if (is.null(inits[["beta"]])) {
      beta_init <- solve(t(X) %*% X) %*% t(X) %*% y_log
    } else{
      beta_init <- inits[["beta"]]
    }
    if (is.null(inits[["sigma2"]])) {
      sigma2_init <- 1
    } else{
      sigma2_init <- inits[["sigma2"]]
    }
    if (is.null(inits[["tau2"]])) {
      tau2_init <- rep(1, s)
    } else{
      tau2_init <- inits[["tau2"]]
    }

      output <-   .Call(`_BayesLN_MCMC_alg`, y_log, X, Zlist, Precmat_Eff,
                       S_beta_pri, l_s, l_t, d_s, d_t, g_s, g_t, s,
                       nsamp, as.numeric(verbose), beta_init, sigma2_init, tau2_init)

      if ("PostPredictive" %in% functional){
        output["yrep"] <- .Call(`_BayesLN_post_pred`, output, Xtilde, Ztilde_list, s, nsamp)
        output$yrep <- exp(output$yrep)
      }


  if (is.null(colnames(X))) {
    colnames(output$beta) <- paste("beta", 0:(dim(X)[2] - 1))
  } else{
    colnames(output$beta) <- colnames(X)
  }
    u_unlist<-output$u[[1]]
    if(s>1){
    for(i in 2:s){
      u_unlist<-cbind(u_unlist,output$u[[i]])
    }
    }
    output$u<-u_unlist
  colnames(output$u) <- colnames(Z)

  out.mat <- cbind(output$tau2, output$sigma2, output$beta, output$u)
  colnames(out.mat)[1:(s + 1)] <- c(paste("tau2", names(Zlist), sep="_"), "sigma2")
  summ1 <-  summary(coda::mcmc(out.mat[seq(burnin + 1, nsamp, by = n_thin), ], start =
                         burnin + 1, thin = n_thin))
  par_summ <- round(cbind(summ1[[1]][, -4], summ1[[2]], "N_eff" = coda::effectiveSize(
                    coda::mcmc(out.mat[seq(burnin + 1, nsamp, by = n_thin), ], start = burnin +
                    1, thin = n_thin))), 3)
  samples <- list(par = coda::mcmc(out.mat[seq(burnin + 1, nsamp, by = n_thin), ], start =
                  burnin + 1, thin = n_thin))
  iter <- paste(summ1$start, ":", nsamp, sep = "")
  ssize <- length(seq(burnin + 1, nsamp, by = n_thin))
  summaries <-
    list(
      Iterations = iter,
      Thinning = n_thin,
      Sample_size = ssize,
      par = par_summ
    )

  if ("Subject" %in% functional) {
    Subj <- matrix(nrow = nsamp, ncol = nrow(Xtilde))
    colnames(Subj) <- rep("null", nrow(Xtilde))
    for (i in 1:nrow(Xtilde)) {
      Subj[, i] <- exp(c(t(matrix(Xtilde[i, ])) %*% t(output$beta)) +
                         c(Ztilde[i, ] %*% t(output$u)) + 0.5 * output$sigma2)
      colnames(Subj)[i] <- paste("Subj", i)
    }
    summ2 <-
      summary(coda::mcmc(Subj[seq(burnin + 1, nsamp, by = n_thin), ], start =
                           burnin + 1, thin = n_thin))
    par_subj <-
      round(cbind(
        matrix(matrix(summ2[[1]], ncol = 4)[, -4], ncol = 3),
        matrix(summ2[[2]], ncol = 5),
        "N_eff" = coda::effectiveSize(coda::mcmc(
          Subj[seq(burnin + 1, nsamp, by = n_thin), ], start = burnin + 1, thin = n_thin
        ))
      ), 3)
    colnames(par_subj) <- colnames(par_summ)
    rownames(par_subj) <- colnames(Subj)
    samples <-
      c(samples, list(subj = coda::mcmc(
        Subj[seq(burnin + 1, nsamp, by = n_thin), ], start = burnin + 1, thin = n_thin
      )))
    summaries <- c(summaries, list(subj = par_subj))
  }

  if ("Marginal" %in% functional) {
    Xtilde_unique <- unique(Xtilde)
    Marginal <- matrix(nrow = nsamp, ncol = nrow(Xtilde_unique))
    colnames(Marginal) <- rep("null", nrow(Xtilde_unique))
    for (i in 1:nrow(Xtilde_unique)) {
      Marginal[, i] <- exp(c(t(matrix(Xtilde_unique[i, ])) %*% t(output$beta)) +
                             0.5 * apply(output$tau2, 1, sum) + 0.5 * output$sigma2)
      colnames(Marginal)[i] <- paste("Marginal", i)
    }
    summ3 <-
      summary(coda::mcmc(Marginal[seq(burnin + 1, nsamp, by = n_thin), ], start =
                           burnin + 1, thin = n_thin))
    par_marg <-
      round(cbind(
        matrix(matrix(summ3[[1]], ncol = 4)[, -4], ncol = 3),
        matrix(summ3[[2]], ncol = 5),
        "N_eff" = coda::effectiveSize(coda::mcmc(
          Marginal[seq(burnin + 1, nsamp, by = n_thin), ], start = burnin + 1, thin = n_thin
        ))
      ), 3)
    colnames(par_marg) <- colnames(par_summ)
    rownames(par_marg) <- colnames(Marginal)
    samples <-
      c(samples, list(marg = coda::mcmc(
        Marginal[seq(burnin + 1, nsamp, by = n_thin), ], start = burnin + 1, thin = n_thin
      )))
    summaries <- c(summaries, list(marg = par_marg))
  }

  if ("PostPredictive" %in% functional) {
    summ4 <-
      summary(coda::mcmc(output$yrep[seq(burnin + 1, nsamp, by = n_thin), ], start =
                           burnin + 1, thin = n_thin))
    par_pred <-
      round(cbind(
        matrix(matrix(summ4[[1]], ncol = 4)[, -4], ncol = 3),
        matrix(summ4[[2]], ncol = 5),
        "N_eff" = coda::effectiveSize(coda::mcmc(
          output$yrep[seq(burnin + 1, nsamp, by = n_thin), ], start = burnin + 1, thin = n_thin
        ))
      ), 3)
    colnames(par_pred) <- colnames(par_summ)
    for (i in 1:nrow(Xtilde)) {
      rownames(par_pred)[i] <- paste("Pred", i)
    }
    samples <-
      c(samples, list(pred = coda::mcmc(
        output$yrep[seq(burnin + 1, nsamp, by = n_thin), ], start = burnin + 1, thin = n_thin
      )))
    summaries <- c(summaries, list(pred = par_pred))
  }

  par_prior <- rbind(par_sigma, matrix(unlist(par_tau), ncol = 3, byrow = T))
  colnames(par_prior)<-c("lambda", "delta", "gamma")
  rownames(par_prior)[1]<-"sigma2"
  rownames(par_prior)[2:(s+1)]<-paste("tau2_",names(Zlist),sep = "")

  return(list(par_prior=par_prior,
              samples = samples,
              summaries = summaries))
}

