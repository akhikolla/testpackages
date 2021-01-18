library(Rcpp)
library(RcppArmadillo)

#' Calculate residuals restricted under H0
#'
#' Given regression \code{model} and \code{clustering}, this function calculates the OLS residuals
#' under the linear null hypothesis, and assigns them to the specified clusters.
#'
#' @param model A regression \code{model}. See \link{example_model} for details.
#' @param clustering A \code{List} that specifies a clustering of indexes 1...n (#datapoints).
#' See \link{example_clustering} for details.
#'
#' @return A \code{List} of the restricted residuals clustered according to \code{clustering}.
#' @examples
#' m = example_model(n=100)
#' cl = list(1:50, 51:100)
#' er = get_clustered_eps(m, cl)
#' stopifnot(length(er) == length(cl))
#' stopifnot(length(er[[1]]) == 50)
get_clustered_eps = function(model, clustering) {
  y = model$y; X = model$X; lam = model$lam; lam0 = model$lam0
  bhat = OLS_c(y, X)
  Q = matrix(lam, ncol=1)
  bhat_r = restricted_OLS_c(y, X, bhat, Q=Q, c=lam0)
  er = y - X %*% bhat_r # restricted residuals.
  stopifnot(all(!is.na(er))) ## no NAs in residuals.
  # return
  lapply(clustering, function(i) er[i])
}

#' One-sided testing
#'
#' Decides to reject or not based on observed test statistic value \code{tobs} and randomization values \code{tvals}.
#'
#' @param tobs The observed value of the test statistic (scalar).
#' @param tvals Vector of randomization values of the test statistic (to compare with \code{tobs}).
#' @param alpha Desired level of the test (between 0 to 1).
#' @param tol Used to check whether \code{tobs} is equal to the 1-\code{alpha} quantile of \code{tvals}.
#' @details
#' The test may randomize to achieve the specified level \code{alpha}
#' when there are very few randomization values.
#' @return Test decision (binary).
#' @seealso Testing Statistical Hypotheses (Ch. 15, Lehman and Romano, 2006)
one_sided_test = function(tobs, tvals, alpha, tol=1e-14) {
  srt = sort(tvals)
  M = length(tvals)
  k = ceiling(M * (1-alpha))
  Tk = srt[k]
  if(abs(tobs - Tk) < tol) {
    # if tobs = Tk
    ax = (M * alpha - sum(tvals > Tk)) / sum(abs(tvals - Tk) < tol)
    return(runif(1) <= ax) ## randomize decision.
  }

  return(tobs > Tk)
}

#' Two-sided testing
#'
#' Decides to reject or not based on observed test statistic value \code{tobs}
#' and randomization values \code{tvals}. The test may randomize to achieve the specified level \code{alpha}
#' when there are very few randomization values.
#'
#' @param tobs The observed value of the test statistic (scalar).
#' @param tvals Vector of randomization values of the test statistic (to compare with \code{tobs}).
#' @param alpha Desired level of the test (between 0 to 1).
#'
#' @return Test decision (binary).
#' @seealso Testing Statistical Hypotheses (Ch. 15, Lehman and Romano, 2006)
two_sided_test = function(tobs, tvals, alpha) {
  m1 = one_sided_test(tobs, tvals, alpha=alpha/2)
  m2 = one_sided_test(-tobs, -tvals, alpha=alpha/2) # only one can be 1.
  return(m1 + m2)
}

#' Calculates p-value or test decision
#'
#' Depending on \code{ret_pval} this function returns either a p-value for the test or the binary decision.
#'
#' @param rtest_out A \code{List} with elements \code{tobs}, \code{tvals} (see \link{one_sided_test} for details.)
#' @param ret_pval A \code{Boolean} indicating whether to return a p-value (TRUE) or not.
#' @param alpha Desired test level (from 0 to 1).
#' @return Binary decision if \code{ret_pval} is TRUE, or the p-value otherwise.
#' @details Returns 1 if the test rejects, 0 otherwise.
out_pval = function(rtest_out, ret_pval, alpha) {
  tobs = rtest_out$tobs
  tvals = c(rtest_out$tvals)

  n_all = length(tvals)
  n_higher = sum(tvals > (tobs + 1e-12))
  n_lower = sum(tvals < (tobs - 1e-12))
  n_equal = n_all - n_lower - n_higher

  p1 = (n_equal + n_higher) / n_all  # P(T >= Tobs)
  p2 = (n_equal + n_lower) / n_all  # P(T <= Tobs)

  pval = min(p1, p2)
  if(ret_pval) return(pval)
  return(two_sided_test(tobs, tvals, alpha = alpha))  # this test is less conservative.
}

#' Example regression model and H0.
#'
#' @param n Number of datapoints.
#' @return List of (y, X, lam, lam0) that corresponds to regression model and null hypothesis:
#' \itemize{
#'   \item y = n-length vector of outcomes
#'   \item X = n x p covariate matrix;
#'   \item lam = p-vector of coefficients
#'   \item lam0 = real number.
#'   }
#'
#' The null we are testing through this specification is
#'
#'         H0: lam' beta = lam[1] * beta[1] + ... + lam[p] * beta[p] = lam0,
#'
#' where beta are the model parameters in the regression, y = X beta + e.
#' By default this example sets p = 2-dim model, lam = (0, 1) and lam0 = 0. In this specification, H0: beta[2] = 0.
#'
#' @examples
#' model = example_model()
#' lm(model$y ~ model$X + 0)
example_model = function(n=100) {
  lam = c(0, 1)
  lam0 = 0
  X = cbind(rep(1, n), rnorm(n))
  eps = rnorm(n)
  beta = c(0, 0)
  y = X %*% beta + eps
  #
  return(list(y=y, X=X, lam=lam, lam0=lam0))
}

#' Checks whether the input \code{model} is valid.
#' @param model A \code{model} object. See \link{example_model} for details.
check_model = function(model) {
  if(!all(model$X[, 1] == 1)) {
    warning("No intercept.")
  }
  with(model, stopifnot(length(y) == nrow(X)))
  with(model, stopifnot(length(lam) == ncol(X)))
  with(model, stopifnot(length(lam0) == 1))
}

#' An example \code{clustering} object. A clustering is a \code{List} that
#' splits indexes 1..#num_datapoints to clusters. Each \code{List} element corresponds to one cluster.
#' The clustering is not necessarily a partition but it usually is.
#'
#' @return A \code{List} for the clustering of indexes 1..#num_datapoints.
example_clustering = function() {
  n = 100
  cl = list(1:30, 31:50, 51:100) # 3 clusters
  return(cl)
}

#' Generic residual randomization test
#'
#' This function tests the specified linear hypothesis in \code{model} assuming the errors are distributionally invariant
#' with respect to stochastic function \code{g_invar}.
#'
#' @param model Regression model and hypothesis. See \link{example_model} for details.
#' @param g_invar Stochastic function that transforms residuals. Accepts n-vector and returns n-vector.
#' @param num_R Number of test statistic values to calculate in the randomization test.
#' @param alpha Nominal test level (between 0 to 1).
#' @param val_type The type of return value.
#' @return
#' If \code{val_type} = "decision" (default) we get the test binary decision (1=REJECT H0).
#'
#' If \code{val_type} = "pval" we get the test p-value.
#'
#' If  \code{val_type} = "full" we get the full test output, i.e., a \code{List} with elements \code{tobs}, \code{tvals},
#' the observed and randomization values of the test statistic, respectively.
#'
#' @details
#' For the regression y = X * beta + e, this function is testing the following linear null hypothesis:
#'
#'         H0: lam' beta = lam[1] * beta[1] + ... + lam[p] * beta[p] = lam0,
#'
#' where y, X, lam, lam0 are specified in \code{model}.
#' The assumption is that the errors, e, have some form of cluster invariance.
#' Specifically:
#'
#' (e_1, e_2, ..., e_n) ~  g_invar(e_1, e_2, ..., e_n),
#'
#' where ~ denotes equality in distribution, and \code{g_invar} is the supplied
#' invariance function.
#'
#' @note
#' There is no guarantee that an arbitrary \code{g_invar} will produce valid tests.
#' The \link{rrtest_clust} function has such guarantees under mild assumptions.
#'
#' @examples
#' model = example_model(n = 100)  # test H0: beta2 = 0 (here, H0 is true)
#' g_invar = function(e) sample(e)   # Assume errors are exchangeable.
#' rrtest(model, g_invar) # same as rrtest_clust(model, "perm")
#'
#' @seealso	Life after bootstrap: residual randomization inference in regression models (Toulis, 2019)
#'
#' \url{https://www.ptoulis.com/residual-randomization}
#'
#' @export
rrtest = function(model, g_invar, num_R=999, alpha=0.05, val_type="decision") {
  check_model(model)
  stopifnot(class(g_invar) == "function")
  stopifnot(class(num_R) == "numeric")

  y = model$y; X = model$X; lam=model$lam; lam0 = model$lam0
  # Tn(u) = (X^T X)^-1 X^T u
  Tn = function(eps) {
    b = fastLm(eps, X + 0)$coefficients
    sum(lam * b)
  }

  # 1. tobs
  bhat = fastLm(y, X + 0)$coefficients
  tobs = (sum(lam * bhat) - lam0) # same as tobs = Tn(y) - lam0
  stopifnot(abs(tobs - Tn(y) + lam0) < 1e-8)

  # 2. e_r : restricted residuals.
  Q = matrix(lam, ncol=1)
  bhat_r = restricted_OLS_c(y, X, bhat, Q=Q, c=lam0)
  stopifnot(abs(sum(lam * bhat_r) - lam0) < 1e-8) # check the restriction condition

  er = y - X %*% bhat_r # restricted residuals.
  tvals = c()
  for(r in 1:num_R) {
    # 1. iid
    er_new = g_invar(er)
    tvals = c(tvals, Tn(er_new))
  }

  out = list(tobs=tobs, tvals=tvals)
  # Return the entire test output?
  if(val_type == "full") return(out)
  # Return pvalue or decision
  ret_pval = as.logical(val_type == "pval")
  out_pval(out, ret_pval, alpha)
}


#' Residual randomization test under cluster invariances
#'
#' This function tests the specified linear hypothesis in \code{model}
#' assuming that the errors have some form of cluster invariance determined by \code{type}
#' within the clusters determined by \code{clustering}.
#'
#' @param model Regression model and hypothesis. See \link{example_model} for details.
#' @param type A \code{character}, either "perm", "sign" or "double".
#' @param clustering A \code{List} that specifies a clustering of datapoint indexes {1, ..., n}.
#' See \link{example_clustering}. If NULL it takes default value according to \code{type} (see Note)
#' @param num_R Number of test statistic values to calculate in the test.
#' @param alpha Nominal test level (between 0 to 1).
#' @param val_type The type of return value.
#' @return
#' If \code{val_type} = "decision" (default) we get the test binary decision (1=REJECT H0).
#'
#' If \code{val_type} = "pval" we get the test p-value.
#'
#' If  \code{val_type} = "full" we get the full test output, i.e., a \code{List} with elements \code{tobs}, \code{tvals},
#' the observed and randomization values of the test statistic, respectively.
#'
#' @details
#' For the regression y = X * beta + e, this function is testing the following linear null hypothesis:
#'
#'         H0: lam' beta = lam[1] * beta[1] + ... + lam[p] * beta[p] = lam0,
#'
#' where y, X, lam, lam0 are specified in \code{model}.
#' The assumption is that the errors, e, have some form of cluster invariance.
#' Specifically:
#' \itemize{
#'   \item If \code{type} = "perm" then the errors are assumed exchangeable within the specified clusters:
#'
#' (e_1, e_2, ..., e_n) ~ cluster_perm(e_1, e_2, ..., e_n),
#'
#' where ~ denotes equality in distribution, and cluster_perm is any random permutation
#' within the clusters defined by \code{clustering}. Internally, the test repeatedly calculates a test statistic
#' by randomly permuting the residuals within clusters.
#'
#'   \item If \code{type} = "sign" then the errors are assumed sign-symmetric within the specified clusters:
#'
#' (e_1, e_2, ..., e_n) ~ cluster_signs(e_1, e_2, ..., e_n),
#'
#' where cluster_signs is a random signs flip of residuals on the cluster level.
#' Internally, the test repeatedly calculates a test statistic by randomly flipping the signs of cluster residuals.
#'
#'   \item If \code{type} = "double" then the errors are assumed both exchangeable and sign symmetric
#'    within the specified clusters:
#'
#' (e_1, e_2, ..., e_n) ~ cluster_signs(cluster_perm(e_1, e_2, ..., e_n)),
#'
#' Internally, the test repeatedly calculates a test statistic by permuting and randomly flipping the signs of residuals on the cluster level.
#' }
#' @note
#' If \code{clustering} is NULL then it will be assigned a default value:
#' \itemize{
#'   \item \code{list(1:n)}if \code{type} = "perm", where n is the number of datapoints;
#'   \item \code{as.list(1:n)} if \code{type} = "sign" or "double".
#'  }
#'
#' As in bootstrap \code{num_R} is usually between 1000-5000.
#'
#' @examples
#' # 1. Validity example
#' set.seed(123)
#' n = 50
#' X = cbind(rep(1, n), 1:n/n)
#' beta = c(0, 0)
#' rej = replicate(200, {
#'   y = X %*% beta  + rt(n, df=5)
#'   model = list(y=y, X=X, lam=c(0, 1), lam0=0)  # H0: beta2 = 0
#'   rrtest_clust(model, "perm")
#' })
#' mean(rej)  # Should be ~ 5% since H0 is true.
#'
#' # 2. Heteroskedastic example
#' set.seed(123)
#' n = 200
#' X = cbind(rep(1, n), 1:n/n)
#' beta = c(-1, 0.2)
#' ind = c(rep(0, 0.9*n), rep(1, .1*n))  # cluster indicator
#' y = X %*% beta + rnorm(n, sd= (1-ind) * 0.1 + ind * 5) # heteroskedastic
#' confint(lm(y ~ X + 0))  # normal OLS does not reject H0: beta2 = 0
#' cl = list(which(ind==0), which(ind==1))
#' model = list(y=y, X=X, lam=c(0, 1), lam0=0)
#'
#' rrtest_clust(model, "sign")  # errors are sign symmetric regardless of cluster.
#' # Cluster sign test does not reject because of noise.
#'
#' rrtest_clust(model, "perm", cl)  # errors are exchangeable within clusters
#' # Cluster permutation test rejects because inference is sharper.
#' @seealso	Life after bootstrap: residual randomization inference in regression models (Toulis, 2019)
#'
#' \url{https://www.ptoulis.com/residual-randomization}
#'
#' @export
rrtest_clust = function(model, type, clustering=NULL,
                        num_R=999, alpha=0.05, val_type="decision") {
  check_model(model)
  n = length(model$y)
  out = NULL
  flag_perm = flag_sign = FALSE
  if(is.null(clustering)) {
    if(type=="perm") {
      clustering = list(1:n)
    } else {
      clustering = as.list(1:n)
    }
  }

  stopifnot(length(clustering) > 0)
  if((length(clustering) < 5) & (type != "perm")) {
    warning("Too few clusters. Should be > 5 when using sign tests.")
  }

  if(type == "perm") {
    flag_perm = TRUE
  } else if(type == "sign") {
    flag_sign = TRUE
  } else if(type == "double") {
    flag_perm = TRUE
    flag_sign = TRUE
  } else {
    stop("Not recognized cluster invariance type.")
  }

  cl_eps_r = get_clustered_eps(model, clustering)
  out = with(model, r_test_c(y, X, lam=lam, lam0 = lam0,
                             cluster_eps_r=cl_eps_r,
                             use_perm = flag_perm, use_sign = flag_sign, num_R = num_R))
  if(val_type=="full") {
    return(out)
  } else {
    ret_pval = as.logical(val_type == "pval")
    return(out_pval(out, ret_pval, alpha))
  }
}

#' Generic residual randomization inference
#
#' This function provides the basis for all other rrinf* functions.
#
#' @param y Vector of outcomes (length n)
#' @param X Covariate matrix (n x p). First column should be ones to include intercept.
#' @param g_or_clust Either \code{clustering} or an invariance function that transforms residuals.
#' @param cover Number from [0, 1] that denotes the confidence interval coverage (e.g., 0.95 denotes 95\%)
#' @param num_R Number of test statistic values to calculate in the randomization test (similar to no. of bootstrap samples).
#' @param control.tinv A \code{List} that determines the test inversion.
#' @return Matrix that includes the confidence interval endpoints, and the interval midpoint estimate.
#' @details
#' This function has similar funtionality as standard \link{confint}.
#' It does so by testing plausible values for each parameter. The plausible values can be controlled as follows.
#' For some parameter beta_i we will test successively
#'
#' H0: beta_i = hat_beta_i - \code{num_se} * se_i
#'
#' ...up to...
#'
#' H0: beta_i = hat_beta_i + \code{num_se} * se_i
#'
#' broken in \code{num_breaks} intervals. Here, hat_beta_i is the OLS estimate of beta_i and se_i is the standard error.
#'
#' The \code{g_or_clust} object should either be (i) a g-invariance function R^n -> R^n; or (ii)
#' a list(type, cl) where type=c("perm", "sign", "double") and cl=\code{clustering} (see \link{example_clustering} for details).
#'
# @seealso	Life after bootstrap: residual randomization inference in regression models (Toulis, 2019)
#'
#' \url{https://www.ptoulis.com/residual-randomization}
rrinfBase = function(y, X, g_or_clust, cover, num_R, control.tinv) {

  # checks
  stopifnot(length(y) == nrow(X))
  stopifnot(cover <= 1 & cover >= 0)
  #
  p = ncol(X) # number of parameters.
  ols = lm(y ~ X + 0)
  bhat = coef(ols)
  se = summary(ols)$coefficients[, 2]  # standard errors
  search  = cbind(bhat, bhat) + control.tinv$num_se * cbind(-se, se)
  colnames(search) = c("left", "right")

  Results = matrix(0, nrow=p, ncol=3)
  colnames(Results) = c("midpoint estimate",
                        sprintf("%.1f%%", 100 * (1-cover) / 2),
                        sprintf("%.1f%%", 100 * (cover + (1 - cover) / 2)))
  for(j in 1:p) {
    lam = rep(0, p)
    lam[j] = 1
    bvals = seq(search[j, 1], search[j, 2], length.out=control.tinv$num_breaks)
    pvals = sapply(bvals, function(b) {
      model = list(y=y, X=X, lam=lam, lam0=b)
      if(class(g_or_clust) == "function") {
        rrtest(model, g_or_clust, num_R=num_R, alpha = 1-cover, val_type = "pval")
      } else if(class(g_or_clust) == "list") {
        rrtest_clust(model,
                     type=g_or_clust$type, clustering=g_or_clust$cl,
                     num_R=num_R, alpha=1-cover, val_type="pval")
      } else {
        stop("Object g_or_cl invalid.")
      }
    })
    ## Determine CI
    w = which(pvals >= (1-cover)/2)   # two-sided p-values
    ci1 = min(bvals[w]); ci2 = max(bvals[w])
    Results[j, ] = c(mean(c(ci1, ci2)), ci1, ci2)
  }
  rownames(Results) = colnames(X)
  Results
}

#' Generic residual randomization confidence intervals
#'
#' This function is a wrapper over \link{rrtest} and gives confidence intervals for all parameters.
#'
#' @param y Vector of outcomes (length n)
#' @param X Covariate matrix (n x p). First column should be ones to include intercept.
#' @param g_invar Function that transforms residuals. Accepts n-vector and returns n-vector.
#' @param cover Number from [0, 1] that denotes the confidence interval coverage (e.g., 0.95 denotes 95\%)
#' @param num_R Number of test statistic values to calculate in the randomization test (similar to no. of bootstrap samples).
#' @param control A \code{List} that constrols the scope of the test inversion.
#' @return Matrix that includes the confidence interval endpoints, and the interval midpoint estimate.
#' @details
#' This function has similar funtionality as standard \link{confint}.
#' It generates confidence intervals by testing plausible values for each parameter.
#' The plausible values are generated as follows.
#' For some parameter beta_i we test successively
#'
#' H0: beta_i = hat_beta_i - \code{num_se} * se_i
#'
#' ...up to...
#'
#' H0: beta_i = hat_beta_i + \code{num_se} * se_i
#'
#' broken in \code{num_breaks} intervals. Here, hat_beta_i is the OLS estimate of beta_i and se_i is the standard error.
#' We then report the minimum and maximum values in this search space which we cannot reject
#' at level \code{alpha}. This forms the desired confidence interval.
#' @note
#' If the confidence interval appears to be a point or is empty, then this means
#' that the nulls we consider are implausible.
#' We can try to improve the search through \code{control.tinv}.
#' For example, we can both increase \code{num_se} to increase the width of search,
#' and increase \code{num_breaks} to make the search space finer.
#'
#' @examples
#' set.seed(123)
#' X = cbind(rep(1, 100), runif(100))
#' beta = c(-1, 1)
#' y = X %*% beta + rnorm(100)
#' g_invar = function(e) sample(e)  # Assume exchangeable errors.
#' M = rrinf(y, X, g_invar, control=list(num_se=4, num_breaks=20))
#' M  # Intervals cover true values
#' @seealso	Life after bootstrap: residual randomization inference in regression models (Toulis, 2019)
#'
#' \url{https://www.ptoulis.com/residual-randomization}
#'
#' @export
rrinf = function(y, X, g_invar, cover=.95, num_R=999,
                 control=list(num_se=6, num_breaks=60)) {

  rrinfBase(y, X, g_or_clust = g_invar, cover = cover,
            num_R = num_R, control.tinv = control)
}

#' Residual randomization inference based on cluster invariances
#'
#' This function is a wrapper over \link{rrtest_clust} and gives confidence intervals for all parameters
#' assuming a particular cluster invariance on the errors.
#'
#' @param y Vector of outcomes (length n)
#' @param X Covariate matrix (n x p). First column should be ones to include intercept.
#' @param type A string, either "perm", "sign" or "double".
#' @param clustering A \code{List} that specifies a clustering of datapoint indexes {1, ..., n}. See \link{example_clustering} for details.
#' @param cover Number from [0, 1] that denotes the confidence interval coverage (e.g., 0.95 denotes 95\%)
#' @param num_R Number of test statistic values to calculate in the randomization test (similar to no. of bootstrap samples).
#' @param control A \code{List} that controls the scope of the test inversion.
#' @return Matrix that includes the OLS estimate, and confidence interval endpoints.
#' @details
#' This function has similar funtionality as standard \link{confint}.
#' It generates confidence intervals by testing plausible values for each parameter.
#' The plausible values are generated as follows.
#' For some parameter beta_i we test successively
#'
#' H0: beta_i = hat_beta_i - \code{num_se} * se_i
#'
#' ...up to...
#'
#' H0: beta_i = hat_beta_i + \code{num_se} * se_i
#'
#' broken in \code{num_breaks} intervals. Here, hat_beta_i is the OLS estimate of beta_i and se_i is the standard error.
#' We then report the minimum and maximum values in this search space which we cannot reject
#' at level \code{alpha}. This forms the desired confidence interval.
#' @note
#' If the confidence interval appears to be a point or is empty, then this means
#' that the nulls we consider are implausible.
#' We can try to improve the search through \code{control.tinv}.
#' For example, we can both increase \code{num_se} to increase the width of search,
#' and increase \code{num_breaks} to make the search space finer.
#'
#' @examples
#' # Heterogeneous example
#' set.seed(123)
#' n = 200
#' X = cbind(rep(1, n), 1:n/n)
#' beta = c(-1, 0.2)
#' ind = c(rep(0, 0.9*n), rep(1, .1*n))  # cluster indicator
#' y = X %*% beta + rnorm(n, sd= (1-ind) * 0.1 + ind * 5) # heteroskedastic
#' confint(lm(y ~ X + 0))  # normal OLS CI is imprecise
#'
#' cl = list(which(ind==0), which(ind==1))  #  define the clustering
#' rrinf_clust(y, X, "perm", cl)  # improved CI through clustered errors
#'
#' @note
#' See \link{rrtest_clust} for a description of \code{type} and \code{clustering}.
#'
#' @seealso	Life after bootstrap: residual randomization inference in regression models (Toulis, 2019)
#'
#' \url{https://www.ptoulis.com/residual-randomization}
#'
#' @export
rrinf_clust = function(y, X, type, clustering=NULL, cover=.95, num_R=999,
                      control=list(num_se=6, num_breaks=60)) {
  g_or_clust = list(type=type, cl=clustering)
  rrinfBase(y, X, g_or_clust, cover = cover, num_R = num_R, control.tinv = control)
}

