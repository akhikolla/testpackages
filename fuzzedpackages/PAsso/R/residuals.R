################################################################################
# S3 generic functions are developed for the class of "clm", "glm", "lrm", "orm",
# "polr", "vgam", or "vglm".
#################################################################################

#' @title Extract Model Residuals
#'
#' @description A generic function to simulate surrogate residuals for cumulative link
#' regression models using the latent method described in Liu and Zhang (2017).
#'
#' It also support the sign-based residuals (Li and Shepherd, 2010), generalized
#' residuals (Franses and Paap, 2001), and deviance residuals for cumulative link
#' regression models.
#'
#' @param object An object of class \code{\link[ordinal]{clm}},
#' \code{\link[stats]{glm}}, \code{\link[rms]{lrm}}, \code{\link[rms]{orm}},
#' \code{\link[MASS]{polr}}, or \code{\link[VGAM]{vglm}}.
#'
#' @param type The type of residuals which should be returned. The alternatives
#' are: "surrogate" (default), "sign", "general", and "deviance". Can be abbreviated.
#' \describe{
#'   \item{\code{surrogate}}{surrogate residuals (Liu and Zhang, 2017);}
#'   \item{\code{sign}}{sign-based residuals;}
#'   \item{\code{general}}{generalized residuals (Franses and Paap, 2001);}
#'   \item{\code{deviance}}{deviance residuals (-2*loglik).}
#' }
#'
#' @param jitter A character string specifying which method to use to generate the
#' surrogate response values. Current options are \code{"latent"} and
#' \code{"uniform"}. Default is \code{"latent"}.
#' \describe{
#'   \item{\code{latent}}{latent approach;}
#'   \item{\code{uniform}}{jittering uniform approach.}
#' }
#'
#' @param jitter.uniform.scale A character string specifying the scale on which to perform
#' the jittering whenever \code{jitter = "uniform"}. Current options are
#' \code{"response"} and \code{"probability"}. Default is \code{"response"}.
#'
#' @param nsim An integer specifying the number of replicates to use.
#' Default is \code{1L} meaning one simulation only of residuals.
#'
#' @param ... Additional optional arguments. (Currently ignored.)
#'
#' @return A numeric vector of class \code{c("numeric", "resids")} containing
#' the simulated surrogate residuals. Additionally, if \code{nsim} > 1,
#' then the result will contain the attributes:
#' \describe{
#'   \item{\code{draws}}{A matrix with \code{nsim} columns, one for each
#'   is a replicate of the surrogate residuals. Note, they correspond
#'   to the original ordering of the data;}
#'   \item{\code{draws_id}}{A matrix  with \code{nsim} columns. Each column
#'   contains the observation number each surrogate residuals corresponds to in
#'   \code{draws}. (This is used for plotting purposes.)}
#' }
#'
#' @note
#' Surrogate response values require sampling from a continuous distribution;
#' consequently, the result will be different with every call to
#' \code{surrogate}. The internal functions used for sampling from truncated
#' distributions are based on modified versions of
#' \code{\link[truncdist]{rtrunc}} and \code{\link[truncdist]{qtrunc}}.
#'
#' For \code{"glm"} objects, only the \code{binomial()} family is supported.
#'
#' @references
#' Liu, Dungang and Zhang, Heping. Residuals and Diagnostics for Ordinal
#' Regression Models: A Surrogate Approach.
#' \emph{Journal of the American Statistical Association}.
#' \url{http://www.tandfonline.com/doi/abs/10.1080/01621459.2017.1292915?journalCode=uasa20}
#'
#' Li, C., & Shepherd, B. E. (2010). Test of association between two ordinal
#' variables while adjusting for covariates. \emph{Journal of the American Statistical
#' Association}, 105(490), 612-620. \url{https://doi.org/10.1198/jasa.2010.tm09386}
#'
#'
#' Franses, Philip Hans, and Richard Paap. Quantitative models in marketing research.
#' Cambridge University Press, 2001.
#' \url{https://pdfs.semanticscholar.org/dad0/820f287a8cf5a4e8039549e35fc111fd86e5.pdf}
#'
#' @importFrom stats residuals
#' @name residuals
#' @method residuals clm
#' @export
#'
#' @examples
#' # Generate data from a quadratic probit model
#' set.seed(101)
#' n <- 2000
#' x <- runif(n, min = -3, max = 6)
#' z <- 10 + 3 * x - 1 * x^2 + rnorm(n)
#' y <- ifelse(z <= 0, yes = 0, no = 1)
#'
#' # Scatterplot matrix
#' pairs(~ x + y + z)
#'
#' # Misspecified mean structure
#' fm1 <- glm(y ~ x, family = binomial(link = "probit"))
#' diagnostic.plot(fm1)
#'
#' # Correctly specified mean structure
#' fm2 <- glm(y ~ x + I(x ^ 2), family = binomial(link = "probit"))
#' diagnostic.plot(fm2)
#'
residuals.clm <- function(object,
                          type = c("surrogate", "sign", "general", "deviance"),
                          jitter = c("latent", "uniform"),
                          jitter.uniform.scale = c("probability", "response"),
                          nsim = 1L, ...) {
  # object=fit;
  # type = "surrogate"
  # jitter = "latent"
  # jitter.uniform.scale = "probability"
  # nsim = 30

  # Sanity check
  if (!inherits(object, c("clm", "glm", "lrm", "orm", "polr", "vglm"))) {
    stop(deparse(substitute(object)), " should be of class \"clm\", \"glm\", ",
         "\"lrm\", \"orm\", \"polr\", \"vgam\", or \"vglm\".")
  }

  # Match arguments
  type <- match.arg(type)
  jitter <- match.arg(jitter)
  jitter.uniform.scale <- match.arg(jitter.uniform.scale)

  # Switch different type of residuals, need to extract "latent" and "jitter" from "surrogate"
  if (type == "surrogate") {
    if (jitter == "uniform") {
      # Issue warning for jittering surrogate method
      message("Jittering uniform is an experimental feature, use at your own risk!")
    }
    gene_method <- jitter

  } else { # When type is "sign", "general", or "deviance", extract "method" for generate_residuals()
    gene_method <- type
  }

  # Generate surrogate response values
  r <- generate_residuals(object, method = gene_method, jitter.uniform.scale = jitter.uniform.scale)

  # Multiple samples
  if (nsim > 1L) {  # multiple draws
    draws <- draws_id <- matrix(nrow = nobs(object), ncol = nsim)
    for(i in seq_len(nsim)) {
      # draws_id[, i] <- sample(nobs(object), replace = TRUE)
      # BUG FIXED: Above original code is not correct! Replicate to get many draws of residuals!
      draws_id[, i] <- seq_along(getResponseValues(object))
      draws[, i] <-
        generate_residuals(object, method = gene_method, jitter.uniform.scale = jitter.uniform.scale,
                           draws_id = draws_id[, i, drop = TRUE])
    }
    attr(r, "draws") <- draws
    attr(r, "draws_id") <- draws_id
  }
  attr(r, "names") <- NULL
  attr(r, "arguments") <- c(type, jitter, jitter.uniform.scale)

  # Return residuals
  class(r) <- c("numeric", "resid")
  r

}


#' @rdname residuals
#' @method residuals lrm
#' @export
residuals.lrm <- residuals.clm



#' @rdname residuals
#' @method residuals orm
#' @export
residuals.orm <- residuals.clm



#' @rdname residuals
#' @method residuals polr
#' @export
residuals.polr <- residuals.clm


#' @rdname residuals
#' @method residuals vglm
#' @export
residuals.vglm <- residuals.clm

#' p_adj_cate
#'
#' @param Z A numerical vector that inputs the latent variable for generating probabilities of adjacent
#' categories regression model.
#'
#' @return A matrix (n by level of respones plus 1) of probabilities of the adjacent categories model.
#'
#' @keywords internal
p_adj_cate <- function(Z){
  k <- ncol(Z)
  p1_pj <- p1_pj_inv<- Z
  ZZ <- 0
  for(j in 1:k){
    ZZ <- ZZ+Z[,j]
    p1_pj[,j] <- exp(ZZ)
    p1_pj_inv[,j] <- 1/p1_pj[,j]
  }
  p1 <- 1/(rowSums(p1_pj_inv)+1)
  pj <- p1_pj
  for(j in 1:k){
    pj[,j] <- p1/p1_pj[,j]
  }
  cbind(1-rowSums(pj), pj)
}


#' generate_residuals_acat
#'
#' @param y A vector inputs the response variable.
#'
#' @param X A data.frame inputs the covariates.
#' @param alpha A vector provides the estimated intercepts of adjacent categories model. If the response
#' has k levels, there should be k+1 numbers in this alpha argument with the k-1 estimated intercepts.
#' The lower bound and upper bound are "-Inf" and "Inf".
#' @param beta A vector provides the estimated coefficients.
#' @param nsim A number to specify the replication of residuals.
#'
#' @return A vector or a matrix (nsim>1) of residuals for the adjacent categories model.
#'
#' @keywords internal
generate_residuals_acat <- function(y, X, alpha, beta, nsim=1){
  # y = y; X = X; alpha = alphas; beta = betas; nsim=1
  # alpha <- matrix(alpha, nrow=1)
  # beta <- matrix(beta, nrow=1)
  n <- length(y)
  z <- sapply(alpha[2:(length(alpha)-1)], function(a) a + tcrossprod(as.matrix(X), beta))
  p_acat <- p_adj_cate(z)
  F_acat <- t(apply(p_acat, 1, cumsum))
  F_acat <- cbind(0, F_acat)
  if(min(y)==0){
    R <- sapply(1:n, function(k) runif(nsim, F_acat[k,y[k]+1], F_acat[k,y[k]+2]))
  }else{
    R <- sapply(1:n, function(k) runif(nsim, F_acat[k,y[k]], F_acat[k,y[k]+1]))
  }
  R - 1/2
}

#' This is a function to deal with the vglm object in S4.
#'
#' @param object An object of class \code{\link[VGAM]{vglm}}.
#'
#' @param type The type of residuals which should be returned. The alternatives
#' are: "surrogate" (default), "sign", "general", and "deviance". Can be abbreviated.
#' \describe{
#'   \item{\code{surrogate}}{surrogate residuals (Liu and Zhang, 2017);}
#'   \item{\code{sign}}{sign-based residuals;}
#'   \item{\code{general}}{generalized residuals (Franses and Paap, 2001);}
#'   \item{\code{deviance}}{deviance residuals (-2*loglik).}
#' }
#' @param jitter A character string specifying which method to use to generate the
#' surrogate response values. Current options are \code{"latent"} and
#' \code{"uniform"}. Default is \code{"latent"}.
#' \describe{
#'   \item{\code{latent}}{latent approach;}
#'   \item{\code{uniform}}{jittering uniform approach.}
#' }
#' @param jitter.uniform.scale A character string specifying the scale on which to perform
#' the jittering whenever \code{jitter = "uniform"}. Current options are
#' \code{"response"} and \code{"probability"}. Default is \code{"response"}.
#' @param nsim An integer specifying the number of replicates to use.
#' Default is \code{1L} meaning one simulation only of residuals.
#' @param ... Additional optional arguments.
#'
#' @return A "resid" object with attributes. It contains a vector or a matrix (nsim>1) of
#' residuals for the adjacent categories model.
#'
#' @export
#' @keywords internal
residualsAcat <- function(object,
                          type = c("surrogate", "sign", "general", "deviance"),
                          jitter = c("latent", "uniform"),
                          jitter.uniform.scale = c("probability", "response"),
                          nsim = 1L, ...)
{
  # object <- fit
  # type = "surrogate"
  # jitter = "latent"
  # jitter.uniform.scale = "probability"
  # nsim = 1L

  # Match arguments -------------
  type <- match.arg(type)
  jitter <- match.arg(jitter)
  jitter.uniform.scale <- match.arg(jitter.uniform.scale)

  # Switch different type of residuals, need to extract "latent" and "jitter" from "surrogate" -------------
  if (type == "surrogate") {
    if (jitter == "uniform") {
      # Issue warning for jittering surrogate method
      message("Jittering uniform is an experimental feature, use at your own risk!")
    }
    gene_method <- jitter
  } else { # When type is "sign", "general", or "deviance", extract "method" for generate_residuals()
    gene_method <- type
  }

  # Generate surrogate response values -------------
  coefs <- coef(object)
  alphas <- matrix(c(-Inf, coefs[1:(ncat(object)-1)], Inf), nrow = 1)
  betas <- matrix(coefs[-c(1:(ncat(object)-1))], nrow = 1)
  # message(alphas, "  ", betas)
  y <- getResponseValues(object)

  # X <- as.matrix(model.frame(object)[,-1])
  X <- as.matrix(object@x[,-1]) # This just work for the "vglm" adjacent categories model!

  r <- generate_residuals_acat(y = y, X = X, alpha = alphas, beta = betas, nsim=1)

  # Multiple samples -------------
  if (nsim > 1L) {  # multiple draws
    draws_id <- matrix(seq(nobs(object)), nrow = nobs(object), ncol = nsim, byrow = F)

    draws <-
      t(generate_residuals_acat(y = y, X = X, alpha = alphas, beta = betas, nsim = nsim))

    attr(r, "draws") <- draws
    attr(r, "draws_id") <- draws_id
  }
  attr(r, "names") <- NULL
  attr(r, "arguments") <- c(type, jitter, jitter.uniform.scale)

  # Return residuals
  class(r) <- c("numeric", "resid")

  return(r)
}



################################################################################
# GLM has its own Generic function for cumulative link models;
# I combine my new features with the recent version of residuals.glm (stats 3.6.3).
# 1. It can generate "surrogate", "sign", "general", and "deviance" residuals.
# 3. It has two approaches to generate surrogate residuals: latent and uniform(jittering).
################################################################################


#' @name residuals
#' @method residuals ord
#'
#' @return
#' @export
residuals.ord <- function (
  object,
  type = c("surrogate", "sign", "general", "deviance",
           "pearson", "working", "response", "partial"),
  jitter = c("latent", "uniform"),
  jitter.uniform.scale = c("probability", "response"),
  nsim = 1L, ...)
{
  type <- match.arg(type)

  # Sanity check: add our new types of residuals to residuals.glm()
  if (type %in% c("surrogate", "sign", "general", "deviance")) {
    res <- residuals.clm(object = object,
                         type = type,
                         jitter = jitter,
                         jitter.uniform.scale = jitter.uniform.scale,
                         nsim = nsim,...)
  } else { # Keep original below: residuals.glm()
    y <- object$y
    r <- object$residuals
    mu <- object$fitted.values
    wts <- object$prior.weights
    switch(type, deviance = , pearson = , response = if (is.null(y)) {
      mu.eta <- object$family$mu.eta
      eta <- object$linear.predictors
      y <- mu + r * mu.eta(eta)
    })
    res <- switch(type, deviance = if (object$df.residual > 0) {
      d.res <- sqrt(pmax((object$family$dev.resids)(y, mu,
                                                    wts), 0))
      ifelse(y > mu, d.res, -d.res)
    } else rep.int(0, length(mu)), pearson = (y - mu) * sqrt(wts)/sqrt(object$family$variance(mu)),
    working = r, response = y - mu, partial = r)
    if (!is.null(object$na.action))
      res <- naresid(object$na.action, res)
    if (type == "partial")
      res <- res + predict(object, type = "terms")
    res
  }
}

#' Extract Residuals from a Partial Association Analysis
#'
#' @param object An object of class \code{PAsso}.
#'
#' @param draw_id A number refers to the i-th draw of residuals.
#'
#' @param ... Additional optional arguments.
#'
#' @return A matrix of class \code{c("matrix", "resids")} containing
#' the simulated surrogate residuals used for the partial association
#' analysis in \code{PAsso}. Additionally, if \code{rep_num} > 1 in \code{PAsso},
#' then the result will contain the attributes:
#' \describe{
#'   \item{\code{draws}}{An array contains all draws of residuals.}
#' }
#'
#'
#' @importFrom stats residuals
#' @name residuals
#' @method residuals PAsso
#'
#' @return
#' @export
#'
#' @examples
#' # Load data
#' data("ANES2016")
#' PAsso_1 <- PAsso(responses = c("PreVote.num", "PID"),
#'                  adjustments = c("income.num", "age", "edu.year"),
#'                  data = ANES2016)
#'
#' # Compute residuals
#' res1 <- residuals(PAsso_1)
#'
residuals.PAsso <- function(object, draw_id=1, ...) {
  if ((draw_id>=1) & (draw_id<=dim(object$rep_SRs)[2])) {
    resids_PAsso <- object$rep_SRs[,draw_id,]
    # resids_PAsso <- PAsso_1$rep_SRs[,1,]
  } else {
    stop("The draw_id is out of bound.")
  }
  attr(resids_PAsso, "draws") <- object$rep_SRs
  attr(resids_PAsso, "arguments") <- attr(object, "arguments")
  # Return residuals
  class(resids_PAsso) <- c("resid", class(resids_PAsso))
  return(resids_PAsso)
}
