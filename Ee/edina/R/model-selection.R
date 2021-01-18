
## Model selection ----

#' Extract the Best Model
#'
#' Extracts the best model from the `auto_*()` search procedure.
#'
#' @param x   An `auto_edina` object
#' @param ic  Information criterion name. Default `"ppp"`.
#' @param ... Not used.
#'
#' @return
#' An `edina` model object corresponding to the smallest value of requested
#' information criterion.
#'
#' @seealso
#' [DIC.edina()],
#' [BIC.edina()],
#' [PPP.edina()]
#'
#' @export
#' @rdname best_model
best_model = function(x, ...) {
    UseMethod("best_model")
}

#' @export
#' @rdname best_model
best_model.auto_edina = function(x, ic = c("ppp", "bic", "dic"), ...) {
    ic = tolower(ic)
    ic = match.arg(ic)

    x$edina_models[[which.min(x$criterion[,ic])]]
}

#' @export
best_model.default = function(x, ...) {
    stop_bad_class(x, "auto_edina")
}


### PPP Criteria ----

#' Posterior Predictive Probabilities (PPPs)
#'
#' Computes posterior predictive probabilities (PPPs) based on the
#' odds ratios for each pair of items.
#'
#' @inheritParams summary.edina
#'
#' @return
#' The PPP value given the specified `alpha` value.
#'
#' @details
#'
#' PPPs that smaller than 0.05 or greater than 0.95 tend to be extreme and
#' evidence of misfit. As a result, this is more of a heuristic metric.
#'
#' @section PPP Computation Procedure:
#'
#' 1. simulate observed responses \eqn{\mathbf Y^{(r)}} using model parameters
#'    from iteration \eqn{r} of the MCMC sampler
#' 2. computing the odds ratio for each pair of items at iteration \eqn{r} as
#'    \deqn{OR^{(r)} = n_{11}^{(r)}n_{00}^{(r)}/\left(n_{10}^{(r)}n_{01}^{(r)}\right)},
#'    where \eqn{n_{11}^{(r)}} is the frequency of ones on both variables at
#'    iteration \eqn{r}, \eqn{n_{10}^{(r)}} is the frequency of ones on the
#'    first item and zeros on the second at iteration \eqn{r}, etc.; and
#' 3. computing PPPs for each item pair as the proportion of generated
#'    \eqn{OR^{(r)}}'s that exceeded elements of the observed odds ratios.
#'
#' @rdname PPP
#' @export
PPP = function(object, ...) {
    UseMethod("PPP")
}

#' @rdname PPP
#' @export
PPP.edina = function(object, alpha = 0.05, ...) {
    or_tested = object$or_tested

    mean(or_tested[upper.tri(or_tested)] < alpha |
             or_tested[upper.tri(or_tested)] > (1-alpha))
}

#' Deviance Information Criterion (DIC)
#'
#' Calculate DIC for EDINA models.
#'
#' @param object An `edina` object
#' @param ...    Not used.
#'
#' @return
#' The DIC value of the given model.
#'
#' @seealso [PPP.edina()], [BIC.edina()]
#'
#' @section DIC Computation Procedure:
#'
#' \eqn{DIC = -2\left({\log p\left( {\mathbf{y}| \mathbf{\hat{\theta}} } \right)  - 2\left( {\log p\left( {\mathbf{y}| \mathbf{\hat{\theta}} } \right) - \frac{1}{N}\sum\limits_{n = 1}^N {\log p\left( {\mathbf{y}|{\mathbf{\theta} _s}} \right)} } \right)} \right)}
#'
#' @rdname DIC
#' @export
DIC = function(object, ...) {
    UseMethod("DIC")
}

#' @rdname DIC
#' @export
DIC.edina = function(object, ...) {

    # LogLikelihood at the Mean of the Posterior Distributioon
    L = object$loglike_pmean

    # Number of parameters in the model
    P = 2 * (L - (1 / object$chain_length * object$loglike_summed))

    -2 * (L - P)
}


## 2 * K - I_K fixed allothers to 1's. Then estimate R* matrix and the J pi*

#' Bayesian Information Criterion (BIC)
#'
#' Calculate BIC for EDINA models.
#'
#' @param object An `edina` object
#' @param ... Not used.
#'
#' @return
#' The BIC value of the given model.
#'
#' @seealso [PPP.edina()], [DIC.edina()]
#'
#' @section BIC Computation Procedure:
#'
#' \eqn{BIC = -2 \log p\left( {\mathbf{y}| \mathbf{\hat{\theta}} } \right) + ((k+2)*j + 2^k)\log(n)}
#'
#' @export
#' @seealso [PPP.edina()], [DIC.edina()]
#' @importFrom stats BIC
BIC.edina = function(object, ...) {
    # K number of attributes and J number of items (count only one slipping/guessing)?
    # -2 * LogLike + ln(n) * (k + j)
    -2*object$loglike_pmean + log(object$n)*((object$k+2)*object$j + 2^object$k)
}
