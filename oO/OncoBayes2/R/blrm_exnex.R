#' @title Bayesian Logistic Regression Model for N-compounds with EXNEX
#'
#' @description Bayesian Logistic Regression Model (BLRM) for N
#'     compounds using EXchangability and NonEXchangability (EXNEX)
#'     modeling.
#'
#' @param formula the model formula describing the linear predictors
#'     of the model. The lhs of the formula is a two-column matrix
#'     which are the number of occured events and the number of times
#'     no event occured. The rhs of the formula defines the linear
#'     predictors for the marginal models for each drug component,
#'     then the interaction model and at last the grouping and
#'     optional stratum factors of the models. These elements of the
#'     formula are separated by a vertical bar. The marginal models
#'     must follow a intercept and slope form while the interaction
#'     model must not include an interaction term. See the examples
#'     below for an example instantiation.
#' @param data optional data frame containing the variables of the
#'     model. If not found in \code{data}, the variables are taken
#'     from \code{environment(formula)}.
#' @param prior_EX_mu_mean_comp,prior_EX_mu_sd_comp Mean and sd for
#'     the prior on the mean parameters \eqn{\boldsymbol\mu_i =
#'     (\mu_{\alpha i}, \mu_{\beta i})} of each component.  Two column
#'     matrix (intercept, log-slope) with one row per component.
#' @param prior_EX_tau_mean_comp,prior_EX_tau_sd_comp Prior mean and
#'     sd for heterogeniety parameter \eqn{\boldsymbol\tau_{si} =
#'     (\tau_{\alpha s i}, \tau_{\beta s i})} of each stratum. If no
#'     differential discounting is required (i.e. if there is only one
#'     stratum \eqn{s = 1}), then it is a two-column matrix
#'     (intercept, log-slope) with one row per component. Otherwise it
#'     is a three-dimensional array whose first dimension indexes the
#'     strata, second dimension indexes the components, and third
#'     dimension of length two for (intercept, log-slope).
#' @param prior_EX_corr_eta_comp Prior LKJ correlation parameter for
#'     each component given as numeric vector. If missing, then a 1 is
#'     assumed corresponding to a marginal uniform prior of the
#'     correlation.
#' @param prior_EX_mu_mean_inter,prior_EX_mu_sd_inter Prior mean and
#'     sd for population mean parameters \eqn{\mu_{\eta k}} of each
#'     interaction parameter. Vector of length equal to the number of
#'     interactions.
#' @param prior_EX_tau_mean_inter,prior_EX_tau_sd_inter Prior mean and
#'     sd for heterogeniety parameter \eqn{\tau_{\eta s k}} of each
#'     stratum. Matrix with one column per interaction and one row per
#'     stratum.
#' @param prior_EX_corr_eta_inter Prior LKJ correlation parameter for
#'     interaction given as numeric. If missing, then a 1 is assumed
#'     corresponding to a marginal uniform prior of the correlations.
#' @param prior_EX_prob_comp Prior probability \eqn{p_{ij}} for
#'     exchangability of each component per group. Matrix with one
#'     column per component and one row per group. Values must lie in
#'     [0-1] range.
#' @param prior_EX_prob_inter Prior probability \eqn{p_{\eta k j}} for
#'     exchangability of each interaction per group. Matrix with one
#'     column per interaction and one row per group. Values must lie
#'     in [0-1] range.
#' @param prior_NEX_mu_mean_comp,prior_NEX_mu_sd_comp Prior mean
#'     \eqn{\boldsymbol m_{ij}} and sd \eqn{\boldsymbol s_{ij} =
#'     \mbox{diag}(\boldsymbol S_{ij})} of each component for
#'     non-exchangable case. Two column matrix (intercept, log-slope)
#'     with one row per component. If missing set to the same prior as
#'     given for the EX part. It is required that the specification be
#'     the same across groups j.
#' @param prior_NEX_mu_mean_inter,prior_NEX_mu_sd_inter Prior mean
#'     \eqn{m_{\eta k j}} and sd \eqn{s_{\eta k j}} for each
#'     interaction parameter for non-exchangable case. Vector of
#'     length equal to the number of interactions. If missing set to
#'     the same prior as given for the EX part.
#' @param prior_is_EXNEX_comp Defines if non-exchangability is
#'     admitted for a given component. Logical vector of length equal
#'     to the number of components. If missing \code{TRUE} is assumed
#'     for all components.
#' @param prior_is_EXNEX_inter Defines if non-exchangability is
#'     admitted for a given interaction parameter. Logical vector of
#'     length equal to the number of interactions. If missing
#'     \code{FALSE} is assumed for all interactions.
#' @param prior_tau_dist Defines the distribution used for
#'     heterogeniety parameters. Choices are 0=fixed to it's mean,
#'     1=log-normal, 2=truncated normal.
#' @param prior_PD Logical flag (defaults to \code{FALSE}) indicating
#'     if to sample the prior predictive distribution instead of
#'     conditioning on the data.
#' @template args-sampling
#' @param verbose Logical flag (defaults to \code{FALSE}) controlling
#'     if additional output like stan progress is reported.
#' @template args-dots-ignored
#' @param digits number of digits to show
#' @param x \code{blrmfit} object to print
#'
#' @details
#'
#' \code{blrm_exnex} is a flexible function for Bayesian meta-analytic modeling of binomial
#' count data. In particular, it is designed to model counts of the number of observed
#' dose limiting toxicities (DLTs) by dose, for guiding dose-escalation studies
#' in Oncology. To accommodate dose escalation over more than one agent, the dose
#' may consist of combinations of study drugs, with any number of treatment components.
#'
#' In the simplest case, the aim is to model the probability \eqn{\pi} that
#' a patient experiences a DLT, by complementing the binomial likelihood with
#' a monotone logistic regression
#'
#' \deqn{\mbox{logit}\,\pi(d) = \log\,\alpha + \beta \, t(d),}
#'
#' where \eqn{\beta > 0}. Most typically, \eqn{d} represents the dose, and \eqn{t(d)}
#' is an appropriate transformation, such as \eqn{t(d) = \log (d \big / d^*)}. A
#' joint prior on \eqn{\boldsymbol \theta = (\log\,\alpha, \log\,\beta)} completes
#' the model and ensures monotonicity \eqn{\beta > 0}.
#'
#' Many extensions are possible. The function supports general combination regimens, and also
#' provides framework for Bayesian meta-analysis of dose-toxicity data from
#' multiple historical and concurrent sources.
#'
#' For an example of a single-agent trial refer to \code{\link{example-single-agent}}.
#'
#' @section Combination of two treatments:
#'
#' For a combination of two treatment components, the basic modeling framework
#' is that the DLT rate \eqn{\pi(d_1,d_2)} is comprised of (1) a "no-interaction"
#' baseline model \eqn{\tilde \pi(d_1,d_2)} driven by the single-agent toxicity of
#' each component, and (2) optional interaction terms \eqn{\gamma(d_1,d_2)}
#' representing synergy or antagonism between the drugs. On the log-odds scale,
#'
#' \deqn{\mbox{logit} \,\pi(d_1,d_2) = \mbox{logit} \, \tilde \pi(d_1,d_2) + \eta \, \gamma(d_1,d_2). }
#'
#' The "no interaction" part \eqn{\tilde \pi(d_1,d_2)} represents the probability
#' of a DLT triggered by either treatment component acting \emph{independently}.
#' That is,
#' \deqn{ \tilde \pi(d_1,d_2) = 1- (1 - \pi_1(d_1))(1 - \pi_2(d_2)). }
#' In simple terms, P(no DLT for combination) = P(no DLT for drug 1) * P(no DLT from drug 2).
#' To complete this part, the treatment components can then be modeled with monotone
#' logistic regressions as before.
#'
#' \deqn{\mbox{logit} \, \pi_i(d_i) = \log\, \alpha_i + \beta_i \, t(d_i),}
#'
#' where \eqn{t(d_i)} is a monotone transformation of the doses, such as
#' \eqn{t(d_i) = \log (d_i \big / d_i^*)}.
#'
#' The inclusion of an interaction term \eqn{\gamma(d_1,d_2)} allows
#' DLT rates above or below the "no-interaction" rate. The magnitude of
#' the interaction term may also be made dependent on the doses (or other covariates)
#' through regression. As an example, one could let
#'
#' \deqn{\gamma(d_1, d_2) = \frac{d_1}{d_1^*} \frac{d_1}{d_2^*}.}
#'
#' The specific functional form is specified in the usual notation for
#' a design matrix. The interaction model must respect the constraint
#' that whenever any dose approaches zero, then the interaction
#' term must vanish as well. Therefore, the interaction model must not
#' include an intercept term which would violate this consistency
#' requirement. A dual combination example can be found in
#' \code{\link{example-combo2}}.
#'
#' @section General combinations:
#'
#' The model is extended to general combination treatments consisting
#' of \eqn{N} components by expressing the probability \eqn{\pi} on
#' the logit scale as
#'
#' \deqn{ \mbox{logit} \, \pi(d_1,\ldots,d_N) = \mbox{logit} \Bigl( 1 - \prod_{i = 1}^N ( 1 - \pi_i(d_i) ) \Bigr) + \sum_{k=1}^K \eta_k \, \gamma_k(d_1,\ldots,d_N), }
#'
#' Multiple drug-drug interactions among the \eqn{N} components are
#' now possible, and are represented through the \eqn{K} interaction
#' terms \eqn{\gamma_k}.
#'
#' Regression models can be again be specified for each \eqn{\pi_i}
#' and \eqn{\gamma_k}, such as
#'
#' \deqn{ \mbox{logit}\, \pi_i(d_i) = \log\, \alpha_i + \beta_i \, t(d_i) }
#'
#' Interactions for some subset \eqn{I(k) \subset \{1,\ldots,N \}} of
#' the treatment components can be modeled with regression as well,
#' for example on products of doses,
#'
#' \deqn{ \gamma_k(d_1,\ldots,d_N) = \prod_{i \in I(k)} \frac{d_i}{d_i^*}.}
#'
#' For example, \eqn{I(k) = \{1,2,3\}} results in the three-way
#' interaction term
#'
#' \deqn{ \frac{d_1}{d_1^*} \frac{d_2}{d_2^*} \frac{d_3}{d_3^*} }
#'
#' for drugs 1, 2, and 3.
#'
#' For a triple combination example please refer to
#' \code{\link{example-combo3}}.
#'
#' @section Meta-analytic framework:
#'
#' Information on the toxicity of a drug may be available from
#' multiple studies or sources. Furthermore, one may wish to stratify
#' observations within a single study (for example into groups of
#' patients corresponding to different geographic regions, or multiple
#' dosing \code{dose_info} corresponding to different schedules).
#'
#' \code{blrm_exnex} provides tools for robust Bayesian hierarchical
#' modeling to jointly model data from multiple sources. An additional
#' index \eqn{j=1, \ldots, J} on the parameters and observations
#' denotes the \eqn{J} groups. The resulting model allows the DLT rate
#' to differ across the groups. The general \eqn{N}-component model
#' becomes
#'
#' \deqn{ \mbox{logit} \, \pi_j(d_1,\ldots,d_N) = \mbox{logit} \Bigl( 1 - \prod_{i = 1}^N ( 1 - \pi_{ij}(d_i) ) \Bigr) + \sum_{k=1}^K \eta_{kj} \, \gamma_{k}(d_1,\ldots,d_N), }
#'
#' for groups \eqn{j = 1,\ldots,J}. The component toxicities
#' \eqn{\pi_{ij}} and interaction terms \eqn{\gamma_{k}} are
#' modelled, as before, through regression. For example,
#' \eqn{\pi_{ij}} could be a logistic regression on \eqn{t(d_i) =
#' \log(d_i/d_i^*)} with intercept and log-slope \eqn{\boldsymbol
#' \theta_{ij}}, and \eqn{\gamma_{k}} regressed with coefficient
#' \eqn{\eta_{kj}} on a product \eqn{\prod_{i\in I(k)} (d_i/d_i^*)}
#' for some subset \eqn{I(k)} of components.
#'
#' Thus, for \eqn{j=1,\ldots,J}, we now have group-specific parameters
#' \eqn{\boldsymbol\theta_{ij} = (\log\, \alpha_{ij}, \log\,
#' \beta_{ij})} and \eqn{\boldsymbol\nu_{j} = (\eta_{1j}, \ldots,
#' \eta_{Kj})} for each component \eqn{i=1,\ldots,N} and interaction
#' \eqn{k=1,\ldots,K}.
#'
#' The structure of the prior on
#' \eqn{(\boldsymbol\theta_{i1},\ldots,\boldsymbol\theta_{iJ})} and
#' \eqn{(\boldsymbol\nu_{1}, \ldots, \boldsymbol\nu_{J})} determines how much
#' information will be shared across groups \eqn{j}. Several modeling
#' choices are available in the function.
#'
#' \itemize{
#'
#' \item \emph{EX (Full exchangeability):} One can assume the
#' parameters are conditionally exchangeable given hyperparameters
#'
#' \deqn{\boldsymbol \theta_{ij} \sim \mbox{N}\bigl( \boldsymbol \mu_{\boldsymbol \theta i}, \boldsymbol \Sigma_{\boldsymbol \theta i} \bigr), }
#'
#' independently across groups \eqn{j = 1,\ldots, J} and treatment
#' components \eqn{i=1,\ldots,N}. The covariance matrix
#' \eqn{\boldsymbol \Sigma_{\boldsymbol \theta i}} captures the patterns of cross-group
#' heterogeneity, and is parametrized with standard deviations
#' \eqn{\boldsymbol \tau_{\boldsymbol\theta i} = (\tau_{\alpha i},
#' \tau_{\beta i})} and the correlation \eqn{\rho_i}. Similarly for the
#' interactions, the fully-exchangeable model is
#'
#' \deqn{\boldsymbol \nu_{j} \sim \mbox{N}\bigl( \boldsymbol \mu_{\boldsymbol \nu}, \boldsymbol \Sigma_{\boldsymbol \nu} \bigr)}
#'
#' for groups \eqn{j = 1,\ldots, J} and interactions
#' \eqn{k=1,\ldots,K}, and the prior on the covariance matrix
#' \eqn{\boldsymbol \Sigma_{\boldsymbol \nu}} captures the amount of
#' heterogeneity expected in the interaction terms a-priori. The
#' covariance is again parametrized with standard deviations
#' \eqn{(\tau_{\eta 1}, \ldots, \tau_{\eta K})} and its correlation matrix.
#'
#' \item \emph{Differential discounting:} For one or more of the groups \eqn{j=1,\ldots,J},
#' larger deviations of \eqn{\boldsymbol\theta_{ij}} may be expected from the mean
#' \eqn{\boldsymbol\mu_i}, or of the interactions \eqn{\eta_{kj}} from the mean
#' \eqn{\mu_{\eta,k}}. Such differential heterogeneity can be modeled by mapping the groups
#' \eqn{j = 1,\ldots,J} to \emph{strata} through \eqn{s_j \in \{1,\ldots,S\}},
#' and modifying the model specification to
#' \deqn{\boldsymbol \theta_{ij} \sim \mbox{N}\bigl( \boldsymbol \mu_{\boldsymbol \theta i}, \boldsymbol \Sigma_{\boldsymbol \theta ij} \bigr), }
#' where
#' \deqn{\boldsymbol \Sigma_{\boldsymbol \theta ij} = \left( \begin{array}{cc}
#' \tau^2_{\alpha s_j i} & \rho_i \tau_{\alpha s_j i} \tau_{\beta s_j i}\\
#' \rho_i \tau_{\alpha s_j i} \tau_{\beta s_j i} & \tau^2_{\beta s_j i}
#' \end{array} \right).}
#' For the interactions, the model becomes
#' \deqn{\boldsymbol \nu_{j} \sim \mbox{N}\bigl( \boldsymbol \mu_{\boldsymbol \nu}, \boldsymbol \Sigma_{\boldsymbol \nu j} \bigr),}
#' where the covariance matrix \eqn{\boldsymbol \Sigma_{\boldsymbol \nu j}} is modelled as stratum specific standard deviations \eqn{(\tau_{\eta 1 s_j}, \ldots, \tau_{\eta K s_j})} and a stratum independent correlation matrix.
#' Each stratum
#' \eqn{s=1,\ldots,S} then corresponds to its own set of standard deviations \eqn{\tau} leading to different discounting per stratum.
#' Independent priors are specified for the component parameters
#' \eqn{\tau_{\alpha s i}} and \eqn{\tau_{\beta s i}} and
#' for the interaction parameters \eqn{\tau_{\eta s k}} for each stratum
#' \eqn{s=1,\ldots,S}. Inference for strata \eqn{s} where the prior is centered
#' on larger values of the \eqn{\tau} parameters will exhibit less shrinkage
#' towards the the means, \eqn{\boldsymbol\mu_{\boldsymbol \theta i}} and \eqn{\boldsymbol \mu_{\boldsymbol \nu}} respectively.
#'
#'
#' \item \emph{EXNEX (Partial exchangeability):} Another mechansim for increasing
#' robustness is to introduce mixture priors for the group-specific parameters,
#' where one mixture component is shared across groups, and the other is group-specific.
#' The result, known as an EXchangeable-NonEXchangeable (EXNEX) type prior, has a form
#'
#' \deqn{\boldsymbol \theta_{ij} \sim p_{\boldsymbol \theta ij}\, \mbox{N}\bigl( \boldsymbol \mu_{\boldsymbol \theta i}, \boldsymbol \Sigma_{\boldsymbol \theta i} \bigr) +(1-p_{\boldsymbol \theta ij})\, \mbox{N}\bigl(\boldsymbol m_{\boldsymbol \theta ij}, \boldsymbol S_{\boldsymbol \theta ij}\bigr)}
#'
#' when applied to the treatment-component parameters, and
#'
#' \deqn{\boldsymbol \nu_{kj} \sim p_{\boldsymbol \nu_{kj}} \,\mbox{N}\bigl(\mu_{\boldsymbol \nu}, \boldsymbol \Sigma_{\boldsymbol \nu}\bigr)_k + (1-p_{\boldsymbol \nu_{kj}})\, \mbox{N}(m_{\boldsymbol \nu_{kj}}, s^2_{\boldsymbol \nu_{kj}})}
#'
#' when applied to the interaction parameters. The \emph{exchangeability weights}
#' \eqn{p_{\boldsymbol \theta ij}} and \eqn{p_{\boldsymbol \nu_{kj}}} are fixed constants in the interval \eqn{[0,1]}
#' that control the degree to which inference for group \eqn{j} is informed
#' by the exchangeable mixture components. Larger values for the
#' weights correspond to greater exchange of information, while smaller values
#' increase robustness in case of outlying observations in individual groups
#' \eqn{j}.
#'
#'
#' }
#'
#'
#' @return The function returns a S3 object of type
#'     \code{blrmfit}.
#'
#'
#' @template ref-mac
#' @template ref-exnex
#' @template ref-critical_aspects
#' @template ref-bayesindustry
#'
#' @template start-example
#' @examples
#' # fit an example model. See documentation for "combo3" example
#' example_model("combo3")
#'
#' # print a summary of the prior
#' prior_summary(blrmfit, digits = 3)
#'
#' # print a summary of the posterior (model parameters)
#' print(blrmfit)
#'
#' # summary of posterior for DLT rate by dose for observed covariate levels
#' summ <- summary(blrmfit, interval_prob = c(0, 0.16, 0.33, 1))
#' print(cbind(hist_combo3, summ))
#'
#' # summary of posterior for DLT rate by dose for new set of covariate levels
#' newdata <- expand.grid(
#'   stratum = "BID", group_id = "Combo",
#'   drug_A = 400, drug_B = 800, drug_C = c(320, 400, 600, 800),
#'   stringsAsFactors = FALSE
#' )
#' summ_pred <- summary(blrmfit, newdata = newdata, interval_prob = c(0, 0.16, 0.33, 1))
#' print(cbind(newdata, summ_pred))
#'
#' # update the model after observing additional data
#' newdata$num_patients <- rep(3, nrow(newdata))
#' newdata$num_toxicities <- c(0, 1, 2, 2)
#' library(dplyr)
#' blrmfit_new <- update(blrmfit,
#'                       data = rbind(hist_combo3, newdata) %>%
#'                                arrange(stratum, group_id))
#'
#' # updated posterior summary
#' summ_upd <- summary(blrmfit_new, newdata = newdata, interval_prob = c(0, 0.16, 0.33, 1))
#' print(cbind(newdata, summ_upd))
#' @template stop-example
#'
#' @export
blrm_exnex <- function(formula,
                       data,
                       prior_EX_mu_mean_comp,
                       prior_EX_mu_sd_comp,
                       prior_EX_tau_mean_comp,
                       prior_EX_tau_sd_comp,
                       prior_EX_corr_eta_comp,
                       prior_EX_mu_mean_inter,
                       prior_EX_mu_sd_inter,
                       prior_EX_tau_mean_inter,
                       prior_EX_tau_sd_inter,
                       prior_EX_corr_eta_inter,
                       prior_is_EXNEX_inter,
                       prior_is_EXNEX_comp,
                       prior_EX_prob_comp,
                       prior_EX_prob_inter,
                       prior_NEX_mu_mean_comp,
                       prior_NEX_mu_sd_comp,
                       prior_NEX_mu_mean_inter,
                       prior_NEX_mu_sd_inter,
                       prior_tau_dist,
                       iter=getOption("OncoBayes2.MC.iter" , 2000),
                       warmup=getOption("OncoBayes2.MC.warmup", 1000),
                       thin=getOption("OncoBayes2.MC.thin", 1),
                       init=getOption("OncoBayes2.MC.init", 0.5),
                       chains=getOption("OncoBayes2.MC.chains", 4),
                       cores=getOption("mc.cores", 1L),
                       control=getOption("OncoBayes2.MC.control", list()),
                       prior_PD=FALSE,
                       verbose=FALSE
                       )
{
    call <- match.call()

    if (missing(data))
        data <- environment(formula)

    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "na.action"), names(mf), 0)
    mf <- mf[c(1, m)]

    f <- Formula::Formula(formula)
    mf[[1]] <- as.name("model.frame")
    ##mf[[1]] <- quote(stats::model.frame)
    mf$formula <- f
    mf <- eval(mf, parent.frame())

    num_rhs_terms <- length(f)[2]
    idx_group_term <- num_rhs_terms
    has_inter <- TRUE
    idx_inter_term <- num_rhs_terms-1

    mt <- attr(mf, "terms")

    for (i in 1:num_rhs_terms) {
        tc <- terms(f, rhs=i)
        nt <- length(tc)-1
        if(nt == 2 && attr(tc, "intercept") == 1 )
            num_comp <- i
        else
            break
    }

    if (num_comp == 1) {
        has_inter <- FALSE
        idx_inter_term <- 0
    }

    ## we only support a single LHS
    assert_that(length(f)[1] == 1)
    ## check that we have an overall intercept for all components
    for (i in 1:num_comp)
        assert_that(attr(terms(f, rhs=i), "intercept") == 1, msg="Intercept must be present for all components.")

    ## the interaction model must not have an intercept
    if(has_inter)
        assert_that(attr(terms(f, rhs=idx_inter_term), "intercept") == 0, msg="No intercept must be present for the interaction model.")

    y <- model.response(mf)
    ## model response must be a two-column matrix (responder,
    ## non-responder)
    assert_matrix(y, ncols=2, any.missing=FALSE)
    nr <- array(y[,2])
    r <- array(y[,1])
    n <- r + nr

    num_obs <- length(r)

    group_index_term <- model.part(f, data = mf, rhs = idx_group_term)

    if(ncol(group_index_term) > 2)
        stop("Grouping factor must have at most two terms (study index and optionally a stratum).")

    if(ncol(group_index_term) == 0)
        stop("Grouping factor must have at least one term (study index).")

    if(ncol(group_index_term) == 2) {
        idx_group_index <- 2
        idx_strata_index <- 1
    } else {
        idx_group_index <- 1
        idx_strata_index <- NA
    }

    group_fct <- group_index_term[,idx_group_index]
    if(!is.factor(group_fct))
        group_fct <- factor(group_index)
    group_index <- as.integer(unclass(group_fct))
    num_groups <- nlevels(group_fct)
    assert_that(NROW(group_index) == num_obs)

    ## per group we must have an assignment to a tau stratum
    if(is.na(idx_strata_index)) {
        strata_fct <- rep(1, num_obs)
    } else {
        strata_fct <- group_index_term[,idx_strata_index]
    }
    if(!is.factor(strata_fct))
        strata_fct <- factor(strata_fct)
    strata_index <- array(as.integer(unclass(strata_fct)))
    num_strata <- nlevels(strata_fct)
    assert_that(NROW(strata_index) == num_obs)

    ## obtain group => stratum map ordered by group id
    group_strata <- data.frame(group_index=1:num_groups) %>%
        left_join(unique(data.frame(group_index=group_index, strata_index=strata_index)), by="group_index")
    ## since we allow groups to be specified through the factors,
    ## there can be groups without a stratum defined in case one level
    ## of the group is not present in the data. These need to be
    ## assigned to a default stratum which is stratum 1 at the moment.
    if(any(is.na(group_strata$strata_index))) {
        group_strata_undef <- which(is.na(group_strata$strata_index))
        message("Info: The group(s) ", paste(levels(group_fct)[group_strata_undef], collapse=", "), " have undefined strata. Assigning first stratum ", levels(strata_fct)[1], ".")
        group_strata$strata_index[is.na(group_strata$strata_index)]  <- 1
    }
    group_index <- array(group_index)

    ## setup design matrices which must have intercept and slope
    X_comp <- list()
    X_comp_cols <- list()
    for (i in 1:num_comp) {
        X_comp <- c(X_comp, list(model.matrix(f, mf, rhs=i)))
        X_comp_cols <- c(X_comp_cols, list(colnames( X_comp[[i]] )))
        assert_matrix(X_comp[[i]], ncols=2, any.missing=FALSE)
    }
    X_comp <- do.call(abind, c(X_comp, list(along=0)))

    if(has_inter)
        X_inter <- model.matrix(f, mf, rhs=idx_inter_term)
    else
        X_inter <- model.matrix(~0, mf)

    num_inter <- ncol(X_inter)
    has_inter <- num_inter > 0

    ## note: most of the consistency checks of the prior are left for
    ## Stan

    assert_matrix(prior_EX_mu_mean_comp, any.missing=FALSE, nrows=num_comp, ncols=2)
    assert_matrix(prior_EX_mu_sd_comp, any.missing=FALSE, nrows=num_comp, ncols=2)

    ## in case a single stratum is used, the user can input a matrix
    if(num_strata == 1) {
        if (is.matrix(prior_EX_tau_mean_comp))
            prior_EX_tau_mean_comp  <- array(prior_EX_tau_mean_comp, c(1, dim(prior_EX_tau_mean_comp)))
        if (is.matrix(prior_EX_tau_sd_comp))
            prior_EX_tau_sd_comp  <- array(prior_EX_tau_sd_comp, c(1, dim(prior_EX_tau_sd_comp)))
    }
    assert_array(prior_EX_tau_mean_comp, any.missing=FALSE, d=3)
    assert_array(prior_EX_tau_sd_comp, any.missing=FALSE, d=3)
    assert_that(all(dim(prior_EX_tau_mean_comp) == c(num_strata, num_comp, 2)),
                msg="prior_EX_tau_mean_comp must have dimensionality of num_strata x num_comp x 2.\nIn case of only one stratum a matrix of num_comp x 2 is sufficient.")
    assert_that(all(dim(prior_EX_tau_sd_comp) == c(num_strata, num_comp, 2)),
                msg="prior_EX_tau_sd_comp must have dimensionality of num_strata x num_comp x 2.\nIn case of only one stratum a matrix of num_comp x 2 is sufficient.")

    if(missing(prior_EX_corr_eta_comp))
        prior_EX_corr_eta_comp <- rep(1.0, times=num_comp)
    assert_numeric(prior_EX_corr_eta_comp, lower=0, finite=TRUE, any.missing=FALSE, len=num_comp)

    if(!has_inter & missing(prior_EX_mu_mean_inter)) {
        prior_EX_mu_mean_inter  <- array(0, dim=0)
    }
    if(!has_inter & missing(prior_EX_mu_sd_inter)) {
        prior_EX_mu_sd_inter  <- array(0, dim=0)
    }
    assert_numeric(prior_EX_mu_mean_inter, any.missing=FALSE, len=num_inter)
    assert_numeric(prior_EX_mu_sd_inter, any.missing=FALSE, len=num_inter, lower=0)

    if(!has_inter & missing(prior_EX_tau_mean_inter)) {
        prior_EX_tau_mean_inter <- matrix(1, nrow=num_strata, ncol=num_inter)
    }
    if(!has_inter & missing(prior_EX_tau_sd_inter)) {
        prior_EX_tau_sd_inter <- matrix(1, nrow=num_strata, ncol=num_inter)
    }
    assert_matrix(prior_EX_tau_mean_inter, any.missing=FALSE, nrows=num_strata, ncols=num_inter)
    assert_matrix(prior_EX_tau_sd_inter, any.missing=FALSE, nrows=num_strata, ncols=num_inter)

    if(!has_inter & missing(prior_EX_prob_inter)) {
        prior_EX_prob_inter <- matrix(1, nrow=num_groups, ncol=num_inter)
    }
    assert_matrix(prior_EX_prob_comp, any.missing=FALSE, nrows=num_groups, ncols=num_comp)
    assert_matrix(prior_EX_prob_inter, any.missing=FALSE, nrows=num_groups, ncols=num_inter)

    if(missing(prior_EX_corr_eta_inter))
        prior_EX_corr_eta_inter <- 1.0
    assert_number(prior_EX_corr_eta_inter, lower=0, finite=TRUE)

    if(missing(prior_NEX_mu_mean_comp))
        prior_NEX_mu_mean_comp <- prior_EX_mu_mean_comp
    if(missing(prior_NEX_mu_sd_comp))
        prior_NEX_mu_sd_comp <- prior_EX_mu_sd_comp

    assert_matrix(prior_NEX_mu_mean_comp, any.missing=FALSE, nrows=num_comp, ncols=2)
    assert_matrix(prior_NEX_mu_sd_comp, any.missing=FALSE, nrows=num_comp, ncols=2)

    if(missing(prior_NEX_mu_mean_inter))
        prior_NEX_mu_mean_inter <- prior_EX_mu_mean_inter
    if(missing(prior_NEX_mu_sd_inter))
        prior_NEX_mu_sd_inter <- prior_EX_mu_sd_inter

    assert_numeric(prior_NEX_mu_mean_inter, any.missing=FALSE, len=num_inter)
    assert_numeric(prior_NEX_mu_sd_inter, any.missing=FALSE, len=num_inter, lower=0)

    if(missing(prior_is_EXNEX_comp))
        prior_is_EXNEX_comp <- rep(TRUE, num_comp)
    if(missing(prior_is_EXNEX_inter))
        prior_is_EXNEX_inter <- rep(FALSE, num_inter)
    assert_logical(prior_is_EXNEX_comp, any.missing=FALSE, len=num_comp)
    assert_logical(prior_is_EXNEX_inter, any.missing=FALSE, len=num_inter)

    assert_number(prior_tau_dist, lower=0, upper=2)

    assert_logical(prior_PD, any.missing=FALSE, len=1)

    stan_data <- list(
        num_obs=num_obs,
        r=r,
        nr=nr,
        num_comp=num_comp,
        num_inter=num_inter,
        ## design matrices
        X_comp=X_comp,
        X_inter=X_inter,
        ## group and strata mapping
        num_groups=num_groups,
        num_strata=num_strata,
        group=group_index,
        stratum=strata_index,
        group_stratum_cid=array(group_strata$strata_index),
        ## priors
        prior_tau_dist = prior_tau_dist,
        prior_EX_prob_comp=prior_EX_prob_comp,
        prior_EX_prob_inter=prior_EX_prob_inter,
        prior_EX_mu_mean_comp=prior_EX_mu_mean_comp,
        prior_EX_mu_sd_comp=prior_EX_mu_sd_comp,
        prior_EX_tau_mean_comp=prior_EX_tau_mean_comp,
        prior_EX_tau_sd_comp=prior_EX_tau_sd_comp,
        prior_EX_corr_eta_comp=array(prior_EX_corr_eta_comp, num_comp),
        prior_EX_mu_mean_inter=array(prior_EX_mu_mean_inter, num_inter),
        prior_EX_mu_sd_inter=array(prior_EX_mu_sd_inter, num_inter),
        prior_EX_tau_mean_inter=prior_EX_tau_mean_inter,
        prior_EX_tau_sd_inter=prior_EX_tau_sd_inter,
        prior_EX_corr_eta_inter=prior_EX_corr_eta_inter,
        prior_NEX_mu_mean_comp=prior_NEX_mu_mean_comp,
        prior_NEX_mu_sd_comp=prior_NEX_mu_sd_comp,
        prior_NEX_mu_mean_inter=array(prior_NEX_mu_mean_inter, num_inter),
        prior_NEX_mu_sd_inter=array(prior_NEX_mu_sd_inter, num_inter),
        prior_is_EXNEX_comp=array(1*prior_is_EXNEX_comp, num_comp),
        prior_is_EXNEX_inter=array(1*prior_is_EXNEX_inter, num_inter),
        prior_PD=1*prior_PD
        )

    control_sampling <- modifyList(list(adapt_delta=0.99, stepsize=0.1), control)

    stan_msg <- capture.output(stanfit <- rstan::sampling(stanmodels$blrm_exnex,
                                                          data=stan_data,
                                                          warmup=warmup,
                                                          iter=iter,
                                                          chains=chains,
                                                          cores=cores,
                                                          thin=thin,
                                                          init=init,
                                                          control=control_sampling,
                                                          algorithm = "NUTS",
                                                          open_progress=FALSE,
                                                          save_warmup=TRUE
                                                          ))
    if(attributes(stanfit)$mode != 0)
        stop("Stan sampler did not run successfully!")

    ## only display Stan messages in verbose mode
    if(verbose) {
        cat(paste(c(stan_msg, ""), collapse="\n"))
    }

    ## label parameters of stanfit object
    labels <- list()
    labels$param_log_beta <- .make_label_factor(c("intercept", "log_slope"))
    labels$param_beta <- .make_label_factor(c("intercept", "slope"))
    labels$component <- .make_label_factor(.abbreviate_label(sapply(X_comp_cols, "[", 2)))
    stanfit <- .label_index(stanfit, "mu_log_beta", labels$component, labels$param_log_beta)
    stanfit <- .label_index(stanfit, "tau_log_beta", strata_fct, labels$component, labels$param_log_beta)
    stanfit <- .label_index(stanfit, "rho_log_beta", labels$component)
    stanfit <- .label_index(stanfit, "beta_group", group_fct, labels$component, labels$param_beta)
    stanfit <- .label_index(stanfit, "beta_EX_prob", group_fct, labels$component)
    stanfit <- .label_index(stanfit, "log_lik_group", group_fct)
    if(has_inter) {
        labels$param_eta <- .make_label_factor(.abbreviate_label(colnames(X_inter)))
        stanfit <- .label_index(stanfit, "eta_group", group_fct, labels$param_eta)
        stanfit <- .label_index(stanfit, "eta_EX_prob", group_fct, labels$param_eta)
        stanfit <- .label_index(stanfit, "mu_eta", labels$param_eta)
        stanfit <- .label_index(stanfit, "tau_eta", strata_fct, labels$param_eta)
        stanfit <- .label_index(stanfit, "Sigma_corr_eta", labels$param_eta, labels$param_eta)
    }

    out <- list(
        call = call,
        group_strata=group_strata,
        standata=stan_data,
        stanfit=stanfit,
        formula = f,
        model = mf,
        terms = mt,
        xlevels = .getXlevels(mt, mf),
        data = data,
        idx_group_term=idx_group_term,
        idx_inter_term=idx_inter_term,
        has_inter=has_inter,
        group_fct=group_fct,
        strata_fct=strata_fct,
        labels=labels
    )
    structure(out, class="blrmfit")
}



#'
#' @describeIn blrm_exnex print function.
#' @template args-prob
#' @method print blrmfit
#' @export
print.blrmfit <- function(x, ..., prob=0.95, digits=2) {
    cat("Bayesian Logistic Regression Model with EXchangeability-NonEXchangeability\n\n")
    cat("Number of observations:", x$standata$num_obs, "\n")
    cat("Number of groups      :", x$standata$num_groups, "\n")
    cat("Number of strata      :", x$standata$num_strata, "\n")
    cat("Number of components  :", x$standata$num_comp, "\n")
    cat("Number of interactions:", x$standata$num_inter, "\n")
    cat("EXNEX components      :", sum(x$standata$prior_is_EXNEX_comp), "\n")
    cat("EXNEX interactions    :", sum(x$standata$prior_is_EXNEX_inter), "\n")

    assert_number(prob, lower=0, upper=1, finite=TRUE)

    probs <- c(0.5-prob/2, 0.5, 0.5+prob/2)

    ## dummy definitions to silence R CMD check
    Stratum <- Group  <- total <- n_total <- NULL

    cat("\nObservations per group:\n")
    ds <- as.data.frame(table(x$group_fct))
    rownames(ds) <- match(ds$Var1, levels(x$group_fct))
    names(ds) <- c("Group", "n")
    totals  <- data.frame(Stratum=x$strata_fct, Group=x$group_fct, total=x$standata$nr+x$standata$r) %>%
        group_by(Stratum, Group) %>%
        summarise(n_total=sum(total)) %>%
        ungroup()
    ds  <- left_join(ds, totals, by="Group")
    ds$Stratum <- levels(x$strata_fct)[x$group_strata$strata_index]
    ds$n_total[is.na(ds$n_total)] <- 0
    print(ds)


    cat("\nGroups per stratum:\n")
    si  <- levels(x$strata_fct)[x$group_strata$strata_index]
    ds <- as.data.frame(table(si), stringsAsFactors = FALSE)
    names(ds) <- c("Stratum", "Groups")
    ds$Stratum <- factor(ds$Stratum, levels=levels(x$strata_fct))
    ds <- ds[order(ds$Stratum, ds$Groups), ]

    totals_stratum  <- totals  %>%
        group_by(Stratum) %>%
        summarise(n_total=sum(n_total))

    ds  <- left_join(ds, totals_stratum, by="Stratum")
    print(ds)

    comp_idx <- function(labels) {
        inter  <- grep("intercept\\]$", labels)
        slope  <- grep("slope\\]$", labels)
        list(inter=inter, slope=slope)
    }
    strip_variable <- function(labels) {
        gsub("^([A-Za-z_]+\\[)(.*)\\]$", "\\2", labels)
    }

    cat("\nComponent posterior:\n")
    cat("Population mean posterior mu_log_beta\n")
    mu_log_beta <- summary(x$stanfit, pars=c("mu_log_beta"), probs=probs)$summary
    rs <- rownames(mu_log_beta)
    idx <- comp_idx(rs)
    rownames(mu_log_beta) <- gsub("^(.*),intercept|,log_slope$", "\\1", strip_variable(rs))
    cat("intercept:\n")
    print(mu_log_beta[idx$inter,], digits=digits)
    cat("log-slope:\n")
    print(mu_log_beta[idx$slope,], digits=digits)

    cat("\nPopulation heterogeniety posterior tau_log_beta\n")
    tau_log_beta <- summary(x$stanfit, pars=c("tau_log_beta"), probs=probs)$summary
    rs <- rownames(tau_log_beta)
    idx <- comp_idx(rs)
    rownames(tau_log_beta) <- gsub("^(.*),intercept|,log_slope$", "\\1", strip_variable(rs))
    cat("intercept:\n")
    print(tau_log_beta[idx$inter,], digits=digits)
    cat("log-slope:\n")
    print(tau_log_beta[idx$slope,], digits=digits)

    cat("\nPopulation correlation posterior rho_log_beta\n")
    rho_log_beta <- summary(x$stanfit, pars=c("rho_log_beta"), probs=probs)$summary
    rownames(rho_log_beta) <- strip_variable(rownames(rho_log_beta))
    print(rho_log_beta, digits=digits)

    if(x$standata$num_inter > 0) {
        cat("\nInteraction model posterior:\n")
        cat("Population mean posterior mu_eta\n")
        mu_eta <- summary(x$stanfit, pars=c("mu_eta"), probs=probs)$summary
        rownames(mu_eta) <- strip_variable(rownames(mu_eta))
        print(mu_eta, digits=digits)

        cat("\nPopulation heterogeniety posterior tau_eta\n")
        tau_eta <- summary(x$stanfit, pars=c("tau_eta"), probs=probs)$summary
        rownames(tau_eta) <- strip_variable(rownames(tau_eta))
        print(tau_eta, digits=digits)

        cat("\nPopulation correlation posterior Sigma_corr_eta\n")
        ## TODO: do not display symmetric values
        Sigma_corr_eta <- summary(x$stanfit, pars=c("Sigma_corr_eta"), probs=probs)$summary
        rownames(Sigma_corr_eta) <- strip_variable(rownames(Sigma_corr_eta))
        print(Sigma_corr_eta, digits=digits)

    } else {
        cat("\nNo interaction model posterior specified.\n")
    }

    invisible(x)
}


## internal -----

#'
#' Utility function to label parameter indices according to factor
#' levels.
#' @param stanfit stan fit which names are being modified
#' @param par parameter selected
#' @param ... must include as many factors as there are indices which
#'     are used in the order given to translate indices to text labels
#'
#' @keywords internal
.label_index <- function(stanfit, par, ...) {
    idx <- grep(paste0("^", par, "\\["), names(stanfit))
    str <- names(stanfit)[idx]
    fct  <- list(...)
    idx_str <- t(sapply(strsplit(gsub("(.*)\\[([0-9,]*)\\]$", "\\2", str), ","), as.numeric))
    if (length(fct) == 1) {
        idx_str  <- matrix(idx_str, ncol=1)
    }
    ni  <- ncol(idx_str)
    colnames(idx_str) <- paste0("idx_", 1:ni)
    idx_str <- as.data.frame(idx_str)
    assert_that(ni == length(fct), msg="Insufficient number of indices specified")
    for(i in 1:ni) {
        f <- fct[[i]]
        key <- data.frame(idx=1:nlevels(f), label=levels(f))
        names(key) <- paste0(names(key), "_", i)
        idx_str <- left_join(idx_str, key, by=paste0("idx_", i))
    }
    labs  <- paste0("label_", 1:ni)
    names(stanfit)[idx] <- paste0(par, "[", do.call(paste, c(idx_str[labs], list(sep=","))), "]")
    stanfit
}

.abbreviate_label <- function(label) {
    minlength <- getOption("OncoBayes2.abbreviate.min", 0)
    if(minlength > 0)
        return(abbreviate(label, minlength=minlength))
    label
}

.make_label_factor <- function(labels) {
    factor(1:length(labels), levels=1:length(labels), labels=labels)
}

#' @method model.matrix blrmfit
#' @export
model.matrix.blrmfit <- function(object, ...) {
  return(model.matrix.default(object, object$data))
}


