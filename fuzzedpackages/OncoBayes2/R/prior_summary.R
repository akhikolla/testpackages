#' Summarise model prior
#'
#' @description
#' Extracts a summary of the prior in a structured data format.
#'
#' @param object \code{blrmfit} (\code{blrm_trial}) object as returned from \code{\link{blrm_exnex}} (\code{\link{blrm_trial}}) analysis
#' @param digits number of digits to show
#' @param ... ignored by the function
#'
#' @details The summary of the prior creates a structured
#'     representation of the specified prior from a
#'     \code{\link{blrm_exnex}} (\code{\link{blrm_trial}}) analysis.
#'
#' @return Returns an analysis specific list, which has it's own
#'     \code{print} function. The returned list contains arrays which
#'     represent the prior in a structured format.
#'
#' @template start-example
#' @examples
#' ## run combo2 analysis which defines blrmfit model object
#' example_model("combo2")
#'
#' prior_summary(blrmfit)
#'
#' prior_sum  <- prior_summary(blrmfit)
#' names(prior_sum)
#'
#' ## the entries of the prior list are labelled arrays
#' dimnames(prior_sum$EX_mu_log_beta)
#'
#' @template stop-example
#'
#' @method prior_summary blrmfit
#' @aliases prior_summary
#' @export
prior_summary.blrmfit <- function(object, digits = 2, ...) {
    standata  <- object$standata
    labels  <- object$labels

    x <- list()
    x$has_inter  <- object$has_inter

    x$is_EXNEX_comp  <- .label_array(standata$prior_is_EXNEX_comp, component=labels$component)

    x$EX_prob_comp  <- .label_array(standata$prior_EX_prob_comp, group=object$group_fct, component=labels$component)
    is_EX_comp_idx  <- which(x$is_EXNEX_comp == 0)
    x$EX_prob_comp[,is_EX_comp_idx] <- 1

    x$EX_mu_log_beta  <- .parse_mu_log_beta(standata$prior_EX_mu_mean_comp, standata$prior_EX_mu_sd_comp, labels)
    x$NEX_mu_log_beta  <- .parse_mu_log_beta(standata$prior_NEX_mu_mean_comp, standata$prior_NEX_mu_sd_comp, labels)

    EX_tau_log_beta_mean <- .label_array(standata$prior_EX_tau_mean_comp, stratum=object$strata_fct, component=labels$component, coefficient=labels$param_log_beta)
    EX_tau_log_beta_sd <- .label_array(standata$prior_EX_tau_sd_comp, stratum=object$strata_fct, component=labels$component, coefficient=labels$param_log_beta)
    x$EX_tau_log_beta <- abind(mean=EX_tau_log_beta_mean, sd=EX_tau_log_beta_sd, along=0)
    names(dimnames(x$EX_tau_log_beta)) <- c("prior", "stratum", "component", "coefficient")
    x$EX_tau_log_beta <- aperm(x$EX_tau_log_beta, c(4,2,3,1))

    x$EX_corr_eta_comp <- .label_array(standata$prior_EX_corr_eta_comp, component=labels$component)

    x$tau_dist <- standata$prior_tau_dist

    if(x$has_inter) {
        x$is_EXNEX_inter  <- .label_array(standata$prior_is_EXNEX_inter, interaction=labels$param_eta)
        x$EX_prob_inter  <- .label_array(standata$prior_EX_prob_inter, group=object$group_fct, interaction=labels$param_eta)
        is_EX_inter_idx  <- which(x$is_EXNEX_inter == 0)
        x$EX_prob_inter[,is_EX_inter_idx] <- 1

        x$EX_mu_eta  <- .parse_mu_eta(standata$prior_EX_mu_mean_inter, standata$prior_EX_mu_sd_inter, labels)
        x$NEX_mu_eta  <- .parse_mu_eta(standata$prior_NEX_mu_mean_inter, standata$prior_NEX_mu_sd_inter, labels)

        EX_tau_eta_mean <- .label_array(standata$prior_EX_tau_mean_inter, stratum=object$strata_fct, interaction=labels$param_eta)
        EX_tau_eta_sd <- .label_array(standata$prior_EX_tau_sd_inter, stratum=object$strata_fct, interaction=labels$param_eta)
        x$EX_tau_eta <- abind(mean=EX_tau_eta_mean, sd=EX_tau_eta_sd, along=0)
        names(dimnames(x$EX_tau_eta)) <- c("prior", "stratum", "interaction")
        x$EX_tau_eta <- aperm(x$EX_tau_eta, c(3,2,1))

        x$EX_corr_eta_inter <- array(standata$prior_EX_corr_eta_inter)
        dimnames(x$EX_corr_eta_inter) <- list("interaction")
    }

    x$num_strata <- standata$num_strata
    x$num_groups <- standata$num_groups

    structure(x, class = "prior_summary.blrmfit",
              model_name = deparse(substitute(object)),
              print_digits = digits)
}

#' @export
#' @method print prior_summary.blrmfit
print.prior_summary.blrmfit <- function(x, digits, ...) {
    cat("Bayesian Logistic Regression Model with EXchangeability-NonEXchangeability\n\n")

    if (missing(digits))
        digits <- attr(x, "print_digits")

    tau_str  <- if(x$tau_dist == 0)
                    "known"
                else if(x$tau_dist == 1)
                    "log-normal"
                else if(x$tau_dist == 2)
                    "half-normal"
                else
                    stop("Unkown tau prior distribution.")

    cat("Mixture configuration\n")
    cat("---------------------\n")
    cat("EXNEX components :", sum(x$is_EXNEX_comp), "\n")
    print(x$is_EXNEX_comp)
    cat("\n")
    if(x$has_inter) {
        cat("EXNEX interactions:", sum(x$is_EXNEX_inter), "\n")
        print(x$is_EXNEX_inter)
        cat("\n")
    }
    cat("Prior probability for exchangeability per group\n")
    print(x$EX_prob_comp)
    cat("\n")
    if(x$has_inter) {
        print(x$EX_prob_inter)
        cat("\n")
    }
    cat("EXchangable hyperparameter priors\n")
    cat("---------------------------------\n")
    cat("Component parameters\n")
    cat("Mean mu_log_beta\n")
    print(ftable(x$EX_mu_log_beta, row.vars=2), digits=digits)
    cat("\n")
    cat(paste0("Heterogeneity tau_log_beta (", tau_str, ")\n"))
    if(x$num_strata > 1) {
        print(ftable(x$EX_tau_log_beta, row.vars=c(2,3)), digits=digits)
    } else {
        print(ftable(x$EX_tau_log_beta, row.vars=c(3)), digits=digits)
    }
    cat("\nCorrelation LKJ\n")
    print(x$EX_corr_eta_comp, digits=digits)

    if(x$has_inter) {
        cat("\nInteraction parameters\n")
        cat("Mean mu_eta\n")
        print(ftable(x$EX_mu_eta, row.vars=1), digits=digits)
        cat("\n")
        cat(paste0("Heterogeneity tau_eta (", tau_str, ")\n"))
        if(x$num_strata > 1) {
            print(ftable(x$EX_tau_eta, row.vars=c(2,1)), digits=digits)
        } else {
            print(ftable(x$EX_tau_eta, row.vars=c(1)), digits=digits)
        }
        cat("\nCorrelation LKJ\n")
        print(x$EX_corr_eta_inter, digits=digits)
    } else {
        cat("\nModel has no interaction parameters.\n")
    }

    cat("\n")
    cat("NonEXchangable priors\n")
    cat("---------------------\n")
    cat("Component parameters\n")
    cat("Mean mu_log_beta\n")
    print(ftable(x$NEX_mu_log_beta, row.vars=2), digits=digits)

    if(x$has_inter) {
        cat("\nInteraction parameters\n")
        cat("Mean mu_eta\n")
        print(ftable(x$NEX_mu_eta, row.vars=1), digits=digits)
    } else {
        cat("\nModel has no interaction parameters.\n")
    }

    invisible(x)
}

#' @method prior_summary blrm_trial
#' @export
prior_summary.blrm_trial <- function(object, ...)
{
    .assert_is_blrm_trial_and_prior_is_set(object)

    x <- list()
    x$prior_summary.blrmfit <- prior_summary(object$blrmfit, ...)

    structure(x, class = "prior_summary.blrm_trial")
}

#' @method print prior_summary.blrm_trial
#' @export
print.prior_summary.blrm_trial <- function(x, ...) {
    print(x$prior_summary.blrmfit, ...)
}

## internal -----

.label_array <- function(data, ...) {
    labs  <- list(...)
    assert_that(length(labs) == length(dim(data)), msg="Number of labels must match dimension of input data.")
    assert_that(all(dim(data) == sapply(labs, nlevels)), msg="Number of factor levels must match array dimensionality.")
    dimnames(data) <- lapply(labs, levels)
    data
}

.parse_mu_log_beta  <- function(mu_mean_comp, mu_sd_comp, labels) {
    mu_log_beta_mean <- .label_array(mu_mean_comp, component=labels$component, coefficient=labels$param_log_beta)
    mu_log_beta_sd <- .label_array(mu_sd_comp, component=labels$component, coefficient=labels$param_log_beta)
    mu_log_beta <- abind(mean=mu_log_beta_mean, sd=mu_log_beta_sd, along=0)
    names(dimnames(mu_log_beta)) <- c("prior", "component", "coefficient")
    mu_log_beta <- aperm(mu_log_beta, c(3,2,1))
    mu_log_beta
}

.parse_mu_eta  <- function(mu_mean_inter, mu_sd_inter, labels) {
    mu_eta_mean <- .label_array(mu_mean_inter, interaction=labels$param_eta)
    mu_eta_sd <- .label_array(mu_sd_inter, interaction=labels$param_eta)
    mu_eta <- abind(mean=mu_eta_mean, sd=mu_eta_sd, along=0)
    mu_eta <- t(mu_eta)
    names(dimnames(mu_eta)) <- c("interaction", "prior")
    mu_eta
}
