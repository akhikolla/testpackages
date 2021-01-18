#' Construct an EDINA object
#'
#' Constructor function to satisfy the EDINA object definition
#'
#' @param coefs          Matrix with means of guessing and slipping coefficients.
#' @param loglike_summed Log-likelihood over the iterations.
#' @param loglike_pmean  Log-likelihood at the mean point
#' @param pis            Estimated latent classes
#' @param avg_q          Averaged Q Matrix over Iterations
#' @param est_q          Estimated Q Matrix
#' @param or_tested      Sample empirical odds ratio compared against simulated
#'                       odds ratio divided by the number of simulations.
#' @param sample_or      Sample empirical odds ratio based on the trial matrix.
#' @inheritParams edina
#' @param n            Number of Observations
#' @param timing       Number of Seconds that have elapsed since run time.
#' @param dataset_name Name of the data set the estimation procedure ran on.
#'
#' @noRd
new_edina = function(coefs, loglike_summed, loglike_pmean,
                     pis, avg_q, est_q, or_tested, sample_or, n, burnin,
                     chain_length,  timing, dataset_name) {

    colnames(coefs) = c("Guessing", "SD(Guessing)", "Slipping", "SD(Slipping)")
    rownames(coefs) = paste0("Item", seq_len(nrow(coefs)) )

    est_q = format_q_matrix(est_q)
    avg_q = format_q_matrix(avg_q)

    structure(list("coefficients"   = coefs,
                   "loglike_summed" = loglike_summed,
                   "loglike_pmean"  = loglike_pmean,
                   "pi_classes"     = pis,
                   "avg_q"          = avg_q,
                   "est_q"          = est_q,
                   "or_tested"      = or_tested,
                   "sample_or"      = sample_or,
                   "n"              = n,
                   "j"              = nrow(est_q),
                   "k"              = ncol(est_q),
                   "burnin"         = burnin,
                   "chain_length"   = chain_length,
                   "timing"         = timing,
                   "dataset_name"   = dataset_name),
              class = "edina")
}


#' Construct a Summary EDINA object
#'
#' Constructor function to satisfy the Summary EDINA object definition
#'
#' @param edina      An `edina` object
#' @param model_fit  Computed model heuristic value
#' @param alpha      The region used in the computation of the heuristic.
#'
#' @noRd
new_edina_summary = function(edina, model_fit, alpha) {

    edina[["model_fit"]] = model_fit
    edina[["alpha"]] = alpha

    class(edina) = c("summary_edina", "edina")

    edina
}


#' EDINA Estimation Routine
#'
#' Performs the Exploratory Deterministic Input, Noise and Gate Model (EDINA)
#' estimation on a given data set with a prespecified `k` value.
#'
#' @param data         Binary responses to assessments in `matrix`
#'                     form with dimensions \eqn{N \times J}{N x J}.
#' @param k            Number of Attribute Levels as a positive `integer`.
#' @param burnin       Number of Observations to discard on the chain.
#' @param chain_length Length of the MCMC chain
#'
#' @return
#' An `edina` object that contains:
#'
#' - `coefficients`: Estimated coefficients of the model fit
#' - `loglike_summed`: Summed log-likelihood
#' - `loglike_pmean`: Mean of log-likelihood
#' - `pi_classes`: Latent classes
#' - `avg_q`: Estimated Averaged Q Matrix
#' - `est_q`: Estimated Dichotomous Q Matrix
#' - `or_tested`: Odds Ratio used in the model selection.
#' - `sample_or`: Odds Ratio for the sample.
#' - `n`: Number of Observations
#' - `j`: Number of Items
#' - `k`: Number of Traits
#' - `burnin`: Amount of iterations to discard
#' - `chain_length`: Amount of iterations to retain.
#' - `timing`: Duration of the run
#' - `dataset_name`: Name of the data set used in estimation.
#'
#' @seealso
#' [auto_edina()],
#' [summary.edina()],
#' [print.edina()]
#'
#' @export
#' @importFrom jjb is_whole
#' @examples
#' if(requireNamespace("simcdm", quietly = TRUE)) {
#'
#' # Set a seed for reproducibility
#' set.seed(1512)
#'
#' # Setup data simulation parameters
#' N = 1    # Number of Examinees / Subjects
#' J = 10   # Number of Items
#' K = 2    # Number of Skills / Attributes
#'
#' # Note:
#' # Sample size and attributes have been reduced to create a minimally
#' # viable example that can be run during CRAN's automatic check.
#' # Please make sure to have a larger sample size...
#'
#' # Assign slipping and guessing values for each item
#' ss = gs = rep(.2, J)
#'
#' # Simulate an identifiable Q matrix
#' Q = simcdm::sim_q_matrix(J, K)
#'
#' # Simulate subject attributes
#' subject_alphas = simcdm::sim_subject_attributes(N, K)
#'
#' # Simulate items under the DINA model
#' items_dina = simcdm::sim_dina_items(subject_alphas, Q, ss, gs)
#'
#' # Compute the edina model
#' edina_model = edina(items_dina, k = K)
#'
#' # Display results
#' edina_model
#'
#' # Provide a summary overview
#' summary(edina_model)
#' }
#'
edina = function(data, k = 3, burnin = 10000, chain_length = 20000){

    stopifnot(is.matrix(data))

    stopifnot(is_whole(k) && length(k) == 1 && k >= 1)

    stopifnot(is_whole(chain_length) && length(chain_length) == 1)

    time_info = system.time({

        edina_model = edina_Gibbs_Q(data, k,
                                    burnin = burnin,
                                    chain_length = chain_length)

    })[1:3]

    dataset_name = deparse(substitute(data))

    new_edina(edina_model$coefficients,
              edina_model$loglike_summed,
              edina_model$loglike_pmean,
              edina_model$pis,
              edina_model$avg_q,
              edina_model$est_q,
              edina_model$or_tested,
              edina_model$sample_or,
              nrow(data),
              burnin, chain_length,
              time_info,
              dataset_name
    )
}

#' Printing out the EDINA Object
#'
#' Custom print method for computing the EDINA.
#'
#' @param x        An `edina` object
#' @param binary   Boolean to indicate whether the _Q_ matrix is shown in
#'                 dichotomous form or in an estimated form.
#' @param ...      Additional methods passed onto the `print.matrix` method.
#'
#' @return
#' None.
#'
#' The function provides a side-effect of displaying the overview of
#' the model estimated.
#'
#' @export
print.edina = function(x, binary = FALSE, ...){
    cat("EDINA model for", x$dataset_name, "with K =", x$k, "\n\n")

    est_mat = cbind(extract_q_matrix(x, binary = binary), x$coefficients)

    print(est_mat, digits = 4, ...)
}

#' Summarize the EDINA Object
#'
#' Provide a more detailed view inside of `edina` model object.
#'
#' @param object An `edina` object
#' @param alpha  Defining region to indicate the level of extremeness
#'               the data must before the model is problematic.
#' @param ...    Not used.
#'
#' @return
#' A summary object that includes everything in the original [edina()] object
#' and:
#'
#' - `model_fit`: Matrix of model fit summary statistics.
#' - `alpha`: Alpha-value used to compute [PPP()]s.
#'
#' @export
summary.edina = function(object, alpha = 0.05, ...) {

    model_fit = matrix(c(object$k, BIC(object), DIC(object),
                         PPP(object, alpha)), ncol = 4)

    colnames(model_fit) = c("k", "bic", "dic", "heuristic")

    new_edina_summary(
        object,
        model_fit = model_fit,
        alpha = alpha
    )

}

#' Printing out the Summary EDINA Object
#'
#' Custom print method for displaying the EDINA model summary information.
#'
#' @param x        A `summary_edina` object
#' @param binary   Boolean to indicate whether the _Q_ matrix is shown in
#'                 dichotomous form or in an estimated form.
#' @param ...      Past onto the `print.data.frame` method.
#'
#' @return
#' None.
#'
#' The function provides a side-effect of displaying the overview of
#' the model estimated.
#'
#' @export
print.summary_edina = function(x, binary = FALSE,  ...) {
    # Rely upon the specification of the `edina` object in the summary class.
    # NextMethod()

    cat("The EDINA model for", x$dataset_name, "with K =", x$k, "\n")

    cat("\nThe model fit is as follows:\n")
    print(as.data.frame(x$model_fit), row.names = FALSE)

    cat("\nThe estimated coefficients for the EDINA model are:\n")
    print(x$coefficients, digits = 4)

    cat("\nThe estimated Q matrix is:\n")
    print(extract_q_matrix(x, binary = binary), digits = 4)
}
