#' Auto EDINA model selection routine
#'
#' Automatically select an appropriate \eqn{K} dimension for a \eqn{Q} matrix
#' under the Exploratory Deterministic Input, Noise And gate (EDINA) Model.
#'
#' @param data          Binary responses to assessments in `matrix`
#'                      form with dimensions \eqn{N \times J}{N x J}.
#' @param k             Number of Attribute Levels as a positive `integer`.
#' @param burnin        Number of Observations to discard on the chain.
#' @param chain_length  Length of the MCMC chain
#'
#' @return
#' An `auto_edina` object that contains:
#'
#' - `edina_models`: A list containing all estimated `edina` model objects.
#' - `criterions`: Information criterions calculated for each model
#' - `k_checked`: Varying `k` dimensions checked.
#' - `j`: Number of Items
#'
#' @seealso
#' [autoplot.auto_edina()],
#' [best_model()],
#' [model_selection_graph()],
#' [parameter_evolution_graph()]
#'
#' @export
#' @examples
#' if(requireNamespace("simcdm", quietly = TRUE)) {
#'
#' # Set a seed for reproducibility
#' set.seed(1512)
#'
#' # Setup data simulation parameters
#' N = 15   # Number of Examinees / Subjects
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
#' \donttest{
#' # Requires at least 15 seconds of execution time.
#' # Three EDINA models will be fit with increasing number of attributes.
#' model_set_edina = auto_edina(items_dina, k = 2:4)
#'
#' # Display results
#' model_set_edina
#'
#' # Retrieve criterion table
#' table = summary(model_set_edina)
#'
#' # Extract "best model"
#' best_model(model_set_edina)
#' }
#' }
#'
auto_edina = function(data, k = 2:4,
                      burnin = 10000, chain_length = 20000) {

    ## Note:
    # Chain length is adjusted in edina

    # Compute the number of _K_ to estimate
    num_k = length(k)
    num_j = ncol(data)

    message("Starting the estimation procedure ... ")

    # Setup storage for EDINA Object
    outobj = vector('list', num_k)

    criterions = matrix(NA, nrow = num_k, ncol = 4)
    colnames(criterions) = c("k", "bic", "dic", "ppp")

    for(i in seq_along(k)) {
        k_idx = k[i]
        message("Working on k = ", k_idx, " ... ")

        modeled_value = edina(data,
                              k = k_idx,
                              burnin = burnin,
                              chain_length = chain_length)

        modeled_value_summary = summary(modeled_value)

        outobj[[i]] = modeled_value_summary

        criterions[i,] = outobj[[i]][["model_fit"]]

        message("Time Elapsed: ",  outobj[[i]][["timing"]][3])

    }

    # Output all EDINA objects
    structure(list("edina_models" = outobj,
                   "criterions" = criterions,
                   "k_checked" = k,
                   "j" = num_j
    )
    , class = "auto_edina" )
}


#' Print method for `auto_edina`
#'
#' Custom print method for displaying the results of the Auto EDINA method.
#'
#' @param x   An `auto_edina` object
#' @param ... Additional values passed onto the `print.data.frame` method.
#'
#' @return
#' None.
#'
#' The function provides a side-effect of displaying the overview of
#' computed results across all models estimated.
#'
#' @export
print.auto_edina = function(x, ...) {
    cat("The results of searching Q-matrices between", min(x$k_checked),
        "and", max(x$k_checked), "...\n")
    print(as.data.frame(x$criterions), digits = 4, row.names = FALSE, ...)
}

#' Summarize `auto_edina` model data
#'
#' Custom method for displaying the results of the `auto_edina`.
#'
#' @param object An `auto_edina` object
#' @param ...    Not used.
#'
#' @return
#' The original `auto_edina` object with an added class of `summary.auto_edina`.
#'
#' @export
summary.auto_edina = function(object, ...) {
   class(object) = c('summary_auto_edina', class(object))
   object
}

#' Print the `auto_edina` model summary
#'
#' Custom method for displaying the results of the `summary(auto_edina)`.
#'
#' @param x   A `summay_auto_edina` object
#' @param ... Additional values passed onto the `print.data.frame` method.
#'
#' @return
#' None.
#'
#' The function provides a side-effect of displaying the overview of
#' computed results across all models estimated.
#'
#' @export
print.summary_auto_edina = function(x, ...) {
    NextMethod()
}

