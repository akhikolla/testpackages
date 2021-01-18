#' Graph the Auto EDINA Object
#'
#' Presents either the fitting of model heuristics or the evolution of parameters
#' on a graph
#'
#' @param object An `auto_edina` object.
#' @param type   Kind of graph to display. Valid types: `"selection"` or `"evolution"`.
#' @param ... Not used.
#'
#' @return
#' A `ggplot2` object.
#'
#' @seealso
#' [auto_edina()],
#' [best_model()],
#' [model_selection_graph()],
#' [parameter_evolution_graph()]
#'
#' @export
#' @importFrom ggplot2 autoplot ggplot geom_line geom_point geom_vline facet_wrap labs aes theme_bw theme element_text
#' @examples
#' if(requireNamespace("simcdm", quietly = TRUE)) {
#'
#' # Set a seed for reproducibility
#' set.seed(1512)
#'
#' # Setup data simulation parameters
#' N = 2    # Number of Examinees / Subjects
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
#' # Visualize results results
#' autoplot(model_set_edina, type = "selection")
#'
#' # Equivalent to:
#' model_selection_graph(model_set_edina)
#'
#' # View model parameters
#' autoplot(model_set_edina, type = "guessing")
#'
#' # Or directly call with:
#' parameter_evolution_graph(model_set_edina, type = "guessing")
#' }
#' }
autoplot.auto_edina = function(object,
                               type = c("selection", "guessing", "slipping", "evolution"),
                               ...) {

    type = tolower(type)
    type = match.arg(type)

    switch(type,
           "selection"  = model_selection_graph(object, ...),
           "guessing"   = parameter_evolution_graph(object, type = type, ...),
           "slipping"   =  parameter_evolution_graph(object, type = type, ...),
           "evolution"  =  parameter_evolution_graph(object, type = type, ...),
           stop('Only the following types are valid: `"selection"`, `"guessing"`, or `"slipping"`')
    )

}

#' View Model Selection Statistics Across Models
#'
#' Displays information about the value of each model information criterion
#' for a given model across the dimensions the Q matrix is estimated.
#'
#' @param x   An `auto_edina` or `auto_errum` object.
#' @param ... Not used
#'
#' @return
#'
#' A `ggplot2` object
#'
#' @seealso
#' [autoplot.auto_edina()]
#'
#' @export
#' @importFrom stats reshape
model_selection_graph = function(x, ...){
    UseMethod("model_selection_graph", x)
}

#' @export
model_selection_graph.auto_edina = function(x, ...) {

    K = ic_value = NULL

    colnames(x$criterions) = toupper(colnames(x$criterions))
    ic_type_names = colnames(x$criterions)[-1]

    df = reshape(as.data.frame(x$criterions),
                 varying   = list(ic_type_names),
                 direction = "long",
                 idvar     = "k",
                 v.names   = "ic_value",
                 timevar   = "ic_type",
                 times     = ic_type_names)

    subset_df = do.call(rbind, by(df, df$ic_type, function(x) x[which.min(x$ic_value), ] ))

    ggplot(df, aes(x = as.factor(K), y = ic_value)) +
        facet_wrap(~ic_type, scales = "free_y") +
        geom_line(aes(group = 1)) +
        geom_point() +
        geom_point(data = subset_df,
                   colour="red", size = 3) +
        labs(title = "Auto EDINA Model Selection",
             y     = "Information Criterion Score",
             x     = "K Dimension of Q Matrix") +
        theme_bw()
}

#' @export
model_selection_graph.default = function(x, ...){
    stop_bad_class(x, "auto_edina")
}

#' View Slipping and Guessing Parameter Changes Across Models
#'
#' Displays the slipping and guessing parameter changes for each model across
#' the dimensions the Q matrix is estimated.
#'
#' @param x   An `auto_edina` or `auto_errum` object.
#' @param ... Not used
#'
#' @return
#'
#' A `ggplot2` object
#'
#' @seealso
#' [autoplot.auto_edina()]
#'
#' @export
parameter_evolution_graph = function(x, ...) {
    UseMethod("parameter_evolution_graph", x)
}

#' @export
parameter_evolution_graph.auto_edina = function(x,
                                                type = c("evolution", "guessing", "slipping"), ...) {

    type = match.arg(type)

    # Globals to quiet CRAN check
    K = ic_value = estimate = param_type = NULL

    J = x$j
    nmodels = length(x$edina_models)

    # Get the length of the string e.g. 200 => 3
    nlen = nchar(J)

    # Potentially add pis class? unlist(m_pi))

    extract_estimates = do.call(rbind, lapply(x$edina_models, `[[`, 1))

    o = data.frame(K          = rep(rep(x$k_checked, each = J), 2),
                   param_name = c(rep(
                       rep(sprintf(paste0("Item %0", nlen, "d"), seq_len(J)),
                           nmodels), 2)
                   ),
                   param_type = c(rep("Guessing", J*nmodels),
                                  rep("Slipping", J*nmodels)
                   ),
                   estimate   = c(extract_estimates[,"Guessing"], extract_estimates[,"Slipping"])

    )

    if(type %in% c("guessing", "slipping")) {
        o = o[grepl(type, o$param_type, ignore.case = TRUE),]
    }

    ggplot(o, aes(x = K, y = estimate, color = param_type)) +
        geom_point() + geom_line() +
        facet_wrap(~param_name) +
        labs(
            title = paste0("Evolution of the DINA Parameters"),
            subtitle = "Over Changes in Q Matrix's K Dimension",
            y = "Parameter Estimate of the DINA Parameters",
            x = "Q Matrix of a given K Dimension",
            color = "Parameter Type"
        ) +
        theme_bw() +
        theme(axis.text.x = element_text(hjust = 1, angle = 75))
}

#' @export
parameter_evolution_graph.default = function(x, ...){
    stop_bad_class(x, "auto_edina")
}
