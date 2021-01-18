#' Graph Q Matrix
#'
#' Provides a heatmap approach to showing the estimated binary or averaged
#' values of the Q Matrix.
#'
#' @param x   Either an `edina`, `auto_edina`, or `q_matrix` object.
#' @param ... Additional parameters not used
#'
#' @return
#' A `ggplot2` object with a heatmap overview of the estimated Q matrix.
#'
#' @export
#' @rdname q_graph
#' @examples
#' q = q_matrix(matrix(c(1, 0, 1, 1, 0, 1), ncol = 3))
#' q_graph(q)
q_graph = function(x, ...) {
    UseMethod("q_graph")
}

#' @inheritParams best_model
#' @export
#' @rdname q_graph
q_graph.auto_edina = function(x, binary = TRUE,
                              ic = c("ppp", "bic", "dic"), ... ) {
    q_graph(best_model(x, ic), binary)
}

#' @param binary   Boolean to indicate if a classified Q (dichotomous by decision rule)
#'                 or an estimate Q (non-dichotomous) or should be shown.
#'                 Default: `TRUE`.
#' @export
#' @rdname q_graph
q_graph.edina = function(x, binary = TRUE, ... ){

    if(binary == TRUE) {
        p = x$est_q
    } else {
        p = x$avg_q
    }

    q_type = if (binary == TRUE) {
        "Binary"
    } else {
        "Average"
    }

    title = paste("Estimated ", q_type, "Q Matrix")

    q_heatmap(p, title)
}

#' @export
#' @rdname q_graph
q_graph.matrix = function(x, ... ){
    q_heatmap(x, title = "Q Matrix")
}

#' @export
#' @rdname q_graph
q_graph.q_matrix = function(x, ... ){
    q_heatmap(x, title = "Q Matrix")
}


#' @importFrom reshape2 melt
#' @importFrom ggplot2 theme_minimal scale_fill_gradient geom_tile
q_heatmap = function(x, title = "Estimated Q Matrix") {

    Trait = Item = Value = NULL

    dgrid = reshape2::melt(x)
    colnames(dgrid) = c("Item", "Trait", "Value")

    ggplot(dgrid, aes(Trait, Item)) +
        geom_tile(aes(fill = Value), color = "white") +
        scale_fill_gradient(low = "white", high = "steelblue") +
        labs(title = title,
             subtitle = paste0("J = ", nrow(x), ", K = ", ncol(x)),
             x = "Items",
             y = "Traits",
             fill = "Estimated Value") +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5),
              plot.subtitle = element_text(hjust = 0.5))
}
