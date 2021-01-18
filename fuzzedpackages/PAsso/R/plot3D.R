#' @title 3-D P-P plot for the inspection of the partial association analysis
#'
#' @description A list of 3-D P-P plots for the inspection of the partial
#' association analysis. Each plot is the 3-D P-P plot from an empirical
#' copula trained from the surrogate residuals of a pair of responses.
#'
#' @param object The input object should be a "PAsso" class that is generated
#' by "PAsso" or "test".
#' @param y1 A string to specify the first response for the 3D plot.
#' @param y2 A string to specify the second response for the 3D plot. If either one of the
#' y1 or y2 is missing. The \code{plot3D} will draw 3D plots for all pairs of responeses.
#' @param ... Additional optional arguments.
#'
#' @details All the plots are based on surrogate residuals generated from \code{"resides"}
#' function in \code{sure}. Graphics are designed based on
#' PAsso and \code{"plotly"}.
#'
#' @return If response y1 or y2 is not specified, a list of \code{"plotly"} objects includes
#' all pairs of responses will be returned (with name "response 1 v.s. response 2" etc.). If
#' responses y1 and y2 are specified, returns a 3D plot as \code{"plotly"} object.
#'
#' @rdname plot3D
#' @importFrom copula C.n pobs
#' @importFrom plotly plotly_build plot_ly add_surface layout
#' @importFrom dplyr %>%
#' @export plot3D
#'
#' @examples
#' # Did not run this to save time
#' # data("ANES2016")
#' # PAsso_3v <- PAsso(responses = c("PreVote.num", "PID", "selfLR"),
#' #                   adjustments = c("income.num", "age", "edu.year"),
#' #                   data = ANES2016)
#'
#' # plot3D(PAsso_3v, y1="PID", y2="selfLR")
#'
plot3D <- function(object, y1, y2, ...) {
  UseMethod("plot3D")
}

#' @rdname plot3D
#' @method plot3D default
#' @export
plot3D.default <- function(object,y1, y2, ...){
  warning(paste("plot3D does not know how to handle object of class ",
                class(object),
                "and can only be used on classes PAsso"))
}

#' @keywords internal
plot3D_one <- function(plot_list, rep_SRs, m, n, plot_titles) {
  resi <- cbind(rep_SRs[,1,m], rep_SRs[,1,n])
  empC <- C.n(pobs(resi), resi)
  empFG1 <- pobs(resi)[,1]*pobs(resi)[,2]

  ## 3-D copula plot
  v1 <- v2 <- seq(0, 1, length.out = 100)
  aa <- matrix(0, length(v1), length(v1))
  for(i in 1:length(v1)){
    for(j in 1:length(v1)){
      aa[i,j] <- C.n(t(as.matrix(c(v1[i], v2[j]))), resi)-v1[i]*v2[j]
    }
  }

  # options(Viewer=NULL)
  # name <- paste("", i_plot, sep = "_")
  plot_list[[plot_titles]] <-
    plotly_build(plot_ly(x = v1, y = v2, z = 12*aa) %>%
                   add_surface() %>%
                   layout(scene = list(xaxis= list(title= "u"),
                                       yaxis= list(title= "v"),
                                       zaxis= list(title= "12(C(u,v)-uv)")),
                          title = paste("3-D P-P Plot: ", plot_titles, sep = "")))
  return(plot_list)
}

#' @param object A PAsso class of object.
#'
#' @param y1 A string to specify the first response for the 3D plot.
#' @param y2 A string to specify the second response for the 3D plot. If either one of the
#' y1 or y2 is missing. The \code{plot3D} will draw 3D plots for all pairs of responses.
#' @param ... Additional optional arguments.
#'
#' @rdname plot3D
#' @method plot3D PAsso
#' @export
plot3D.PAsso <- function(
  object,
  y1, y2, ...
) {
  # object = PAsso_2; y1 = "selfLR"; y2 = "PID"

  if (!inherits(object, "PAsso")) stop("Input object must be 'PAsso' class.")

  rep_SRs <- object$rep_SRs
  resp_name <- attr(object, "responses")
  n_resp <- length(resp_name)

  if (missing(y1) | missing(y2)) { # If no input for the responses pair, draw all.
    all_pairs <- apply(expand.grid(resp_name, resp_name), 1, paste, collapse=" v.s. ")
    mat_pairs <- matrix(all_pairs, nrow = n_resp, byrow = F)
    plot_titles <- mat_pairs[upper.tri(mat_pairs)]

    plot_list <- list()
    m <- 1; n <- 2

    for (i_plot in 1:(n_resp*(n_resp-1)/2)) {
      plot_list <- plot3D_one(plot_list = plot_list, rep_SRs = rep_SRs,
                              m = m, n = n,
                              plot_titles = plot_titles[i_plot])
      # Update response ---------------------------------------------------------
      if (n == n_resp) { m <- m + 1; n <- m + 1 } else n <- n + 1
    }
  } else {
    m <- which(y1 == resp_name)
    n <- which(y2 == resp_name)

    i_plot <- plot_titles <- apply(expand.grid(y1, y2), 1, paste, collapse=" v.s. ")
    plot_list <- list()

    plot_list <- plot3D_one(plot_list = plot_list, rep_SRs = rep_SRs,
                            m = m, n = n,
                            plot_titles = plot_titles)
    plot_list <- plot_list[[1]]
  }

  return(plot_list)
}
