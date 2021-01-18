#' @include allClasses.R
NULL

#' Obtaining the cutpoints and / or regroupments of a discretization.
#' @name cutpoints
#' @rdname cutpoints-method
#' @aliases cutpoints,glmdisc-method
#' @param glmdisc The trained glmdisc S4 object.
#' @description This defines the method to provide the cutpoints of a trained glmdisc.
#' @examples
#' # Simulation of a discretized logit model
#' set.seed(1)
#' x <- matrix(runif(300), nrow = 100, ncol = 3)
#' cuts <- seq(0, 1, length.out = 4)
#' xd <- apply(x, 2, function(col) as.numeric(cut(col, cuts)))
#' theta <- t(matrix(c(0, 0, 0, 2, 2, 2, -2, -2, -2), ncol = 3, nrow = 3))
#' log_odd <- rowSums(t(sapply(seq_along(xd[, 1]), function(row_id) {
#'   sapply(
#'     seq_along(xd[row_id, ]),
#'     function(element) theta[xd[row_id, element], element]
#'   )
#' })))
#' y <- rbinom(100, 1, 1 / (1 + exp(-log_odd)))
#'
#' sem_disc <- glmdisc(x, y,
#'   iter = 50, m_start = 4, test = FALSE,
#'   validation = FALSE, criterion = "aic"
#' )
#' cutpoints(sem_disc)
methods::setMethod("cutpoints", methods::signature(object = "glmdisc"), function(object) {
  cutpoints <- list()

  for (j in which(object@parameters$types_data == "factor")) {
    cutpoints[[j]] <- apply(prop.table(object@best.disc$bestLinkFunction[[j]], 2), 2, which.max)
  }

  for (j in which(object@parameters$types_data == "numeric")) {
    data_disc <- data.frame(disc = object@disc.data[, j], cont = object@cont.data[, j], stringsAsFactors = TRUE)

    cut1 <- stats::aggregate(cont ~ disc, data = data_disc, min)
    cut2 <- stats::aggregate(cont ~ disc, data = data_disc, max)

    cut1 <- cut1[order(cut1$cont), ]
    cut2 <- cut2[order(cut2$cont), ]

    cut1 <- cut1[-1, ]
    cut2 <- cut2[-nrow(cut2), ]

    cutpoints[[j]] <- rowMeans(cbind(cut1$cont, cut2$cont))
  }

  names(cutpoints) <- colnames(object@disc.data)[-ncol(object@disc.data)]

  return(cutpoints)
})
