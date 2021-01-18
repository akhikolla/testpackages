# #' Discretize a dataset using a trained SEM discretization scheme
# #'
# #' This function discretizes a user-specified dataset using a pre-trained SEM discretization scheme.
# #' @param link A multinomial logit model.
# #' @param df The dataframe, containing the same variables as the one used to train the discretization scheme, to be discretized.
# #' @keywords discretization, predict
# #' @importFrom stats predict
# #' @export
# #' @references
# #' Asad Hasan, Wang Zhiyu and Alireza S. Mahani (2015). mnlogit: Multinomial Logit Model. R package version 1.2.4. \url{https://CRAN.R-project.org/package=mnlogit}
# #' @examples
# #' set.seed(1)
# #' x = matrix(runif(100), nrow = 100, ncol = 1)
# #' cuts = seq(0,1,length.out = 4)
# #' xd = as.numeric(cut(x,cuts))
# #'
# #' long_dataset <- data.frame(e = as.vector(sapply(xd,function(var) (seq(1:3)[seq(1:3)]==var))),
# #' x = as.vector(sapply(x[,1], function(var) rep(var,3))),
# #' names = as.character(as.vector(rep(seq(1:3)[seq(1:3)],length(x)))),stringsAsFactors=FALSE)
# #' link <- mnlogit::mnlogit(e ~ 1 | x | 1, data=long_dataset, choiceVar = "names")
# #' discretize_link(link,as.data.frame(x))


discretize_link <- function(link, df, m_start) {
  n <- nrow(df)
  d <- ncol(df)
  types_data <- sapply(df[1, ], class)
  # Cas complets
  continu_complete_case <- !is.na(df)
  if (sum(!(types_data %in% c("numeric", "factor"))) > 0) {
    stop(simpleError("Unsupported data types. Columns of predictors must be numeric or factor."))
  }
  emap <- array(0, c(n, d))

  if (d > 1) {
    for (j in (1:d)) {
      if (types_data[j] == "numeric") {
        if (sum(is.na(link[[j]])) == 0) {
          t <- predict(link[[j]], newdata = data.frame(x = df[continu_complete_case[, j], j], stringsAsFactors = TRUE), type = "probs")

          if (is.vector(t)) {
            t <- cbind(1 - t, t)
            colnames(t) <- c("1", "2")
          }

          if (sum(!continu_complete_case[, j]) > 0) {
            t_bis <- matrix(NA, nrow = nrow(df), ncol = ncol(t) + 1)
            t_bis[continu_complete_case[, j], 1:ncol(t)] <- t
            t_bis[continu_complete_case[, j], ncol(t) + 1] <- 0
            t_bis[!continu_complete_case[, j], ] <- t(matrix(c(rep(0, ncol(t)), 1), nrow = ncol(t) + 1, ncol = sum(!continu_complete_case[, j])))
            colnames(t_bis) <- c(colnames(t), m_start + 1)
            t <- t_bis
          }
        } else {
          t <- matrix(1, nrow = n, ncol = 1)
          colnames(t) <- "1"
        }
      } else {
        if (!is.null(link[[j]])) {
          t <- prop.table.robust(t(sapply(df[, j], function(row) link[[j]][, row])), 1)
        } else {
          t <- matrix(1, nrow = n, ncol = 1)
          colnames(t) <- "1"
        }
      }
      emap[, j] <- unlist(apply(t, 1, function(p) names(which.max(p))))
    }
  } else if (types_data == "numeric") {
    t <- predict(link[[1]], newdata = data.frame(x = df, stringsAsFactors = TRUE), type = "probs")
    emap[, 1] <- apply(t, 1, function(p) names(which.max(p)))
  } else {
    t <- prop.table(t(apply(df, 1, function(row) link[[1]][, row])), 1)
    emap[, 1] <- apply(t, 1, function(p) names(which.max(p)))
  }

  return(emap)
}
