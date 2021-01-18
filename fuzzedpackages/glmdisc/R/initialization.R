first_argument_checks <- function(criterion, labels, predictors, validation) {
  if (!criterion %in% c("gini", "aic", "bic")) {
    stop(simpleError("criterion must be gini, aic or bic"))
  }
  if (!(any(class(predictors) %in% c("matrix", "array", "data.frame")))) {
    stop(simpleError(paste0("`predictors` must be an array, a matrix or a dataframe, you provided a "), class(predictors)))
  }
  if (!length(labels) == length(predictors[, 1])) {
    stop(simpleError("labels and predictors must be of same length"))
  }
  if ((criterion == "gini") & (validation == FALSE)) {
    warning("Using Gini index on training set might yield an overfitted model.")
  }

  if ((criterion %in% c("aic", "bic")) & (validation == TRUE)) {
    warning("No need to penalize the log-likelihood when a validation set is used. Using log-likelihood instead of AIC/BIC.")
  }
}


second_argument_checks <- function(types_data, d, interact) {
  if (sum(!(types_data %in% c("numeric", "factor"))) > 0) {
    stop(simpleError("Unsupported data types. Columns of predictors must be numeric or factor."))
  }
  if (d < 2 && interact) {
    stop(simpleError("Cannot perform interaction screening with less than 2 features, set `interact = FALSE`"))
  }
}


initialize_e_emap <- function(n, d, types_data, continu_complete_case, m_start, predictors) {
  e <- emap <- array(0, c(n, d))
  for (j in which(types_data == "numeric")) {
    e[continu_complete_case[, j], j] <- emap[continu_complete_case[, j], j] <- as.factor(sample(1:m_start, sum(continu_complete_case[, j]), replace = TRUE))
    e[!continu_complete_case[, j], j] <- emap[!continu_complete_case[, j], j] <- m_start + 1
  }
  for (j in which(types_data == "factor")) {
    # For categorical features, if m_start is above the number of original levels, then the quantized version is initialized at the original number of levels
    if (m_start > nlevels(predictors[, j])) {
      e[, j] <- emap[, j] <- as.factor(sample(1:nlevels(predictors[, j]), n, replace = TRUE))
    } else {
      # Otherwise proceed as for continuous features
      e[, j] <- emap[, j] <- as.factor(sample(1:m_start, n, replace = TRUE))
    }
  }
  return(list(e, emap))
}


initialize_interaction <- function(d, predictors, continu_complete_case, labels, ensemble) {
  delta <- matrix(sample(0:1, d^2, replace = TRUE, prob = c(0.5, 0.5)), nrow = d, ncol = d)
  delta[lower.tri(delta)] <- 0
  diag(delta) <- 0
  xPrincipal <- paste0("V", 1:d)

  p_delta <- matrix(0, nrow = d, ncol = d)
  p_delta[lower.tri(p_delta)] <- 0
  diag(p_delta) <- 0

  # The proposal distribution of Metropolis-Hastings, p_delta, is calculated as follows
  for (j in 1:(d - 1)) {
    for (k in (j + 1):d) {
      if (nlevels(factor(predictors[continu_complete_case[, j] & continu_complete_case[, k] & ensemble[[1]], j])) == 1) stop(paste0("Some levels are scarce such that feature ", j, " only has 1 level in train. Try with validation = F"))
      if (nlevels(factor(predictors[continu_complete_case[, j] & continu_complete_case[, k] & ensemble[[1]], k])) == 1) stop(paste0("Some levels are scarce such that feature ", k, " only has 1 level in train. Try with validation = F"))
      sans_inter <- stats::glm(labels ~ X1 + X2, family = stats::binomial(link = "logit"), data = data.frame(labels = labels[continu_complete_case[, j] & continu_complete_case[, k] & ensemble[[1]]], X1 = predictors[continu_complete_case[, j] & continu_complete_case[, k] & ensemble[[1]], j], X2 = predictors[continu_complete_case[, j] & continu_complete_case[, k] & ensemble[[1]], k], stringsAsFactors = TRUE))
      avec_inter <- stats::glm(labels ~ X1 + X2 + X1:X2, family = stats::binomial(link = "logit"), data = data.frame(labels = labels[continu_complete_case[, j] & continu_complete_case[, k] & ensemble[[1]]], X1 = predictors[continu_complete_case[, j] & continu_complete_case[, k] & ensemble[[1]], j], X2 = predictors[continu_complete_case[, j] & continu_complete_case[, k] & ensemble[[1]], k], stringsAsFactors = TRUE))
      p_delta[j, k] <- 1 / (1 + exp(-sans_inter$deviance - log(sum(ensemble[[1]])) * length(sans_inter$coefficients) + avec_inter$deviance + log(sum(ensemble[[1]])) * length(avec_inter$coefficients)))
    }
  }

  return(list(delta, p_delta, xPrincipal))
}
