#### Function check ####
popCheck = function(
  y_S,
  X_S,
  X_Sa,
  Z_Sa,
  Z_U
) {
  if (nrow(X_S) != nrow(as.matrix(y_S))) {
    stop('y_S and X_S have different number of rows.')
  }

  if (nrow(X_Sa) != nrow(Z_Sa)) {
    stop('X_Sa and Z_Sa have different number of rows.')
  }

  if (ncol(X_S) != ncol(X_Sa)) {
    stop('X_S and X_Sa have different number of columns.')
  }

  if (ncol(Z_Sa) != ncol(Z_U)) {
    stop('Z_Sa and Z_U have different number of columns.')
  }
}
