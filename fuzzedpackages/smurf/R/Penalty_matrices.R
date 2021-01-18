###############################################
#
# Penalty matrices
#
###############################################


# Penalty matrix for Lasso penalized predictor
#
# n.par: Number of parameters for the predictor
.pen.mat.lasso <- function(n.par) {
  
  return(diag(n.par))
}


# Penalty matrix for Group Lasso penalized predictor
#
# n.par: Number of parameters for the predictor
.pen.mat.grouplasso <- function(n.par) {
  
  return(diag(n.par))
}


# Penalty matrix for Fused Lasso penalized predictor
#
# n.par: Number of parameters for the predictor
# refcat: Numeric indicating which of the levels is the reference category, default is first one
# when no reference category is present it is 0
.pen.mat.flasso <- function(n.par, refcat = 1) {
  
  matr <- .pen_aux(1:(n.par + (refcat > 0) * 1), 1)

  if (refcat > 0) {
    # Remove row which corresponds to reference category (first row by default)
    return(t(matr[-refcat, ]))
    
  } else {
    return(t(matr))
  }
  
}


# Penalty matrix for Generalized Fused Lasso penalized predictor
#
# n.par: Number of parameters for the predictor
# refcat: Logical which indicates if a reference category is present
.pen.mat.gflasso <- function(n.par, refcat = TRUE) {
  
  # No reference category
  if (!refcat) {
    n.par <- n.par - 1
  }
  
  matr <- .pen_aux(1:(n.par+1), n.par)
  
  if (refcat) {
    # Remove first row which corresponds to reference category
    return(t(matr[-1, ]))
    
  } else {
    return(t(matr))
  }
}


# Penalty matrix for 2D Fused Lasso penalty
#
# n.par.row: Number of parameters for the column predictor
# n.par.col: Number of parameters for the row predictor
.pen.mat.2dflasso <- function(n.par.row, n.par.col) {
  
  total.columns <- n.par.row * n.par.col
  total.rows <- 3 * 1 + 3 * (n.par.row -2 + n.par.col - 2) + 2 * (2 + (n.par.col-2) * (n.par.row-2)) + 
    1 * (n.par.row-2+n.par.col-2) + 0 * 1
  matr <- matrix(0, nrow = total.rows, ncol = total.columns)
  
  i.row <- 1
  
  for (j in 1:total.columns) {
    
    if ((j-1) < n.par.row | (j-1) %% n.par.row == 0 ) {
      matr[i.row, j] <- 1
      i.row <- i.row + 1
    }
    if (j %% n.par.row != 0 & j <= n.par.row * n.par.col - n.par.row) {
      matr[i.row, j] <- -1
      matr[i.row, j+1] <-  1
      i.row <- i.row + 1
      matr[i.row, j] <- -1
      matr[i.row, j + n.par.row] <-  1
      i.row <- i.row + 1
    }
    if (j %% n.par.row == 0 & j != n.par.row * n.par.col) {
      matr[i.row, j] <- -1
      matr[i.row, j + n.par.row] <-  1
      i.row <- i.row + 1
    }
    
    if (j > n.par.row * n.par.col - n.par.row & j != n.par.row * n.par.col) {
      matr[i.row,j] <- -1
      matr[i.row,j+1] <-  1
      i.row <- i.row + 1
    }
  }
  
  return(matr)
}


# Compute penalty matrix for Graph-Guided Fused Lasso penalized predictor based on adjacency matrix
#
# adj.matrix: Adjacency matrix of the graph
# lev.names: Names of the levels of the predictor
# refcat: Numeric indicating which of the levels is the reference category, default is first one
# when no reference category is present it is 0
.pen.mat.ggflasso <- function(adj.matrix, lev.names = NULL, refcat = 1) {
  
  # Check if square matrix
  if (nrow(adj.matrix) != ncol(adj.matrix)) {
    stop("An adjacency matrix needs to be square.")
  }
  
  # Check if symmetric
  if (!isSymmetric(adj.matrix)) {
    stop("An adjacency matrix needs to be symmetric (including row and column names).")
  }
  
  # Convert to sparse matrix if not done already
  if (!(class(adj.matrix)[1] %in% c("dgCMatrix", "dgTMatrix"))) {
    adj.matrix <- as(adj.matrix, "sparseMatrix")
  }
  
  # Check if all elements are zero or one
  if (any(adj.matrix@x != 1)) {
    stop("All elements of an adjacency matrix need to be zero or one.")
  }
  
  # lev.names is not NULL, check if all levels are present
  if (!is.null(lev.names)) {
    
    if (is.null(colnames(adj.matrix))) {
      stop("An adjacency matrix needs to have column names.")
    }
    
    if (is.null(rownames(adj.matrix))) {
      stop("An adjacency matrix needs to have row names.")
    }
    
    # Check if all levels are present
    if (!all(as.character(lev.names) %in% as.character(colnames(adj.matrix)))) {
      stop("Not all levels are present in the column names of the adjacency matrix.")
    }
    
  }
  
  if (refcat > 1) {
    # Change reference category by moving columns and rows in adjacency matrix
    
    adj.matrix <- adj.matrix[, c(refcat, (1:ncol(adj.matrix))[-refcat])]
    adj.matrix <- adj.matrix[c(refcat, (1:nrow(adj.matrix))[-refcat]), ]
    
    refcat <- 1
  }
  
  # Number of rows for penalty matrix, this is equal to total number of unique edges
  total <- sum(triu(adj.matrix))
  
  # Penalty matrix
  pen.mat <- matrix(0, total, ncol(adj.matrix) - refcat)
  # Change column names
  if (refcat) {
    colnames(pen.mat) <- colnames(adj.matrix)[-1]
    
  } else {
    colnames(pen.mat) <- colnames(adj.matrix)
  }
  
  
  
  # Counter for row of penalty matrix
  ind <- 0L
  
  for (j in 1:ncol(pen.mat)) {
    # Count number of edges connecting level j+1 with levels below j+1
    # First level is reference category, hence j + 1 if refcat = TRUE; j otherwise
    total.j <- sum(adj.matrix[1:j, j + refcat])
    # Indices of elements that are one
    ind.j <- which(adj.matrix[1:j, j + refcat] == 1)
    
    # More than one edge with level below j+1
    if (total.j > 0) {
      
      for (l in 1:total.j) {
        # 1 for level j+1
        pen.mat[ind + l, j] <- 1
        
        # If not reference category when refcat=TRUE, otherwise always TRUE
        if (ind.j[l] > refcat) {
          # Column minus -1 because of reference category
          pen.mat[ind + l, ind.j[l] - refcat] <- -1 
        }
        
      }
    }
    # Increase counter
    ind <- ind + total.j
  }
  return(pen.mat)
}


# When m=1, a k by k-1 matrix, where k=length(p),
# with -1 on the diagonal and 1 on the first diagonal below this
# This corresponds to all possible (-1,1) pairs (looking by the columns)
#
# When m>=1, matrix where all possibilities with -1 above 1 with distance at most m-1
# are considered
# Note that m is always set to min(k-1, m)
.pen_aux <- function(p, m = 1, ...) {
  
  k <- length(p)
  m <- min(k - 1, m)
  if (k > 1) {
    if (k > 2) {
      mat <- matrix(nrow = k, ncol = 0)
      for (j in 1:m) {
        mat2 <- diag(1, ncol = k, nrow = k)
        for (i in 1:(k - j)) {
          mat2[p[i], p[i + j]] <- -1
        }
        mat2 <- mat2[, -p[1:j]]
        mat <- cbind(mat, mat2)
      }
    } else {
      mat <- matrix(c(-1, 1), ncol = 1)
    }
  } else {
    mat <- matrix(ncol = 0, nrow = 1)
  }
  colnames(mat) <- NULL
  return(mat)
}



# Compute eigenvalue decomposition of t(pen.mat[[j]]) %*% pen.mat[[j]] for all penalty types
# except "none", "lasso" and "grouplasso"
#
# pen.mat: List with (weighted) penalty matrix per predictor
# pen.cov: List with penalty type per predictor
.pen.mat.eig <- function(pen.mat, pen.cov) {
  
  pen.mat.aux <- vector("list", length(pen.cov))
  ind <- which(pen.cov %in% c("flasso", "gflasso", "2dflasso", "ggflasso"))
  if (length(ind > 0)) {
    
    for (j in ind) {
      # Return NULL if error
      tmp <- tryCatch(eigen(t(pen.mat[[j]]) %*% pen.mat[[j]]), error = function(e) NULL)
      
      if (!is.null(tmp)) {
        # Get eigenvectors and -values if eigen did not give an error
        pen.mat.aux[[j]]$Q <- tmp$vectors
        pen.mat.aux[[j]]$eigval <- tmp$values
        
      } else {
        # eigen gave an error, use slower ADMM version (check happens when calling C++ code)
        pen.mat.aux[[j]]$Q <- as.matrix(0)
        pen.mat.aux[[j]]$eigval <- 0
      }
      
    }
  }
  return(pen.mat.aux)
}