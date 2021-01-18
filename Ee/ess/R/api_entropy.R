na_tab <- function(df, a) {
  ct <- table(df[, a, drop = FALSE])
  names(dimnames(ct)) <- a # Needed for the onedimensional separators
  ct
}

joint_entropy <- function(df) {
  ## if( class(df) == "character" ) stop( "From entropy function: df is not a data.frame!" )
  x  <- na_tab(df, colnames(df))
  Nx <- sum(x)
  entropy_table <- apply(x, seq_along(dim(x)), function(y) {
    ifelse(y == 0 , 0, y/Nx * log(y/Nx) )
  })
  -sum(entropy_table)
}

joint_entropy2 <- function(df) {
  A  <- apply(df, 1, paste0, collapse = "")
  x  <- table(A)
  Nx <- sum(x)
  -sum(x/Nx * log(x/Nx))
}

#' Joint Entropy
#' 
#' @description Calculates the joint entropy over discrete variables in \code{df}
#' 
#' @param df data.frame
#' @param thres A threshold mechanism for choosing between two different ways of calculating the entropy. Can Speed up the procedure with the "correct" value.
#' @return A number representing the entropy of the variables in \code{df}.
#' @examples
#' entropy(derma[1:100, 1:3])
#' 
#' @export
entropy <- function(df, thres = 5) {
  if (ncol(df) <= thres) return(joint_entropy(df))
  else return(joint_entropy2(df))
}

entropy_difference <- function(e, S, df, mem, thres = 5) {
  # FIX: Replace ifelse with if() ... else ...
  # FIX: For large problems, delete some entropies?
  #      - maybe as for memory accessible and
  #      - clean mem if it reaches 50% of this or so?
  v <- unlist(es_to_vs(e))

  H_S <- 0L        
  if (neq_empt_chr(S)) {
    S_ <- sort_(S)
    if (exists(S_, envir = mem, inherits = FALSE)) {
      H_S <- mem[[S_]]
    } else {
      H_S  <- entropy(df[S], thres)
      mem[[S_]] <- H_S
    }
  }

  H_S_x <- 0L        
  Sx <- sort_(c(S, v[1]))
  if (exists(Sx, envir = mem, inherits = FALSE)) {
    H_S_x <- mem[[Sx]]
  } else {
    H_S_x  <- entropy(df[c(S, v[1])], thres)
    mem[[Sx]] <- H_S_x 
  }
  
  H_S_y <- 0L
  Sy <- sort_(c(S, v[2]))
  if (exists(Sy, envir = mem, inherits = FALSE)) {
    H_S_y <- mem[[Sy]]
  } else {
    H_S_y  <- entropy(df[c(S, v[2])], thres)
    mem[[Sy]] <- H_S_y 
  }
  
  H_S_xy <- 0L
  Sxy <- sort_(c(S, v))
  if (exists(Sxy, envir = mem, inherits = FALSE)) {
    H_S_xy <- mem[[Sxy]]
  } else {
    H_S_xy  <- entropy(df[c(S, v)], thres)
    mem[[Sxy]] <- H_S_xy  
  }
  
  H_S_x_S_y <- H_S_x + H_S_y
  # Test needed to avoid < 0 due to floating point errors
  edge_ent <- ifelse(isTRUE(all.equal(H_S_x_S_y, H_S_xy)), 0L,  H_S_x_S_y - H_S_xy - H_S)
  return(list(ent = edge_ent, mem = mem ))
}

