#' Stepwise efficient forward selection in decomposable graphical models
#' @description Stepwise efficient forward selection in decomposable graphical models
#' @param x A \code{fwd} object
#' @param df data.frame
#' @param q Penalty term in the stopping criterion  (\code{0} = AIC and \code{1} = BIC)
#' @param thres A threshold mechanism for choosing between two different ways of calculating the entropy. Can Speed up the procedure with the "correct" value.
#' @details A \code{fwd} object can be created using the \code{gengraph} constructor with \code{type = "fwd"}
#' @return A \code{fwd} object; a subclass of \code{gengraph}) used for forward selection.
#' @references \url{https://arxiv.org/abs/1301.2267}, \url{https://doi.org/10.1109/ictai.2004.100}
#' @examples
#'
#' d <- derma[, 10:25]
#'
#' g <- gengraph(d, type = "fwd")
#' s <- walk(g, d)
#' print(s)
#' plot(s)
#' adj_lst(s)
#' adj_mat(s)
#'
#' @seealso \code{\link{fit_graph}}, \code{\link{walk.bwd}}, \code{\link{gengraph}}
#' @export
walk.fwd <- function(x, df, q = 0.5, thres = 5) {
  ## -----------------------------------------------------------------------------
  ##                    STORE CURRENT INFORMATION
  ## -----------------------------------------------------------------------------
  vab   <- unlist(strsplit(x$e, "\\|"))
  va    <- vab[1]
  vb    <- vab[2]
  Ca    <- x$MSI[[attr(x$e, "idx")]]$C1
  Cb    <- x$MSI[[attr(x$e, "idx")]]$C2
  Sab   <- x$MSI[[attr(x$e, "idx")]]$S
  Cab   <- c(Sab, va, vb)

  # G_prime           <- x
  G_prime_adj       <- x$G_adj
  G_prime_A         <- x$G_A

  # Adding the new edge (va, vb)
  G_prime_adj[[va]] <- c(G_prime_adj[[va]], vb)
  G_prime_adj[[vb]] <- c(G_prime_adj[[vb]], va)
  G_prime_A[va, vb] <- 1L 
  G_prime_A[vb, va] <- 1L
  CG_prime     <- x$CG    
  CG_prime_A   <- x$CG_A
  G_dbl_prime  <- subgraph(Sab, x$G_A) # make_G_dbl_prime(Sab, x$G_A)
  msi_prime    <- x$MSI
        
  ## Vertices connected to a and b in G_dbl_prime
  G_dbl_prime_lst <- as_adj_lst(G_dbl_prime)
  cta <- dfs(G_dbl_prime_lst, va)
  ctb <- dfs(G_dbl_prime_lst, vb)

  ## -----------------------------------------------------------------------------
  ##       INSERTING Cab BETWEEN Ca AND Cb IN CG_prime AND REMOVE (Ca, Cb)
  ## -----------------------------------------------------------------------------
  ## Add Cab to CG_prime and add the edges (Ca, Cab) and (Cb, Cab) to CG_prime_A
  ins          <- vector("numeric", nrow(CG_prime_A))
  ins[attr(x$e, "ins")] <- 1L
  CG_prime_A   <- rbind(CG_prime_A, ins)
  CG_prime_A   <- cbind(CG_prime_A, c(ins, 0L)) # 0L since no loops
  rownames(CG_prime_A) <- NULL
  colnames(CG_prime_A) <- NULL
  CG_prime[[length(CG_prime) + 1L]]   <- Cab

  ## Delete (Ca, Cb) in CG_prime_A
  CG_prime_A[matrix(attr(x$e, "ins"), 1)]      <- 0L
  CG_prime_A[matrix(rev(attr(x$e, "ins")), 1)] <- 0L

  ## -----------------------------------------------------------------------------
  ##                  DELETING EDGES FROM CG IN CG_prime
  ## -----------------------------------------------------------------------------
  TVL  <- vector("list", 0L) # Temporary Vertex List (see Altmueller)
  # msi corresponds to CG
  Sabs <- .map_lgl(x$MSI,  function(s) setequal(s$S, Sab))
  prone_to_deletion <- x$MSI[Sabs]
  MSab <- .map_lgl(prone_to_deletion, function(z) { # See Altmueller
    es <- names(z$e)
    if (x$e %in% es) {
      return(TRUE)
    } else if(rev_es(x$e) %in% es) { # Bottleneck
      return(TRUE)
    } else {
      return(FALSE)
    }
  })

  etd <- edges_to_delete(prone_to_deletion, TVL, MSab, Cab, Sab, cta, ctb) 
  delete_edges <- etd$del
  TVL <- etd$TVL

  if (neq_empt_lst(delete_edges)) {
    delete_idx <- lapply(delete_edges, function(z) {
      sapply(seq_along(CG_prime), function(k) {
        y <- CG_prime[[k]]
        ifelse(setequal(z$C1, y) || setequal(z$C2, y), k, NA)
      })
    })
    for (i in delete_idx) {
      k <- stats::na.omit(i)
      CG_prime_A[k[1], k[2]] <- 0L
      CG_prime_A[k[2], k[1]] <- 0L
    }
    ## Step 3.3.2 in Jordan
    msi_prime[which(msi_prime %in% delete_edges)] <- NULL
  }

  ## -----------------------------------------------------------------------------
  ##                       ADDING EDGES TO CG_prime
  ## -----------------------------------------------------------------------------
  ## Prime clique indicies
  Cp_idx     <- CG_prime_A[attr(x$e, "ins"),]
  C_prime_Ca <- which(Cp_idx[1, , drop = TRUE] == 1L)
  C_prime_Ca <- C_prime_Ca[-length(C_prime_Ca)]
  C_prime_Cb <- which(Cp_idx[2, , drop = TRUE] == 1L) 
  C_prime_Cb <- C_prime_Cb[-length(C_prime_Cb)]

  add_a   <- which_Cp_from_Cx_to_Cab(CG_prime, C_prime_Ca, Ca, va, Cab, Sab, ctb, TVL)
  add     <- c(add_a$add, add_a$add_tvl)
  add_b   <- which_Cp_from_Cx_to_Cab(CG_prime, C_prime_Cb, Cb, vb, Cab, Sab, cta, TVL)
  add     <- unique(c(add, add_b$add, add_b$add_tvl))

  if (neq_empt_num(add)) {
    CG_prime_A[add, length(CG_prime)] <- 1L
    CG_prime_A[length(CG_prime), add] <- 1L
  }

  ## Needed to update msi_prime
  C_primes <- CG_prime[add]
  
  ## -----------------------------------------------------------------------------
  ##                       DELETE Ca AND Cb IF IN Cab
  ## -----------------------------------------------------------------------------
  Ca_in_Cab <- all(Ca %in% Cab)
  Cb_in_Cab <- all(Cb %in% Cab)
  Ca_Cb_idx <- attr(x$e, "ins")
  CG_Ca_idx <- attr(x$e, "ins")[1]
  CG_Cb_idx <- attr(x$e, "ins")[2]

  if (Ca_in_Cab || Cb_in_Cab) {
    if (Ca_in_Cab && Cb_in_Cab) {
      CG_prime     <- CG_prime[-Ca_Cb_idx]
      CG_prime_A   <- CG_prime_A[-Ca_Cb_idx, -Ca_Cb_idx]
      msi_prime    <- msi_prime[!is_Ca_or_Cb(msi_prime, Ca, Cb)]
    }
    else if (Ca_in_Cab) {
      CG_prime     <- CG_prime[-CG_Ca_idx]
      CG_prime_A   <- CG_prime_A[-CG_Ca_idx, -CG_Ca_idx]
      msi_prime  <- msi_prime[!is_Cx(msi_prime, Ca)]
      ## We need to update C_primes to include Ca or Cb if they are not in Cab!!!
      C_primes <- c(C_primes, list(Cb)) # Since Cb not in Cab then
    } else {
      CG_prime     <- CG_prime[-CG_Cb_idx]
      CG_prime_A   <- CG_prime_A[-CG_Cb_idx, -CG_Cb_idx]
      msi_prime    <- msi_prime[!is_Cx(msi_prime, Cb)]
      ## We need to update C_primes to include Ca or Cb if they are not in Cab!!!
      C_primes  <- c(C_primes, list(Ca)) # Since Ca not in Cab then
    }
  } else {
    msi_prime <- msi_prime[!is_Ca_and_Cb(msi_prime, Ca, Cb)]
    C_primes <- c(C_primes, list(Ca), list(Cb)) # Since Ca and Cb not in Cab then
  }

  ## ---------------------------------------------------------
  ##                CALCULATE NEW ENTROPIES
  ## ---------------------------------------------------------
  ue        <- update_edges_from_C_primes_to_Cab(df, C_primes, Cab, va, vb, x$MEM, x$LV, q, thres)
  msi_prime <- c(msi_prime, ue$msi)
  x$MEM     <- ue$mem
  ## ---------------------------------------------------------
  ##                   RETURN THE GRAPH
  ## ---------------------------------------------------------
  x$G_adj <- G_prime_adj
  x$G_A   <- G_prime_A
  x$CG    <- CG_prime
  x$CG_A  <- CG_prime_A
  x$MSI   <- msi_prime
  if (!neq_empt_lst(msi_prime)) {
    # If the graph is complete
    x$e <- new_edge(d_qic = NULL)
    return(x)
  } 
  x$e <- find_new_edge(msi_prime, CG_prime)
  return(x)
}
