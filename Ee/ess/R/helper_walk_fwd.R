fwd_init <- function(x, df, q = 0.5) {
  # x : fwd object
  M       <- nrow(df)
  penalty <- log(M)*q + (1 - q)*2
  nodes <- colnames(df)
  n     <- length(nodes)
  pairs <- utils::combn(nodes, 2,  simplify = FALSE)
  for (j in 1:n) x$MEM[[nodes[j]]] <- entropy(df[nodes[j]])
  x$MSI <- lapply(seq_along(pairs), function(p) {
    v      <- pairs[[p]]
    edge_v <- sort_(v)
    x$MEM[[edge_v]] <<- entropy(df[v])
    ed       <- x$MEM[[v[1]]] + x$MEM[[v[2]]] - x$MEM[[edge_v]]
    dev      <- 2 * M * ed
    d_parms  <- -prod(x$LV[v] - 1)
    d_qic    <- dev + penalty * d_parms
    if (d_qic >= attr(x$e, "d_qic")) {
      x$e <<- new_edge(edge_v, d_qic, p)
    }
    list(C1 = v[1], C2 = v[2], S = character(0L), e = structure(ed, names = edge_v))
  })
  # If the null graph is the best, d_qic = 0 and x$e = character(0)
  if (neq_empt_chr(as.vector(x$e))) attr(x$e, "ins") <- match(es_to_vs(x$e)[[1]], x$CG)
  return(x)
}

find_new_edge <- function(msi_prime, CG_prime) {
  max_es   <- sapply(msi_prime, function(z) z$e[which.max(z$e)])
  max_idx  <- which.max(max_es)
  max_e    <- names(max_idx)
  max_atom <- msi_prime[[max_idx]]
  max_ins  <- unlist(sapply(seq_along(CG_prime), function(x) {
    cond <- setequal(CG_prime[[x]], max_atom$C1) || setequal(CG_prime[[x]], max_atom$C2)
    if( cond ) return(x)
  }))
  new_edge(max_e, unname(max_es[max_idx]), max_idx, max_ins)
}

# m : msi object
# Cx: clique (vector of characters) 
is_Cx <- function(m, Cx) {
  .map_lgl(m, function(x) setequal(x$C1, Cx) || setequal(x$C2, Cx))
} 

is_Ca_or_Cb <- function(m, x, y) {
  .map_lgl(m, function(z) {
    is_CaCb <- setequal(z$C1, x) || setequal(z$C2, y)
    is_CbCa <- setequal(z$C1, y) || setequal(z$C2, x)
    is_CaCb || is_CbCa
  })
}

is_Ca_and_Cb <- function(m, x, y) {
  .map_lgl(m, function(z) {
    is_CaCb <- setequal(z$C1, x) && setequal(z$C2, y)
    is_CbCa <- setequal(z$C1, y) && setequal(z$C2, x)
    is_CaCb || is_CbCa
  })
}

edges_to_delete <- function(prone_to_deletion, TVL, MSab, Cab, Sab, cta, ctb) {
  ## TVL  : Temporary vertex list
  ## MSab : logical indicating which prone_to_deletion (those with C1 \cap C2 = Sab) that has eab = (va, vb)
  etd <- .map_lgl(seq_along(MSab), function(k) {
    delete <- FALSE
    x <- prone_to_deletion[[k]]
    if (MSab[k]) {
      if (!all(x$C1 %in% Cab) && !all(x$C2 %in% Cab)) {
        delete <- TRUE
        TVL <- c(TVL, x$C1, x$C2)
      }
    }
    else {
      C1_minus_Sab <- setdiff(x$C1, Sab)
      C2_minus_Sab <- setdiff(x$C2, Sab)
      va_in_C1_Sab <- all(C1_minus_Sab %in% cta)
      va_in_C2_Sab <- all(C2_minus_Sab %in% cta)
      vb_in_C1_Sab <- all(C1_minus_Sab %in% ctb)
      vb_in_C2_Sab <- all(C2_minus_Sab %in% ctb)
      delete <- FALSE
      if ((va_in_C1_Sab && vb_in_C2_Sab) || (va_in_C2_Sab && vb_in_C1_Sab)) {
        delete <- TRUE
      }
    }
    delete
  })
  list(del = Filter(neq_null, prone_to_deletion[etd]), TVL = TVL)
}

which_Cp_from_Cx_to_Cab <- function(CG_prime, C_prime_Cx, Cx, vx, Cab, Sab,  cty, TVL) {
  add <- vector("numeric", 0L)
  for (k in C_prime_Cx) {
    Cp    <- CG_prime[[k]]
    Sp    <- intersect(Cp, Cx)
    Sn    <- intersect(Cp, Cab)
    Sab_x <- c(Sab, vx)
    if (neq_empt_chr(Sab)) {
      if (all(Sp %in% Sab_x)) {
        add <- c(add, k)
      }
      if (setequal(Sn, Sab_x) && !(any(setdiff(Cp, Sab_x) %in% cty))) {
        add <- c(add, k)
      }
    }
    else {
      if (vx %in% Sn || !neq_empt_chr(Sp)) {
        add <- c(add, k)
      }
    }
  }
  add_tvl <- which(CG_prime %in% TVL) # For all C' in TVL add (C', Cab) to CG_prime
  list(add = unique(add), add_tvl = unique(add_tvl)) 
}

update_edges_from_C_primes_to_Cab <- function(df, Cps, Cab, va, vb, mem, LV, q = 0.5, thres = 5) {
  # Cps : C_primes
  M       <- nrow(df)
  penalty <- log(M)*q + (1 - q)*2
  sep <- lapply(Cps, function(Cp) {
    Sp          <- intersect(Cp, Cab)
    eligs_Cab   <- setdiff(Cab, Cp)
    eligs_Cp    <- setdiff(Cp, Cab)
    eligs       <- expand.grid(eligs_Cp, eligs_Cab)
    eligs       <- apply(eligs, 1, paste, collapse = "|")
    eligs_names <- eligs
    eligs <- sapply(eligs, function(e) {
      el <- entropy_difference(e, Sp, df, mem)
      mem <<- el$mem
      el$ent
    })
    names(eligs) <- eligs_names
    dev <- 2 * M * eligs
    d_parms <- .map_dbl(es_to_vs(names(eligs)), function(x) {
      -prod(LV[x] - 1) * prod(LV[Sp])
    })
    d_qic <- dev + penalty * d_parms
    list(C1 = Cp, C2 = Cab, S = Sp, e = d_qic)
  })
  return(list(msi = Filter(neq_null, sep) , mem = mem))
}
