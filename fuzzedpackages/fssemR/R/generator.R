## generate simulated network structure for FSSEM algorithm
##' @title randomFSSEMdata
##' @param n  number of observations
##' @param p  number of genes
##' @param k  number of eQTLs
##' @param sparse  ratio of edges / gene_number
##' @param df   ratio of differential edges among two network
##' @param sigma2 noise variance of error
##' @param u  variance of bias in SEM model.
##' @param type type of generated network, can be selected as DG, ER, Scale-free network
##' @param dag  network is directed-acyclic or not. Default TRUE
##' @param coef Range of absolute value of coefficients in simulated network matrices. Default (0.2, 0.4), or (0.5, 1)
##' @param nhub If you select to generate ER network, nhub is the number of pre-defined hub node number. Default 2
##' @return list of generated data
##' \describe{
##' \item{Data}{ List of observed, Xs, Ys, Sk }
##' \item{Vars}{ List of model, Bs, Fs, mu, n, p, k }
##' }
##' @importFrom igraph graph_from_adjacency_matrix is.dag
##' @importFrom mvtnorm rmvnorm
##' @importFrom qtl sim.map sim.cross pull.geno find.markerpos sim.geno makeqtl
##' @export
randomFSSEMdata = function(n, p, k, sparse = 0.1, df = 0.2, sigma2 = 0.01, u = 5, type = c("DG", "ER"), dag = TRUE,
                           coef = c(0.2, 0.4), nhub = 2) {
  mincoef = coef[1]
  maxcoef = coef[2]
  type = match.arg(type)
  DG = function() {
    B = vector("list", 2)
    B[[1]] = matrix(0, nrow = p, ncol = p)
    d = p * p
    ne = rbinom(1, d, sparse / (p - 1))  # edges
    niter1 = 0
    while (sum(B[[1]]) <= ne & niter1 < 2 * d) {
      ix = runif(1, min = 1, max = d)
      B[[1]][ix] = TRUE
      if (dag) {
        graph = graph_from_adjacency_matrix(B[[1]])
        B[[1]][ix] = is.dag(graph)
      }
      niter1 = niter1 + 1
    }
    B[[2]] = B[[1]]
    nonzero = which(B[[1]] != 0)
    nonedge = which(B[[1]] == 0)
    nd = ceiling(ne * df)
    rmed = rbinom(1, nd, 0.5)
    while (sum(abs(B[[1]] - B[[2]])) < rmed) {
      ix = sample(nonzero, 1)
      B[[2]][ix] = FALSE
    }
    niter2 = 0
    while (sum(B[[2]]) <= ne & sum(B[[2]] != B[[1]]) < nd & niter2 < 2 * d) {
      ix = sample(nonedge, 1)
      B[[2]][ix] = TRUE
      if (dag) {
        graph = graph_from_adjacency_matrix(B[[2]])
        B[[2]][ix] = is.dag(graph)
      }
      niter2 = niter2 + 1
    }
    ei = which(B[[1]] & B[[2]])
    B[[1]][ei] = B[[2]][ei] = runif(length(ei), min = mincoef, max = maxcoef) * sample(c(-1, 1), length(ei), replace = T)
    di = which(B[[1]] & !B[[2]])
    B[[1]][di] = runif(length(di), min = mincoef, max = maxcoef) * sample(c(-1, 1), length(di), replace = T)
    di = which(!B[[1]] & B[[2]])
    B[[2]][di] = runif(length(di), min = mincoef, max = maxcoef) * sample(c(-1, 1), length(di), replace = T)
    # check
    if (niter1 < d & niter2 < d & any(B[[1]] != B[[2]])) {
      if (!dag) {
        if (abs(det(diag(p) - B[[1]])) > 1e-2 & abs(det(diag(p) - B[[2]])) > 1e-2) {
          B
        } else {
          NULL
        }
      } else {
        B
      }
    } else {
      NULL
    }
  }
  ER = function() {
    B = vector("list", 2)
    B[[1]] = matrix(0, nrow = p, ncol = p)
    d = p * p
    ne = rbinom(1, d, sparse / (p - 1))
    trace = sapply(1:p, function(x){ (x - 1) * p + x })
    niter1 = 1
    while (sum(B[[1]]) < ne & niter1 < d) {
      ix = sample(setdiff(seq(1, d), trace), 1)
      B[[1]][ix] = TRUE
      if (dag) {
        graph = graph_from_adjacency_matrix(B[[1]])
        B[[1]][ix] = is.dag(graph)
      }
      niter1 = niter1 + 1
    }
    ## hub node
    regulon = sample(1:p, size = nhub)
    nsize = max(p * sparse / p * 10, 0.03 * (p - 1))
    for (g in regulon) {
      B[[1]][, g] = FALSE
      B[[1]][g, ] = FALSE
      while (TRUE) {
        ix = sample(setdiff(seq(1, p), g), as.integer(nsize))
        B[[1]][ix, g] = 1
        if (dag) {
          graph = graph_from_adjacency_matrix(B[[1]])
          B[[1]][ix, g] = is.dag(graph)
        }
        if (any(B[[1]][, g] != 0)) {
          break
        }
      }
    }
    ni = which(B[[1]] != 0)
    B[[1]][ni] = runif(length(ni), min = mincoef, max = maxcoef) * sample(c(-1, 1), length(ni), replace = T)
    B[[2]] = B[[1]]
    pnode = sample(setdiff(1:p, regulon), size = sparse * (p - 1) * df - nhub, replace = FALSE)
    pnode = c(pnode, regulon)
    for (j in pnode) {
      xi = which(B[[1]][,j] != 0)
      B[[2]][xi, j] = B[[2]][xi, j] + runif(length(xi), min = mincoef, max = maxcoef) *
        sample(c(-1, 1), length(xi), replace = T)
    }
    B[[2]][abs(B[[2]]) < 0.1] = 0
    # check
    if (any(B[[1]] != B[[2]])) {
      if (!dag) {
        if (abs(det(diag(p) - B[[1]])) > 1e-2 & abs(det(diag(p) - B[[2]])) > 1e-2) {
          B
        } else {
          NULL
        }
      } else {
        B
      }
    } else {
      NULL
    }
  }
  ## generate eQTL
  eQTLs = function(n, p, k, len = 100) {
    e = k / p
    X = round(2 * matrix(runif(n * k), nrow = k)) + 1
    Qtlmap = sim.map(
      len = rep(len, p),
      n.mar = e,
      eq.spacing = FALSE,
      include.x = FALSE
    )
    Cross   = sim.cross(Qtlmap, n.ind = n, type = "f2")
    m       = lapply(Qtlmap, names)
    F = matrix(0, nrow = p, ncol = k)
    Sk     = lapply(1:p, function(i) {
      s = seq(0, e - 1) * p + i
      F[i, s] <<- 1
      data = t(X[s,])
      colnames(data) = colnames(Cross$geno[[i]]$data)
      Cross$geno[[i]]$data <<- data
      s
    })
    Cross   = sim.geno(Cross)
    eQTLs = lapply(1:p, function(i) {
      pos = find.markerpos(Cross, m[[i]])
      makeqtl(Cross, chr = pos[, "chr"], pos = pos[, "pos"])
    })
    list(F = F, X = X, Sk = Sk, eQTLs = eQTLs, Cross = Cross, marker = m)
  }

  B = NULL
  while (is.null(B)) {
    if (type == "DG") {
      B = DG()
    } else {
      B = ER()
    }
  }
  QTL = eQTLs(n, p, k)
  F = QTL$F
  X = QTL$X
  Sk = QTL$Sk
  E = lapply(1:2, function(i) {
    sqrt(sigma2) * t(rmvnorm(n, mean = rep(0, p), sigma = diag(p)))
  })
  mu = tcrossprod(rnorm(p, 0, u), rep(1, n))
  Y = lapply(1:2, function(i) {
    solve(diag(p) - B[[i]]) %*% (F %*% X + mu + E[[i]])
  })
  names = paste0("g", seq(1, p))
  names(QTL$marker) = names
  names(QTL$eQTLs)  = names
  Cross = lapply(1:2, function(i) {
    cross = QTL$Cross
    cross$pheno = as.data.frame(t(Y[[i]]))
    colnames(cross$pheno) = names
    cross
  })
  list(
    Data = list(
      X = X, Y = Y, Sk = Sk
    ),
    Vars = list(
      B = lapply(B, function(x){Matrix(x, sparse = T)}),
      F = F,
      mu = mu, n = n, p = p, k = k
    ),
    QTL = list(
      Cross  = Cross,
      eQTLs  = QTL$eQTLs,
      Pheno  = names,
      marker = QTL$marker
    )
  )
}






## generate simulated network structure for FSSEM algorithm
##' @title randomFSSEMdata2
##' @param n  number of observations. Vector for unbalance observations
##' @param p  number of genes
##' @param k  number of eQTLs
##' @param sparse  ratio of edges / gene_number
##' @param df   ratio of differential edges among two network
##' @param sigma2 noise variance of error
##' @param u  variance of bias in SEM model.
##' @param type type of generated network, can be selected as DG, ER, Scale-free network
##' @param dag  network is directed-acyclic or not. Default TRUE
##' @param coef Range of absolute value of coefficients in simulated network matrices. Default (0.2, 0.4), or (0.5, 1)
##' @param nhub If you select to generate ER network, nhub is the number of pre-defined hub node number. Default 2
##' @return list of generated data
##' \describe{
##' \item{Data}{ List of observed, Xs, Ys, Sk }
##' \item{Vars}{ List of model, Bs, Fs, mu, n, p, k }
##' }
##' @importFrom igraph graph_from_adjacency_matrix is.dag
##' @importFrom mvtnorm rmvnorm
##' @importFrom qtl sim.map sim.cross pull.geno find.markerpos sim.geno makeqtl
##' @export
randomFSSEMdata2 = function(n, p, k, sparse = 0.1, df = 0.2, sigma2 = 0.01, u = 5, type = c("DG", "ER"), dag = TRUE,
                           coef = c(0.2, 0.4), nhub = 2) {
  mincoef = coef[1]
  maxcoef = coef[2]
  type = match.arg(type)
  DG = function() {
    B = vector("list", 2)
    B[[1]] = matrix(0, nrow = p, ncol = p)
    d = p * p
    ne = rbinom(1, d, sparse / (p - 1))  # edges
    niter1 = 0
    while (sum(B[[1]]) <= ne & niter1 < 2 * d) {
      ix = runif(1, min = 1, max = d)
      B[[1]][ix] = TRUE
      if (dag) {
        graph = graph_from_adjacency_matrix(B[[1]])
        B[[1]][ix] = is.dag(graph)
      }
      niter1 = niter1 + 1
    }
    B[[2]] = B[[1]]
    nonzero = which(B[[1]] != 0)
    nonedge = which(B[[1]] == 0)
    nd = ceiling(ne * df)
    rmed = rbinom(1, nd, 0.5)
    while (sum(abs(B[[1]] - B[[2]])) < rmed) {
      ix = sample(nonzero, 1)
      B[[2]][ix] = FALSE
    }
    niter2 = 0
    while (sum(B[[2]]) <= ne & sum(B[[2]] != B[[1]]) < nd & niter2 < 2 * d) {
      ix = sample(nonedge, 1)
      B[[2]][ix] = TRUE
      if (dag) {
        graph = graph_from_adjacency_matrix(B[[2]])
        B[[2]][ix] = is.dag(graph)
      }
      niter2 = niter2 + 1
    }
    ei = which(B[[1]] & B[[2]])
    B[[1]][ei] = B[[2]][ei] = runif(length(ei), min = mincoef, max = maxcoef) * sample(c(-1, 1), length(ei), replace = T)
    di = which(B[[1]] & !B[[2]])
    B[[1]][di] = runif(length(di), min = mincoef, max = maxcoef) * sample(c(-1, 1), length(di), replace = T)
    di = which(!B[[1]] & B[[2]])
    B[[2]][di] = runif(length(di), min = mincoef, max = maxcoef) * sample(c(-1, 1), length(di), replace = T)
    # check
    if (niter1 < d & niter2 < d & any(B[[1]] != B[[2]])) {
      if (!dag) {
        if (abs(det(diag(p) - B[[1]])) > 1e-2 & abs(det(diag(p) - B[[2]])) > 1e-2) {
          B
        } else {
          NULL
        }
      } else {
        B
      }
    } else {
      NULL
    }
  }
  ER = function() {
    B = vector("list", 2)
    B[[1]] = matrix(0, nrow = p, ncol = p)
    d = p * p
    ne = rbinom(1, d, sparse / (p - 1))
    trace = sapply(1:p, function(x){ (x - 1) * p + x })
    niter1 = 1
    while (sum(B[[1]]) < ne & niter1 < d) {
      ix = sample(setdiff(seq(1, d), trace), 1)
      B[[1]][ix] = TRUE
      if (dag) {
        graph = graph_from_adjacency_matrix(B[[1]])
        B[[1]][ix] = is.dag(graph)
      }
      niter1 = niter1 + 1
    }
    ## hub node
    regulon = sample(1:p, size = nhub)
    nsize = max(p * sparse / p * 10, 0.03 * (p - 1))
    for (g in regulon) {
      B[[1]][, g] = FALSE
      B[[1]][g, ] = FALSE
      while (TRUE) {
        ix = sample(setdiff(seq(1, p), g), as.integer(nsize))
        B[[1]][ix, g] = 1
        if (dag) {
          graph = graph_from_adjacency_matrix(B[[1]])
          B[[1]][ix, g] = is.dag(graph)
        }
        if (any(B[[1]][, g] != 0)) {
          break
        }
      }
    }
    ni = which(B[[1]] != 0)
    B[[1]][ni] = runif(length(ni), min = mincoef, max = maxcoef) * sample(c(-1, 1), length(ni), replace = T)
    B[[2]] = B[[1]]
    pnode = sample(setdiff(1:p, regulon), size = sparse * (p - 1) * df - nhub, replace = FALSE)
    pnode = c(pnode, regulon)
    for (j in pnode) {
      xi = which(B[[1]][,j] != 0)
      B[[2]][xi, j] = B[[2]][xi, j] + runif(length(xi), min = mincoef, max = maxcoef) *
        sample(c(-1, 1), length(xi), replace = T)
    }
    B[[2]][abs(B[[2]]) < 0.1] = 0
    # check
    if (any(B[[1]] != B[[2]])) {
      if (!dag) {
        if (abs(det(diag(p) - B[[1]])) > 1e-2 & abs(det(diag(p) - B[[2]])) > 1e-2) {
          B
        } else {
          NULL
        }
      } else {
        B
      }
    } else {
      NULL
    }
  }
  ## generate eQTL
  eQTLs = function(n, p, k, len = 100) {
    e = k / p
    X = lapply(1:2, function(i) {
      round(2 * matrix(runif(n[i] * k), nrow = k)) + 1
    })
    Qtlmap = sim.map(
      len = rep(len, p),
      n.mar = e,
      eq.spacing = FALSE,
      include.x = FALSE
    )
    Cross   = lapply(1:2, function(i) {
      sim.cross(Qtlmap, n.ind = n[i], type = "f2")
    })
    m = lapply(Qtlmap, names)
    F = matrix(0, nrow = p, ncol = k)
    Sk = lapply(1:p, function(i) {
      s = seq(0, e - 1) * p + i
      F[i, s] <<- 1
      for(j in 1:2) {
        data = t(X[[j]][s,])
        colnames(data) = colnames(Cross[[j]]$geno[[i]]$data)
        Cross[[j]]$geno[[i]]$data <<- data
      }
      s
    })
    Cross   = lapply(1:2, function(i) { sim.geno(Cross[[i]]) })
    eQTLs = lapply(1:p, function(i) {
      pos = find.markerpos(Cross[[1]], m[[i]])
      makeqtl(Cross[[1]], chr = pos[, "chr"], pos = pos[, "pos"])
    })
    list(F = F, X = X, Sk = Sk, eQTLs = eQTLs, Cross = Cross, marker = m)
  }

  B = NULL
  while (is.null(B)) {
    if (type == "DG") {
      B = DG()
    } else {
      B = ER()
    }
  }
  QTL = eQTLs(n, p, k)
  F = QTL$F
  X = QTL$X
  Sk = QTL$Sk
  E = lapply(1:2, function(i) {
    sqrt(sigma2) * t(rmvnorm(n[i], mean = rep(0, p), sigma = diag(p)))
  })
  mu = lapply(1:2, function(i) {
    tcrossprod(rnorm(p, 0, u), rep(1, n[i]))
  })
  Y = lapply(1:2, function(i) {
    solve(diag(p) - B[[i]]) %*% (F %*% X[[i]] + mu[[i]] + E[[i]])
  })
  names = paste0("g", seq(1, p))
  names(QTL$marker) = names
  names(QTL$eQTLs)  = names
  Cross = lapply(1:2, function(i) {
    cross = QTL$Cross[[i]]
    cross$pheno = as.data.frame(t(Y[[i]]))
    colnames(cross$pheno) = names
    cross
  })
  list(
    Data = list(
      X = X, Y = Y, Sk = Sk
    ),
    Vars = list(
      B = lapply(B, function(x){Matrix(x, sparse = T)}),
      F = F,
      mu = mu, n = n, p = p, k = k
    ),
    QTL = list(
      Cross  = Cross,
      eQTLs  = QTL$eQTLs,
      Pheno  = names,
      marker = QTL$marker
    )
  )
}
