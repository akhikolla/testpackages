## generate random graph(DAG) of G genes and their correspond differential network
## B_1 & B_2
## @param Ng gene number nodes
## @param e Expected number of edges per node
## @param d Expected ratio of differential edges per node (0.1)
## @param dag DAG or not
require(igraph)
require(Matrix)
require(glmnet)
getrandDAG = function(Ng,
                      e,
                      dag = TRUE,
                      d = 0.1,
                      Bmin = 0.5,
                      Bmax = 1,
                      maxit = Ng * Ng) {
  B1 = matrix(0,
              nrow = Ng,
              ncol = Ng)
  Nc = Ng * Ng
  Ne = rbinom(1, Nc, e / (Ng - 1))
  ## iteration mark
  iter1 = 0
  while (sum(B1) < Ne & iter1 < maxit) {
    edge     = runif(1, min = 1, max = Nc)
    B1[edge] = TRUE
    if (dag) {
      g        = graph_from_adjacency_matrix(B1)
      B1[edge] = is.dag(g)
    }
    iter1 = iter1 + 1
  }
  B2  = B1
  nn  = which(B1 != 0)
  nz  = which(B1 == 0)
  Nd  = ceiling(Ne * d)
  Ndf = rbinom(1, Nd, 0.5)
  while (sum(abs(B1 - B2)) < Ndf) {
    edge     = sample(nn, 1)
    B2[edge] = FALSE
  }
  iter2 = 0
  while (sum(B2) < Ne & iter2 < maxit) {
    edge     = sample(nz, 1)
    B2[edge] = TRUE
    if (dag) {
      g        = graph_from_adjacency_matrix(B2)
      B2[edge] = is.dag(g)
    }
    iter2 = iter2 + 1
  }
  ne  = which(B1 & B2)
  n1  = which(B1 & !(B2))
  n2  = which(!(B1) & B2)
  B1[ne] = B2[ne] = runif(length(ne), min = Bmin, max = Bmax) * sample(c(-1, 1), length(ne), replace = T)
  B1[n1] = runif(length(n1), min = Bmin, max = Bmax) * sample(c(-1, 1), length(n1), replace = T)
  B2[n2] = runif(length(n2), min = Bmin, max = Bmax) * sample(c(-1, 1), length(n2), replace = T)

  if(iter1 < maxit & iter2 < maxit & any(B1 != B2)) {
    if(!dag) {
      detIB1 = det(diag(Ng) - B1)
      detIB2 = det(diag(Ng) - B2)
      if(abs(detIB1) > 1e-6 & abs(detIB2) > 1e-6){
        list(B1 = B1, B2 = B2)
      } else {
        NULL
      }
    } else {
      list(B1 = B1, B2 = B2)
    }
  } else {
    NULL
  }
}


#' @param N   number of sample
#' @param Ng  number of gene
#' @param k   number of eQTL
getrandeQTL = function(N, Ng, Nk) {
  step = Nk / Ng
  X = round(2 * matrix(runif(Nk * N), nrow = Nk)) + 1
  G = matrix(0,
             nrow = Ng,
             ncol = Nk)
  ix = lapply(1:Ng,
              function(i) {
                s = seq(0, step - 1) * Ng + i
                G[i, s] <<- 1
                s
              })
  list(G = G, X = X, sk = ix)
}

## randNetinit
## randomly generate regulatory network for fixed seed
randNetinit = function(Ng = 10,
                       Nk = 10,
                       r = 0.3,
                       d = 0.1,
                       ...) {
  B = getrandDAG(Ng,
                 e = Ng * r,
                 dag = dag,
                 d = d,
                 ...)
  while (is.null(B)) {
    B = getrandDAG(Ng,
                   e = Ng * r,
                   dag = dag,
                   d = d,
                   ...)
  }
  B
}

require(mvtnorm)
getrandsem = function(N = 200,
                      Ng = 10,
                      Nk = 10,
                      r = 0.3,
                      d = 0.1,
                      dag = TRUE,
                      sigma = 0.1,
                      B = NULL,
                      ...) {
  if (is.null(B)) {
    B = getrandDAG(Ng,
                   e = Ng * r,
                   dag = dag,
                   d = d,
                   ...)
    while (is.null(B)) {
      B = getrandDAG(Ng,
                     e = Ng * r,
                     dag = dag,
                     d = d,
                     ...)
    }
  }
  Q = getrandeQTL(N, Ng, Nk)
  F = Q[[1]]
  X = Q[[2]]
  sk = Q[[3]]
  E1 = sigma * t(rmvnorm(N, mean = rep(0, Ng), sigma = diag(Ng)))
  E2 = sigma * t(rmvnorm(N, mean = rep(0, Ng), sigma = diag(Ng)))
  Y1 = solve(diag(Ng) - B[[1]]) %*% (F %*% X + E1)
  Y2 = solve(diag(Ng) - B[[2]]) %*% (F %*% X + E2)
  list(
    obs = list(
      Y1 = Y1,
      Y2 = Y2,
      X = X,
      sk = sk
    ),
    var = list(
      B1 = Matrix(B[[1]], sparse = T),
      B2 = Matrix(B[[2]], sparse = T),
      F = Matrix(F, sparse = T),
      N = N,
      Ng = Ng,
      Nk = Nk
    )
  )
}

## data = getrandsem(N = 200, Ng = 30, Nk = 90, r = 0.1, d = 0.1, sigma = 1, dag = TRUE)
## datn = getrandsem(N = 200, Ng = 30, Nk = 90, r = 0.1, d = 0.1, sigma = 1, dag = FALSE)

## ultility funcitons
center = function(X) {
  apply(X, 1, function(x) {
    x - mean(x)
  })
}

submatX = function(data) {
  ## submatrix for X on eQTL
  sk = data$obs$sk
  X  = data$obs$X
  lapply(sk, function(ix) {
    X[ix, , drop = F]
  })
}

## X(X^TX)^{-1}X^T
projection = function(X) {
  X %*% solve(crossprod(X)) %*% t(X)
}


## use QR more slowly
projection.QR = function(X) {
  qr = qr.default(X, LAPACK = TRUE)
  Q  = qr.qy(qr, diag(1, nrow = nrow(qr$qr), ncol = qr$rank))
  tcrossprod(Q)
}

## centeralized Y (gene expression) and X (eQTL quantitive)
centralize = function(X, Y) {
  meanX = lapply(X, rowMeans)
  meanY = rowMeans(Y)
  X = lapply(X, center)
  Y = center(Y)
  list(X = X,
       Y = Y,
       muX = meanX,
       muY = meanY)
}

## ridge regression for estimate sigma2 in gene expression
#' @example
#' X = submatX(data)
#' Y = data$obs$Y1
#' B = constrained_L2reg(X, Y, rho = 0.1)
constrained_L2reg = function(X, Y, rho) {
  # gene number(M) & sample number(N)
  std = centralize(X, Y)
  X = std$X
  Y = std$Y
  meanX = std$muX
  meanY = std$muY
  M = ncol(Y)
  N = nrow(Y)
  B = Matrix(0,
             nrow = M,
             ncol = M,
             sparse = T)
  f = list()
  err = 0
  for (i in 1:M) {
    Xi = X[[i]]                   ## n x sk
    Pi = diag(N) - projection(Xi) ## n x n
    yi = Y[, i, drop = F]         ## n x 1
    Yi = Y[,-i, drop = F]         ## n x (p-1)
    ## (Y^TPY + rho)^{-1}Y^TPy
    bi = solve(crossprod(Yi, Pi %*% Yi) + rho * diag(M - 1)) %*% t(Yi) %*% Pi %*% yi
    ## bi = glmnet(Pi %*% Yi, Pi %*% yi, alpha = 0, lambda = rho)[["beta"]][, 1]
    B[i, -i] = bi
    f[[i]] = solve(crossprod(Xi)) %*% t(Xi) %*% (yi - Yi %*% bi)
    err = err + crossprod(yi - Yi %*% bi - Xi %*% f[[i]])
  }
  sigma2 = err / (M * N - 1)
  mu     = (diag(M) - B) %*% meanY - sapply(1:M, function(i) {
    meanX[[i]] %*% f[[i]]
  })
  list(
    B = as.matrix(B),
    F = f,
    sigma2 = sigma2,
    mu = mu
  )
}

## cross-validation on ridge regression to estimate sigma2
#' @param nrho number of L2 penalty's coefficient
#' @param ncv  number of cross-validation
#' @example
#' sigma2 = getsigma2_L2reg(X, Y, nrho = 15, ncv = 5)
getsigma2_L2reg = function(X, Y, nrho = 10, ncv = 5) {
  rho_factors = 10 ** (seq(-6, 2, length.out = nrho))
  N = ncol(Y)
  M = nrow(Y)
  cv.err  = matrix(0, nrow = nrho, ncol = ncv)
  cv.fold = sample(seq(1, ncv), size = N, replace = T)
  irho    = 1
  for (rho in rho_factors) {
    for (cv in 1:ncv) {
      ytrain = Y[, cv.fold != cv]
      xtrain = lapply(X, function(x) {
        x[, cv.fold != cv, drop = F]
      })
      ytest  = Y[, cv.fold == cv]
      xtest  = lapply(X, function(x) {
        x[, cv.fold == cv, drop = F]
      })
      fit    = constrained_L2reg(xtrain, ytrain, rho)
      ftest  = lapply(1:M, function(i) {
        crossprod(fit$F[[i]], xtest[[i]])
      })
      ftest  = do.call(rbind, ftest)
      cv.err[irho, cv] = norm((diag(M) - fit$B) %*% ytest - ftest - fit$mu, type = "f") ^
        2
    }
    irho = irho + 1
  }
  cv.mean = rowMeans(cv.err)
  rho.min = rho_factors[which.min(cv.mean)]
  fit = constrained_L2reg(X, Y, rho.min)
  list(rho.opt = rho.min, sigma2.opt = fit$sigma2)
}


##---------------------------------
# utility functions for SML-lasso #
##---------------------------------
obj_cwiseSML = function(N, a0, a1, a2, lambda, w, sigma2) {
  function(x) {
    -N / 2 * sigma2 * log((a0 - x) ^ 2) - a1 * x + 1 / 2 * a2 * x ^ 2 + lambda * w * abs(x)
  }
}

## a0 = det(I - B) / cij + Bij => c = 1
grad_cwiseSML = function(N, a0, a1, a2, lambda, w, sigma2) {
  ## t := conditions of Bij
  ## t  = {1  | Bij > 0}
  ## t  = {-1 | Bij < 0}
  ## t  = {0  | Bij = 0}
  function(t) {
    list(
      a = -a2,
      b = a1 + a2 * a0 - lambda * w * t ,
      c = N * sigma2 + (lambda * w * t - a1) * a0
    )
  }
}

# ax^2 + bx + c = 0
poly2_solver = function(a, b, c) {
  r = b ^ 2 - 4 * a * c
  if (r < 0) {
    list(n = 0, x = NULL)
  } else if (r == 0) {
    list(n = 1, x = c(-b / (2 * a)))
  } else {
    list(n = 2, x = c((-b - sqrt(r)) / (2 * a), (-b + sqrt(r)) / (2 * a)))
  }
}


## solve SML problem by component-wise update
## component-wise --> row-wise update Bij
#' @param sigma2 extimated from constrained_L2reg
#' @param B B0 initialization
#' @param f F initialization
#' @example
#' X = submatX(data)
#' Y = data$obs$Y1
#' sigma2 = getsigma2_L2reg(X, Y, nrho = 20, ncv = 5)
#' param.init = constrained_L2reg(X, Y, rho = sigma2$rho.opt)
#' param.opt  = sparse_maximum_likehood_cwise(B = param.init$B, f = param.init$F, Y = Y, X = X, sigma2 = param.init$sigma2[1], N = data$var$N, Ng = data$var$Ng, lambda = 15, maxit = 100)
#' B = param.init$B; f=param.init$F; Y = Y;  X = X; sigma2 = param.init$sigma2; N = data$var$N; Ng = data$var$Ng; Nk = data$var$Nk; lambda = 0.1
sparse_maximum_likehood_cwise = function(B,
                                         f,
                                         Y,
                                         X,
                                         sigma2,
                                         N,
                                         Ng,
                                         lambda,
                                         weighted = TRUE,
                                         maxit = 100,
                                         verbose = 2) {
  ## data centralization
  std = centralize(X, Y)
  X = std$X
  Y = std$Y
  meanX = std$muX
  meanY = std$muY
  ## update for eQTL coeffs
  f0 = list()
  f1 = list()
  for (i in 1:Ng) {
    Xi = X[[i]]                   # n x sk
    yi = Y[, i, drop = F]         # n x 1 (for specific gene i)
    Pi = solve(crossprod(Xi)) %*% t(Xi)
    f0[[i]] = Pi %*% yi
    f1[[i]] = Pi %*% Y            # f[[i]] = f0[[i]] - f1[[i]] %*% B[i,]
  }
  ## update for gnet coeffs
  niter  = 1
  ImB    = diag(Ng) - B
  IBinv  = solve(ImB)
  wB     = 1 / abs(B)
  while (niter <= maxit) {
    B.prev = B
    f.prev = f
    for (i in 1:Ng) {
      ## IBinv i column and j row -> c_{ij}
      ## c_{ij} / det(I - B) = (I - B)^{-1}_{j, i}
      ci = IBinv[, i]
      dbi = vector("numeric", Ng)
      for (j in 1:Ng) {
        ## update B[i, j] for i != j
        if (i != j) {
          bij.prev = bij = B[i, j]
          wij      = if (weighted) {
            wB[i, j]
          } else {
            1
          }
          mij      = ci[j]
          ## i-th row of B
          bi       = ImB[i, ]
          bi[j]    = 0
          ## Yej is the j-th column of Y
          Yej      = Y[, j, drop = F]
          a1       = crossprod(Y %*% bi - X[[i]] %*% f[[i]], Yej)
          a2       = crossprod(Yej)
          ## if mij == 0, cij == 0
          if (mij == 0) {
            if (abs(a1) > lambda * wij) {
              bij = sign(a1) * (abs(a1) - lambda * wij) / a2
            } else {
              bij = 0
            }
          } else {
            a0     = 1 / mij + bij.prev
            cond   = list(1,-1, 0)  # bij condition
            obj    = obj_cwiseSML(N, a0, a1, a2, lambda, wij, sigma2)
            grad   = grad_cwiseSML(N, a0, a1, a2, lambda, wij, sigma2)
            params = lapply(cond, function(t) {
              x = if (t != 0) {
                tmp  = do.call(poly2_solver, grad(t))
                tmp$x[tmp$x * t > 0]
              } else {
                0
              }
              list(x = x)
            })
            params = unlist(params)
            objval = sapply(params, obj)
            mix    = which.min(objval)
            bij    = params[mix]
          }
          dbij   = bij.prev - bij
          dbi[j] = dbij
          B[i, j] = bij
          ci     = ci / (1 + dbij * mij)
          ##IBinv  = IBinv - IBinv[,i,drop = F] %*% IBinv[j,,drop = F] / (1/dbij + mij)
          ImB    = diag(Ng) - B
        }
      }  ## for(j in 1:Ng)
      ## (ImB + ei^T %*% 1 %*% dbi)^{-1}
      IBinv = IBinv - IBinv[, i, drop = F] %*% dbi %*% IBinv / (1 + dbi %*% ImB[, i, drop =
                                                                                  F])[1]
      f[[i]] = f0[[i]] - f1[[i]] %*% B[i, ]
    } ## for(i in 1:Ng)
    Berr = norm(B - B.prev, type = "f") / (1 + norm(B.prev, type = "f"))
    Ferr = sum(sapply(1:Ng, function(i) {
      norm(f[[i]] - f.prev[[i]], type = "f")
    })) / (1 + sum(sapply(1:Ng, function(i) {
      norm(f.prev[[i]], type = "f")
    })))
    err  = Berr + Ferr
    if (verbose >= 2) {
      cat(sprintf("SML: iteration = %d, error = %f\n", niter, err))
    }
    niter = niter + 1
    if (err < 1e-4 || niter > maxit) {
      mu     = (diag(Ng) - B) %*% meanY - sapply(1:Ng, function(i) {
        meanX[[i]] %*% f[[i]]
      })
      B = Matrix(B, sparse = T)
      break
    }
  } ## while(niter <= maxit)
  list(
    B = B,
    f = f,
    mu = mu,
    niter = niter,
    err = err
  )
}

##---------------------------------
# utility functions for SML-CD    #
##---------------------------------
## objective function for one condition
#' @param B gene network coefficients
#' @param f eQTL coefficients
#' @param Y gene expression (centralized)
#' @param X eQTL quantitive (centralized)
#' @param sigma2 estimated sigma2 in logliklihood
#' @param N number of sample
#' @param lambda lambda coefficient in weighted lasso
#' @param weighted weighted lasso
#' @description argmin -N / 2 * log(det(I - B))^2 + 1 / (2 * sigma2) * ||(I - B) %*% Y - F %*% X - mu||_F^2 + lambda * abs(weight * B)
sparse_likelihood = function(B,
                             f,
                             Y,
                             X,
                             mu,
                             sigma2,
                             N,
                             lambda,
                             weight,
                             detIB,
                             type = c("lang", "prim", "err")) {
  logdet = -N / 2 * log(detIB ^ 2)
  IBerr2 = 0
  for (i in 1:length(f)) {
    err = Y[i, ] - B[i,-i] %*% Y[-i, ] - crossprod(f[[i]], X[[i]]) - mu[i]
    IBerr2 = IBerr2 + sum(err ^ 2)
  }
  if (match.arg(type) == "lang") {
    logdet + IBerr2 / (2 * sigma2) + lambda * sum(abs(weight * B))
  } else if (match.arg(type) == "prim") {
    logdet + IBerr2 / (2 * sigma2)
  } else {
    IBerr2
  }
}

#' @description gradient row-wise
#' @param c vector c = ci / det(I-B); (ng-1) x 1
grad_rwise_SML = function(N, c, Yp, Hy, sigma2) {
  function(x) {
    N * c + (Yp %*% x - t(Hy)) / sigma2
  }
}


#' @description lipshitz moduli row-wise
#' @param o solve((I-B)[-i,] %*% t((I-B)[-i,]))
#' @param gs (I-B)[-i,-i]
#' @param si gi[,i]
#' @param Yp
lips_rwise_SML = function(N,
                          o,
                          gs,
                          si,
                          c2i,
                          detIBi,
                          maxEigenYp,
                          sigma2,
                          Ng) {
  ogs = o %*% gs
  ImO = diag(Ng - 1) - crossprod(gs, ogs)
  sOg = crossprod(si, ogs)
  c   = 1 - crossprod(si, o %*% si)
  lambda = 1e-6
  x   = -1 * tcrossprod(chol2inv(ImO + diag(Ng - 1) * lambda), sOg)
  L   = N * c2i / (crossprod(x, ImO %*% x) + 2 * sOg %*% x + c) / (detIBi + 1e-6) + maxEigenYp / sigma2
  while(L < 0) {
    lambda = lambda * 10
    x   = -1 * tcrossprod(chol2inv(ImO + diag(Ng - 1) * lambda), sOg)
    L   = N * c2i / (crossprod(x, ImO %*% x) + 2 * sOg %*% x + c) / (detIBi + 1e-6) + maxEigenYp / sigma2
  }
  L
}

#' @description proximal operator for lasso
#' argmin lambda * |w * x| + c / 2 * |x - u|_2^2
prox_lasso = function(lambda, c, u, w) {
  pmax(u - lambda * w / c, 0) + pmin(u + lambda * w / c, 0)
}

## solve SML problem by coordinate descent
## row-wise --> row-wise update B[i,]
#' @param sigma2 extimated from constrained_L2reg
#' @param B B0 initialization (Derived from ridge regression)
#' @param f F initialization
#' @example
#' X = submatX(data)
#' Y = data$obs$Y1
#' sigma2 = getsigma2_L2reg(X, Y, nrho = 20, ncv = 5)
#' param.init = constrained_L2reg(X, Y, rho = sigma2$rho.opt)
#' param.opt1  = sparse_maximum_likehood_bcd(B = param.init$B, f = param.init$F, Y = Y, X = X, sigma2 = param.init$sigma2[1], N = data$var$N, Ng = data$var$Ng, lambda = 15, maxit = 1000, rho = sigma2$rho.opt)
sparse_maximum_likehood_bcd = function(B,
                                       f,
                                       Y,
                                       X,
                                       sigma2,
                                       N,
                                       Ng,
                                       lambda,
                                       weighted = TRUE,
                                       maxit = 100,
                                       verbose = 2) {
  ## data centralization
  std = centralize(X, Y)
  X = std$X
  Y = std$Y
  meanX = std$muX
  meanY = std$muY
  ## update for eQTL row-wise (specific for genes)
  f0 = list()
  f1 = list()
  Yp = list()
  Hy = list()
  Yp.maxEigen = list()
  for (i in 1:Ng) {
    Xi = X[[i]]                    # n x sk
    yi = Y[, i, drop = F]          # n x 1 (for specific gene i)
    Yi = Y[, -i]                   # n x (ng-1) (for specific gene i)
    Pi = solve(crossprod(Xi)) %*% t(Xi)
    f0[[i]] = Pi %*% yi            # f0 \in sk x 1; f1 \in sk x (ng - 1)
    f1[[i]] = Pi %*% Yi            # f[[i]] = f0[[i]] - f1[[i]] %*% bi | bi = B[i,-i]
    Hi = diag(N) - Xi %*% Pi       # n x n projection matrix
    Yp[[i]] = t(Yi) %*% Hi %*% Yi  # (ng-1) x (ng-1)
    Hy[[i]] = t(yi) %*% Hi %*% Yi  # 1 x (ng-1)
    ## maximized eigen-value for Yp
    Yp.maxEigen[[i]] = eigen(Yp[[i]])$values[1]
  }
  ## update for gnet row-wise
  niter = 1
  ImB   = diag(Ng) - B
  IBinv = solve(ImB)
  detIB = det(ImB)
  wB    = 1 / abs(B)
  while (niter <= maxit) {
    B.prev = B
    f.prev = f
    for (i in 1:Ng) {
      ## -N*sigma2*log(det(I-B)^2) + \sum_{i=1}^{Ng} bi^T %*% Yp %*% bi - Hy %*% bi
      ## ci / det(I-B); (ng-1) x 1
      ci = IBinv[-i, i, drop = F]
      bi = t(B[i,-i, drop = F])
      grad = grad_rwise_SML(N, ci, Yp[[i]], Hy[[i]], sigma2[1])
      ## Lipschitz for row-i
      oi  = solve(tcrossprod(ImB[-i,]))
      deti = det(tcrossprod(ImB[-i,]))
      gii = ImB[-i, -i]
      si  = ImB[-i, i, drop = F]
      c2i = sum((ci * detIB) ^ 2)
      gi  = grad(bi)
      Li  = lips_rwise_SML(N, oi, gii, si, c2i, deti, Yp.maxEigen[[i]], sigma2, Ng)
      ## proximal operator for lasso
      ## argmin(lambda * w * |x| + Li/2||x - ui||_2^2)
      ui  = bi - gi / Li[1]
      wBi = wB[i,-i]
      B[i, -i] = prox_lasso(lambda, Li[1], ui[, 1], wBi)
      dbi = B.prev[i,] - B[i,]
      ImB = diag(Ng) - B
      ## update det(I-B) & (I-B)^{-1}
      detIB  = (ImB[i,] %*% IBinv[, i, drop = F])[1] * detIB
      # IBinv  = IBinv - IBinv[, i, drop = F] %*% dbi %*% IBinv / (1 + dbi %*% ImB[, i, drop = F])[1]
      IBinv  = solve(ImB)
      f[[i]] = f0[[i]] - f1[[i]] %*% B[i, -i]
    }
    Berr = norm(B - B.prev, type = "f") / norm(B.prev, type = "f")
    Ferr = sum(sapply(1:Ng, function(i) {
      norm(f[[i]] - f.prev[[i]], type = "f")
    })) / sum(sapply(1:Ng, function(i) {
      norm(f.prev[[i]], type = "f")
    }))
    err  = Berr + Ferr
    if (verbose >= 2) {
      cat(sprintf("SML: iteration = %d, error = %f\n", niter, err))
    }
    niter = niter + 1
    if (err < 1e-4 || niter > maxit) {
      mu     = (diag(Ng) - B) %*% meanY - sapply(1:Ng, function(i) {
        meanX[[i]] %*% f[[i]]
      })
      B = Matrix(B, sparse = T)
      break
    }
  } ## while (niter <= maxit)
  list(
    B = B,
    f = f,
    mu = mu,
    niter = niter,
    err = err,
    detIB = detIB
  )
}

## inertial PALM
# utility functions
inertial_pars = function(opts = c("cont", "lin"),
                         init = 0) {
  switch(
    opts,
    "cont" = function(k) {
      init
    },
    "lin"  = function(k) {
      (k - 1) / (k + 2)
    }
  )
}

## solve SML problem by coordinate descent
## row-wise --> row-wise update B[i,]
#' @param sigma2 extimated from constrained_L2reg
#' @param B B0 initialization (Derived from ridge regression)
#' @param f F initialization
#' @param rho stable for lipschitz calculation(deprecated)
#' @example
#' X = submatX(data)
#' Y = data$obs$Y1
#' sigma2 = getsigma2_L2reg(X, Y, nrho = 20, ncv = 5)
#' param.init = constrained_L2reg(X, Y, rho = sigma2$rho.opt)
#' param.opt2  = sparse_maximum_likehood_iPALM(B = param.init$B, f = param.init$F, Y = Y, X = X, sigma2 = param.init$sigma2[1], N = data$var$N, Ng = data$var$Ng, lambda = 15, maxit = 200)
sparse_maximum_likehood_iPALM = function(B,
                                         f,
                                         Y,
                                         X,
                                         sigma2,
                                         N,
                                         Ng,
                                         lambda,
                                         weighted = TRUE,
                                         inertial = inertial_pars("lin"),
                                         maxit = 100,
                                         verbose = 2,
                                         threshold = 1e-4) {
  ## data centralization
  std = centralize(X, Y)
  X = std$X
  Y = std$Y
  meanX = std$muX
  meanY = std$muY
  ## update for eQTL row-wise (specific for genes)
  f0 = list()
  f1 = list()
  Yp = list()
  Hy = list()
  Yp.maxEigen = list()
  for (i in 1:Ng) {
    Xi = X[[i]]                    # n x sk
    yi = Y[, i, drop = F]          # n x 1 (for specific gene i)
    Yi = Y[, -i]                   # n x (ng-1) (for specific gene i)
    Pi = solve(crossprod(Xi)) %*% t(Xi)
    f0[[i]] = Pi %*% yi            # f0 \in sk x 1; f1 \in sk x (ng - 1)
    f1[[i]] = Pi %*% Yi            # f[[i]] = f0[[i]] - f1[[i]] %*% bi | bi = B[i,-i]
    Hi = diag(N) - Xi %*% Pi       # n x n projection matrix
    Yp[[i]] = t(Yi) %*% Hi %*% Yi  # (ng-1) x (ng-1)
    Hy[[i]] = t(yi) %*% Hi %*% Yi  # 1 x (ng-1)
    ## maximized eigen-value for Yp
    Yp.maxEigen[[i]] = eigen(Yp[[i]])$values[1]
  }
  ## update for gnet row-wise
  niter   = 1
  ImB     = diag(Ng) - B
  IBinv   = solve(ImB)
  detIB   = det(ImB)
  wB      = 1 / abs(B)
  B.prevs = list(B, B)
  while (niter <= maxit) {
    inert.pars = inertial(niter)
    B.inert = B.prevs[[2]] + inert.pars * (B.prevs[[2]] - B.prevs[[1]])
    f.prev  = f
    for (i in 1:Ng) {
      ## -N*sigma2*log(det(I-B)^2) + \sum_{i=1}^{Ng} bi^T %*% Yp %*% bi - Hy %*% bi
      ## ci / det(I-B); (ng-1) x 1
      ci = IBinv[-i, i, drop = F]
      bi = t(B.inert[i,-i, drop = F])
      grad = grad_rwise_SML(N, ci, Yp[[i]], Hy[[i]], sigma2[1])
      ## Lipschitz for row-i
      oi  = solve(tcrossprod(ImB[-i,]))
      deti = det(tcrossprod(ImB[-i,]))
      gii = ImB[-i, -i]
      si  = ImB[-i, i, drop = F]
      c2i = sum((ci * detIB) ^ 2)
      gi  = grad(bi)
      Li  = lips_rwise_SML(N, oi, gii, si, c2i, deti, Yp.maxEigen[[i]], sigma2, Ng)
      ## proximal operator for lasso
      ## argmin(lambda * w * |x| + Li/2||x - ui||_2^2)
      ui  = bi - gi / Li[1]
      wBi = wB[i,-i]
      B[i, -i] = prox_lasso(lambda, Li[1], ui[, 1], wBi)
      dbi = B.prevs[[2]][i,] - B[i,]
      ImB = diag(Ng) - B
      ## update det(I-B) & (I-B)^{-1}
      detIB  = (ImB[i,] %*% IBinv[, i, drop = F])[1] * detIB
      IBinv  = IBinv - IBinv[, i, drop = F] %*% dbi %*% IBinv / (1 + dbi %*% IBinv[, i, drop = F])[1]
      f[[i]] = f0[[i]] - f1[[i]] %*% B[i, -i]
    }
    sigma2 = sigma2_sml(X, Y, B, f, Ng, N)
    Berr = norm(B - B.prevs[[2]], type = "f") / (1 + norm(B.prevs[[2]], type = "f"))
    Ferr = sum(sapply(1:Ng, function(i) {
      norm(f[[i]] - f.prev[[i]], type = "f")
    })) / (1 + sum(sapply(1:Ng, function(i) {
      norm(f.prev[[i]], type = "f")
    })))
    err  = Berr + Ferr
    if (verbose >= 2) {
      cat(sprintf("SML: iteration = %d, error = %f\n", niter, err))
    }
    niter = niter + 1
    B.prevs = list(B.prevs[[2]], B)
    if (err < threshold || niter > maxit) {
      mu     = (diag(Ng) - B) %*% meanY - sapply(1:Ng, function(i) {
        meanX[[i]] %*% f[[i]]
      })
      break
    }
  } ## while (niter <= maxit)
  list(
    B = B,
    f = f,
    mu = mu,
    sigma2 = sigma2,
    niter = niter,
    err = err,
    detIB = detIB
  )
}


###### fused lasso ########
## centeralized Ys (gene expression) and Xs (eQTL quantitive)
#' @example
#' Ys = data$obs[c("Y1", "Y2")]
#' Xs = data$obs$X
centralize_mult = function(Xs, Ys) {
  meanX = lapply(Xs, rowMeans)
  meanY = lapply(Ys, rowMeans)
  Xs = lapply(Xs, center)
  Ys = lapply(Ys, center)
  list(X = Xs,
       Y = Ys,
       muX = meanX,
       muY = meanY)
}

## ridge regression for estimate sigma2 initialization in gene expression
#' @param M number of gene
#' @param N number of sample
#' @example
#' M = data$var$Ng
#' N = data$var$N
#' B = constrained_L2reg_multi(Xs, Ys, sigma2$rho.opt, M, N)
constrained_L2reg_multi = function(X, Ys, rho, M, N) {
  K = length(Ys)
  B = list()
  F = list()
  mu = list()
  err = 0
  df  = 0
  for (i in 1:K) {
    fit = constrained_L2reg(X, Ys[[i]], rho)
    B[[i]]  = as.matrix(fit$B)
    F[[i]]  = fit$F
    mu[[i]] = fit$mu
    err = err + fit$sigma2 * (N * M - 1)
    df  = df + (N * M - 1)
  }
  sigma2 = err / df
  list(
    B = B,
    F = F,
    sigma2 = sigma2,
    mu = mu
  )
}

## cross-validation on ridge regression to estimate sigma2
#' @param nrho number of L2 penalty's coefficient
#' @param ncv  number of cross-validation
#' @example
#' Ys = data$obs[c("Y1", "Y2")]
#' Xs = submatX(data)
#' M  = data$var$Ng
#' N  = data$var$N
#' sigma2 = getsigma2_L2reg_multi(Xs, Ys, nrho = 20, M = M, N = N)
getsigma2_L2reg_multi = function(X,
                                 Ys,
                                 nrho = 10,
                                 ncv = 5,
                                 M,
                                 N) {
  rho_factors = 10 ** (seq(-6, 2, length.out = nrho))
  cv.err  = matrix(0, nrow = nrho, ncol = ncv)
  cv.fold = sample(seq(1, ncv), size = N, replace = T)
  irho    = 1
  for (rho in rho_factors) {
    for (cv in 1:ncv) {
      ytrain = lapply(Ys, function(y) {
        y[, cv.fold != cv, drop = F]
      })
      xtrain = lapply(X, function(x) {
        x[, cv.fold != cv, drop = F]
      })
      ytest  = lapply(Ys, function(y) {
        y[, cv.fold == cv, drop = F]
      })
      xtest  = lapply(X, function(x) {
        x[, cv.fold == cv, drop = F]
      })
      Ntrain = sum(cv.fold != cv)
      fit    = constrained_L2reg_multi(xtrain, ytrain, rho, M, Ntrain)
      for (k in 1:length(Ys)) {
        ftest  = lapply(1:M, function(i) {
          crossprod(fit$F[[k]][[i]], xtest[[i]])
        })
        ftest  = do.call(rbind, ftest)
        cv.err[irho, cv] = cv.err[irho, cv] + norm((diag(M) - fit$B[[k]]) %*% ytest[[k]] - ftest - fit$mu[[k]], type = "f")
      }
    }
    irho = irho + 1
  }
  cv.mean = rowMeans(cv.err)
  rho.min = rho_factors[which.min(cv.mean)]
  fit = constrained_L2reg_multi(X, Ys, rho.min, M, N)
  list(
    rho.opt = rho.min,
    sigma2.opt = fit$sigma2[1],
    cv.ram = list(rho = rho_factors, cvm = cv.mean)
  )
}


##---------------------------------------
# utility functions for SML-fused_lasso #
##---------------------------------------
#' @param lambda lasso penalty
#' @param rho    fused lasso penalty
#' @param w      weight for lasso term
#' @param r      weight for fused lasso term
#' @param c      ci / det(B)
#' @param b      bij[1,...,K]
obj_multiSML = function(N, c, b, a1, a2, lambda, rho, w, r, sigma2) {
  a0 = 1 / c + b
  if (c[1] == 0 & c[2] == 0) {
    function(x, y) {
      -a1[1] * x - a1[2] * y +
        1 / 2 * (a2[1] * x ^ 2 + a2[2] * y ^ 2) +
        lambda * (w[1] * abs(x) +  w[2] * abs(y)) +
        rho * r * (abs(x - y))
    }
  } else if (c[1] == 0 & c[2] != 0) {
    function(x, y) {
      sigma2 * (-N[2] / 2 * log((a0[2] - y) ^ 2)) -
        a1[1] * x - a1[2] * y +
        1 / 2 * (a2[1] * x ^ 2 + a2[2] * y ^ 2) +
        lambda * (w[1] * abs(x) +  w[2] * abs(y))  +
        rho * r * (abs(x - y))
    }
  } else if (c[1] != 0 & c[2] == 0) {
    function(x, y) {
      sigma2 * (-N[1] / 2 * log((a0[1] - x) ^ 2)) -
        a1[1] * x - a1[2] * y +
        1 / 2 * (a2[1] * x ^ 2 + a2[2] * y ^ 2) +
        lambda * (w[1] * abs(x) +  w[2] * abs(y))  +
        rho * r * (abs(x - y))
    }
  } else {
    function(x, y) {
      sigma2 * (-N[1] / 2 * log((a0[1] - x) ^ 2) - N[2] / 2 * log((a0[2] - y) ^
                                                                    2)) -
        a1[1] * x - a1[2] * y +
        1 / 2 * (a2[1] * x ^ 2 + a2[2] * y ^ 2) +
        lambda * (w[1] * abs(x) +  w[2] * abs(y)) +
        rho * r * (abs(x - y))
    }
  }
}

grad_multiSML = function(N, c, b, a1, a2, lambda, rho, w, r, sigma2) {
  a0 = 1 / c + b
  if (c[1] == 0 & c[2] == 0) {
    function(t) {
      xt  = t[1]
      yt  = t[2]
      dxy = t[3]
      if (dxy != 0) {
        x = if (xt != 0) {
          structure((a1[1] - lambda * w[1] * xt - rho * r * dxy) / a2[1], class = "value")
        } else {
          structure(0, class = "value")
        }
        y = if (yt != 0) {
          structure((a1[2] - lambda * w[2] * yt + rho * r * dxy) / a2[2], class = "value")
        } else {
          structure(0, class = "value")
        }
        list(x = x, y = y)
      } else {
        tau = xt
        wxy = w[1] + w[2]
        xy = if (tau != 0) {
          structure((a1[1] + a1[2] - lambda * wxy * tau) / (a2[1] + a2[2]), class = "value")
        } else {
          structure(0, class = "value")
        }
        list(xy = xy)
      }
    }
  } else if (c[1] != 0 & c[2] == 0) {
    function(t) {
      xt  = t[1]
      yt  = t[2]
      dxy = t[3]
      if (dxy != 0) {
        x = if (xt != 0) {
          structure(
            list(
              a = -a2[1],
              b = a1[1] + a2[1] * a0[1] - lambda * w[1] * xt - rho * r * dxy,
              c = N[1] * sigma2 + (lambda * w[1] * xt + rho * r * dxy - a1[1]) * a0[1]
            ),
            class = "grad2"
          )
        } else {
          structure(0, class = "value")
        }
        y = if (yt != 0) {
          structure((a1[2] - lambda * w[2] * yt + rho * r * dxy) / a2[2], class = "value")
        } else {
          structure(0, class = "value")
        }
        list(x = x, y = y)
      } else {
        tau = xt
        wxy = w[1] + w[2]
        xy = if (tau != 0) {
          structure(list(
            a = -(a2[1] + a2[2]),
            b = (a1[1] + a1[2]) + (a2[1] + a2[2]) * a0[1] - lambda * wxy * tau,
            c = N[1] * sigma2 + (lambda * wxy * tau - a1[1] - a1[2])  * a0[1]
          ),
          class = "grad2")
        } else {
          structure(0, class = "value")
        }
        list(xy = xy)
      }
    }
  } else if (c[1] == 0 & c[2] != 0) {
    function(t) {
      xt  = t[1]
      yt  = t[2]
      dxy = t[3]
      if (dxy != 0) {
        x = if (xt != 0) {
          structure((a1[1] - lambda * w[1] * xt - rho * r * dxy) / a2[1], class = "value")
        } else {
          structure(0, class = "value")
        }
        y = if (yt != 0) {
          structure(
            list(
              a = -a2[2],
              b = a1[2] + a2[2] * a0[2] - lambda * w[2] * yt + rho * r * dxy,
              c = N[2] * sigma2 + (lambda * w[2] * yt - rho * r * dxy - a1[2]) * a0[2]
            ),
            class = "grad2"
          )
        } else {
          structure(0, class = "value")
        }
        list(x = x, y = y)
      } else {
        tau = xt
        wxy = w[1] + w[2]
        xy = if (tau != 0) {
          structure(list(
            a = -(a2[1] + a2[2]),
            b = (a1[1] + a1[2]) + (a2[1] + a2[2]) * a0[2] - lambda * wxy * tau,
            c = N[1] * sigma2 + (lambda * wxy * tau - a1[1] - a1[2])  * a0[2]
          ),
          class = "grad2")
        } else {
          structure(0, class = "value")
        }
        list(xy = xy)
      }
    }
  } else {
    function(t) {
      xt  = t[1]
      yt  = t[2]
      dxy = t[3]
      if (dxy != 0) {
        x = if (xt != 0) {
          structure(
            list(
              a = -a2[1],
              b = a1[1] + a2[1] * a0[1] - lambda * w[1] * xt - rho * r * dxy,
              c = N[1] * sigma2 + (lambda * w[1] * xt + rho * r * dxy - a1[1]) * a0[1]
            ),
            class = "grad2"
          )
        } else {
          structure(0, class = "value")
        }
        y = if (yt != 0) {
          structure(
            list(
              a = -a2[2],
              b = a1[2] + a2[2] * a0[2] - lambda * w[2] * yt + rho * r * dxy,
              c = N[2] * sigma2 + (lambda * w[2] * yt - rho * r * dxy - a1[2]) * a0[2]
            ),
            class = "grad2"
          )
        } else {
          structure(0, class = "value")
        }
        list(x = x, y = y)
      } else {
        tau = xt
        wxy = w[1] + w[2]
        xy = if (tau != 0) {
          structure(
            list(
              a = a2[1] + a2[2],
              b = lambda * wxy * tau - (a1[1] + a1[2]) - (a2[1] + a2[2]) * (a0[1] + a0[2]),
              c = (a1[1] + a1[2] - lambda * wxy * tau) * (a0[1] + a0[2]) + (a2[1] + a2[2]) * a0[1] * a0[2] - (N[1] + N[2]) * sigma2,
              d = (N[2] * a0[1] + N[1] * a0[2]) * sigma2 + (lambda * wxy * tau - (a1[1] + a1[2])) * a0[1] * a0[2]
            ),
            class = "grad3"
          )
        } else {
          structure(0, class = "value")
        }
        list(xy = xy)
      }
    }
  }
}

poly3_solver = function(a, b, c, d, eps = 1e-6) {
  r = solver3P_(b / a, c / a, d / a)
  r$x = r$x * (abs(r$x) > eps)
  r
}

# x^3 + ax^2 + bx + c = 0
solver3P_ = function(a, b, c) {
  a2 = a ^ 2
  q  = (a2 - 3 * b) / 9
  r  = (a * (2 * a2 - 9 * b) + 27 * c) / 54
  r2 = r ^ 2
  q3 = q ^ 3
  if (r2 <= q3) {
    t = r / sqrt(q3)
    if (t < -1) {
      t = -1
    }
    if (t > 1) {
      t = 1
    }
    t = acos(t)
    a = a / 3
    q = -2 * sqrt(q)
    list(n = 3, x = c(q * cos(t / 3) - a, q * cos((t + 2 * pi) / 3) - a, q * cos((t - 2 * pi) / 3) - a))
  } else {
    4
    A = -(abs(r) + sqrt(r2 - q3)) ^ (1 / 3)
    if (r < 0) {
      A = -A
    }
    B = if (A == 0) {
      0
    } else {
      q / A
    }

    a = a / 3
    Re = -(A + B) / 2 - a
    Im = sqrt(3) / 2 * (A - B)
    if (abs(Re) <= 1e-6) {
      Im = Re
      list(n = 2, x = c(A + B - a, Re))
    } else {
      list(n = 1, x = c(A + B - a))
    }
  }
}

### call function for grad list
grad_solver = function(g, t) {
  gsolver_ = function(g) {
    if (class(g) == "value") {
      list(n = 1, x = 0)
    } else if (class(g) == "grad2") {
      do.call(poly2_solver, g)
    } else {
      do.call(poly3_solver, g)
    }
  }
  xt  = t[1]
  yt  = t[2]
  dxy = t[3]
  res = list()
  if (dxy != 0) {
    cand.x = gsolver_(g[["x"]])
    cand.x = cand.x[["x"]][sign(cand.x[["x"]]) == xt]
    cand.y = gsolver_(g[["y"]])
    cand.y = cand.y[["x"]][sign(cand.y[["x"]]) == yt]
    i = 1
    for (x in cand.x) {
      for (y in cand.y) {
        if (sign(x - y) == dxy) {
          res[[i]] = list(x = x, y = y)
          i = i + 1
        }
      }
    }
  } else {
    cand.xy = gsolver_(g[["xy"]])
    cand.xy = cand.xy[["x"]][sign(cand.xy[["x"]]) == xt]
    i = 1
    for (xy in cand.xy) {
      res[[i]] = list(x = xy, y = xy)
      i = i + 1
    }
  }
  res
}


## solve SML problem by component-wise update
#' @param lambda lasso penalty
#' @param rho    fused lasso penalty
#' @example
#' M = data$var$Ng
#' N = data$var$N
#' Ys = data$obs[c("Y1", "Y2")]
#' Xs = submatX(data)
#' sigma2 = getsigma2_L2reg_multi(Xs, Ys, nrho = 20, M = M, N = N)
#' params.init = constrained_L2reg_multi(Xs, Ys, sigma2$rho.opt, M, N)
#' params.opt = multiSML_cwise(Bs = params.init$B, fs = params.init$F, Ys = Ys, Xs = Xs,
#'                             sigma2 = params.init$sigma2[1], Ng = data$var$Ng,
#'                             lambda = 0.1, rho = 0.1, maxit = 1000)
multiSML_cwise = function(Bs,
                          fs,
                          Ys,
                          Xs,
                          sigma2,
                          Ng,
                          lambda,
                          rho,
                          weighted = TRUE,
                          wBs = inverse(Bs),
                          rB  = flinv(Bs),
                          maxit = 100,
                          threshold = 1e-4,
                          verbose = 2) {
  std = centralize_mult(Xs, Ys)
  Xs  = std$X
  Ys  = std$Y
  meanXs = std$muX
  meanYs = std$muY
  K   = length(Ys)  # number of conditions = 2
  if (verbose == 2) {
    cat(sprintf("conditions must be restricted to 2.. K = %d\n", K))
  }
  ## update for eQTL coeffs
  f0 = vector("list", K)
  f1 = vector("list", K)
  for (i in 1:Ng) {
    Xi = Xs[[i]]
    Pi = solve(crossprod(Xi)) %*% t(Xi)
    for (k in 1:K) {
      yi_k = Ys[[k]][, i, drop = F]     # n x 1 for gene i
      f0[[k]][[i]] = Pi %*% yi_k
      f1[[k]][[i]] = Pi %*% Ys[[k]]   # f = f0 - f1 %*% B[i,]
    }
  }
  ## update for gnet coeffs
  niter  = 1
  Ns     = sapply(Ys, nrow)
  ImBs   = lapply(Bs, function(B) {
    diag(Ng) - B
  })
  IBsinv = lapply(ImBs, solve)
  while (niter <= maxit) {
    Bs.prev = Bs
    fs.prev = fs
    for (i in 1:Ng) {
      ci  = lapply(IBsinv, function(IBi) {
        IBi[, i]
      })
      dbi = lapply(1:K, function(k) {
        vector("numeric", Ng)
      })
      for (j in 1:Ng) {
        ## update B[[k]][i,j] for i != j
        if (i != j) {
          bij.prev = bij = sapply(Bs, function(B)
            (B[i, j]))
          wij      = sapply(wBs, function(w) {
            if (weighted) {
              w[i, j]
            } else {
              1
            }
          })
          rij     = if (weighted) {
            rB[i, j]
          } else {
            1
          }
          mij     = sapply(ci, function(c) {
            c[j]
          })
          bi      = lapply(ImBs, function(ImB) {
            bi_k = ImB[i, ]
            bi_k[j] = 0
            bi_k
          })
          ## j-th column of Ys
          Yej    = lapply(Ys, function(Y) {
            Y[, j, drop = F]
          })
          a1     = sapply(1:K, function(k) {
            crossprod(Ys[[k]] %*% bi[[k]] - Xs[[i]] %*% fs[[k]][[i]], Yej[[k]])
          })
          a2     = sapply(1:K, function(k) {
            crossprod(Yej[[k]])
          })
          ## a0 = 1/mij + bij.prev
          cond   = list(
            c(1, 1, 1),
            c(1,-1, 1),
            c(1, 0, 1),
            c(-1,-1, 1),
            c(0,-1, 1),
            c(1, 1,-1),
            c(0, 1,-1),
            c(-1,-1,-1),
            c(-1, 0,-1),
            c(-1, 1,-1),
            c(1, 1, 0),
            c(-1,-1, 0),
            c(0, 0, 0)
          )
          obj    = obj_multiSML(Ns, mij, bij.prev, a1, a2, lambda, rho, wij, rij, sigma2)
          grad   = grad_multiSML(Ns, mij, bij.prev, a1, a2, lambda, rho, wij, rij, sigma2)
          params = list()
          for (t in cond) {
            cand.grad = grad(t)
            params    = c(params, grad_solver(cand.grad, t))
          }
          objval = sapply(params, function(args) {
            do.call(obj, args)
          })
          mix    = which.min(objval)
          bij    = unlist(params[[mix]])
          dbij = bij.prev - bij
          for (k in 1:K) {
            dbi[[k]][j]  = dbij[k]
            Bs[[k]][i, j] = bij[k]
            ci[[k]]      = ci[[k]] / (1 + dbij[k] * mij[k])
            ImBs[[k]]    = diag(Ng) - Bs[[k]]
          }
        }
      } ## for(j in 1:Ng)
      ## (ImB + ei^T %*% dbi)^{-1}
      for (k in 1:K) {
        ## IBsinv[[k]] = IBsinv[[k]] - IBsinv[[k]][, i, drop = F] %*% dbi[[k]] %*% IBsinv[[k]] / (1 + dbi[[k]] %*% IBsinv[[k]][, i, drop = F])[1]
        IBsinv[[k]] = solve(ImBs[[k]])
        fs[[k]][[i]]  = f0[[k]][[i]] - f1[[k]][[i]] %*% Bs[[k]][i,]
      }
    } ## for(i in 1:Ng)
    Berr = sum(sapply(1:K, function(k) {
      norm(Bs[[k]] - Bs.prev[[k]], type = "f") / norm(Bs.prev[[k]], type = "f")
    }))
    Ferr = sum(sapply(1:K, function(k) {
      sum(sapply(1:Ng, function(i) {
        norm(fs[[k]][[i]] - fs.prev[[k]][[i]], type = "f")
      })) / sum(sapply(1:Ng, function(i) {
        norm(fs.prev[[k]][[i]], type = "f")
      }))
    }))
    err = Berr + Ferr
    cat(sprintf("SML: iteration = %d, error = %f\n", niter, err))
    niter = niter + 1
    if (err < threshold || niter > maxit || is.nan(err)) {
      mu = lapply(1:K, function(k) {
        (diag(Ng) - Bs[[k]]) %*% meanYs[[k]] - sapply(1:Ng, function(i) {
          meanXs[[i]] %*% fs[[k]][[i]]
        })
      })
      Bs = lapply(Bs, Matrix, sparse = T)
      break
    }
  } ## while(niter <= maxit)
  list(
    B = Bs,
    f = fs,
    mu = mu,
    niter = niter,
    err = err
  )
}


## ultility function for multiple condition
#' @param Bs multiple gene regulatory networks (list)
#' @param fs multiple eQTL coefficient (list)
#' @param Ys multiple gene expression matrix (list)
#' @param Xs multiple eQTLs quantitive (list)
#' @param Ng number of genes
#' @param lambda lambda coefficient in lasso penalty term
#' @param rho rho coefficient in fused penalty term
#' @param type type of likelihood:
#'             o  objective function
#'             c  cross validation function
#'             e  independent error function
SML_error = function(Xs, Ys, Bs, fs, Ng, Ns, K) {
  std = centralize_mult(Xs, Ys)
  X  = lapply(std$X, t)
  Y  = lapply(std$Y, t)
  err = 0
  for (k in 1:K) {
    for (i in 1:Ng) {
      Xi = X[[i]]                     # sk x N
      bi = Bs[[k]][i, -i, drop = F]   # (ng-1) x 1
      yi = Y[[k]][i, , drop = F]      # 1 x N
      Yi = Y[[k]][-i, , drop = F]     # (ng-1) x N
      fi = fs[[k]][[i]]               # sk x 1
      err = err + tcrossprod(yi - bi %*% Yi - crossprod(fi, Xi))
    }
  }
  err
}

SML_logLik = function(Xs, Ys, Bs, fs, Ng, Ns, K, detIBs, sigma2) {
  std = centralize_mult(Xs, Ys)
  X  = lapply(std$X, t)
  Y  = lapply(std$Y, t)
  Ls  = 0
  err = 0
  for (k in 1:K) {
    Ls = Ls - Ns[k] / 2 * log(detIBs[k] ^ 2)
    for (i in 1:Ng) {
      Xi = X[[i]]                     # sk x N
      bi = Bs[[k]][i, -i, drop = F]   # (ng-1) x 1
      yi = Y[[k]][i, , drop = F]      # 1 x N
      Yi = Y[[k]][-i, , drop = F]     # (ng-1) x N
      fi = fs[[k]][[i]]               # sk x 1
      err = err + tcrossprod(yi - bi %*% Yi - crossprod(fi, Xi))
    }
  }
  Ls + err / (2 * sigma2) + Ng * sum(Ns) / 2 * log(sigma2)
}

#' @description proximal operator for fused lasso
#' argmin lambda * |w1 * x| + lambda * |w2 * y| + rho * r * |y - x| + c/2 (|x - u1|_2^2 + |y - u_2|_2^2)
#' @param lambda lasso parameter
#' @param rho    fused lasso parameter
#' @note FLSA algorithm is used in this step
prox_flsa = function(lambda, rho, c, us, ws, r) {
  ## lambda = 0
  du = us[[1]] - us[[2]]
  eq = (abs(du) <= 2 * rho * r / c)
  df = 1 - eq
  rho = min(rho, 1e16)
  x = list((us[[1]] + us[[2]]) / 2 * eq + (us[[1]] - sign(du) * rho * r / c) * df,
           (us[[1]] + us[[2]]) / 2 * eq + (us[[2]] + sign(du) * rho * r / c) * df)
  x = lapply(x, as.numeric)
  lapply(1:2, function(i) {
    xe = pmax(x[[i]] - lambda * (ws[[1]] + ws[[2]]) / 2 / c, 0) + pmin(x[[i]] + lambda * (ws[[1]] + ws[[2]]) / 2 / c, 0)
    xd = pmax(x[[i]] - lambda * ws[[i]] / c, 0) + pmin(x[[i]] + lambda * ws[[i]] / c, 0)
    xe * eq + xd * df
  })
}


## sigma2 estimation from SEM logLikihood function
sigma2_sem = function(Xs, Ys, B, f, Ng, Ns, K) {
  X = lapply(Xs, t)
  Y = lapply(Ys, t)
  err = 0
  for (k in 1:K) {
    for (i in 1:Ng) {
      Xi = X[[i]]                     # sk x N
      bi = B[[k]][i, -i, drop = F]     # (ng-1) x 1
      yi = Y[[k]][i, , drop = F]      # 1 x N
      Yi = Y[[k]][-i, , drop = F]     # (ng-1) x N
      fi = f[[k]][[i]]                # sk x 1
      err = err + tcrossprod(yi - bi %*% Yi - crossprod(fi, Xi))
    }
  }
  sigma2 = err / (Ng * sum(Ns))
  sigma2[1]
}

## solve SML problem by block coordinate descent
#' @param lambda lambda hyper-parameter for lasso term
#' @param rho    fused lasso hyper-parameter for fused lasso term
#' @param gamma  invertible matrix stablize parameter gamma
#' M = data$var$Ng
#' N = data$var$N
#' Ys = data$obs[c("Y1", "Y2")]
#' Xs = submatX(data)
#' sigma2 = getsigma2_L2reg_multi(Xs, Ys, nrho = 20, M = M, N = N)
#' params.init = constrained_L2reg_multi(Xs, Ys, sigma2$rho.opt, M, N)
#' params.opt2 = multiSML_bcd(Bs = params.init$B, fs = params.init$F, Ys = Ys, Xs = Xs,
#'                            sigma2 = params.init$sigma2[1], Ng = data$var$Ng,
#'                            lambda = 0.1, rho = 0.1, maxit = 500, gamma = sigma2$rho.opt)
multiSML_bcd = function(Bs,
                        fs,
                        Ys,
                        Xs,
                        sigma2,
                        Ng,
                        lambda,
                        rho,
                        weighted = TRUE,
                        maxit = 100,
                        threshold = 1e-3,
                        verbose = 2) {
  std = centralize_mult(Xs, Ys)
  Xs  = std$X
  Ys  = std$Y
  meanXs = std$muX
  meanYs = std$muY
  K      = length(Ys)      # number of conditions = 2
  if (verbose == 2) {
    cat(sprintf("conditions must be restricted to 2.. K = %d\n", K))
  }
  ## update for eQTL row-wise (specific for genes)
  f0 = vector("list", K)
  f1 = vector("list", K)
  Yp = vector("list", K)
  Hy = vector("list", K)
  Yp.maxEigen = vector("list", K)
  Ns = sapply(Ys, nrow)
  for (i in 1:Ng) {
    Xi = Xs[[i]]
    Pi = solve(crossprod(Xi)) %*% t(Xi)
    for (k in 1:K) {
      # specific condition
      yi_k = Ys[[k]][, i, drop = F]              # n[k] x 1 for gene i
      Yi_k = Ys[[k]][, -i]                       # n[k] x (ng-1) (for specific gene i)
      f0[[k]][[i]] = Pi %*% yi_k
      f1[[k]][[i]] = Pi %*% Yi_k                 # f[[k]][[i]] = f0[[k]][[i]] - f1[[k]][[i]] %*% bi^k | bi^k = B[[k]][i,-i]
      Hi_k = diag(Ns[k]) - Xi %*% Pi             # n[k] x n[k] projection matrix
      Yp[[k]][[i]] = t(Yi_k) %*% Hi_k %*% Yi_k   # (ng - 1) x (ng - 1)
      Hy[[k]][[i]] = t(yi_k) %*% Hi_k %*% Yi_k   # 1 x (ng - 1)
      ## maximized eigen-value for Yp
      Yp.maxEigen[[k]][[i]] = eigen(Yp[[k]][[i]])$values[1]
    }
  }
  ## update for gnet row-wise
  niter  = 1
  ImBs   = lapply(Bs, function(B) {
    diag(Ng) - B
  })
  detIBs = sapply(ImBs, det)
  IBsinv = lapply(ImBs, solve)
  wBs    = lapply(Bs, function(B) {
    1 / abs(B)
  })
  rB     = 1 / abs(Bs[[1]] - Bs[[2]])
  Ls     = logLik(detIBs, Bs, wBs, rB, lambda, rho, K) + Ng * sum(Ns) / 2 * log(sigma2)
  while (niter <= maxit) {
    Bs.prev = Bs
    fs.prev = fs
    Ls.prev = Ls
    for (i in 1:Ng) {
      ## sum(-Ns[k] * sigma2 * log(det(I-B[[k]])^2)) + \sum_{i=1}^{Ng} bi^T %*% Yp %*% bi - Hy %*% bi)
      ci = lapply(IBsinv, function(IBi) {
        # ci / det(I - B); from (I - B)^{-1}
        IBi[-i, i, drop = F]
      })
      bi = lapply(Bs, function(B) {
        t(B[i, -i, drop = F])
      })
      gi = lapply(1:K, function(k) {
        grad = grad_rwise_SML(Ns[k], ci[[k]], Yp[[k]][[i]], Hy[[k]][[i]], sigma2[1])
        grad(bi[[k]])
      })
      ## Lipschitz moduli for row-i
      Lis  = sapply(1:K, function(k) {
        gtg  = tcrossprod(ImBs[[k]][-i,])
        oi   = chol2inv(gtg)
        deti = det(gtg)
        gii  = ImBs[[k]][-i, -i]
        si   = ImBs[[k]][-i, i, drop = F]
        c2i  = sum((ci[[k]] * detIBs[k]) ^ 2)
        lips_rwise_SML(Ns[k],
                       oi,
                       gii,
                       si,
                       c2i,
                       deti,
                       Yp.maxEigen[[k]][[i]],
                       sigma2,
                       Ng)[1]
      })
      Li   = max(Lis)
      ui   = lapply(1:K, function(k) {
        bi[[k]] - gi[[k]] / Li
      })
      wBi  = lapply(wBs, function(wB) {
        wB[i, -i]
      })
      rBi  = rB[i, -i]
      xi   = prox_flsa(lambda, rho, Li, ui, wBi, rBi)
      for (k in 1:K) {
        Bs[[k]][i, -i] = xi[[k]]
        dbi = Bs.prev[[k]][i,] - Bs[[k]][i,]
        IBsinv[[k]] = IBsinv[[k]] - IBsinv[[k]][, i, drop = F] %*% dbi %*% IBsinv[[k]] / (1 + dbi %*% IBsinv[[k]][, i, drop = F])[1]
        ImBs[[k]]   = diag(Ng) - Bs[[k]]
        detIBs[k]   = (ImBs[[k]][i,] %*% IBsinv[[k]][, i, drop = F])[1] * detIBs[k]
        ## IBsinv[[k]] = solve(ImBs[[k]])
        fs[[k]][[i]] = f0[[k]][[i]] - f1[[k]][[i]] %*% Bs[[k]][i, -i]
      }
    } # row-wise update
    Berr = sum(sapply(1:K, function(k) {
      norm(Bs[[k]] - Bs.prev[[k]], type = "f") / norm(Bs.prev[[k]], type = "f")
    }))
    Ferr = sum(sapply(1:K, function(k) {
      sum(sapply(1:Ng, function(i) {
        norm(fs[[k]][[i]] - fs.prev[[k]][[i]], type = "f")
      })) / sum(sapply(1:Ng, function(i) {
        norm(fs.prev[[k]][[i]], type = "f")
      }))
    }))
    err = Berr + Ferr
    sigma2 = sigma2_sem(Xs, Ys, Bs, fs, Ng, Ns, K)
    Ls     = logLik(detIBs, Bs, wBs, rB, lambda, rho, K) + Ng * sum(Ns) / 2 * log(sigma2)
    Lerr   = abs(Ls.prev - Ls) / (1 + abs(Ls.prev))
    if (verbose >= 2) {
      cat(
        sprintf(
          "SML: iteration = %d, error = %f, logLik = %f, sigma2 = %f\n",
          niter,
          err,
          Ls,
          sigma2
        )
      )
    }
    niter = niter + 1
    if ((err <= threshold & Lerr <= threshold) || niter > maxit) {
      mu = lapply(1:K, function(k) {
        (diag(Ng) - Bs[[k]]) %*% meanYs[[k]] - sapply(1:Ng, function(i) {
          meanXs[[i]] %*% fs[[k]][[i]]
        })
      })
      sigma2 = sigma2_sem(Xs, Ys, Bs, fs, Ng, Ns, K)
      Bs = lapply(Bs, Matrix, sparse = T)
      break
    }
  } # while(niter <= maxit)
  list(
    B = Bs,
    f = fs,
    mu = mu,
    sigma2 = sigma2,
    niter = niter,
    err = err,
    detIB = detIBs
  )
}

## stepwise update sigma2
err2_sem = function(Xs, Ys, B, f, Ng, Ns, K) {
  X = lapply(Xs, t)
  Y = lapply(Ys, t)
  err2 = 0
  for (k in 1:K) {
    for (i in 1:M) {
      Xi = X[[i]]                     # sk x N
      bi = B[[k]][i, -i, drop = F]     # (ng-1) x 1
      yi = Y[[k]][i, , drop = F]      # 1 x N
      Yi = Y[[k]][-i, , drop = F]     # (ng-1) x N
      fi = f[[k]][[i]]                # sk x 1
      err2 = err2 + tcrossprod(yi - bi %*% Yi - crossprod(fi, Xi))
    }
  }
  err2
}


## ultility function for objective function
logLik = function(detIBs, Bs, ws, r, lambda, rho, Ns, K) {
  Ls     = 0
  l1norm = 0
  rho    = min(rho, 1e16)
  lfnorm = rho * r * abs(Bs[[2]] - Bs[[1]])
  for (k in 1:K) {
    l1norm = l1norm + lambda * (ws[[k]] * abs(Bs[[k]]))
    Ls = Ls - Ns[k] / 2 * log(detIBs[k] ^ 2)
  }
  diag(l1norm) = 0
  diag(lfnorm) = 0
  Ls + sum(l1norm) + sum(lfnorm)
}

## inverse
inverse = function(Bs) {
  lapply(Bs, function(B) {
    1 / abs(B)
  })
}

## invone
invone = function(Bs) {
  lapply(Bs, function(B) {
    w = matrix(1, nrow = nrow(B), ncol = ncol(B))
    diag(w) = Inf
    w
  })
}

## flinv
flinv = function(Bs) {
  1 / abs(Bs[[1]] - Bs[[2]])
}

## flone
flone = function(Bs) {
  w = matrix(1, nrow = nrow(Bs[[1]]), ncol = ncol(Bs[[2]]))
  diag(w) = Inf
  w
}

## solve SML problem by block coordinate descent by backtracking inert-PALM
#' @param lambda lambda hyper-parameter for lasso term
#' @param rho    fused lasso hyper-parameter for fused lasso term
#' @param gamma  invertible matrix stablize parameter gamma
#' params.opt4 = multiSML_iPALM(Bs = params.init$B, fs = params.init$F, Ys = Ys, Xs = Xs,
#'                            sigma2 = params.init$sigma2[1], Ng = data$var$Ng,
#'                            lambda = 0.1, rho = 0.1, maxit = 500)
multiSML_iPALM = function(Bs,
                          fs,
                          Ys,
                          Xs,
                          sigma2,
                          Ng,
                          lambda,
                          rho,
                          wBs = inverse(Bs),
                          rB  = flinv(Bs),
                          maxit = 100,
                          acc = TRUE,
                          inertial = inertial_pars("lin"),
                          threshold = 1e-3,
                          sparse = FALSE,
                          use.strict = TRUE,
                          verbose = 2) {
  std = centralize_mult(Xs, Ys)
  Xs  = std$X
  Ys  = std$Y
  meanXs = std$muX
  meanYs = std$muY
  K      = length(Ys)      # number of conditions = 2
  if (verbose == 2) {
    cat(sprintf("conditions must be restricted to 2.. K = %d\n", K))
  }
  ## update for eQTL row-wise (specific for genes)
  f0 = vector("list", K)
  f1 = vector("list", K)
  Yp = vector("list", K)
  Hy = vector("list", K)
  Yp.maxEigen = vector("list", K)
  Ns = sapply(Ys, nrow)
  for (i in 1:Ng) {
    Xi = Xs[[i]]
    Pi = solve(crossprod(Xi)) %*% t(Xi)
    for (k in 1:K) {
      # specific condition
      yi_k = Ys[[k]][, i, drop = F]              # n[k] x 1 for gene i
      Yi_k = Ys[[k]][, -i]                       # n[k] x (ng-1) (for specific gene i)
      f0[[k]][[i]] = Pi %*% yi_k
      f1[[k]][[i]] = Pi %*% Yi_k                 # f[[k]][[i]] = f0[[k]][[i]] - f1[[k]][[i]] %*% bi^k | bi^k = B[[k]][i,-i]
      Hi_k = diag(Ns[k]) - Xi %*% Pi             # n[k] x n[k] projection matrix
      Yp[[k]][[i]] = t(Yi_k) %*% Hi_k %*% Yi_k   # (ng - 1) x (ng - 1)
      Hy[[k]][[i]] = t(yi_k) %*% Hi_k %*% Yi_k   # 1 x (ng - 1)
      ## maximized eigen-value for Yp
      Yp.maxEigen[[k]][[i]] = eigen(Yp[[k]][[i]])$values[1]
    }
  }
  ## update for gnet row-wise
  niter  = 1
  ImBs   = lapply(Bs, function(B) {
    diag(Ng) - B
  })
  detIBs = sapply(ImBs, det)
  IBsinv = lapply(ImBs, solve)
  Ls     = logLik(detIBs, Bs, wBs, rB, lambda, rho, Ns, K) + N * sum(Ns) / 2 * log(sigma2)
  ## history
  Bs.prevs = list(Bs, Bs)
  inert    = acc
  while (niter <= maxit) {
    inert.pars = inertial(niter)
    Bs.inert = if (inert) {
      lapply(1:K, function(k) {
        Bs.prevs[[2]][[k]] + inert.pars * (Bs.prevs[[2]][[k]] - Bs.prevs[[1]][[k]])
      })
    } else {
      Bs
    }
    fs.prev = fs
    Ls.prev = Ls
    for (i in 1:Ng) {
      ## sum(-Ns[k] * sigma2 * log(det(I-B[[k]])^2)) + \sum_{i=1}^{Ng} bi^T %*% Yp %*% bi - Hy %*% bi)
      ci = lapply(IBsinv, function(IBi) {
        # ci / det(I - B); from (I - B)^{-1}
        IBi[-i, i, drop = F]
      })
      bi = lapply(Bs.inert, function(B.inert) {
        t(B.inert[i, -i, drop = F])
      })
      gi = lapply(1:K, function(k) {
        grad = grad_rwise_SML(Ns[k], ci[[k]], Yp[[k]][[i]], Hy[[k]][[i]], sigma2[1])
        grad(bi[[k]])
      })
      ## Lipschitz moduli for row-i
      Lis  = sapply(1:K, function(k) {
        gtg  = tcrossprod(ImBs[[k]][-i,])
        oi   = chol2inv(gtg)
        deti = crossprod(IBsinv[[k]][,i])[1] * (detIBs[k]^2)
        gii  = ImBs[[k]][-i, -i]
        si   = ImBs[[k]][-i, i, drop = F]
        c2i  = sum((ci[[k]] * detIBs[k]) ^ 2)
        Li   = lips_rwise_SML(Ns[k],
                       oi,
                       gii,
                       si,
                       c2i,
                       deti,
                       Yp.maxEigen[[k]][[i]],
                       sigma2,
                       Ng)[1]
        Li
      })
      Li   = max(Lis)
      Li   = (1 + 2 * inert.pars) * Li / (2 * (1 - inert.pars))
      detZero = TRUE
      cl      = 1
      while(detZero) {
        ui   = lapply(1:K, function(k) {
          bi[[k]] - gi[[k]] / (cl * Li)
        })
        wBi  = lapply(wBs, function(wB) {
          wB[i, -i]
        })
        rBi  = rB[i, -i]
        xi   = prox_flsa(lambda, rho, Li, ui, wBi, rBi)
        dIBu  = sapply(1:K, function(k) {
          IBsinv[[k]][i, i] - (t(xi[[k]]) %*% IBsinv[[k]][-i, i, drop = F])[1]
        })
        cl = cl * 2
        detZero = any(dIBu == 0)
      }
      for (k in 1:K) {
        Bs[[k]][i, -i] = xi[[k]]
        ImBs[[k]]      = diag(Ng) - Bs[[k]]
        detIBs[k]      = (ImBs[[k]][i,] %*% IBsinv[[k]][, i, drop = F])[1] * detIBs[k]
        # detIBs[k]      = det(ImBs[[k]])
        dbi            = Bs.prevs[[2]][[k]][i,] - Bs[[k]][i,]
        # IBsinv[[k]]    = solve(ImBs[[k]])
        IBsinv[[k]]    = IBsinv[[k]] - IBsinv[[k]][, i, drop = F] %*% dbi %*% IBsinv[[k]] / (1 + dbi %*% IBsinv[[k]][, i, drop = F])[1]
        fs[[k]][[i]]   = f0[[k]][[i]] - f1[[k]][[i]] %*% Bs[[k]][i, -i]
      }
      # sigma2 = sigma2_sem(Xs, Ys, Bs, fs, Ng, Ns, K)
    } # row-wise update
    Berr = sum(sapply(1:K, function(k) {
      norm(Bs[[k]] - Bs.prevs[[2]][[k]], type = "f") / (1 + norm(Bs.prevs[[2]][[k]], type = "f"))
    }))
    Ferr = sum(sapply(1:K, function(k) {
      sum(sapply(1:Ng, function(i) {
        norm(fs[[k]][[i]] - fs.prev[[k]][[i]], type = "f")
      })) / (1 + sum(sapply(1:Ng, function(i) {
        norm(fs.prev[[k]][[i]], type = "f")
      })))
    }))
    err = Berr + Ferr
    sigma2 = sigma2_sem(Xs, Ys, Bs, fs, Ng, Ns, K)
    Ls     = logLik(detIBs, Bs, wBs, rB, lambda, rho, Ns, K) + Ng * sum(Ns) / 2 * log(sigma2)
    Lerr   = abs(Ls.prev - Ls) / (1 + abs(Ls.prev))
    # inert  = ifelse(Lerr < 1e-6, FALSE, acc)
    if (verbose >= 2) {
      cat(
        sprintf(
          "SML: iteration = %d, error = %f, logLik = %f, sigma2 = %f, inert=%s\n",
          niter,
          err,
          Ls,
          sigma2,
          inert
        )
      )
    }
    niter = niter + 1
    Bs.prevs = list(Bs.prevs[[2]], Bs)
    opt.cond = if (use.strict) {
      (err < threshold && Lerr < threshold)
    } else {
      (err < threshold || Lerr < threshold)
    }
    if (opt.cond || niter > maxit || is.nan(err)) {
      mu = lapply(1:K, function(k) {
        (diag(Ng) - Bs[[k]]) %*% meanYs[[k]] - sapply(1:Ng, function(i) {
          meanXs[[i]] %*% fs[[k]][[i]]
        })
      })
      sigma2 = sigma2_sem(Xs, Ys, Bs, fs, Ng, Ns, K)
      if (sparse) {
        Bs = lapply(Bs, Matrix, sparse = T)
      }
      break
    }
  } # while(niter <= maxit)
  list(
    B = Bs,
    f = fs,
    mu = mu,
    sigma2 = sigma2,
    niter = niter,
    err = err,
    detIB = detIBs
  )
}


## cross-validation and EBIC for hyper-parameter tuning
## lambda max can be estimated
#' @description get_lambda.max
#' @example
#' lamax = get_lambda.max(params.init$B, Ys, Xs, Ng)
get_lambda.max = function(Bs, Ys, Xs, Ng, weighted = TRUE) {
  std  = centralize_mult(Xs, Ys)
  Xs   = std$X               ## N x sk
  Ys   = std$Y               ## N x p
  K    = length(Ys)
  Ns   = sapply(Ys, nrow)
  R    = vector("list", K)   ## Ng
  w    = if(weighted) {
    inverse(Bs)
  } else {
    invone(Bs)
  }
  for (k in 1:K) {
    R[[k]] = matrix(0, nrow = Ng, ncol = Ns[k])   ## Ng x N
  }
  for (i in 1:Ng) {
    Xi = Xs[[i]]
    Pi = solve(crossprod(Xi)) %*% t(Xi)
    for (k in 1:K) {
      yi = Ys[[k]][, i, drop = F]   # n x 1
      fi = solve(crossprod(Xi)) %*% t(Xi) %*% yi
      Xf = Xi %*% fi               # n x 1
      R[[k]][i,] = yi - Xf
    }
  }
  err = 0
  for (k in 1:K) {
    err = err + norm(R[[k]], type = "f") ^ 2
  }
  sigma2 = err / (Ng * sum(Ns))
  Ry  = vector("list", K)   ## Ng
  for (k in 1:K) {
    Ry[[k]] = R[[k]] %*% Ys[[k]]
    Ry[[k]] = abs(Ry[[k]] / sigma2) / w[[k]]
  }
  max(sapply(Ry, max))
}


## cross-validation and EBIC for hyper-parameter tuning
## rho max can be estimated, rho is the fused lasso
## regularized hyper parameter
#' @description get_rho.max
#' @example
#' rhomax = get_rho.max(params.init$B, params.init$F, Ys, Xs, params.init$sigma2[1], data$var$Ng)
get_rho.max = function(Bs, fs, Ys, Xs, sigma2, Ng, weighted = TRUE) {
  if(weighted) {
    wBs = inverse(Bs)
    rB  = flinv(Bs)
  } else {
    wBs = invone(Bs)
    rB  = flone(Bs)
  }
  params.rho = multiSML_iPALM(
    Bs,
    fs,
    Ys,
    Xs,
    sigma2,
    Ng,
    lambda = 0,
    rho = Inf,
    wBs = wBs,
    rB  = rB,
    maxit = 2000,
    threshold = 1e-4,
    use.strict = F,
    sparse = T,
    verbose = 1
  )
  weight.rho = if(weighted) {
    flinv(Bs)
  } else {
    flone(Bs)
  }
  Bs = params.rho$B[[1]]
  fs = params.rho$f
  sigma2 = params.rho$sigma2
  std = centralize_mult(Xs, Ys)
  Xs = std$X           ## Ng (n x sk)
  Ys = std$Y           ## n x ng
  ## multiple
  K   = length(Ys)
  Ns  = sapply(Ys, nrow)
  Bc  = vector("list", K)
  YY  = vector("list", K)
  FX  = vector("list", K)
  Dx  = vector("list", K)
  for (k in 1:K) {
    Bc[[k]] = -Ns[k] * t(solve(diag(Ng) - Bs))
    YY[[k]] = crossprod(Ys[[k]])   ## Y %*% t(Y)
    FX[[k]] = matrix(0, nrow = Ng, ncol = Ns[k])  ## F %*% X (p x k x k x n = p x n)
    for (i in 1:Ng) {
      FX[[k]][i, ] = as.numeric(Xs[[i]] %*% fs[[k]][[i]])
    }
    Dx[[k]] = abs(Bc[[k]] + ((diag(Ng) - Bs) %*% YY[[k]] - FX[[k]] %*% Ys[[k]]) / sigma2) / weight.rho
    diag(Dx[[k]]) = -Inf
  }
  ## Dxy = abs((diag(Ng) - Bs) %*% (YY[[2]] - YY[[1]]) - (FX[[2]] %*% Ys[[2]] - FX[[1]] %*% Ys[[1]])) / sigma2 / 2 / weight.rho
  ## diag(Dxy) = -Inf
  max(c(max(Dx[[1]]), max(Dx[[2]])))
  ## max(Dxy)
}


## cross validation for hyper-parameter tuning
#' @description 5-fold cross-validation
#' @param dyn dynamic updated rho.max by given lambda
#' @example
#' cv.params = cv_multiSML(Bs = params.init$B, fs = params.init$F, Ys = Ys, Xs = Xs, sigma2 = params.init$sigma2[1], Ng = data$var$Ng, nlambda = 20, nrho = 20)
cv_multiSML = function(Bs,
                       fs,
                       Ys,
                       Xs,
                       sigma2,
                       Ng,
                       nlambda = 20,
                       nrho = 20,
                       threshold = 1e-4,
                       verbose = 1) {
  lambda.max = get_lambda.max(Bs, Ys, Xs, Ng)
  lambda.factors = 10 ^ seq(0, -4, length.out = nlambda) * lambda.max
  wBs = inverse(Bs)
  rB  = flinv(Bs)
  rho.max = get_rho.max(Bs, fs, Ys, Xs, sigma2, Ng)
  rho.factors = 10 ^ seq(0, -4, length.out = nrho) * rho.max
  ncv = 5
  Ns  = sapply(Ys, ncol)
  K   = length(Ys)
  Ytrain  = vector("list", ncv)
  Xtrain  = vector("list", ncv)
  Ytest   = vector("list", ncv)
  Xtest   = vector("list", ncv)
  cv.fold = sample(seq(1, ncv), size = Ns[1], replace = T)
  for (i in 1:ncv) {
    Ytrain[[i]] = lapply(Ys, function(y) {
      y[, cv.fold != i, drop = F]
    })
    Xtrain[[i]] = lapply(Xs, function(x) {
      x[, cv.fold != i, drop = F]
    })
    Ytest[[i]]  = lapply(Ys, function(y) {
      y[, cv.fold == i, drop = F]
    })
    Xtest[[i]] = lapply(Xs, function(x) {
      x[, cv.fold == i, drop = F]
    })
  }
  cverrs = vector("list", nrho * nlambda)
  cvlls  = vector("list", nrho * nlambda)
  hyper.params = list()
  for (cv in 1:ncv) {
    params.opt = list()
    Nc = sapply(Ytest[[cv]], ncol)
    ix = 1
    for (rho in rho.factors) {
      for (lambda in lambda.factors) {
        cat(sprintf("lambda = %4f, rho = %4f, foldid = %d\n", lambda, rho, cv))
        if (ix %% nlambda == 1) {
          params.opt[[ix]] = multiSML_iPALM(
            Bs,
            fs,
            Ytrain[[cv]],
            Xtrain[[cv]],
            sigma2,
            Ng,
            lambda,
            rho,
            wBs,
            rB,
            maxit = 1000,
            threshold = threshold,
            acc = TRUE,
            sparse = FALSE,
            use.strict = FALSE,
            verbose = verbose
          )
        } else {
          params.opt[[ix]] = multiSML_iPALM(
            params.opt[[ix - 1]]$B,
            params.opt[[ix - 1]]$f,
            Ytrain[[cv]],
            Xtrain[[cv]],
            params.opt[[ix - 1]]$sigma2,
            Ng,
            lambda,
            rho,
            wBs,
            rB,
            maxit = 1000,
            threshold = threshold,
            acc = TRUE,
            sparse = FALSE,
            use.strict = FALSE,
            verbose = verbose
          )
        }
        loglik = SML_logLik(
            Xtest[[cv]],
            Ytest[[cv]],
            params.opt[[ix]]$B,
            params.opt[[ix]]$f,
            Ng,
            Nc,
            K,
            params.opt[[ix]]$detIB,
            params.opt[[ix]]$sigma2
          )[1]
        err = SML_error(Xtest[[cv]],
                    Ytest[[cv]],
                    params.opt[[ix]]$B,
                    params.opt[[ix]]$f,
                    Ng,
                    Nc,
                    K)[1]
        cverrs[[ix]] = c(cverrs[[ix]], err)
        cvlls[[ix]]  = c(cvlls[[ix]], loglik)
        if (cv == 1) {
          hyper.params[[ix]] = c(lambda, rho)
        }
        ix = ix + 1
      }
    }
  }
  list(opt.hyperparams = hyper.params,
       cverrs = cverrs,
       loglik = cvlls)
}



## ultility funciton
cvsurface = function(cvparams, type = c("err", "loglik")) {
  cvm = if(type == "err") {
    cvparams$cverrs
  } else {
    cvparams$loglik
  }
  cvfuns = data.frame(
    lambda = sapply(cvparams$opt.hyperparams, `[`, 1),
    rho    = sapply(cvparams$opt.hyperparams, `[`, 2),
    cvmean = sapply(cvm, mean),
    cvsd   = sapply(cvm, sd)
  )
  lambda = sort(unique(cvfuns$lambda))
  rho    = sort(unique(cvfuns$rho))
  cvmean = matrix(nrow = length(lambda), ncol = length(rho))
  cvfuns$lambda = sapply(cvfuns$lambda, function(l) {
    which(lambda == l)
  })
  cvfuns$rho    = sapply(cvfuns$rho, function(r) {
    which(rho == r)
  })
  apply(cvfuns, 1, function(x) {
    cvmean[x[1], x[2]] <<- x[3]
  })
  require(plotly)
  if(type == "err") {
    surface = plot_ly(x = log10(lambda),
                    y = log10(rho),
                    z = log10(cvmean)) %>% add_surface()
  } else {
    surface = plot_ly(x = log10(lambda),
                    y = log10(rho),
                    z = cvmean) %>% add_surface()
  }

  list(
    lambda = lambda,
    rho = rho,
    cvm = cvmean,
    surf = surface
  )
}

## ultility functions
## pick lambda
optimLambda_cv = function(cvparams, type = c("err", "loglik"), se = TRUE, fused.sparse = TRUE) {
  cvm = if(type == "err") {
    cvparams$cverrs
  } else {
    cvparams$loglik
  }
  cvfuns = data.frame(
    lambda = sapply(cvparams$opt.hyperparams, `[`, 1),
    rho    = sapply(cvparams$opt.hyperparams, `[`, 2),
    cvmean = sapply(cvm, mean),
    cvsd   = sapply(cvm, sd)
  )
  cv.min     = which.min(cvfuns$cvmean)
  cv.1se     = cvfuns[cv.min, 3] + cvfuns[cv.min, 4]
  cvfun.1se  = cvfuns[cvfuns$cvmean <= cv.1se, c(1, 2, 3)]
  lambda = cvfun.1se[, 1]
  rho    = cvfun.1se[, 2]
  if(fused.sparse) {
    rho.1se = max(rho)
    lambda.1se = lambda[which.min(cvfun.1se[cvfun.1se$rho == rho.1se, 3])]
  } else {
    lambda.1se = max(lambda)
    # rho.1se    = min(rho[lambda == lambda.1se])
    rho.1se = rho[which.min(cvfun.1se[cvfun.1se$lambda == lambda.1se, 3])]
  }

  if (se) {
    list(lambda = lambda.1se, rho = rho.1se)
  } else {
    list(lambda = cvfuns[cv.min, 1], rho = cvfuns[cv.min, 2])
  }
}

## optimLambda_eBIC
optimLambda_eBIC = function(eBICparams) {
  eBICfuns = data.frame(
    lambda = sapply(eBICparams$opt.hyperparams, `[`, 1),
    rho    = sapply(eBICparams$opt.hyperparams, `[`, 2),
    eBIC = eBICparams$eBIC
  )
  eBIC.min     = which.min(eBICfuns$eBIC)
  lambda.min   = eBICfuns[eBIC.min, 1]
  rho.min      = eBICfuns[eBIC.min, 2]
  list(lambda = lambda.min, rho = rho.min)
}

## optimGamma_eBIC
optimGamma_eBIC = function(BICparams) {
  optim = vector("list", 11)
  i = 1
  eBIC  = vector("numeric", 11)
  for (gamma in seq(0, 1, by = 0.1)) {
    eBICs    = BICparams$BIC + gamma * BICparams$extend
    eBICparams = list(opt.hyperparams = BICparams$opt.hyperparams,
                      eBIC = eBICs)
    optim[[i]] = optimLambda_eBIC(eBICparams)
    eBIC[i]  = min(eBICs)
    i = i + 1
  }
  list(
    gamma = seq(0, 1, by = 0.1),
    optimLambda = optim,
    eBIC = eBIC
  )
}


#' extended BIC
#' @param gamma [0,1]
eBIC_multSML = function(Bs,
                        fs,
                        Ys,
                        Xs,
                        sigma2,
                        Ng,
                        Nk,
                        nlambda = 20,
                        nrho = 20,
                        verbose = 1) {
  lambda.max = get_lambda.max(Bs, Ys, Xs, Ng)
  lambda.factors = 10 ^ seq(0, -4, length.out = nlambda) * lambda.max
  wBs = inverse(Bs)
  rB  = flinv(Bs)
  rho.max = get_rho.max(Bs, fs, Ys, Xs, sigma2, Ng)
  rho.factors = 10 ^ seq(0, -4, length.out = nrho) * rho.max
  ncv = 5
  Ns  = sapply(Ys, ncol)
  K   = length(Ys)
  hyper.params = vector("list", nrho * nlambda)
  params.prev  = NULL
  BIC      = vector("numeric", nrho * nlambda)
  extend   = vector("numeric", nrho * nlambda)
  ix  = 1
  for (rho in rho.factors) {
    for (lambda in lambda.factors) {
      cat(sprintf("lambda = %f, rho = %f\n", lambda, rho))
      if (ix %% nlambda == 1) {
        params.opt = multiSML_iPALM(
          Bs,
          fs,
          Ys,
          Xs,
          sigma2,
          Ng,
          lambda,
          rho,
          wBs,
          rB,
          maxit = 1000,
          threshold = 1e-4,
          acc = TRUE,
          sparse = FALSE,
          verbose = verbose
        )
      } else {
        params.opt = multiSML_iPALM(
          params.prev$B,
          params.prev$f,
          Ys,
          Xs,
          params.prev$sigma2,
          Ng,
          lambda,
          rho,
          wBs,
          rB,
          maxit = 1000,
          threshold = 1e-4,
          acc = TRUE,
          sparse = FALSE,
          verbose = verbose
        )
      }
      logLik = SML_logLik(
        Xs,
        Ys,
        params.opt$B,
        params.opt$f,
        Ng,
        Ns,
        K,
        params.opt$detIB,
        params.opt$sigma2
      )[1]
      df = sum(params.opt$B[[1]] != 0) +
        sum(params.opt$B[[2]] != 0) -
        sum(params.opt$B[[2]] == params.opt$B[[1]] &
        params.opt$B[[1]] != 0) + 2 * Nk
      BIC[ix] = 2 * logLik + df * log(sum(Ns))
      extend[ix] = 2 * log(choose(2 * (Ng ^ 2 - Ng + Nk), df))
      hyper.params[[ix]] = c(lambda, rho)
      params.prev = params.opt
      ix = ix + 1
    }
  }
  list(opt.hyperparams = hyper.params,
       BIC = BIC,
       extend = extend)
}

## ultility
## ultility funciton
eBICsurface = function(eBICparams, gamma = 0) {
  ebicfuns = data.frame(
    lambda = sapply(eBICparams$opt.hyperparams, `[`, 1),
    rho    = sapply(eBICparams$opt.hyperparams, `[`, 2),
    eBIC   = eBICparams$BIC + gamma * eBICparams$extend
  )
  lambda = sort(unique(ebicfuns$lambda))
  rho    = sort(unique(ebicfuns$rho))
  eBIC   = matrix(nrow = length(lambda), ncol = length(rho))
  ebicfuns$lambda = sapply(ebicfuns$lambda, function(l) {
    which(lambda == l)
  })
  ebicfuns$rho    = sapply(ebicfuns$rho, function(r) {
    which(rho == r)
  })
  apply(ebicfuns, 1, function(x) {
    eBIC[x[1], x[2]] <<- x[3]
  })
  require(plotly)
  surface = plot_ly(x = log(lambda),
                    y = log(rho),
                    z = log(eBIC)) %>% add_surface()
  list(
    lambda = lambda,
    rho = rho,
    ebics = eBIC,
    surf = surface
  )
}

## eQTL and gene regulatory network are all different with each other under
## different conditions

#' @description getdiffeQTL
#' @param N   number of sample
#' @param Ng  number of gene
#' @param k   number of eQTL
#' @param d   differential ratio = 0.1
getdiffeQTL = function(N, Ng, Nk, d = 0.1) {
  step = Nk / Ng
  X = vector("list", 2)
  X[[1]] = round(2 * matrix(runif(Nk * N), nrow = Nk)) + 1
  X[[2]] = apply(X[[1]], 2, function(x) {
    Nd = Nk * d
    dx = sample(1:Nk, Nd, replace = F)
    x[dx] = round(2 * runif(Nd)) + 1
    x
  })
  G = matrix(0,
             nrow = Ng,
             ncol = Nk)
  ix = lapply(1:Ng,
              function(i) {
                s = seq(0, step - 1) * Ng + i
                G[i, s] <<- 1
                s
              })
  list(G = G, X = X, sk = ix)
}


#' @description getdiffsem
#' @details eQTL measurement for different condition are generated with proportional difference.
#' @param f difference proportion of each gene's eQTL measurement, such as SNP
getdiffsem = function(N = 200,
                      Ng = 10,
                      Nk = 10,
                      r = 0.3,
                      d = 0.1,
                      f = 0.1,
                      dag = TRUE,
                      sigma = 0.1) {
  B = getrandDAG(Ng, e = Ng * r, dag = dag, d = d)
  Q = getdiffeQTL(N, Ng, Nk, f)
  F = Q[[1]]
  X = Q[[2]]
  sk = Q[[3]]
  E1 = sigma * t(rmvnorm(N, mean = rep(0, Ng), sigma = diag(Ng)))
  E2 = sigma * t(rmvnorm(N, mean = rep(0, Ng), sigma = diag(Ng)))
  Y1 = solve(diag(Ng) - B[[1]]) %*% (F %*% X[[1]] + E1)
  Y2 = solve(diag(Ng) - B[[2]]) %*% (F %*% X[[2]] + E2)

  list(
    obs = list(
      Y1 = Y1,
      Y2 = Y2,
      X1 = X[[1]],
      X2 = X[[2]],
      sk = sk
    ),
    var = list(
      B1 = Matrix(B[[1]], sparse = T),
      B2 = Matrix(B[[2]], sparse = T),
      F = Matrix(F, sparse = T),
      N = N,
      Ng = Ng,
      Nk = Nk
    )
  )
}


## data = getdiffsem(N = 200, Ng = 30, Nk = 90, r = 0.1, d = 0.1, f = 0.1, sigma = 1, dag = TRUE)
#' @description build sumbatrix for eQTLs observation for subset of corresponding genes
submatXs = function(data) {
  sk = data$obs$sk
  Xs  = list(X1 = data$obs$X1, X2 = data$obs$X2)
  lapply(Xs, function(X) {
    lapply(sk, function(ix) {
      X[ix, , drop = F]
    })
  })
}

## centralized for multiple Ys and Xs
#' @description generalized centralization of Ys and Xs
#'              Xs -> n x sk
#'              Ys -> n x ng
#' @example
#' Ys = data$obs[c("Y1", "Y2")]
#' Xs = submatXs(data)
centralize_gen = function(Xs, Ys) {
  meanX = lapply(Xs, function(X) {
    lapply(X, rowMeans)
  })
  meanY = lapply(Ys, rowMeans)
  Xs = lapply(Xs, function(X) {
    lapply(X, center)
  })
  Ys = lapply(Ys, center)
  list(X = Xs,
       Y = Ys,
       muX = meanX,
       muY = meanY)
}


## ridge regression for estimate sigma2 initialization
## on different gene expressionand different eQTLs
#' @param M number of gene
#' @param N number of sample
#' @example
#' M = data$var$Ng
#' N = data$var$N
#' params.init = constrained_L2reg_gen(Xs, Ys, sigma2$rho.opt, M, N)
constrained_L2reg_gen = function(Xs, Ys, rho, M, N) {
  K = length(Ys)
  B = list()
  F = list()
  mu = list()
  err = 0
  df  = 0
  for (i in 1:K) {
    fit = constrained_L2reg(Xs[[i]], Ys[[i]], rho)
    B[[i]]  = as.matrix(fit$B)
    F[[i]]  = fit$F
    mu[[i]] = fit$mu
    err = err + fit$sigma2 * (N[i] * M - 1)
    df  = df + (N[i] * M - 1)
  }
  sigma2 = err / df
  list(
    B = B,
    F = F,
    sigma2 = sigma2,
    mu = mu
  )
}


## generalized cross-validation on ridge regression to estimate sigma2
## on different gene expressionand different eQTLs
#' @param nrho number of L2 penalty's coefficient
#' @param ncv  number of cross-validation
#' @example
#' Ys = data$obs[c("Y1", "Y2")]
#' Xs = submatXs(data)
#' M  = data$var$Ng
#' N  = data$var$N
#' sigma2 = getsigma2_L2reg_gen(Xs, Ys, nrho = 20, M = M, N = N)
getsigma2_L2reg_gen = function(Xs,
                               Ys,
                               nrho = 10,
                               ncv = 5,
                               M,
                               N) {
  rho_factors = 10 ** (seq(-6, 2, length.out = nrho))
  cv.err  = matrix(0, nrow = nrho, ncol = ncv)
  cv.fold = list()
  cv.fold[[1]] = sample(seq(1, ncv), size = N[1], replace = T)
  cv.fold[[2]] = sample(seq(1, ncv), size = N[2], replace = T)
  irho    = 1
  for (rho in rho_factors) {
    for (cv in 1:ncv) {
      ytrain = lapply(1:2, function(ix) {
        Ys[[ix]][, cv.fold[[ix]] != cv, drop = F]
      })
      xtrain = lapply(1:2, function(ix) {
        lapply(Xs[[ix]], function(x) {
          x[, cv.fold[[ix]] != cv, drop = F]
        })
      })
      ytest  = lapply(1:2, function(ix) {
        Ys[[ix]][, cv.fold[[ix]] == cv, drop = F]
      })
      xtest  = lapply(1:2, function(ix) {
        lapply(Xs[[ix]], function(x) {
          x[, cv.fold[[ix]] == cv, drop = F]
        })
      })
      Ntrain = sapply(cv.fold, function(f){ sum(f != cv) })
      fit    = constrained_L2reg_gen(xtrain, ytrain, rho, M, Ntrain)
      for (k in 1:length(Ys)) {
        ftest  = lapply(1:M, function(i) {
          crossprod(fit$F[[k]][[i]], xtest[[k]][[i]])
        })
        ftest  = do.call(rbind, ftest)
        cv.err[irho, cv] = cv.err[irho, cv] + norm((diag(M) - fit$B[[k]]) %*% ytest[[k]] - ftest - fit$mu[[k]], type = "f")
      }
    }
    irho = irho + 1
  }
  cv.mean = rowMeans(cv.err)
  rho.min = rho_factors[which.min(cv.mean)]
  fit = constrained_L2reg_gen(Xs, Ys, rho.min, M, N)
  list(
    rho.opt = rho.min,
    sigma2.opt = fit$sigma2[1],
    cv.ram = list(rho = rho_factors, cvm = cv.mean)
  )
}


## solve SML problem by component-wise update with generalized configuration
## different gene expression and different eQTLs
#' @param lambda lasso penalty
#' @param rho    fused lasso penalty
#' @example
#' M = data$var$Ng
#' N = data$var$N
#' Ys = data$obs[c("Y1", "Y2")]
#' Xs = submatX(data)
#' sigma2 = getsigma2_L2reg_multi(Xs, Ys, nrho = 20, M = M, N = N)
#' params.init = constrained_L2reg_multi(Xs, Ys, sigma2$rho.opt, M, N)
#' params.opt = genSML_cwise(Bs = params.opt4$B, fs = params.opt4$f, Ys = Ys, Xs = Xs,
#'                           sigma2 = params.init$sigma2[1], Ng = data$var$Ng,
#'                           wBs = inverse(params.init$B), rB = flinv(params.init$B),
#'                           lambda = 10, rho = 40, maxit = 1000)
genSML_cwise = function(Bs,
                        fs,
                        Ys,
                        Xs,
                        sigma2,
                        Ng,
                        lambda,
                        rho,
                        weighted = TRUE,
                        wBs = inverse(Bs),
                        rB  = flinv(Bs),
                        maxit = 100,
                        threshold = 1e-4,
                        verbose = 2) {
  std = centralize_gen(Xs, Ys)
  Xs  = std$X
  Ys  = std$Y
  meanXs = std$muX
  meanYs = std$muY
  K   = length(Ys)  # number of conditions = 2
  if (verbose == 2) {
    cat(sprintf("conditions must be restricted to 2.. K = %d\n", K))
  }
  ## update for eQTL coeffs
  f0 = vector("list", K)
  f1 = vector("list", K)
  for (i in 1:Ng) {
    for (k in 1:K) {
      Xi = Xs[[k]][[i]]
      Pi = solve(crossprod(Xi)) %*% t(Xi)
      yi_k = Ys[[k]][, i, drop = F]     # n x 1 for gene i
      f0[[k]][[i]] = Pi %*% yi_k
      f1[[k]][[i]] = Pi %*% Ys[[k]]     # f = f0 - f1 %*% B[i,]
    }
  }
  ## update for gnet coeffs
  niter  = 1
  Ns     = sapply(Ys, nrow)
  ImBs   = lapply(Bs, function(B) {
    diag(Ng) - B
  })
  IBsinv = lapply(ImBs, solve)
  while (niter <= maxit) {
    Bs.prev = Bs
    fs.prev = fs
    for (i in 1:Ng) {
      ci  = lapply(IBsinv, function(IBi) {
        IBi[, i]
      })
      dbi = lapply(1:K, function(k) {
        vector("numeric", Ng)
      })
      for (j in 1:Ng) {
        ## update B[[k]][i,j] for i != j
        if (i != j) {
          bij.prev = bij = sapply(Bs, function(B)
            (B[i, j]))
          wij      = sapply(wBs, function(w) {
            if (weighted) {
              w[i, j]
            } else {
              1
            }
          })
          rij     = if (weighted) {
            rB[i, j]
          } else {
            1
          }
          mij     = sapply(ci, function(c) {
            c[j]
          })
          bi      = lapply(ImBs, function(ImB) {
            bi_k = ImB[i,]
            bi_k[j] = 0
            bi_k
          })
          ## j-th column of Ys
          Yej    = lapply(Ys, function(Y) {
            Y[, j, drop = F]
          })
          a1     = sapply(1:K, function(k) {
            crossprod(Ys[[k]] %*% bi[[k]] - Xs[[k]][[i]] %*% fs[[k]][[i]], Yej[[k]])
          })
          a2     = sapply(1:K, function(k) {
            crossprod(Yej[[k]])
          })
          ## a0 = 1/mij + bij.prev
          cond   = list(
            c(1, 1, 1),
            c(1, -1, 1),
            c(1, 0, 1),
            c(-1, -1, 1),
            c(0, -1, 1),
            c(1, 1, -1),
            c(0, 1, -1),
            c(-1, -1, -1),
            c(-1, 0, -1),
            c(-1, 1, -1),
            c(1, 1, 0),
            c(-1, -1, 0),
            c(0, 0, 0)
          )
          obj    = obj_multiSML(Ns, mij, bij.prev, a1, a2, lambda, rho, wij, rij, sigma2)
          grad   = grad_multiSML(Ns, mij, bij.prev, a1, a2, lambda, rho, wij, rij, sigma2)
          params = list()
          for (t in cond) {
            cand.grad = grad(t)
            params    = c(params, grad_solver(cand.grad, t))
          }
          objval = sapply(params, function(args) {
            do.call(obj, args)
          })
          mix    = which.min(objval)
          bij    = unlist(params[[mix]])
          dbij = bij.prev - bij
          for (k in 1:K) {
            dbi[[k]][j]  = dbij[k]
            Bs[[k]][i, j] = bij[k]
            ci[[k]]      = ci[[k]] / (1 + dbij[k] * mij[k])
            ImBs[[k]]    = diag(Ng) - Bs[[k]]
          }
        }
      } ## for(j in 1:Ng)
      ## (ImB + ei^T %*% dbi)^{-1}
      for (k in 1:K) {
        ## IBsinv[[k]] = IBsinv[[k]] - IBsinv[[k]][, i, drop = F] %*% dbi[[k]] %*% IBsinv[[k]] / (1 + dbi[[k]] %*% IBsinv[[k]][, i, drop = F])[1]
        IBsinv[[k]]   = solve(ImBs[[k]])
        fs[[k]][[i]]  = f0[[k]][[i]] - f1[[k]][[i]] %*% Bs[[k]][i, ]
      }
    } ## for(i in 1:Ng)
    Berr = sum(sapply(1:K, function(k) {
      norm(Bs[[k]] - Bs.prev[[k]], type = "f") / norm(Bs.prev[[k]], type = "f")
    }))
    Ferr = sum(sapply(1:K, function(k) {
      sum(sapply(1:Ng, function(i) {
        norm(fs[[k]][[i]] - fs.prev[[k]][[i]], type = "f")
      })) / sum(sapply(1:Ng, function(i) {
        norm(fs.prev[[k]][[i]], type = "f")
      }))
    }))
    err = Berr + Ferr
    cat(sprintf("SML: iteration = %d, error = %f\n", niter, err))
    niter = niter + 1
    if (err < threshold || niter > maxit || is.nan(err)) {
      mu = lapply(1:K, function(k) {
        (diag(Ng) - Bs[[k]]) %*% meanYs[[k]] - sapply(1:Ng, function(i) {
          meanXs[[k]][[i]] %*% fs[[k]][[i]]
        })
      })
      Bs = lapply(Bs, Matrix, sparse = T)
      break
    }
  } ## while(niter <= maxit)
  list(
    B = Bs,
    f = fs,
    mu = mu,
    niter = niter,
    err = err
  )
}

sigma2_gen = function(Xs, Ys, B, f, Ng, Ns, K) {
  X = lapply(Xs, function(X) {
    lapply(X, t)
  })
  Y = lapply(Ys, t)
  err = 0
  for (k in 1:K) {
    for (i in 1:Ng) {
      Xi = X[[k]][[i]]                 # sk x N
      bi = B[[k]][i,-i, drop = F]     # (ng-1) x 1
      yi = Y[[k]][i, , drop = F]       # 1 x N
      Yi = Y[[k]][-i, , drop = F]      # (ng-1) x N
      fi = f[[k]][[i]]                 # sk x 1
      err = err + tcrossprod(yi - bi %*% Yi - crossprod(fi, Xi))
    }
  }
  sigma2 = err / (Ng * sum(Ns))
  sigma2[1]
}


## solve SML problem by block coordinate descent by backtracking inert-PALM on
## different gene expression and different eQTLs
#' @param lambda lambda hyper-parameter for lasso term
#' @param rho    fused lasso hyper-parameter for fused lasso term
#' @param gamma  invertible matrix stablize parameter gamma
#' params.init = constrained_L2reg_gen(Xs, Ys, sigma2$rho.opt, M, N)
#' params.opt = genSML_iPALM(Bs = params.init$B, fs = params.init$F, Ys = Ys, Xs = Xs,
#'                            sigma2 = params.init$sigma2[1], Ng = data$var$Ng,
#'                            lambda = 0.1, rho = 0.1, maxit = 500)
genSML_iPALM = function(Bs,
                        fs,
                        Ys,
                        Xs,
                        sigma2,
                        Ng,
                        lambda,
                        rho,
                        wBs = inverse(Bs),
                        rB  = flinv(Bs),
                        maxit = 100,
                        acc = TRUE,
                        inertial = inertial_pars("lin"),
                        threshold = 1e-3,
                        sparse = FALSE,
                        use.strict = TRUE,
                        verbose = 2) {
  std = centralize_gen(Xs, Ys)
  Xs  = std$X
  Ys  = std$Y
  meanXs = std$muX
  meanYs = std$muY
  K      = length(Ys)      # number of conditions = 2
  if (verbose == 2) {
    cat(sprintf("conditions must be restricted to 2.. K = %d\n", K))
  }
  ## update for eQTL row-wise (specific for genes)
  f0 = vector("list", K)
  f1 = vector("list", K)
  Yp = vector("list", K)
  Hy = vector("list", K)
  Yp.maxEigen = vector("list", K)
  Ns = sapply(Ys, nrow)
  for (i in 1:Ng) {
    for (k in 1:K) {
      Xi = Xs[[k]][[i]]
      Pi = solve(crossprod(Xi)) %*% t(Xi)
      # specific condition
      yi_k = Ys[[k]][, i, drop = F]              # n[k] x 1 for gene i
      Yi_k = Ys[[k]][,-i]                       # n[k] x (ng-1) (for specific gene i)
      f0[[k]][[i]] = Pi %*% yi_k
      f1[[k]][[i]] = Pi %*% Yi_k                 # f[[k]][[i]] = f0[[k]][[i]] - f1[[k]][[i]] %*% bi^k | bi^k = B[[k]][i,-i]
      Hi_k = diag(Ns[k]) - Xi %*% Pi             # n[k] x n[k] projection matrix
      Yp[[k]][[i]] = t(Yi_k) %*% Hi_k %*% Yi_k   # (ng - 1) x (ng - 1)
      Hy[[k]][[i]] = t(yi_k) %*% Hi_k %*% Yi_k   # 1 x (ng - 1)
      ## maximized eigen-value for Yp
      Yp.maxEigen[[k]][[i]] = eigen(Yp[[k]][[i]])$values[1]
    }
  }
  ## update for gnet row-wise
  niter  = 1
  ImBs   = lapply(Bs, function(B) {
    diag(Ng) - B
  })
  detIBs = sapply(ImBs, det)
  IBsinv = lapply(ImBs, solve)
  Ls     = logLik(detIBs, Bs, wBs, rB, lambda, rho, Ns, K) + N * sum(Ns) / 2 * log(sigma2)
  ## history
  Bs.prevs = list(Bs, Bs)
  inert    = acc
  while (niter <= maxit) {
    inert.pars = inertial(niter)
    Bs.inert = if (inert) {
      lapply(1:K, function(k) {
        Bs.prevs[[2]][[k]] + inert.pars * (Bs.prevs[[2]][[k]] - Bs.prevs[[1]][[k]])
      })
    } else {
      Bs
    }
    fs.prev = fs
    Ls.prev = Ls
    for (i in 1:Ng) {
      ## sum(-Ns[k] * sigma2 * log(det(I-B[[k]])^2)) + \sum_{i=1}^{Ng} bi^T %*% Yp %*% bi - Hy %*% bi)
      ci = lapply(IBsinv, function(IBi) {
        # ci / det(I - B); from (I - B)^{-1}
        IBi[-i, i, drop = F]
      })
      bi = lapply(Bs.inert, function(B.inert) {
        t(B.inert[i,-i, drop = F])
      })
      gi = lapply(1:K, function(k) {
        grad = grad_rwise_SML(Ns[k], ci[[k]], Yp[[k]][[i]], Hy[[k]][[i]], sigma2[1])
        grad(bi[[k]])
      })
      ## Lipschitz moduli for row-i
      Lis  = sapply(1:K, function(k) {
        gtg  = tcrossprod(ImBs[[k]][-i, ])
        oi   = chol2inv(gtg)
        deti = det(gtg)
        gii  = ImBs[[k]][-i,-i]
        si   = ImBs[[k]][-i, i, drop = F]
        c2i  = sum((ci[[k]] * detIBs[k]) ^ 2)
        lips_rwise_SML(Ns[k],
                       oi,
                       gii,
                       si,
                       c2i,
                       deti,
                       Yp.maxEigen[[k]][[i]],
                       sigma2,
                       Ng)[1]
      })
      Li   = max(Lis)
      Li   = (1 + 2 * inert.pars) * Li / (2 * (1 - inert.pars))
      ui   = lapply(1:K, function(k) {
        bi[[k]] - gi[[k]] / Li
      })
      wBi  = lapply(wBs, function(wB) {
        wB[i,-i]
      })
      rBi  = rB[i,-i]
      xi   = prox_flsa(lambda, rho, Li, ui, wBi, rBi)
      for (k in 1:K) {
        Bs[[k]][i,-i] = xi[[k]]
        ImBs[[k]]     = diag(Ng) - Bs[[k]]
        detIBs[k]     = (ImBs[[k]][i, ] %*% IBsinv[[k]][, i, drop = F])[1] * detIBs[k]
        dbi           = Bs.prevs[[2]][[k]][i, ] - Bs[[k]][i, ]
        IBsinv[[k]]   = IBsinv[[k]] - IBsinv[[k]][, i, drop = F] %*% dbi %*% IBsinv[[k]] / (1 + dbi %*% IBsinv[[k]][, i, drop = F])[1]
        fs[[k]][[i]]  = f0[[k]][[i]] - f1[[k]][[i]] %*% Bs[[k]][i,-i]
      }
    } # row-wise update
    Berr = sum(sapply(1:K, function(k) {
      norm(Bs[[k]] - Bs.prevs[[2]][[k]], type = "f") / (1 + norm(Bs.prevs[[2]][[k]], type = "f"))
    }))
    Ferr = sum(sapply(1:K, function(k) {
      sum(sapply(1:Ng, function(i) {
        norm(fs[[k]][[i]] - fs.prev[[k]][[i]], type = "f")
      })) / (1 + sum(sapply(1:Ng, function(i) {
        norm(fs.prev[[k]][[i]], type = "f")
      })))
    }))
    err = Berr + Ferr
    sigma2 = sigma2_gen(Xs, Ys, Bs, fs, Ng, Ns, K)
    Ls     = logLik(detIBs, Bs, wBs, rB, lambda, rho, Ns, K) + Ng * sum(Ns) / 2 * log(sigma2)
    Lerr   = abs(Ls.prev - Ls) / (1 + abs(Ls.prev))
    # inert  = ifelse(Lerr < 1e-8, FALSE, acc)
    if (verbose >= 2) {
      cat(
        sprintf(
          "SML: iteration = %d, error = %f, logLik = %f, sigma2 = %f, inert = %s\n",
          niter,
          err,
          Ls,
          sigma2,
          inert
        )
      )
    }
    niter = niter + 1
    Bs.prevs = list(Bs.prevs[[2]], Bs)
    opt.cond = if (use.strict) {
      (err < threshold && Lerr < threshold)
    } else {
      (err < threshold || Lerr < threshold)
    }
    if (opt.cond || niter > maxit || is.nan(err)) {
      mu = lapply(1:K, function(k) {
        (diag(Ng) - Bs[[k]]) %*% meanYs[[k]] - sapply(1:Ng, function(i) {
          meanXs[[k]][[i]] %*% fs[[k]][[i]]
        })
      })
      sigma2 = sigma2_gen(Xs, Ys, Bs, fs, Ng, Ns, K)
      if (sparse) {
        Bs = lapply(Bs, Matrix, sparse = T)
      }
      break
    }
  } # while(niter <= maxit)
  list(
    B = Bs,
    f = fs,
    mu = mu,
    sigma2 = sigma2,
    niter = niter,
    err = err,
    detIB = detIBs
  )
}

######################
# comparable method
######################
## cross validation for hyper-parameter tuning
#' @description 5-fold cross-validation
#' @param dyn dynamic updated rho.max by given lambda
#' @example
#' cv.params = cv_multiSML(Bs = params.init$B, fs = params.init$F, Ys = Ys, Xs = Xs, sigma2 = params.init$sigma2[1], Ng = data$var$Ng, nlambda = 20, nrho = 20)
cv_genSML = function(Bs,
                     fs,
                     Ys,
                     Xs,
                     sigma2,
                     Ng,
                     nlambda = 20,
                     nrho = 20,
                     weighted = TRUE,
                     threshold = 1e-4,
                     use.strict = FALSE,
                     verbose = 1) {
  lambda.max = gen_lambda.max(Bs, Ys, Xs, Ng, weighted)
  lambda.factors = 10 ^ seq(0,-5, length.out = nlambda) * lambda.max
  if(weighted) {
    wBs = inverse(Bs)
    rB  = flinv(Bs)
  } else{
    wBs = invone(Bs)
    rB  = flone(Bs)
  }

  rho.max = gen_rho.max(Bs, fs, Ys, Xs, sigma2, Ng, weighted)
  rho.factors = 10 ^ seq(0,-5, length.out = nrho) * rho.max
  ncv = 5
  Ns  = sapply(Ys, ncol)
  K   = length(Ys)
  Ytrain  = vector("list", ncv)
  Xtrain  = vector("list", ncv)
  Ytest   = vector("list", ncv)
  Xtest   = vector("list", ncv)
  cv.fold = list()
  cv.fold[[1]] = sample(seq(1, ncv), size = Ns[1], replace = T)
  cv.fold[[2]] = sample(seq(1, ncv), size = Ns[2], replace = T)
  for (i in 1:ncv) {
    Ytrain[[i]] = lapply(1:2, function(ix) {
      Ys[[ix]][, cv.fold[[ix]] != i, drop = F]
    })
    Xtrain[[i]] = lapply(1:2, function(ix) {
      lapply(Xs[[ix]], function(x) {
          x[, cv.fold[[ix]] != i, drop = F]
      })
    })
    Ytest[[i]]  = lapply(1:2, function(ix) {
      Ys[[ix]][, cv.fold[[ix]] == i, drop = F]
    })
    Xtest[[i]] = lapply(1:2, function(ix) {
      lapply(Xs[[ix]], function(x) {
          x[, cv.fold[[ix]] == i, drop = F]
      })
    })
  }
  cverrs = vector("list", nrho * nlambda)
  cvlls  = vector("list", nrho * nlambda)
  hyper.params = list()
  for (cv in 1:ncv) {
    params.opt = list()
    Nc = sapply(Ytest[[cv]], ncol)
    ix = 1
    for (rho in rho.factors) {
      for (lambda in lambda.factors) {
        cat(sprintf("lambda = %4f, rho = %4f, foldid = %d\n", lambda, rho, cv))
        if (ix %% nlambda == 1) {
          params.opt[[ix]] = genSML_iPALM(
            Bs,
            fs,
            Ytrain[[cv]],
            Xtrain[[cv]],
            sigma2,
            Ng,
            lambda,
            rho,
            wBs,
            rB,
            maxit = 1000,
            threshold = threshold,
            acc = TRUE,
            sparse = FALSE,
            use.strict = use.strict,
            verbose = verbose
          )
        } else {
          params.opt[[ix]] = genSML_iPALM(
            params.opt[[ix - 1]]$B,
            params.opt[[ix - 1]]$f,
            Ytrain[[cv]],
            Xtrain[[cv]],
            params.opt[[ix - 1]]$sigma2,
            Ng,
            lambda,
            rho,
            wBs,
            rB,
            maxit = 1000,
            threshold = threshold,
            acc = TRUE,
            sparse = FALSE,
            use.strict = use.strict,
            verbose = verbose
          )
        }
        loglik = gen_logLik(
          Xtest[[cv]],
          Ytest[[cv]],
          params.opt[[ix]]$B,
          params.opt[[ix]]$f,
          Ng,
          Nc,
          K,
          params.opt[[ix]]$detIB,
          params.opt[[ix]]$sigma2
        )[1]
        err = gen_error(Xtest[[cv]],
                        Ytest[[cv]],
                        params.opt[[ix]]$B,
                        params.opt[[ix]]$f,
                        Ng,
                        Nc,
                        K)[1]
        cverrs[[ix]] = c(cverrs[[ix]], err)
        cvlls[[ix]]  = c(cvlls[[ix]], loglik)
        if (cv == 1) {
          hyper.params[[ix]] = c(lambda, rho)
        }
        ix = ix + 1
      }
    }
  }
  list(opt.hyperparams = hyper.params,
       cverrs = cverrs,
       loglik = rbind, cvlls)
}

gen_lambda.max = function(Bs, Ys, Xs, Ng, weighted = TRUE) {
  std  = centralize_gen(Xs, Ys)
  Xs   = std$X               ## N x sk
  Ys   = std$Y               ## N x p
  K    = length(Ys)
  Ns   = sapply(Ys, nrow)
  R    = vector("list", K)   ## Ng
  w    = if(weighted) {
    inverse(Bs)
  } else {
    invone(Bs)
  }
  for (k in 1:K) {
    R[[k]] = matrix(0, nrow = Ng, ncol = Ns[k])   ## Ng x N
  }
  for (i in 1:Ng) {
    for (k in 1:K) {
      Xi = Xs[[k]][[i]]
      Pi = solve(crossprod(Xi)) %*% t(Xi)
      yi = Ys[[k]][, i, drop = F]   # n x 1
      fi = solve(crossprod(Xi)) %*% t(Xi) %*% yi
      Xf = Xi %*% fi               # n x 1
      R[[k]][i,] = yi - Xf
    }
  }
  err = 0
  for (k in 1:K) {
    err = err + norm(R[[k]], type = "f") ^ 2
  }
  sigma2 = err / (Ng * sum(Ns))
  Ry  = vector("list", K)   ## Ng
  for (k in 1:K) {
    Ry[[k]] = R[[k]] %*% Ys[[k]]
    Ry[[k]] = abs(Ry[[k]] / sigma2 - Ns[[k]]) / w[[k]]
  }
  max(sapply(Ry, max))
}

## cross-validation and EBIC for hyper-parameter tuning
## rho max can be estimated, rho is the fused lasso
## regularized hyper parameter
#' @description get_rho.max
#' @example
#' rhomax = get_rho.max(params.init$B, params.init$F, Ys, Xs, params.init$sigma2[1], data$var$Ng)
gen_rho.max = function(Bs, fs, Ys, Xs, sigma2, Ng, weighted = TRUE) {
  if(weighted) {
    wBs = inverse(Bs)
    rB  = flinv(Bs)
  } else {
    wBs = invone(Bs)
    rB  = flone(Bs)
  }
  params.rho = genSML_iPALM(
    Bs,
    fs,
    Ys,
    Xs,
    sigma2,
    Ng,
    lambda = 0,
    rho = Inf,
    wBs = wBs,
    rB = rB,
    maxit = 2000,
    threshold = 1e-4,
    use.strict = F,
    sparse = T,
    verbose = 1
  )
  weight.rho = flinv(Bs)
  Bs = params.rho$B[[1]]
  fs = params.rho$f
  sigma2 = params.rho$sigma2
  std = centralize_gen(Xs, Ys)
  Xs = std$X           ## Ng (n x sk)
  Ys = std$Y           ## n x ng
  ## multiple
  K   = length(Ys)
  Ns  = sapply(Ys, nrow)
  Bc  = vector("list", K)
  YY  = vector("list", K)
  FX  = vector("list", K)
  Dx  = vector("list", K)
  for (k in 1:K) {
    Bc[[k]] = -Ns[k] * t(solve(diag(Ng) - Bs))
    YY[[k]] = crossprod(Ys[[k]])   ## Y %*% t(Y)
    FX[[k]] = matrix(0, nrow = Ng, ncol = Ns[k])  ## F %*% X (p x k x k x n = p x n)
    for (i in 1:Ng) {
      FX[[k]][i, ] = as.numeric(Xs[[k]][[i]] %*% fs[[k]][[i]])
    }
    Dx[[k]] = abs(Bc[[k]] + ((diag(Ng) - Bs) %*% YY[[k]] - FX[[k]] %*% Ys[[k]]) / sigma2 / weight.rho)
    diag(Dx[[k]]) = -Inf
  }
  ## Dxy = abs((diag(Ng) - Bs) %*% (YY[[2]] - YY[[1]]) - (FX[[2]] %*% Ys[[2]] - FX[[1]] %*% Ys[[1]])) / sigma2 / 2 / weight.rho
  ## diag(Dxy) = -Inf
  max(c(max(Dx[[1]]), max(Dx[[2]])))
  ## max(Dxy)
}

gen_logLik = function(Xs, Ys, Bs, fs, Ng, Ns, K, detIBs, sigma2) {
  std = centralize_gen(Xs, Ys)
  X  = lapply(1:K, function(k) {
    lapply(std$X[[k]], t)
  })
  Y  = lapply(std$Y, t)
  Ls  = 0
  err = 0
  for (k in 1:K) {
    Ls = Ls - Ns[k] / 2 * log(detIBs[k] ^ 2)
    for (i in 1:Ng) {
      Xi = X[[k]][[i]]                     # sk x N
      bi = Bs[[k]][i, -i, drop = F]   # (ng-1) x 1
      yi = Y[[k]][i, , drop = F]      # 1 x N
      Yi = Y[[k]][-i, , drop = F]     # (ng-1) x N
      fi = fs[[k]][[i]]               # sk x 1
      err = err + tcrossprod(yi - bi %*% Yi - crossprod(fi, Xi))
    }
  }
  Ls + err / (2 * sigma2) + Ng * sum(Ns) / 2 * log(sigma2)
}

gen_error = function(Xs, Ys, Bs, fs, Ng, Ns, K) {
  std = centralize_gen(Xs, Ys)
  X  = lapply(1:K, function(k) {
    lapply(std$X[[k]], t)
  })
  Y  = lapply(std$Y, t)
  err = 0
  for (k in 1:K) {
    for (i in 1:Ng) {
      Xi = X[[k]][[i]]                # sk x N
      bi = Bs[[k]][i, -i, drop = F]   # (ng-1) x 1
      yi = Y[[k]][i, , drop = F]      # 1 x N
      Yi = Y[[k]][-i, , drop = F]     # (ng-1) x N
      fi = fs[[k]][[i]]               # sk x 1
      err = err + tcrossprod(yi - bi %*% Yi - crossprod(fi, Xi))
    }
  }
  err
}



## single SML for fused lasso method
#' @description SML method for single problem and fused lasso
cv_SMLasso = function(Bs,
                      fs,
                      Ys,
                      Xs,
                      sigma2,
                      Ng,
                      nlambda = 20,
                      threshold = 1e-4) {
  lambda.max = get_lambda.max(Bs, Ys, Xs, Ng)
  lambda.factors = 10 ^ seq(0,-4, length.out = nlambda) * lambda.max
  wBs = inverse(Bs)
  ncv = 5
  Ns  = sapply(Ys, ncol)
  K   = length(Ys)
  Ytrain  = vector("list", ncv)
  Xtrain  = vector("list", ncv)
  Ytest   = vector("list", ncv)
  Xtest   = vector("list", ncv)
  cv.fold = sample(seq(1, ncv), size = Ns[1], replace = T)
  for (i in 1:ncv) {
    Ytrain[[i]] = lapply(Ys, function(y) {
      y[, cv.fold != i, drop = F]
    })
    Xtrain[[i]] = lapply(Xs, function(x) {
      x[, cv.fold != i, drop = F]
    })
    Ytest[[i]]  = lapply(Ys, function(y) {
      y[, cv.fold == i, drop = F]
    })
    Xtest[[i]] = lapply(Xs, function(x) {
      x[, cv.fold == i, drop = F]
    })
  }
  cverrs = vector("list", nlambda)
  hyper.params = NULL
  for (cv in 1:ncv) {
    params.opt = list()
    Nt = sapply(Ytrain[[cv]], ncol)
    ix = 1
    for (lambda in lambda.factors) {
      cat(sprintf("lambda = %4f, foldid = %d\n", lambda, cv))
      params.opt[[ix]] = vector("list", 2)
      params.opt[[ix]][[1]] = sparse_maximum_likehood_iPALM(
        B = Bs[[1]],
        f = fs[[1]],
        Y = Ytrain[[cv]][[1]],
        X = Xtrain[[cv]],
        sigma2 = sigma2[1],
        N = Nt[1],
        Ng = Ng,
        lambda = lambda,
        maxit = 50,
        verbose = 1,
        threshold = 1e-3
      )
      params.opt[[ix]][[2]] = sparse_maximum_likehood_iPALM(
        B = Bs[[2]],
        f = fs[[2]],
        Y = Ytrain[[cv]][[2]],
        X = Xtrain[[cv]],
        sigma2 = sigma2[1],
        N = Nt[2],
        Ng = Ng,
        lambda = lambda,
        maxit = 50,
        verbose = 1,
        threshold = 1e-3
      )
      Nc = sapply(Ytest[[cv]], ncol)
      err = SML_error(
        Xtest[[cv]],
        Ytest[[cv]],
        list(params.opt[[ix]][[1]]$B, params.opt[[ix]][[2]]$B),
        list(params.opt[[ix]][[1]]$f, params.opt[[ix]][[2]]$f),
        Ng,
        Nc,
        K
      )[1]
      cverrs[[ix]] = c(cverrs[[ix]], err)
      if (cv == 1) {
        hyper.params = c(hyper.params, lambda)
      }
      ix = ix + 1
    }
  }
  list(opt.hyperparams = hyper.params,
       cverrs = do.call(rbind, cverrs))
}

optLasso_cv = function(cvparams, se = TRUE) {
  cvfuns = data.frame(
    lambda = cvparams$opt.hyperparams,
    cvmean = apply(cvparams$cverrs, 1, mean),
    cvsd   = apply(cvparams$cverrs, 1, sd)
  )
  cv.min     = which.min(cvfuns$cvmean)
  cv.1se     = cvfuns[cv.min, 2] + cvfuns[cv.min, 3]
  cvfun.1se  = cvfuns[cvfuns$cvmean <= cv.1se, c(1, 2, 3)]
  lambda = cvfun.1se[, 1]
  lambda.1se = max(lambda)
  if (se) {
    lambda.1se
  } else {
    cvfuns[cv.min, 1]
  }
}


## stability selection
ssSML_iPALM = function(Bs,
                       fs,
                       Ys,
                       Xs,
                       sigma2,
                       Ng,
                       lambda,
                       rho,
                       wBs = inverse(Bs),
                       rB  = flinv(Bs),
                       maxit = 100,
                       acc = TRUE,
                       inertial = inertial_pars("lin"),
                       threshold = 1e-3,
                       sparse = FALSE,
                       use.strict = TRUE,
                       Nbootstrap = 100,
                       Nsample = 0.75,
                       verbose = 2)  {
  N = ncol(Ys[[1]])
  ss.fold = lapply(1:Nbootstrap,
                   function(n) {
                     sample(seq(1, N), ceiling(N * Nsample), replace = F)
                   })
  ss.fit = list()
  N.ss = vector("list", 2)
  D.ss = NULL
  for (i in 1:Nbootstrap) {
    Yss = lapply(Ys, function(Y){Y[, ss.fold[[i]]]})
    Xss = list()
    for(k in 1:length(Xs)) {
      Xss[[k]] = lapply(Xs[[k]], function(X){X[, ss.fold[[k]], drop = F]})
    }
    ss.fit[[i]] = genSML_iPALM(
      Bs = params.init$B,
      fs = params.init$F,
      Ys = Yss,
      Xs = Xss,
      sigma2 = params.init$sigma2[1],
      Ng = data$var$Ng,
      lambda = cvlambda.opt$lambda,
      rho = cvlambda.opt$rho,
      maxit = maxit,
      threshold = threshold,
      use.strict = use.strict,
      acc = acc,
      sparse = sparse
    )
    err2abs = ss.fit[[i]]$err
    if(is.null(N.ss[[1]])) {
      N.ss[[1]] = ifelse(abs(ss.fit[[i]]$B[[1]]) > err2abs, 1, 0)
    } else {
      N.ss[[1]] = N.ss[[1]] + ifelse(abs(ss.fit[[i]]$B[[1]]) > err2abs, 1, 0)
    }
    if(is.null(N.ss[[2]])) {
      N.ss[[2]] = ifelse(abs(ss.fit[[i]]$B[[2]]) > err2abs, 1, 0)
    } else {
      N.ss[[2]] = N.ss[[2]] + ifelse(abs(ss.fit[[i]]$B[[2]]) > err2abs, 1, 0)
    }
    nonzero = c(as.numeric(ss.fit[[i]]$B[[1]]), as.numeric(ss.fit[[i]]$B[[2]]))
    nonzero = nonzero[nonzero != 0]
    thresh.2 = sort(abs(nonzero))[round(0.2 * length(nonzero))+1]
    if(is.null(D.ss)) {
      D.ss = ifelse(abs(ss.fit[[i]]$B[[1]] - ss.fit[[i]]$B[[2]]) > pmin(abs(ss.fit[[i]]$B[[1]]), abs(ss.fit[[i]]$B[[2]])) &
                    abs(ss.fit[[i]]$B[[1]] - ss.fit[[i]]$B[[2]]) > thresh.2, 1, 0)
    } else {
      D.ss = D.ss + ifelse(abs(ss.fit[[i]]$B[[1]] - ss.fit[[i]]$B[[2]]) > pmin(abs(ss.fit[[i]]$B[[1]]), abs(ss.fit[[i]]$B[[2]])) &
                    abs(ss.fit[[i]]$B[[1]] - ss.fit[[i]]$B[[2]]) > thresh.2, 1, 0)
    }
  }
  list(fit = ss.fit, Ns = N.ss, Ds = D.ss)
}

################################################################# adaptive hyper-parameter version ##########################
## cross-validation for adaptive hyper-parameter                                                                          ###
## version 2, added Jun 25                                                                                                ###
##                                                                                                                        ###
#############################################################################################################################
## Shrinkage fused lasso regularizer parameter
#--------------------------------------------------------------------
## cross-validation on adaptive rho(fused lasso params) of lambda
## same function names but different scheme
#--------------------------------------------------------------------
get_rho.max = function(Bs, fs, Ys, Xs, lambda, sigma2, Ng, weighted = TRUE) {
  params.rho = multiSML_iPALM(
    Bs,
    fs,
    Ys,
    Xs,
    sigma2,
    Ng,
    lambda = lambda,
    rho = Inf,
    maxit = 2000,
    threshold = 1e-4,
    use.strict = F,
    sparse = T,
    verbose = 1
  )
  weight.rho = flinv(Bs)
  weight.lambda = inverse(Bs)[[1]]
  Bs = params.rho$B[[1]]
  fs = params.rho$f
  sigma2 = params.rho$sigma2
  std = centralize_mult(Xs, Ys)
  Xs = std$X           ## Ng (n x sk)
  Ys = std$Y           ## n x ng
  ## multiple
  K   = length(Ys)
  Ns  = sapply(Ys, nrow)
  Bc  = vector("list", K)
  YY  = vector("list", K)
  FX  = vector("list", K)
  Dx  = vector("list", K)
  for (k in 1:K) {
    Bc[[k]] = -Ns[k] * t(solve(diag(Ng) - Bs))
    YY[[k]] = crossprod(Ys[[k]])   ## Y %*% t(Y)
    FX[[k]] = matrix(0, nrow = Ng, ncol = Ns[k])  ## F %*% X (p x k x k x n = p x n)
    for (i in 1:Ng) {
      FX[[k]][i, ] = as.numeric(Xs[[i]] %*% fs[[k]][[i]])
    }
    Dx[[k]] = abs(Bc[[k]] + ((diag(Ng) - Bs) %*% YY[[k]] - FX[[k]] %*% Ys[[k]]) / sigma2 + lambda * weight.lambda * sign(Bs)) / weight.rho
    diag(Dx[[k]]) = -Inf
  }
  max(c(max(Dx[[1]]), max(Dx[[2]])))
}

### debug switch
cv_multiSML = function(Bs,
                       fs,
                       Ys,
                       Xs,
                       sigma2,
                       Ng,
                       nlambda = 20,
                       nrho = 20,
                       threshold = 1e-4,
                       logLik = TRUE,
                       verbose = 1) {
  lambda.max = get_lambda.max(Bs, Ys, Xs, Ng)
  lambda.factors = 10 ^ seq(0,-4, length.out = nlambda) * lambda.max
  wBs = inverse(Bs)
  rB  = flinv(Bs)
  ncv = 5
  Ns  = sapply(Ys, ncol)
  K   = length(Ys)
  Ytrain  = vector("list", ncv)
  Xtrain  = vector("list", ncv)
  Ytest   = vector("list", ncv)
  Xtest   = vector("list", ncv)
  cv.fold = sample(seq(1, ncv), size = Ns[1], replace = T)
  for (i in 1:ncv) {
    Ytrain[[i]] = lapply(Ys, function(y) {
      y[, cv.fold != i, drop = F]
    })
    Xtrain[[i]] = lapply(Xs, function(x) {
      x[, cv.fold != i, drop = F]
    })
    Ytest[[i]]  = lapply(Ys, function(y) {
      y[, cv.fold == i, drop = F]
    })
    Xtest[[i]] = lapply(Xs, function(x) {
      x[, cv.fold == i, drop = F]
    })
  }
  cverrs = vector("list", nrho * nlambda)
  cvlls  = vector("list", nrho * nlambda)
  params = NULL
  rho.factors = list()
  il = 1
  for (lambda in lambda.factors) {
    rho.max = get_rho.max(Bs, fs, Ys, Xs, lambda, sigma2, Ng)
    rho.factors[[il]] = 10 ^ seq(0,-4, length.out = nrho) * rho.max
    params = c(params, lapply(rho.factors[[il]], function(rho) {
      c(lambda, rho)
    }))
    il = il + 1
  }
  for (cv in 1:ncv) {
    params.opt = list()
    Nc = sapply(Ytest[[cv]], ncol)
    ix = 1
    il = 1
    for (lambda in lambda.factors) {
      for (rho in rho.factors[[il]]) {
        cat(sprintf("lambda = %4f, rho = %4f, foldid = %d\n", lambda, rho, cv))
        if (ix %% nrho == 1) {
          params.opt[[ix]] = multiSML_iPALM(
            Bs,
            fs,
            Ytrain[[cv]],
            Xtrain[[cv]],
            sigma2,
            Ng,
            lambda,
            rho,
            wBs,
            rB,
            maxit = 1000,
            threshold = threshold,
            acc = TRUE,
            sparse = FALSE,
            use.strict = FALSE,
            verbose = verbose
          )
        } else {
          params.opt[[ix]] = multiSML_iPALM(
            params.opt[[ix - 1]]$B,
            params.opt[[ix - 1]]$f,
            Ytrain[[cv]],
            Xtrain[[cv]],
            params.opt[[ix - 1]]$sigma2,
            Ng,
            lambda,
            rho,
            wBs,
            rB,
            maxit = 1000,
            threshold = threshold,
            acc = TRUE,
            sparse = FALSE,
            use.strict = FALSE,
            verbose = verbose
          )
        }
        loglik = SML_logLik(
            Xtest[[cv]],
            Ytest[[cv]],
            params.opt[[ix]]$B,
            params.opt[[ix]]$f,
            Ng,
            Nc,
            K,
            params.opt[[ix]]$detIB,
            params.opt[[ix]]$sigma2
          )[1]
        err = SML_error(Xtest[[cv]],
                    Ytest[[cv]],
                    params.opt[[ix]]$B,
                    params.opt[[ix]]$f,
                    Ng,
                    Nc,
                    K)[1]
        cverrs[[ix]] = c(cverrs[[ix]], err)
        cvlls[[ix]]  = c(cvlls[[ix]], loglik)
        ix = ix + 1
      } ## rho
      il = il + 1
    } ## lambda
  } ## ncv
  list(opt.hyperparams = do.call(rbind, params),
       cverrs = do.call(rbind, cverrs),
       loglik = do.call(rbind, cvlls))
}

## ultility functions
## pick lambda
optimLambda_cv = function(cvparams, type = c("err", "loglik"), se = TRUE, fused.sparse = TRUE) {
  cvm = if(type == "err") {
    cvparams$cverrs
  } else {
    cvparams$loglik
  }
  cvfuns = data.frame(
    lambda = cvparams$opt.hyperparams[, 1],
    rho    = cvparams$opt.hyperparams[, 2],
    cvmean = apply(cvm, 1, mean),
    cvsd   = apply(cvm, 1, sd)
  )
  cv.min   = which.min(cvfuns$cvmean)
  cvlambda = split(cvfuns, cvfuns$lambda)
  cvms = data.frame(
    lambda = as.numeric(names(cvlambda)),
    cvmean = sapply(cvlambda, function(x){mean(x$cvmean)}),
    cvsd   = sapply(cvlambda, function(x){sum(x$cvsd)/sqrt(length(x))})
  )
  lambda.1se = min_Lambda(cvms$lambda, cvms$cvmean, cvms$cvsd)
  rhos       = cvlambda[[as.character(lambda.1se)]]
  rho.1se    = min_Lambda(rhos$rho, rhos$cvmean, rhos$cvsd)
  rho.min    = rhos$rho[which.min(rhos$cvmean)]
  # ggplot(cvlambda, aes(x = lambda, y = cvmean)) +
  #  geom_errorbar(aes(ymin = cvmean - cvsd, ymax = cvmean + cvsd)) +
  #  geom_line() + geom_point() + scale_x_log10() + xlab(expression(log(lambda))) +
  #  ylab("cv-mean") + ggtitle("cvmean ~ lambda")
  if (se) {
    list(lambda = lambda.1se, rho = rho.1se)
  } else {
    list(lambda = cvfuns[cv.min, 1], rho = cvfuns[cv.min, 2])
  }
}

min_Lambda = function(lambda, cvmean, cvsd) {
  cvmin = min(cvmean, na.rm = TRUE)
  ixmin = cvmean <= cvmin
  lambda.min = max(lambda[ixmin], na.rm = TRUE)
  ixmin = match(lambda.min, lambda)
  ix1se = cvmean <= (cvmean + cvsd)[ixmin]
  max(lambda[ix1se], na.rm = TRUE)
}



################## adaptive version of genSML
######################
# comparable method
######################
## cross validation for hyper-parameter tuning
#' @description 5-fold cross-validation
#' @param dyn dynamic updated rho.max by given lambda
#' @example
#' cv.params = cv_multiSML(Bs = params.init$B, fs = params.init$F, Ys = Ys, Xs = Xs, sigma2 = params.init$sigma2[1], Ng = data$var$Ng, nlambda = 20, nrho = 20)
cv_genSML = function(Bs,
                     fs,
                     Ys,
                     Xs,
                     sigma2,
                     Ng,
                     nlambda = 20,
                     nrho = 20,
                     weighted = TRUE,
                     threshold = 1e-4,
                     use.strict = FALSE,
                     verbose = 1) {
  lambda.max = gen_lambda.max(Bs, Ys, Xs, Ng, weighted)
  lambda.factors = 10 ^ seq(0,-4, length.out = nlambda) * lambda.max
  if(weighted) {
    wBs = inverse(Bs)
    rB  = flinv(Bs)
  } else{
    wBs = invone(Bs)
    rB  = flone(Bs)
  }
  ncv = 5
  Ns  = sapply(Ys, ncol)
  K   = length(Ys)
  Ytrain  = vector("list", ncv)
  Xtrain  = vector("list", ncv)
  Ytest   = vector("list", ncv)
  Xtest   = vector("list", ncv)
  cv.fold = list()
  cv.fold[[1]] = sample(seq(1, ncv), size = Ns[1], replace = T)
  cv.fold[[2]] = sample(seq(1, ncv), size = Ns[2], replace = T)
  for (i in 1:ncv) {
    Ytrain[[i]] = lapply(1:2, function(ix) {
      Ys[[ix]][, cv.fold[[ix]] != i, drop = F]
    })
    Xtrain[[i]] = lapply(1:2, function(ix) {
      lapply(Xs[[ix]], function(x) {
          x[, cv.fold[[ix]] != i, drop = F]
      })
    })
    Ytest[[i]]  = lapply(1:2, function(ix) {
      Ys[[ix]][, cv.fold[[ix]] == i, drop = F]
    })
    Xtest[[i]] = lapply(1:2, function(ix) {
      lapply(Xs[[ix]], function(x) {
          x[, cv.fold[[ix]] == i, drop = F]
      })
    })
  }
  cverrs = vector("list", nrho * nlambda)
  cvlls  = vector("list", nrho * nlambda)
  hyper.params = list()
  rho.factors  = list()
  il = 1
  for (lambda in lambda.factors) {
    rho.max = gen_rho.max(Bs, fs, Ys, Xs, lambda, sigma2, Ng)
    rho.factors[[il]] = 10 ^ seq(0,-4, length.out = nrho) * rho.max
    params = c(params, lapply(rho.factors[[il]], function(rho) {
      c(lambda, rho)
    }))
    il = il + 1
  }

  for (cv in 1:ncv) {
    params.opt = list()
    Nc = sapply(Ytest[[cv]], ncol)
    ix = 1
    il = 1
    for (lambda in lambda.factors) {
      for (rho in rho.factors[[il]]) {
        cat(sprintf("lambda = %4f, rho = %4f, foldid = %d\n", lambda, rho, cv))
        if (ix %% nrho == 1) {
          params.opt[[ix]] = genSML_iPALM(
            Bs,
            fs,
            Ytrain[[cv]],
            Xtrain[[cv]],
            sigma2,
            Ng,
            lambda,
            rho,
            wBs,
            rB,
            maxit = 1000,
            threshold = threshold,
            acc = TRUE,
            sparse = FALSE,
            use.strict = use.strict,
            verbose = verbose
          )
        } else {
          params.opt[[ix]] = genSML_iPALM(
            params.opt[[ix - 1]]$B,
            params.opt[[ix - 1]]$f,
            Ytrain[[cv]],
            Xtrain[[cv]],
            params.opt[[ix - 1]]$sigma2,
            Ng,
            lambda,
            rho,
            wBs,
            rB,
            maxit = 1000,
            threshold = threshold,
            acc = TRUE,
            sparse = FALSE,
            use.strict = use.strict,
            verbose = verbose
          )
        }
        loglik = gen_logLik(
          Xtest[[cv]],
          Ytest[[cv]],
          params.opt[[ix]]$B,
          params.opt[[ix]]$f,
          Ng,
          Nc,
          K,
          params.opt[[ix]]$detIB,
          params.opt[[ix]]$sigma2
        )[1]
        err = gen_error(Xtest[[cv]],
                        Ytest[[cv]],
                        params.opt[[ix]]$B,
                        params.opt[[ix]]$f,
                        Ng,
                        Nc,
                        K)[1]
        cverrs[[ix]] = c(cverrs[[ix]], err)
        cvlls[[ix]]  = c(cvlls[[ix]], loglik)
        if (cv == 1) {
          hyper.params[[ix]] = c(lambda, rho)
        }
        ix = ix + 1
      } ## rho
      il = il + 1
    }
  }
  list(opt.hyperparams = hyper.params,
       cverrs = cverrs,
       loglik = cvlls)
}

## cross-validation and EBIC for hyper-parameter tuning
## rho max can be estimated, rho is the fused lasso
## regularized hyper parameter
#' @description get_rho.max
#' @example
#' rhomax = get_rho.max(params.init$B, params.init$F, Ys, Xs, params.init$sigma2[1], data$var$Ng)
gen_rho.max = function(Bs, fs, Ys, Xs, lambda, sigma2, Ng, weighted = TRUE) {
  if(weighted) {
    wBs = inverse(Bs)
    rB  = flinv(Bs)
  } else {
    wBs = invone(Bs)
    rB  = flone(Bs)
  }
  params.rho = genSML_iPALM(
    Bs,
    fs,
    Ys,
    Xs,
    sigma2,
    Ng,
    lambda = lambda,
    rho = Inf,
    wBs = wBs,
    rB = rB,
    maxit = 2000,
    threshold = 1e-4,
    use.strict = F,
    sparse = T,
    verbose = 1
  )
  weight.rho = flinv(Bs)
  weight.lambda = inverse(Bs)[[1]]
  Bs = params.rho$B[[1]]
  fs = params.rho$f
  sigma2 = params.rho$sigma2
  std = centralize_gen(Xs, Ys)
  Xs = std$X           ## Ng (n x sk)
  Ys = std$Y           ## n x ng
  ## multiple
  K   = length(Ys)
  Ns  = sapply(Ys, nrow)
  Bc  = vector("list", K)
  YY  = vector("list", K)
  FX  = vector("list", K)
  Dx  = vector("list", K)
  for (k in 1:K) {
    Bc[[k]] = -Ns[k] * t(solve(diag(Ng) - Bs))
    YY[[k]] = crossprod(Ys[[k]])   ## Y %*% t(Y)
    FX[[k]] = matrix(0, nrow = Ng, ncol = Ns[k])  ## F %*% X (p x k x k x n = p x n)
    for (i in 1:Ng) {
      FX[[k]][i, ] = as.numeric(Xs[[k]][[i]] %*% fs[[k]][[i]])
    }
    Dx[[k]] = abs(Bc[[k]] + ((diag(Ng) - Bs) %*% YY[[k]] - FX[[k]] %*% Ys[[k]]) / sigma2 + lambda * weight.lambda * sign(Bs)) / weight.rho
    diag(Dx[[k]]) = -Inf
  }
  ## Dxy = abs((diag(Ng) - Bs) %*% (YY[[2]] - YY[[1]]) - (FX[[2]] %*% Ys[[2]] - FX[[1]] %*% Ys[[1]])) / sigma2 / 2 / weight.rho
  ## diag(Dxy) = -Inf
  max(c(max(Dx[[1]]), max(Dx[[2]])))
  ## max(Dxy)
}

############# new version
## ultility functions
## pick lambda
optimLambda_cv = function(cvparams, type = c("err", "loglik"), se = TRUE, fused.sparse = TRUE) {
  cvm = if(type == "err") {
    cvparams$cverrs
  } else {
    cvparams$loglik
  }
  cvfuns = data.frame(
    lambda = cvparams$opt.hyperparams[, 1],
    rho    = cvparams$opt.hyperparams[, 2],
    cvmean = apply(cvm, 1, mean),
    cvsd   = apply(cvm, 1, sd)
  )
  cv.min     = which.min(cvfuns$cvmean)
  cv.1se     = cvfuns[cv.min, 3] + cvfuns[cv.min, 4]
  cvfun.1se  = cvfuns[cvfuns$cvmean <= cv.1se, c(1, 2, 3)]
  lambda = cvfun.1se[, 1]
  rho    = cvfun.1se[, 2]
  if(fused.sparse) {
    rho.1se = max(rho)
    lambda.1se = lambda[which.min(cvfun.1se[cvfun.1se$rho == rho.1se, 3])]
  } else {
    lambda.1se = max(lambda)
    # rho.1se    = min(rho[lambda == lambda.1se])
    rho.1se = rho[which.min(cvfun.1se[cvfun.1se$lambda == lambda.1se, 3])]
  }

  if (se) {
    list(lambda = lambda.1se, rho = rho.1se)
  } else {
    list(lambda = cvfuns[cv.min, 1], rho = cvfuns[cv.min, 2])
  }
}

#############################################
## new version of genSML
#############################################
## solve SML problem by block coordinate descent by backtracking inert-PALM on
## different gene expression and different eQTLs
#' @param lambda lambda hyper-parameter for lasso term
#' @param rho    fused lasso hyper-parameter for fused lasso term
#' @param gamma  invertible matrix stablize parameter gamma
#' params.init = constrained_L2reg_gen(Xs, Ys, sigma2$rho.opt, M, N)
#' params.opt = genSML_iPALM(Bs = params.init$B, fs = params.init$F, Ys = Ys, Xs = Xs,
#'                            sigma2 = params.init$sigma2[1], Ng = data$var$Ng,
#'                            lambda = 0.1, rho = 0.1, maxit = 500)
genSML_iPALM = function(Bs,
                        fs,
                        Ys,
                        Xs,
                        sigma2,
                        Ng,
                        lambda,
                        rho,
                        wBs = inverse(Bs),
                        rB  = flinv(Bs),
                        maxit = 100,
                        acc = TRUE,
                        inertial = inertial_pars("lin"),
                        threshold = 1e-3,
                        sparse = FALSE,
                        use.strict = TRUE,
                        verbose = 2) {
  std = centralize_gen(Xs, Ys)
  Xs  = std$X
  Ys  = std$Y
  meanXs = std$muX
  meanYs = std$muY
  K      = length(Ys)      # number of conditions = 2
  if (verbose == 2) {
    cat(sprintf("conditions must be restricted to 2.. K = %d\n", K))
  }
  ## update for eQTL row-wise (specific for genes)
  f0 = vector("list", K)
  f1 = vector("list", K)
  Yp = vector("list", K)
  Hy = vector("list", K)
  Yp.maxEigen = vector("list", K)
  Ns = sapply(Ys, nrow)
  for (i in 1:Ng) {
    for (k in 1:K) {
      Xi = Xs[[k]][[i]]
      Pi = solve(crossprod(Xi)) %*% t(Xi)
      # specific condition
      yi_k = Ys[[k]][, i, drop = F]              # n[k] x 1 for gene i
      Yi_k = Ys[[k]][,-i]                       # n[k] x (ng-1) (for specific gene i)
      f0[[k]][[i]] = Pi %*% yi_k
      f1[[k]][[i]] = Pi %*% Yi_k                 # f[[k]][[i]] = f0[[k]][[i]] - f1[[k]][[i]] %*% bi^k | bi^k = B[[k]][i,-i]
      Hi_k = diag(Ns[k]) - Xi %*% Pi             # n[k] x n[k] projection matrix
      Yp[[k]][[i]] = t(Yi_k) %*% Hi_k %*% Yi_k   # (ng - 1) x (ng - 1)
      Hy[[k]][[i]] = t(yi_k) %*% Hi_k %*% Yi_k   # 1 x (ng - 1)
      ## maximized eigen-value for Yp
      Yp.maxEigen[[k]][[i]] = eigen(Yp[[k]][[i]])$values[1]
    }
  }
  ## update for gnet row-wise
  niter  = 1
  ImBs   = lapply(Bs, function(B) {
    diag(Ng) - B
  })
  detIBs = sapply(ImBs, det)
  IBsinv = lapply(ImBs, solve)
  Ls     = logLik(detIBs, Bs, wBs, rB, lambda, rho, Ns, K) + N * sum(Ns) / 2 * log(sigma2)
  ## history
  Bs.prevs = list(Bs, Bs)
  inert    = acc
  while (niter <= maxit) {
    inert.pars = inertial(niter)
    Bs.inert = if (inert) {
      lapply(1:K, function(k) {
        Bs.prevs[[2]][[k]] + inert.pars * (Bs.prevs[[2]][[k]] - Bs.prevs[[1]][[k]])
      })
    } else {
      Bs
    }
    fs.prev = fs
    Ls.prev = Ls
    for (i in 1:Ng) {
      ## sum(-Ns[k] * sigma2 * log(det(I-B[[k]])^2)) + \sum_{i=1}^{Ng} bi^T %*% Yp %*% bi - Hy %*% bi)
      ci = lapply(IBsinv, function(IBi) {
        # ci / det(I - B); from (I - B)^{-1}
        IBi[-i, i, drop = F]
      })
      bi = lapply(Bs.inert, function(B.inert) {
        t(B.inert[i,-i, drop = F])
      })
      gi = lapply(1:K, function(k) {
        grad = grad_rwise_SML(Ns[k], ci[[k]], Yp[[k]][[i]], Hy[[k]][[i]], sigma2[1])
        grad(bi[[k]])
      })
      ## Lipschitz moduli for row-i
      Lis  = sapply(1:K, function(k) {
        gtg  = tcrossprod(ImBs[[k]][-i, ])
        oi   = chol2inv(gtg)
        deti = det(gtg)
        gii  = ImBs[[k]][-i,-i]
        si   = ImBs[[k]][-i, i, drop = F]
        c2i  = sum((ci[[k]] * detIBs[k]) ^ 2)
        lips_rwise_SML(Ns[k],
                       oi,
                       gii,
                       si,
                       c2i,
                       deti,
                       Yp.maxEigen[[k]][[i]],
                       sigma2,
                       Ng)[1]
      })
      Li   = max(Lis)
      Li   = (1 + 2 * inert.pars) * Li / (2 * (1 - inert.pars))
      detZero = TRUE
      cl   = 1
      while (detZero) {
        ui   = lapply(1:K, function(k) {
          bi[[k]] - gi[[k]] / Li
        })
        wBi  = lapply(wBs, function(wB) {
          wB[i, -i]
        })
        rBi  = rB[i, -i]
        xi   = prox_flsa(lambda, rho, Li, ui, wBi, rBi)
        dIBu  = sapply(1:K, function(k) {
          IBsinv[[k]][i, i] - (t(xi[[k]]) %*% IBsinv[[k]][-i, i, drop = F])[1]
        })
        cl = cl * 2
        detZero = any(dIBu == 0)
      }
      for (k in 1:K) {
        Bs[[k]][i,-i] = xi[[k]]
        ImBs[[k]]     = diag(Ng) - Bs[[k]]
        detIBs[k]     = (ImBs[[k]][i, ] %*% IBsinv[[k]][, i, drop = F])[1] * detIBs[k]
        dbi           = Bs.prevs[[2]][[k]][i, ] - Bs[[k]][i, ]
        IBsinv[[k]]   = IBsinv[[k]] - IBsinv[[k]][, i, drop = F] %*% dbi %*% IBsinv[[k]] / (1 + dbi %*% IBsinv[[k]][, i, drop = F])[1]
        fs[[k]][[i]]  = f0[[k]][[i]] - f1[[k]][[i]] %*% Bs[[k]][i,-i]
      }
    } # row-wise update
    Berr = sum(sapply(1:K, function(k) {
      norm(Bs[[k]] - Bs.prevs[[2]][[k]], type = "f") / (1 + norm(Bs.prevs[[2]][[k]], type = "f"))
    }))
    Ferr = sum(sapply(1:K, function(k) {
      sum(sapply(1:Ng, function(i) {
        norm(fs[[k]][[i]] - fs.prev[[k]][[i]], type = "f")
      })) / (1 + sum(sapply(1:Ng, function(i) {
        norm(fs.prev[[k]][[i]], type = "f")
      })))
    }))
    err = Berr + Ferr
    sigma2 = sigma2_gen(Xs, Ys, Bs, fs, Ng, Ns, K)
    Ls     = logLik(detIBs, Bs, wBs, rB, lambda, rho, Ns, K) + Ng * sum(Ns) / 2 * log(sigma2)
    Lerr   = abs(Ls.prev - Ls) / (1 + abs(Ls.prev))
    # inert  = ifelse(Lerr < 1e-8, FALSE, acc)
    if (verbose >= 2) {
      cat(
        sprintf(
          "SML: iteration = %d, error = %f, logLik = %f, sigma2 = %f, inert = %s\n",
          niter,
          err,
          Ls,
          sigma2,
          inert
        )
      )
    }
    niter = niter + 1
    Bs.prevs = list(Bs.prevs[[2]], Bs)
    opt.cond = if (use.strict) {
      (err < threshold && Lerr < threshold)
    } else {
      (err < threshold || Lerr < threshold)
    }
    if (opt.cond || niter > maxit || is.nan(err)) {
      mu = lapply(1:K, function(k) {
        (diag(Ng) - Bs[[k]]) %*% meanYs[[k]] - sapply(1:Ng, function(i) {
          meanXs[[k]][[i]] %*% fs[[k]][[i]]
        })
      })
      sigma2 = sigma2_gen(Xs, Ys, Bs, fs, Ng, Ns, K)
      if (sparse) {
        Bs = lapply(Bs, Matrix, sparse = T)
      }
      break
    }
  } # while(niter <= maxit)
  list(
    B = Bs,
    f = fs,
    mu = mu,
    sigma2 = sigma2,
    niter = niter,
    err = err,
    detIB = detIBs
  )
}


############# new version
## ultility functions
## pick lambda
optimLambda_cv1 = function(cvparams,
                          type = c("err", "loglik"),
                          se = TRUE,
                          fused.sparse = TRUE) {
  cvm = if (type == "err") {
    cvparams$cverrs
  } else {
    cvparams$loglik
  }
  cvfuns = data.frame(
    lambda = sapply(cvparams$opt.hyperparams, `[`, 1),
    rho    = sapply(cvparams$opt.hyperparams, `[`, 2),
    cvmean = sapply(cvm, mean),
    cvsd   = sapply(cvm, sd)
  )
  cv.min     = which.min(cvfuns$cvmean)
  cv.1se     = cvfuns[cv.min, 3] + cvfuns[cv.min, 4]
  cvfun.1se  = cvfuns[cvfuns$cvmean <= cv.1se, c(1, 2, 3)]
  lambda = cvfun.1se[, 1]
  rho    = cvfun.1se[, 2]
  if (fused.sparse) {
    rho.1se = max(rho)
    lambda.1se = lambda[which.min(cvfun.1se[cvfun.1se$rho == rho.1se, 3])]
  } else {
    lambda.1se = max(lambda)
    # rho.1se    = min(rho[lambda == lambda.1se])
    rho.1se = rho[which.min(cvfun.1se[cvfun.1se$lambda == lambda.1se, 3])]
  }

  if (se) {
    list(lambda = lambda.1se, rho = rho.1se)
  } else {
    list(lambda = cvfuns[cv.min, 1], rho = cvfuns[cv.min, 2])
  }
}

#### select subset of hyper-parameter with in 1 sd
############# new version
## ultility functions
## pick lambda
subLambda_ss = function(cvparams,
                        type = c("err", "loglik"),
                        ntop = 10) {
  cvm = if (type == "err") {
    cvparams$cverrs
  } else {
    cvparams$loglik
  }
  cvfuns = data.frame(
    lambda = sapply(cvparams$opt.hyperparams, `[`, 1),
    rho    = sapply(cvparams$opt.hyperparams, `[`, 2),
    cvmean = sapply(cvm, mean),
    cvsd   = sapply(cvm, sd)
  )
  cv.min     = which.min(cvfuns$cvmean)
  cv.1se     = cvfuns[cv.min, 3] + cvfuns[cv.min, 4]
  cvfun.1se  = cvfuns[cvfuns$cvmean < cv.1se, c(1, 2, 3, 4)]
  cvfun.1se[order(cvfun.1se[,3], decreasing = F)[1:10],c(1,2)]
}

####################################################
## stability selection on a class of parameters
##################
##' @title ss_fssem
## stability selection
ss_fssem = function(Bs,
                    fs,
                    Ys,
                    Xs,
                    sigma2,
                    Ng,
                    params = NULL,
                    wBs = inverse(Bs),
                    rB  = flinv(Bs),
                    maxit = 100,
                    acc = TRUE,
                    inertial = inertial_pars("lin"),
                    threshold = 1e-3,
                    sparse = FALSE,
                    use.strict = TRUE,
                    Nbootstrap = 100,
                    Nsample = 0.75,
                    verbose = 2)  {
  N = ncol(Ys[[1]])
  ss.fold = vector("list", Nbootstrap)
  i = 1
  while(i <= Nbootstrap) {
    subs = sample(seq(1, N), ceiling(N * Nsample), replace = F)
    ss.fold[[i]] = sort(subs)
    ss.fold[[i+1]] = setdiff(seq(1, N), subs)
    i = i + 2
  }
  ss.fit = list()
  N.ss = vector("list", 2)
  D.ss = NULL
  for (j in 1:nrow(params)) {
    lambda = params[j, 1]
    rho    = params[j, 2]
    for (i in 1:Nbootstrap) {
      Yss = lapply(Ys, function(Y) {
        Y[, ss.fold[[i]]]
      })
      Xss = list()
      for (k in 1:length(Xs)) {
        Xss[[k]] = lapply(Xs[[k]], function(X) {
          X[, ss.fold[[k]], drop = F]
        })
      }
      ss.fit[[i]] = genSML_iPALM(
        Bs = params.init$B,
        fs = params.init$F,
        Ys = Yss,
        Xs = Xss,
        sigma2 = params.init$sigma2[1],
        Ng = data$var$Ng,
        lambda = lambda,
        rho = rho,
        maxit = maxit,
        threshold = threshold,
        use.strict = use.strict,
        acc = acc,
        sparse = sparse
      )
      err2abs = ss.fit[[i]]$err
      if (is.null(N.ss[[1]])) {
        N.ss[[1]] = ifelse(abs(ss.fit[[i]]$B[[1]]) > err2abs, 1, 0)
      } else {
        N.ss[[1]] = N.ss[[1]] + ifelse(abs(ss.fit[[i]]$B[[1]]) > err2abs, 1, 0)
      }
      if (is.null(N.ss[[2]])) {
        N.ss[[2]] = ifelse(abs(ss.fit[[i]]$B[[2]]) > err2abs, 1, 0)
      } else {
        N.ss[[2]] = N.ss[[2]] + ifelse(abs(ss.fit[[i]]$B[[2]]) > err2abs, 1, 0)
      }
      nonzero = c(as.numeric(ss.fit[[i]]$B[[1]]), as.numeric(ss.fit[[i]]$B[[2]]))
      nonzero = nonzero[nonzero != 0]
      thresh.2 = sort(abs(nonzero))[round(0.2 * length(nonzero)) + 1]
      if (is.null(D.ss)) {
        D.ss = ifelse(
          abs(ss.fit[[i]]$B[[1]] - ss.fit[[i]]$B[[2]]) > pmin(abs(ss.fit[[i]]$B[[1]]), abs(ss.fit[[i]]$B[[2]])) &
            abs(ss.fit[[i]]$B[[1]] - ss.fit[[i]]$B[[2]]) > thresh.2,
          1,
          0
        )
      } else {
        D.ss = D.ss + ifelse(
          abs(ss.fit[[i]]$B[[1]] - ss.fit[[i]]$B[[2]]) > pmin(abs(ss.fit[[i]]$B[[1]]), abs(ss.fit[[i]]$B[[2]])) &
            abs(ss.fit[[i]]$B[[1]] - ss.fit[[i]]$B[[2]]) > thresh.2,
          1,
          0
        )
      }
    }
  }
  list(N = N.ss, D = D.ss)
}


bic_SMLasso = function(Bs,
                      fs,
                      Ys,
                      Xs,
                      sigma2,
                      Ng,
                      nlambda = 20,
                      threshold = 1e-3) {
  lambda.max = get_lambda.max(Bs, Ys, Xs, Ng)
  lambda.factors = 10 ^ seq(0,-3, length.out = nlambda) * lambda.max
  wBs = inverse(Bs)
  ncv = 5
  Ns  = sapply(Ys, ncol)
  K   = length(Ys)
  allfit = list()
  berr = vector("list", nlambda)
  ix = 1
  for (lambda in lambda.factors) {
    cat(sprintf("lambda = %4f\n", lambda))
    df = 0
    fit = vector("list", 2)
    for (i in 1:2) {
      fit[[i]] = sparse_maximum_likehood_iPALM(
        B = Bs[[i]],
        f = fs[[i]],
        Y = Ys[[i]],
        X = Xs,
        sigma2 = sigma2,
        N = Ns[i],
        Ng = Ng,
        lambda = lambda,
        maxit = 30,
        verbose = 1,
        threshold = 1e-3
      )
    }
    BIC = SML_BIC(Xs, Ys, fit, Ng, Ns, 2)
    berr[[ix]] = c(lambda, BIC)
    allfit[[ix]] = fit
    ix = ix + 1
  }
  berr = do.call(rbind, berr)
  BICmin = which.min(berr[,2])
  list(lambda = berr[BICmin,1], fit = allfit[[BICmin]])
}

SML_BIC = function(Xs, Ys, fit, Ng, Ns, K) {
  std = centralize_mult(Xs, Ys)
  X  = lapply(std$X, t)
  Y  = lapply(std$Y, t)
  err = 0
  BIC = 0
  Bs = lapply(fit, function(x){x$B})
  fs = lapply(fit, function(x){x$f})
  for (k in 1:K) {
    nll = - Ns[k] / 2 * log(fit[[k]]$detIB**2)
    for (i in 1:Ng) {
      Xi = X[[i]]                     # sk x N
      bi = Bs[[k]][i, -i, drop = F]   # (ng-1) x 1
      yi = Y[[k]][i, , drop = F]      # 1 x N
      Yi = Y[[k]][-i, , drop = F]     # (ng-1) x N
      fi = fs[[k]][[i]]               # sk x 1
      err = err + tcrossprod(yi - bi %*% Yi - crossprod(fi, Xi))
    }
    nll = nll + err / (2 * fit[[k]]$sigma2) + Ng * Ns[k] / 2 * log(2 * pi * fit[[k]]$sigma2)
    BIC = BIC + (2 * nll + sum(Bs[[k]] != 0) * log(Ns[k])) 
  }
  as.numeric(BIC[1,1])
}

sigma2_sml = function(X, Y, B, f, Ng, N) {
  X = lapply(X, t)
  Y = t(Y)
  err = 0
  for (i in 1:Ng) {
    Xi = X[[i]]                     # sk x N
    bi = B[i,-i, drop = F]     # (ng-1) x 1
    yi = Y[i, , drop = F]      # 1 x N
    Yi = Y[-i, , drop = F]     # (ng-1) x N
    fi = f[[i]]                # sk x 1
    err = err + tcrossprod(yi - bi %*% Yi - crossprod(fi, Xi))
  }
  sigma2 = err / (Ng * N)
  sigma2[1]
}




