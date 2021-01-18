vp_cp <- function(Z_mat, B.mat, control) {
  #throughout, Z_mat is the data frame containing test score data
  #and B.mat holds the binary outcome-data and student links
  if(control$cpp.benchmark){
    cpp.rmstep.times=c()
    R.rmstep.times=c()
  }
  fam.binom <- control$outcome.family
  persistence <- control$persistence
  # H.eta is the inverse of \tilde{v} from equation 18 of Karl, Yang, Lohr (2013),
  #Efficient Maximum Likelihood Estimation of multiple membership linear mixed
  #models, with an application to educational value-added assesments, Computational
  #Statistics and Data Analysis. 
  H.eta <- function(G.inv, Z, R.full.inv) {
    h.eta <- G.inv
    h.y <- crossprod(Z, R.full.inv) %*% Z
    as(symmpart(h.eta + h.y), "sparseMatrix")
  }
  #if x is a matrix, ltriangle extracts the lower triangle of x as a vector.
  #if x is a vector, ltriangle(x) builds a symmetric matrix with x as the
  #lower triangle.
  ltriangle <- function(x) {
    if (!is.null(dim(x)[2])) {
      resA <- as.vector(suppressMessages(x[lower.tri(x, diag = TRUE)]))
      resA
    } else {
      nx <- length(x)
      d <- 0.5 * (-1 + sqrt(1 + 8 * nx))
      resB <- .symDiagonal(d)
      resB[lower.tri(resB, diag = TRUE)] <- x
      if (nx > 1) {
        resB <- resB + t(resB) - diag(diag(resB))
      }
      as(resB, "sparseMatrix")
    }
  }
  #If G is a matrix, then reduce.G extracts the unique parametres from G (G is the 
  #block diagonal matrix from Equation (9) of Karl, Yang, Lohr (2013)
  #If G is a vector, then reduce.G builds the full (sparse, block diagonal) matrix
  #G from the list of uniqu parameters.
  reduce.G <- function(G, nyear.score, nteacher) {
    if (!is.null(dim(G)[2])) {
      temp_mat <- G
      index1 <- 0
      resA <- c(NULL)
      for (j in 1:nyear.score) {
        temp_mat_j <- ltriangle(as.matrix(suppressMessages(temp_mat[(index1 + 
                                                                       1):(index1 + 2), (index1 + 1):(index1 + 2)])))
        resA <- c(resA, temp_mat_j)
        index1 <- index1 + 2 * nteacher[j]
      }
      if (control$school.effects) {
        temp_mat_j <- ltriangle(as.matrix(suppressMessages(temp_mat[(index1 + 
                                                                       1):(index1 + 2), (index1 + 1):(index1 + 2)])))
        resA <- c(resA, temp_mat_j)
      }
      resA
    } else {
      resB <- Matrix(0, 0, 0)
      for (j in 1:nyear.score) {
        resB <- bdiag(resB, suppressMessages(kronecker(Diagonal(nteacher[j]), 
                                                       ltriangle(G[(3 * j - 2):(3 * j)]))))
      }
      if (control$school.effects) {
        resB <- bdiag(resB, suppressMessages(kronecker(Diagonal(nschool_effects), 
                                                       ltriangle(G[(3 * (j + 1) - 2):(3 * (j + 1))]))))
      }
      rm(j)
      resB
    }
  }
  #update.eta returns both H.inv (which is \tilde{v} from equation (18) of
  #Karl, Yang, Lohr (2013)) and eta (which is eta from equation (19).
  update.eta <- function(X, Y, Z, R.full.inv, ybetas, G, cons.logLik, Ny, nstudent, 
                         n_eta) {
    G.chol <- chol(G)
    G.inv <- chol2inv(G.chol)
    H <- H.eta(G.inv = G.inv, Z = Z, R.full.inv = R.full.inv)
    chol.H <- chol(H)
    H.inv <- chol2inv(chol.H)
    c.temp <- crossprod(X, R.full.inv) %*% Z
    c.1 <- rbind(crossprod(X, R.full.inv) %*% X, t(c.temp))
    c.2 <- rbind(c.temp, H)
    C_inv <- cbind(c.1, c.2)
    chol.C_inv <- chol(forceSymmetric(symmpart(C_inv)))
    cs <- chol2inv(chol.C_inv)
    C12 <- as.matrix(cs[1:n_ybeta,(n_ybeta+1):ncol(cs)])
    C.mat <- cs[-c(1:n_ybeta), -c(1:n_ybeta)]
    betacov<-as.matrix(cs[c(1:n_ybeta), c(1:n_ybeta)])
    if (control$REML) {
      var.eta <- C.mat
    } else {
      var.eta <- H.inv
    }
    res <- var.eta
    sign.mult <- function(det) {
      det$modulus * det$sign
    }
    eta <- H.inv %*% as.vector(crossprod(Z, R.full.inv) %*% (Y - X %*% ybetas))
    log.p.eta <- -(n_eta/2) * log(2 * pi) - sign.mult(determinant(G.chol)) - 
      0.5 * crossprod(eta, G.inv) %*% eta
    
    log.p.y <- -(Ny/2) * log(2 * pi) + sign.mult(determinant(chol(R.full.inv))) - 
      0.5 * crossprod(Y - X %*% ybetas - Z %*% eta, R.full.inv) %*% (Y - X %*% 
                                                                       ybetas - Z %*% eta)
    if (control$REML) {
      attr(res, "likelihood") <- as.vector(cons.logLik + log.p.eta + log.p.y - 
                                             sign.mult(determinant(chol.C_inv)))
                                            
    } else {
      attr(res, "likelihood") <- as.vector(cons.logLik + log.p.eta + log.p.y - 
                                             sign.mult(determinant(chol.H)))
    }
    attr(res, "eta") <- eta
    attr(res,"betacov")<-betacov
    attr(res,"C12")<- C12
    attr(res, "h.inv") <- H.inv
    res
  }
  #Given a vector, thetas, of current parameter values, Scor() returns a vector
  #of the first-derivative of the likelihood function at that vector
  #This function is not used to produce the parameter estimates, only
  #to calculate the Hessian (matrix of second derivatives). After a few EM steps
  #to improve the initial value, it would be possible to use the Score function in a
  #Newton-Rhaphson algorithm to obtain the parameter estimates. We do not do that,
  #however, since that routine is unstable in the presence of potentially strong
  #correlations in the G matrix
  Score <- function(thetas, eta = eta.hat, ybetas, X, Y, Z, n_ybeta, n_eta, sqrt.W, 
                    inv.sqrt.W, nstudent, nteacher, Kg, cons.logLik, con = control, mis.list, 
                    pattern.parmlist2, pattern.count, pattern.length, pattern.Rtemplate, pattern.diag, 
                    pattern.key, Ny, pattern.sum = pattern.sum, persistence, P, alpha.diag, nalpha, 
                    alpha) {
    n_Rparm <- nyear.pseudo * (nyear.pseudo + 1)/2 - 1
    G <- thetas[(n_Rparm + 1):(n_Rparm + ncol(G.mat))]
    G <- reduce.G(G = G, nyear.score = nyear.score, nteacher = nteacher)
    if (persistence == "VP") {
      alpha.parm <- thetas[(n_Rparm + ncol(G.mat) + 1):length(thetas)]
      alpha[!((1:nalpha) %in% alpha.diag)] <- alpha.parm
      if (!control$school.effects) {
        Z <- Matrix(0, nrow(Z_mat), nteach_effects)
      } else {
        Z <- Matrix(0, nrow(Z_mat), nteach_effects + nschool_effects)
      }
      for (i in 1:nalpha) {
        comp <- which(tril(ltriangle(1:nalpha)) == i, arr.ind = TRUE)
        Z <- Z + alpha[i] * P[[comp[1]]][[comp[2]]]
      }
      
      
      if (control$school.effects) {
        colnames(Z) <- eta_effects
        Z <- Z + Z.school.only
      }
      Z <- drop0(Matrix(Z))
      # BZ structure
      if (!control$school.effects) {
        Z.expand <- Matrix(0, nrow(Z_mat), 2 * nteach_effects)
        Z.expand[, seq(1, 2 * nteach_effects, by = 2)] <- Z
      } else {
        Z.expand <- Matrix(0, nrow(Z_mat), 2 * nteach_effects + 2 * nschool_effects)
        Z.expand[, seq(1, 2 * nteach_effects + 2 * nschool_effects, by = 2)] <- Z
      }
      colnames(Z) <- eta_effects
      J.Z <- rbind(Z.expand, B.Z.expand)
      J.Z <- J.Z[order(J.mat.original$student, J.mat.original$year), ]
      colnames(J.Z) <- interleave(colnames(Z), colnames(B.Z))
      for (p in unique(J.mat$pat)) {
        J.Z.p[[p]] <- J.Z[pat[[p]], , drop = FALSE]
      }
      Z <- J.Z
    }
    R_i <- ltriangle(c(as.vector(thetas[1:n_Rparm]), 1))
    R_i.parm <- c(as.vector(thetas[1:n_Rparm]), 1)
    LRI <- length(R_i.parm)
    R_i.parm.constrained <- R_i.parm[-LRI]
    if (length(mis.list) > 0) {
      R <- symmpart(suppressMessages(kronecker(suppressMessages(Diagonal(nstudent)), 
                                               R_i)[-mis.list, -mis.list]))
      rinv <- Matrix(0, 0, 0)
      R_i.inv <- chol2inv(chol(R_i))
      for (i in 1:nstudent) {
        if (!any(mis.list %in% seq((i - 1) * nyear.pseudo + 1, i * nyear.pseudo))) {
          rinv <- bdiag(rinv, R_i.inv)
        } else {
          inv.indx <- which(!(seq((i - 1) * nyear.pseudo + 1, i * nyear.pseudo) %in% 
                                mis.list))
          rinv <- bdiag(rinv, chol2inv(chol(R_i[inv.indx, inv.indx])))
        }
      }
      R_inv <- rinv
    } else {
      R <- symmpart(suppressMessages(kronecker(suppressMessages(Diagonal(nstudent)), 
                                               R_i)))
      R_inv <- symmpart(suppressMessages(kronecker(suppressMessages(Diagonal(nstudent)), 
                                                   chol2inv(chol(R_i)))))
    }
    R.full <- drop0(symmpart(sqrt.W %*% R %*% sqrt.W))
    R.full.inv <- drop0(symmpart(inv.sqrt.W %*% R_inv %*% inv.sqrt.W))
    new.eta <- update.eta(X = X, Y = Y, Z = Z, R.full.inv = R.full.inv, ybetas = ybetas, 
                          G = G, cons.logLik = cons.logLik, Ny = Ny, nstudent = nstudent, n_eta = n_eta)
    eta.hat <- attr(new.eta, "eta")
    C12 <- attr(new.eta, "C12")
    betacov <- matrix(attr(new.eta, "betacov"),nrow=length(ybetas))
    var.eta.hat <- new.eta
    temp_mat <- var.eta.hat + tcrossprod(eta.hat, eta.hat)
    temp_mat_R <- attr(new.eta, "h.inv") + tcrossprod(eta.hat, eta.hat)
    #### 
    pattern.sum <- list()         
    for (p in unique(J.mat$pat)) {
      if(control$REML){
        pattern.sum[[p]] <- Matrix(REML_Rm(invsqrtW_ = as.matrix(diag(inv.sqrt.W)), betacov_ = betacov, C12_ = C12,
                                          JYp_ = as.matrix(J.Y.p[[p]]), loopsize_ = pattern.count[[p]]/pattern.length[[p]], 
                                          patternlength_ = pattern.length[[p]], rownumber_ = as.matrix(J.Y.p.rownumber[[p]]), 
                                          ybetas_ = as.matrix(ybetas), etahat_ = as.matrix(eta.hat), tempmatR_ = as.matrix(temp_mat), 
                                          JXpi_ = as.matrix(J.X.p[[p]]@i), JXpp_ = as.matrix(J.X.p[[p]]@p), 
                                          JXpx_ = as.matrix(J.X.p[[p]]@x), JXpdim_ = as.matrix(J.X.p[[p]]@Dim), 
                                          JZpi_ = as.matrix(J.Z.p[[p]]@i), JZpp_ = as.matrix(J.Z.p[[p]]@p), 
                                          JZpx_ = as.matrix(J.Z.p[[p]]@x), JZpdim_ = as.matrix(J.Z.p[[p]]@Dim)))
        
      }else{
        
        pattern.sum[[p]] <- Matrix(R_mstep2(invsqrtW_ = as.matrix(diag(inv.sqrt.W)), 
                                            JYp_ = as.matrix(J.Y.p[[p]]), loopsize_ = pattern.count[[p]]/pattern.length[[p]], 
                                            patternlength_ = pattern.length[[p]], rownumber_ = as.matrix(J.Y.p.rownumber[[p]]), 
                                            ybetas_ = as.matrix(ybetas), etahat_ = as.matrix(eta.hat), tempmatR_ = as.matrix(temp_mat_R), 
                                            JXpi_ = as.matrix(J.X.p[[p]]@i), JXpp_ = as.matrix(J.X.p[[p]]@p), 
                                            JXpx_ = as.matrix(J.X.p[[p]]@x), JXpdim_ = as.matrix(J.X.p[[p]]@Dim), 
                                            JZpi_ = as.matrix(J.Z.p[[p]]@i), JZpp_ = as.matrix(J.Z.p[[p]]@p), 
                                            JZpx_ = as.matrix(J.Z.p[[p]]@x), JZpdim_ = as.matrix(J.Z.p[[p]]@Dim)))
        
      }
    }
    #### 
    score.R <- -pattern.f.score.constrained(R_i.parm.constrained = R_i.parm.constrained, 
                                            nyear = nyear.pseudo, pattern.parmlist2 = pattern.parmlist2, pattern.count = pattern.count, 
                                            pattern.length = pattern.length, pattern.Rtemplate = pattern.Rtemplate, 
                                            pattern.diag = pattern.diag, pattern.key = pattern.key, pattern.sum = pattern.sum)
    G.originial <- G
    gam_t <- list()
    sv_gam_t <- list()
    index1 <- 0
    for (j in 1:nyear.score) {
      gam_t[[j]] <- matrix(0, Kg[j], Kg[j])
      G.originial_j <- G.originial[(index1 + 1):(index1 + nteacher[j] * Kg[j]), 
                                   (index1 + 1):(index1 + nteacher[j] * Kg[j])]
      gam_t[[j]] <- G.originial_j[1:Kg[j], 1:Kg[j]]
      sv_gam_t[[j]] <- chol2inv(chol(gam_t[[j]]))
      index1 <- index1 + nteacher[j] * Kg[j]
    }
    if (control$school.effects) {
      gam_t_school <- matrix(0, 2, 2)
      G.originial_j <- G.originial[(index1 + 1):(index1 + nschool_effects * 
                                                   2), (index1 + 1):(index1 + nschool_effects * 2)]
      gam_t_school <- G.originial_j[1:2, 1:2]
      sv_gam_t_school <- chol2inv(chol(gam_t_school))
      index1 <- index1 + nschool_effects * 2
    }
    rm(j)
    gam_t_sc <- list()
    index1 <- 0
    score.G <- Matrix(0, 0, 0)
    for (j in 1:nyear.score) {
      gam_t_sc[[j]] <- matrix(0, Kg[j], Kg[j])
      temp_mat_j <- temp_mat[(index1 + 1):(index1 + nteacher[j] * Kg[j]), (index1 + 
                                                                             1):(index1 + nteacher[j] * Kg[j])]
      index2 <- c(1)
      for (k in 1:nteacher[j]) {
        gam_t_sc[[j]] <- gam_t_sc[[j]] + temp_mat_j[(index2):(index2 + Kg[j] - 
                                                                1), (index2):(index2 + Kg[j] - 1)]
        index2 <- index2 + Kg[j]
      }
      index1 <- index1 + nteacher[j] * Kg[j]
      der <- -0.5 * (nteacher[j] * sv_gam_t[[j]] - sv_gam_t[[j]] %*% gam_t_sc[[j]] %*% 
                       sv_gam_t[[j]])
      if (is.numeric(drop(sv_gam_t[[j]]))) {
        score.eta.t <- der
      } else {
        score.eta.t <- 2 * der - diag(diag(der))
      }
      for (k in 1:nteacher[j]) {
        score.G <- bdiag(score.G, score.eta.t)
      }
    }
    
    if (control$school.effects) {
      gam_t_school <- matrix(0, 2, 2)
      temp_mat_j <- temp_mat[(index1 + 1):(index1 + nschool_effects * 2), (index1 + 
                                                                             1):(index1 + nschool_effects * 2)]
      index2 <- c(1)
      for (k in 1:nschool_effects) {
        gam_t_school <- gam_t_school + temp_mat_j[(index2):(index2 + 2 - 
                                                              1), (index2):(index2 + 2 - 1)]
        index2 <- index2 + 2
      }
      index1 <- index1 + nschool_effects * 2
      der <- -0.5 * (nschool_effects * sv_gam_t_school - sv_gam_t_school %*% 
                       gam_t_school %*% sv_gam_t_school)
      if (is.numeric(drop(sv_gam_t_school))) {
        score.eta.t <- der
      } else {
        score.eta.t <- 2 * der - diag(diag(der))
      }
      for (k in 1:nschool_effects) {
        score.G <- bdiag(score.G, score.eta.t)
      }
    }
    rm(j, k)
    if (persistence == "CP" | persistence == "ZP") {
      -c(score.R, reduce.G(G = score.G, nyear.score = nyear.score, nteacher = nteacher))
    } else if (persistence == "VP") {
      alpha.parm <- alpha[!((1:nalpha) %in% alpha.diag)]
      score.a <- alpha.score(alpha.parm, alpha = alpha, temp_mat = temp_mat_R[seq(1, 
                                                                                  n_eta, 2), seq(1, n_eta, 2)], nalpha = nalpha, alpha.diag = alpha.diag, 
                             P = P, R_inv = R.full.inv[J.mat$response == "score", J.mat$response == 
                                                         "score"], eta.hat = eta.hat[seq(1, n_eta, 2)], ybetas = ybetas, 
                             X = X[J.mat$response == "score", ], Y = Y[J.mat$response == "score"])
      -c(score.R, reduce.G(G = score.G, nyear.score = nyear.score, nteacher = nteacher), 
         score.a)
    }
  }
  #hessian.f() is used to calculate the Hessian for the purpose of calculating
  #standard errors. It is not used to obtain the parameter estimates.
  hessian.f <- function() {
    lgLik.hist <- lgLik
    lgLik <- lgLik[it]
    Hessian <- NA
    std_errors <- c(rep(NA, length(thetas)))
    if (control$hessian == TRUE) {
      if (control$verbose) 
        # cat('Calculating the Hessian...\n')
        flush.console()
      if (control$hes.method == "richardson") {
        Hessian <- ltriangle(ltriangle(jacobian(Score, thetas, eta = eta.hat, 
                                                ybetas = ybetas, X = J.X, Y = J.Y, Z = J.Z, sqrt.W = sqrt.W, inv.sqrt.W = inv.sqrt.W, 
                                                pattern.sum = pattern.sum, con = control, n_ybeta = n_ybeta, n_eta = n_eta, 
                                                nstudent = nstudent, nteacher = nteacher, Kg = Kg, cons.logLik = cons.logLik, 
                                                mis.list = mis.list, pattern.parmlist2 = pattern.parmlist2, pattern.count = pattern.count, 
                                                pattern.length = pattern.length, pattern.Rtemplate = pattern.Rtemplate, 
                                                pattern.diag = pattern.diag, pattern.key = pattern.key, Ny = Ny, 
                                                persistence = persistence, P = P, alpha.diag = alpha.diag, nalpha = nalpha, 
                                                alpha = alpha)))
      } else {
        Hessian <- ltriangle(ltriangle(jacobian(Score, thetas, method = "simple", 
                                                eta = eta.hat, ybetas = ybetas, X = J.X, Y = J.Y, Z = J.Z, sqrt.W = sqrt.W, 
                                                inv.sqrt.W = inv.sqrt.W, pattern.sum = pattern.sum, con = control, 
                                                n_ybeta = n_ybeta, n_eta = n_eta, nstudent = nstudent, nteacher = nteacher, 
                                                Kg = Kg, cons.logLik = cons.logLik, mis.list = mis.list, pattern.parmlist2 = pattern.parmlist2, 
                                                pattern.count = pattern.count, pattern.length = pattern.length, 
                                                pattern.Rtemplate = pattern.Rtemplate, pattern.diag = pattern.diag, 
                                                pattern.key = pattern.key, Ny = Ny, persistence = persistence, 
                                                P = P, alpha.diag = alpha.diag, nalpha = nalpha, alpha = alpha)))
      }
      std_errors <- try(c(sqrt(diag(solve(Hessian)))), silent = TRUE)
      hes.warn <- FALSE
      if (any(eigen(Hessian)$values <= 0)) {
        if (control$verbose) 
          cat("Warning: Hessian not PD", "\n")
        std_errors <- c(rep(NA, length(thetas)))
        hes.warn <- TRUE
      }
    }
    chol2inv(chol(Hessian))
  }
  #update.ybeta is only used to calculate the initial fixed effects vector
  #given inital values for R and eta
  #a different calculation is used throughout the EM algorithm
  #The souce of this function is Equation (15) of Karl, Yang, Lohr (2013)
  update.ybeta <- function(X, Y, Z, R_inv, eta.hat) {
    A.ybeta <- crossprod(X, R_inv) %*% X
    B.ybeta <- crossprod(X, R_inv) %*% (Y - Z %*% eta.hat)
    as.vector(solve(A.ybeta, B.ybeta))
  }
  #bin2dec and dec2 bin are used to move between the columns of the
  #base2 calculations shown in Table 1 of Karl, Yang, Lohr (2013)
  bin2dec <- function(s) sum(s * 2^(rev(seq_along(s)) - 1))
  dec2bin <- function(s) {
    L <- length(s)
    maxs <- max(s)
    digits <- floor(logb(maxs, base = 2)) + 1
    res <- array(NA, dim = c(L, digits))
    for (i in 1:digits) {
      res[, digits - i + 1] <- (s%%2)
      s <- (s%/%2)
    }
    if (L == 1) 
      res[1, ] else res
  }
  #alpha.score gives the score vector of the persistence parameters
  #See section 3.4 of Karl, Yang, Lohr (2013)
  alpha.score <- function(alpha.parm, alpha, temp_mat_R, nalpha, alpha.diag, P, 
                          R_inv, eta.hat, ybetas, X, Y) {
    score.a <- c(NULL)
    alpha.s <- alpha
    alpha.s[!((1:nalpha) %in% alpha.diag)] <- alpha.parm
    alpha.set <- (1:nalpha)[!((1:nalpha) %in% alpha.diag)]
    if (!control$school.effects) {
      Z.a <- Matrix(0, nrow(Z_mat), nteach_effects)
    } else {
      Z.a <- Matrix(0, nrow(Z_mat), nteach_effects + nschool_effects)
    }
    for (i in 1:nalpha) {
      comp <- which(tril(ltriangle(1:nalpha)) == i, arr.ind = TRUE)
      Z.a <- Z.a + alpha.s[i] * P[[comp[1]]][[comp[2]]]
    }
    if (control$school.effects) {
      colnames(Z.a) <- eta_effects
      Z.a <- Z.a + Z.school.only
    }
    for (i in alpha.set) {
      comp <- which(tril(ltriangle(1:nalpha)) == i, arr.ind = TRUE)
      score.a <- c(score.a, as.numeric((t(Y) - t(ybetas) %*% t(X)) %*% R_inv %*% 
                                         P[[comp[1]]][[comp[2]]] %*% eta.hat - sum(diag(crossprod(Z.a, R_inv %*% 
                                                                                                    P[[comp[1]]][[comp[2]]] %*% temp_mat_R)))))
    }
    score.a
  }
  #pattern.f.score.constrained calculates the matrices shown in Equation 16
  #of Karl, Yang, Lohr (2013). The "constrained" refers to the fact that
  #the bottom right element of R is fixed at 1.
  pattern.f.score.constrained <- function(R_i.parm.constrained, nyear, pattern.parmlist2, 
                                          pattern.count, pattern.length, pattern.Rtemplate, pattern.diag, pattern.key, 
                                          pattern.sum) {
    R_i <- as.matrix(ltriangle(c(as.vector(R_i.parm.constrained), 1)))
    pattern.score <- numeric(nyear/2 * (nyear + 1) - 1)
    for (p in nonempty.patterns) {
      pattern.y <- solve(pattern.f.R(R_i, p, nyear, pattern.key))
      YSY <- pattern.y %*% pattern.sum[[p]] %*% pattern.y
      PCL <- pattern.countoverlength[[p]]
      for (r.parm in 1:(nyear/2 * (nyear + 1) - 1)) {
        if (is.null(pat.coord <- pat.coord.guide[[r.parm]][[p]])) 
          next
        pattern.yc <- pattern.y[pat.coord]
        pattern.score[r.parm] <- pattern.score[r.parm] - (PCL * pattern.yc - 
                                                            YSY[pat.coord])
      }
    }
    pattern.score[1:(nyear/2 * (nyear + 1) - 1) %in% pattern.diag] <- 0.5 * pattern.score[1:(nyear/2 * 
                                                                                               (nyear + 1) - 1) %in% pattern.diag]
    -pattern.score
  }
  #pattern.f.R takes a vector of unique parameters and build it into a
  #symmetric matrix of appropriate size. Even though R is block diagonal,
  #not all of the blocks are the same due to missing data. For example,
  #if the year-2 observation is missing, then the second row and column
  #need to be deleted from the particular R-block.
  pattern.f.R <- function(R, p, nyear, pattern.key) {
    R[pattern.key[p, ] * (1:nyear), pattern.key[p, ] * (1:nyear), drop = FALSE]
  }
  #mark the time
  ptm.total <- proc.time()[3]
  #make sure year is numeric
  Z_mat$year <- as.numeric(Z_mat$year)
  #student.delete.list records any students that appear in the data
  #set that do not have any scores recorded
  student.delete.list <- c(NULL)
  for (s in unique(Z_mat$student)) {
    if (sum(!is.na(Z_mat[Z_mat$student == s, "y"])) == 0) 
      student.delete.list <- c(student.delete.list, s)
  }
  #any students in student.delete.list are deleted
  if (length(student.delete.list) > 0) {
    Z_mat <- Z_mat[!Z_mat$student %in% student.delete.list, ]
    B.mat <- B.mat[!B.mat$student %in% student.delete.list, ]
  }
  #how many years of score data are there?
  nyear.score <- length(unique(Z_mat$year))
  #In the program, the binary observations are just treated as the last
  #year of the data set. For example, if there are 3 years of score
  #data, then "year 4" will actually record the binary observations
  #(Think about how the R matrix is structured to see why this makes
  #sense)
  B.year <- nyear.score + 1
  #B.mat holds the binary outcome observations and student links
  B.mat$year <- B.year
  B.mat$mis <- 0
  B.mat$y <- B.mat$r
  B.mat$r <- NULL
  B.mat <- B.mat[order(B.mat$student, decreasing = TRUE), ]
  #Z_mat holds the test score data
  #Z_mat$mis will hold indicators for whether or not the observation missing
  #These will be needed to figure out which subset of the R-blocks are needed
  #for each student (see pattern.f.r, etc).
  Z_mat$mis <- rep(0, dim(Z_mat)[1])
  student_list <- unique(Z_mat$student)
  Z_mat[is.na(Z_mat$y), ]$mis <- rep(1, dim(Z_mat[is.na(Z_mat$y), ])[1])
  #the following loop "fills in the gaps" of missing student scores,
  #making sure that each student has a row in each year
  for (g in 1:nyear.score) {
    mis_stu <- student_list[!(student_list %in% unique(Z_mat[Z_mat$year == g, 
                                                             ]$student))]
    le <- length(mis_stu)
    if (le > 0) {
      temp.exp <- Z_mat[1:le, ]
      temp.exp$year <- rep(g, le)
      temp.exp$mis <- rep(1, le)
      temp.exp$student <- mis_stu
      temp.exp$teacher <- rep(NA, le)
      temp.exp$y <- rep(NA, le)
      Z_mat <- rbind(Z_mat, temp.exp)
    }
  }
  rm(g, le)
  #Z_mat.full will be the current Z_mat matrix, while
  #any missing observations will be deleted from Z_mat
  Z_mat.full <- Z_mat
  Z_mat <- Z_mat[!is.na(Z_mat$y), ]
  Z_mat.full <- Z_mat.full[which((Z_mat.full$student %in% Z_mat$student)), ]
  #create missing indicators for the binary observations
  B.mat[is.na(B.mat$y), ]$mis <- rep(1, dim(B.mat[is.na(B.mat$y), ])[1])
  #the following lines and loop "fill in the blanks" to make sure that
  #each student has a row in B.mat, the data frame of binary outcomes
  mis_stu <- student_list[!(student_list %in% unique(B.mat[B.mat$year == B.year, 
                                                           ]$student))]
  le <- length(mis_stu)
  if (le > 0) {
    temp.exp <- B.mat[1:le, ]
    temp.exp$year <- rep(B.year, le)
    temp.exp$mis <- rep(1, le)
    temp.exp$student <- mis_stu
    temp.exp$y <- rep(NA, le)
    B.mat <- rbind(B.mat, temp.exp)
  }
  rm(le)
  B.mat.full <- B.mat
  B.mat <- B.mat[!is.na(B.mat$y), ]
  nstudent <- length(unique(Z_mat$student))
  #Kg is the dimension of the square matrix of each block of the G matrix
  #E.g. year 1 has a teacher-score variance component and a teacher-outcome 
  #variance compoenent (and an implied covariance). Same for each year.
  Kg <- rep(2, nyear.score)
  Z_mat <- Z_mat[order(Z_mat$year, Z_mat$teacher), ]
  Z_mat.full <- Z_mat.full[order(Z_mat.full$year, Z_mat.full$teacher), ]    
  na_list <- grep("^NA", Z_mat$teacher)
  #teachyearcomb tracks the unique teachers in the data set.
  #if a teacher teaches multiple years in the data set, those are treated
  #as different effects and count as different "teachers"
  if (length(na_list) > 0) {
    teachyearcomb <- unique(cbind(Z_mat[-na_list, ]$year, Z_mat[-na_list, ]$teacher))
  } else {
    teachyearcomb <- unique(cbind(Z_mat$year, Z_mat$teacher))
  }
  Z_mat <- Z_mat[order(Z_mat$student, Z_mat$year, Z_mat$teacher), , drop = FALSE]
  Z_mat.full <- Z_mat.full[order(Z_mat.full$student, Z_mat.full$year, Z_mat.full$teacher), 
                           ]
  nteach_effects <- dim(teachyearcomb)[1]
  teacheffZ_mat <- Matrix(0, nrow = nrow(Z_mat), ncol = nteach_effects)
  t_effects <- rep(NA, nteach_effects)
  indx <- 1
  eblup.tracker <- matrix(0, 0, 3)
  #dP is used as a building block of indicators for the persistence parameters
  dP <- list()
  for (i in 1:nyear.score) {
    dP[[i]] <- list()
    for (j in 1:i) dP[[i]][[j]] <- matrix(0, 0, 2)
  }
  #nalpha is the number of persistence parameters
  nalpha <- nyear.score/2 * (nyear.score + 1)
  if (persistence == "ZP") {
    alpha <- ltriangle(diag(nyear.score))
  } else if (persistence == "CP" | persistence == "VP") {
    #initial values of persistence parameters are all 1
    alpha <- rep(1, nalpha)
  }
  #alpha_key is a lower trianglular matrix that indexes the alpha
  #it tells us where in the matrix of alphas to find each particular alpha
  #E.g. the alpha for the effect of year1 teachers on year2 scores appears in
  #column 1, row 2
  alpha_key <- tril(ltriangle(1:nalpha))
  #which alpha are on the diagonal? E.g. with two years of observations, these
  #these would be alpha 1 and alpha 3. But these are all assumed to equal 1
  alpha.diag <- diag(alpha_key)
  #alpha.parm are the unconstrained alpha (remove the diagonal elements equal
  #to 1
  alpha.parm <- alpha[!((1:nalpha) %in% alpha.diag)]
  #populate dP with persistence information
  for (k in 1:nrow(teachyearcomb)) {
    student_subset <- Z_mat.full$student[Z_mat.full$year == teachyearcomb[k, 
                                                                          1] & Z_mat.full$teacher == teachyearcomb[k, 2]]
    # t_effects[k] <- paste(teachyearcomb[k, 1], '_', teachyearcomb[k, 2], sep = '')
    t_effects[k] <- paste(teachyearcomb[k, 2], sep = "")
    eblup.tracker <- rbind(eblup.tracker, c(teachyearcomb[k, ], teachyearcomb[k, 
                                                                              1]))
    for (yr in as.numeric(teachyearcomb[k, 1]):nyear.score) {
      if (sum(is.element(Z_mat$student, student_subset) & Z_mat$year == yr) != 
          0) {
        q1 <- (1:nrow(Z_mat))[is.element(Z_mat$student, student_subset) & 
                                Z_mat$year == yr & !is.na(Z_mat$y)]
        q2 <- rep(k, length(q1))
        dP[[yr]][[as.numeric(teachyearcomb[k, 1])]] <- rbind(dP[[yr]][[as.numeric(teachyearcomb[k, 
                                                                                                1])]], cbind(q1, q2))
      }
    }
  }
  #school effects calculates an extra variance component for school-level
  #variation. One new instance of this component will be added to G for each
  #school. This information needs to be righ-concatonated to the Z matrix
  if (control$school.effects) {
    school.start <- nteach_effects + 1
    z.school <- t(fac2sparse(Z_mat$schoolID))
    # school.length<-ncol(z.school)
    Z.school.only <- cbind(Matrix(0, nrow(Z_mat), nteach_effects), z.school)
    nschool_effects <- ncol(z.school)
  }
  #P holds indicator matrices (built from dP) to show which persistence
  #parameters are modeled against each row
  P <- list()
  for (i in 1:nyear.score) {
    P[[i]] <- list()
    for (j in 1:i) {
      P[[i]][[j]] <- as(sparseMatrix(i = dP[[i]][[j]][, 1], j = dP[[i]][[j]][, 
                                                                             2], dims = c(nrow(Z_mat), nteach_effects)), "dgCMatrix")
      if (control$school.effects) 
        P[[i]][[j]] <- cbind(P[[i]][[j]], Matrix(0, nrow(Z_mat), nschool_effects))
    }
  }
  if (!control$school.effects) {
    Z <- Matrix(0, nrow(Z_mat), nteach_effects)
  } else {
    Z <- Matrix(0, nrow(Z_mat), nteach_effects + nschool_effects)
  }
  #build persistence information into Z. See the discussion in section
  #3.4 of Karl, Yang, Lohr 2013
  for (i in 1:nalpha) {
    comp <- which(tril(ltriangle(1:nalpha)) == i, arr.ind = TRUE)
    Z <- Z + alpha[i] * P[[comp[1]]][[comp[2]]]
  }
  if (control$school.effects) {
    eta_effects <- colnames(Z) <- colnames(Z.school.only) <- c(t_effects, colnames(z.school))
    colnames(Z) <- eta_effects
    Z <- Z + Z.school.only
  } else {
    eta_effects <- colnames(Z) <- t_effects
  }
  #temp_Z  is used to build the portion of the Z matrix (B.Z, below)
  #corresponding to teacher random effects on the binary outcomes
  temp_Z <- Z[order(Z_mat$year, Z_mat$student, decreasing = TRUE), ]
  temp_Z_mat <- Z_mat[order(Z_mat$year, Z_mat$student, decreasing = TRUE), ]
  temp_Z <- temp_Z[!(duplicated(temp_Z_mat$student)), ]
  temp_Z_mat <- temp_Z_mat[!(duplicated(temp_Z_mat$student)), ]
  temp_Z <- temp_Z[temp_Z_mat$student %in% B.mat$student, ]
  temp_Z_mat <- temp_Z_mat[temp_Z_mat$student %in% B.mat$student, ]
  temp_Z_mat <- temp_Z_mat[, "student", drop = FALSE]
  B.Z <- temp_Z
  colnames(B.Z) <- paste(colnames(temp_Z), "_outcome", sep = "")
  B.mat <- merge(temp_Z_mat, B.mat, by = "student", sort = FALSE)
  B.mat$response <- "outcome"
  Z_mat$response <- "score"
  B.mat.full$response <- "outcome"
  Z_mat.full$response <- "score"
  #J.mat is the "joint" maxtix, combining score and outcome data
  J.mat <- J.mat.original <- rbind(Z_mat[, c("student", "year", "mis", "y", "response")], 
                                   B.mat[, c("student", "year", "mis", "y", "response")])
  J.mat.full <- J.mat.full.original <- rbind(Z_mat.full[, c("student", "year", 
                                                            "mis", "y", "response")], B.mat.full[, c("student", "year", "mis", "y", "response")])
  #nyear.pseudo is nyear+1. The "pseudo" refers to the fact that we treat
  #the binary observations as a pseudo additional-year. This is not related
  #to the "pseudo-likelihood" routine that is used by this program to estimate
  #GLMM    
  nyear.pseudo <- length(unique(J.mat$year))
  #create X for the score submodel
  X_mat <- sparse.model.matrix(control$score.fixed.effects, Z_mat, drop.unused.levels = TRUE)
  X_mat <- X_mat[, !(colSums(abs(X_mat)) == 0), drop = FALSE]
  if (rankMatrix(X_mat, method = "qrLINPACK")[1] != dim(X_mat)[2]) {
    cat("WARNING: Fixed-effects design matrix for test scores is not full-rank", 
        "\n")
    flush.console()
    ANSWER <- readline("Continue anyway? (Y/N)")
    if (substr(ANSWER, 1, 1) == "n" | substr(ANSWER, 1, 1) == "N") 
      stop("WARNING: Fixed-effects design matrix not full-rank")
  }
  #create the X matrix for the binary submodel
  B_X_mat <- sparse.model.matrix(control$outcome.fixed.effects, B.mat, drop.unused.levels = TRUE)
  B_X_mat <- B_X_mat[, !(colSums(abs(B_X_mat)) == 0), drop = FALSE]
  if (rankMatrix(B_X_mat, method = "qrLINPACK")[1] != dim(B_X_mat)[2]) {
    cat("WARNING: Fixed-effects design matrix for outcomes is not full-rank", 
        "\n")
    flush.console()
    ANSWER <- readline("Continue anyway? (Y/N)")
    if (substr(ANSWER, 1, 1) == "n" | substr(ANSWER, 1, 1) == "N") 
      stop("WARNING: Fixed-effects design matrix not full-rank")
  }
  #combine the fixed effects design matrices for both models into a single
  #matrix. 
  J.X <- bdiag(X_mat, B_X_mat)
  colnames(J.X) <- c(paste(colnames(X_mat), "_score", sep = ""), paste(colnames(B_X_mat), 
                                                                       "_outcome", sep = ""))
  #The ".expand" matrices are used to place the elements of Z and B.Z into
  #appropriate places in matrices of equal dimension, so that they may
  #be interleaved (Teacher-1 outcome effect needs to be next to Teacher-1
  #score effect so that G can be block-diagonal
  if (!control$school.effects) {
    Z.expand <- Matrix(0, nrow(Z_mat), 2 * nteach_effects)
    Z.expand[, seq(1, 2 * nteach_effects, by = 2)] <- Z
  } else {
    Z.expand <- Matrix(0, nrow(Z_mat), 2 * nteach_effects + 2 * nschool_effects)
    Z.expand[, seq(1, 2 * nteach_effects + 2 * nschool_effects, by = 2)] <- Z
  }
  if (!control$school.effects) {
    B.Z.expand <- Matrix(0, nrow(B.mat), 2 * nteach_effects)
    B.Z.expand[, seq(2, 2 * nteach_effects, by = 2)] <- B.Z
  } else {
    B.Z.expand <- Matrix(0, nrow(B.mat), 2 * nteach_effects + 2 * nschool_effects)
    B.Z.expand[, seq(2, 2 * nteach_effects + 2 * nschool_effects, by = 2)] <- B.Z
  }
  #J.Z will be the master Z matrix going forward. We are modelling
  #J.Y = J.X * ybetas + J.Z * eta.hat
  J.Z <- rbind(Z.expand, B.Z.expand)
  J.Z <- J.Z[order(J.mat$student, J.mat$year), ]
  J.X <- J.X[order(J.mat$student, J.mat$year), ]
  interleave <- function(v1, v2) {
    ord1 <- 2 * (1:length(v1)) - 1
    ord2 <- 2 * (1:length(v2))
    c(v1, v2)[order(c(ord1, ord2))]
  }
  colnames(J.Z) <- interleave(colnames(Z), colnames(B.Z))
  J.mat.full <- J.mat.full[order(J.mat.full$student, J.mat.full$year), ]
  J.mat <- J.mat[order(J.mat$student, J.mat$year), ]
  mis.list <- which(J.mat.full$mis == 1)
  #nteacher is a vector of length equal to the number of years of score data
  #that records how many teachers taught in each year
  nteacher <- as.vector(tapply(teachyearcomb[, 2], teachyearcomb[, 1], length))
  n_eta <- ncol(J.Z)
  #ybeta would be perhaps more approriately named "jbeta" since it contains
  #both the score and the outcome fixed effects
  n_ybeta <- ncol(J.X)
  J.Y <- as.vector(J.mat$y)
  #indicator for if observation present
  J.mat.full$r <- 1 - J.mat.full$mis
  #pattern.student: each student has a row of indicators showing which year present
  #see Table1 of Karl, Yang, Lohr (2013)
  pattern.student <- matrix(J.mat.full$r, nstudent, nyear.pseudo, byrow = TRUE)
  #for each row of J.mat, what missingness pattern the student in that row
  #belongs to
  #all of the pattern variables below are used to track the attendence
  #patterns discussed in Section 4 of Karl, Yang, Lohr (2013)
  J.mat.full$pat <- rep(apply(pattern.student, 1, bin2dec), each = nyear.pseudo)
  J.mat$pat <- J.mat.full[J.mat.full$r == 1, ]$pat
  pat <- list()
  pattern.count <- list()
  pattern.length <- list()
  pattern.countoverlength <- list()
  J.Y.p <- list()
  J.X.p <- list()
  J.Z.p <- list()
  J.Y.p.rownumber <- list()
  J.mat$rownumber <- 1:nrow(J.mat)
  pattern.key <- dec2bin(1:(2^nyear.pseudo - 1))
  for (p in unique(J.mat$pat)) {
    pat[[p]] <- which(J.mat$pat == p)
    J.X.p[[p]] <- J.X[pat[[p]], , drop = FALSE]
    J.Z.p[[p]] <- J.Z[pat[[p]], , drop = FALSE]
    J.Y.p[[p]] <- J.Y[pat[[p]]]
    J.Y.p.rownumber[[p]] <- J.mat$rownumber[pat[[p]]]
    pattern.count[[p]] <- length(J.Y.p[[p]])
    pattern.length[[p]] <- sum(pattern.key[p, ])
    pattern.countoverlength[[p]] <- pattern.count[[p]]/pattern.length[[p]]
  }
  pattern.yguide <- list()
  for (g in 1:nyear.pseudo) {
    pattern.yguide[[g]] <- which(pattern.key[, g] == 1)
  }
  pattern.Rtemplate <- ltriangle(1:(nyear.pseudo/2 * (nyear.pseudo + 1)))
  pattern.diag <- diag(pattern.Rtemplate)
  pattern.parmlist1 <- list()
  pattern.parmlist2 <- list()
  for (p in unique(J.mat$pat)) {
    pattern.parmlist1[[p]] <- sort(unique(as.vector(pattern.f.R(pattern.Rtemplate, 
                                                                p, nyear.pseudo, pattern.key))))
  }
  for (r.parm in 1:(nyear.pseudo/2 * (nyear.pseudo + 1))) {
    pattern.parmlist2[[r.parm]] <- which(sapply(pattern.parmlist1, f <- function(x) r.parm %in% 
                                                  x))
  }
  #eta.hat is eta from Karl, Yang, Lohr (2013) (The vector of random effects)
  eta.hat <- numeric(n_eta)
  #var.eta.hat is the quantity defined in Equation (18) of
  #Karl, Yang, Lohr (2013)
  var.eta.hat <- Matrix(0, n_eta, n_eta)
  #R.temp.comp are the initial estimates for the diagonal elements of R
  R.temp.comp <- rep(1, nyear.pseudo)
  for (g in 1:nyear.score) {
    R.temp.comp[g] <- var(J.mat[J.mat$year == g, ]$y)/4
  }
  #R_i is the block of R corresponding to a student with no missing data.
  #R is built from the blocks of R with appropriate rows and columns deleted
  #where observations are missing
  R_i <- diag(R.temp.comp)
  if (length(mis.list) > 0) {
    R <- symmpart(suppressMessages(kronecker(suppressMessages(Diagonal(nstudent)), 
                                             R_i)[-mis.list, -mis.list]))
    rinv <- Matrix(0, 0, 0)
    R_i.inv <- chol2inv(chol(R_i))
    for (i in 1:nstudent) {
      if (!any(mis.list %in% seq((i - 1) * nyear.pseudo + 1, i * nyear.pseudo))) {
        rinv <- bdiag(rinv, R_i.inv)
      } else {
        inv.indx <- which(!(seq((i - 1) * nyear.pseudo + 1, i * nyear.pseudo) %in% 
                              mis.list))
        rinv <- bdiag(rinv, chol2inv(chol(R_i[inv.indx, inv.indx])))
      }
    }
    R_inv <- rinv
  } else {
    R <- symmpart(suppressMessages(kronecker(suppressMessages(Diagonal(nstudent)), 
                                             R_i)))
    R_inv <- symmpart(suppressMessages(kronecker(suppressMessages(Diagonal(nstudent)), 
                                                 chol2inv(chol(R_i)))))
  }
  #the only use of update.ybeta() in the program. This calculates the initial
  #values of the fixed effects (for both score and response submodels).
  #calculation for score model fixed effects
  ybetas <- update.ybeta(J.X, J.Y, J.Z, R_inv, eta.hat)
  #calculation for binary model fixed effects
  ybetas[(ncol(X_mat) + 1):n_ybeta] <- glm.fit(x = B_X_mat, y = B.mat$y, family = binomial("probit"))[[1]]
  names(ybetas) <- colnames(J.X)
  #initial value for G. A digaonal matrix with 100 on the diagonal
  G <- 100 * .symDiagonal(n_eta)
  #if(control$independent.responses)   G <- 0.01 * .symDiagonal(n_eta)
  #number of observations (both score and binary).
  Ny <- sum(J.mat$mis == 0)
  #pat.coord.guide tells where to find each of the R parameters in each of the
  #possible blocks of R (across all possible missing data patterns).
  pat.coord.guide <- list()
  for (r.parm in 1:(nyear.pseudo/2 * (nyear.pseudo + 1) - 1)) {
    pat.coord.guide[[r.parm]] <- list()
    for (p in pattern.parmlist2[[r.parm]]) {
      pat.coord.guide[[r.parm]][[p]] <- which(tril(pattern.f.R(pattern.Rtemplate, 
                                                               p, nyear.pseudo, pattern.key)) == r.parm)
    }
  }
  nonempty.patterns <- NULL
  for (p in 1:length(pat.coord.guide[[1]])) {
    if (is.null(retrieve.parm <- pattern.parmlist1[[p]])) 
      next
    nonempty.patterns <- c(nonempty.patterns, p)
  }
  #the loglikelihood function. Do not use this for likelihood ratio tests, 
  #which are not valid for pseudo-likelihood estimation
  if (control$REML) {
    cons.logLik <- 0.5 * (n_eta + n_ybeta) * log(2 * pi)
  } else {
    cons.logLik <- 0.5 * (n_eta) * log(2 * pi)
  }
  J.mat$fit <- J.mat$mu <- J.mat$mu.eta.val <- J.mat$nu <- J.mat$sqrt.w <- NA
  #J.mat$fit is the vector of estimated responses (and psuedoresponses)
  J.mat[J.mat$response == "outcome", ]$fit <- as.vector(J.X[J.mat$response == "outcome", 
                                                            ] %*% ybetas + J.Z[J.mat$response == "outcome", ] %*% eta.hat)
  #J.mat$mu is defined on page 236 of Wolfinger and O'Connell (1993)
  J.mat[J.mat$response == "outcome", ]$mu <- as.vector(fam.binom$linkinv(J.mat[J.mat$response == 
                                                                                 "outcome", ]$fit))
  #J.mat$mu.eta.val is the value of the first derivative of the link function
  #at each predicted point. This is g' on page 236 of Wolfinger and O'Connell (1993)
  J.mat[J.mat$response == "outcome", ]$mu.eta.val <- fam.binom$mu.eta(J.mat[J.mat$response == 
                                                                              "outcome", ]$fit)
  #J.mat$nu is nu from Step 2 in Section 3 of Wolfinger and O'Connell (1993)
  J.mat[J.mat$response == "outcome", ]$nu <- as.vector(J.mat[J.mat$response == 
                                                               "outcome", ]$fit) + (J.mat[J.mat$response == "outcome", ]$y - J.mat[J.mat$response == 
                                                                                                                                     "outcome", ]$mu)/J.mat[J.mat$response == "outcome", ]$mu.eta.val
  #In Wolfinger and O'Connell (1993), diag(J.mat$sqrt.w) is called R_{\mu}
  J.mat[J.mat$response == "outcome", ]$sqrt.w <- sqrt(fam.binom$variance(J.mat[J.mat$response == 
                                                                                 "outcome", ]$mu)/J.mat[J.mat$response == "outcome", ]$mu.eta.val^2)
  #no transformation needed for the normal-distributed scores. The
  #pseudoresponses are the responses
  J.mat[J.mat$response == "score", ]$sqrt.w <- 1
  J.mat[J.mat$response == "score", ]$nu <- J.mat[J.mat$response == "score", ]$y
  sqrt.W <- Diagonal(x = J.mat$sqrt.w)
  inv.sqrt.W <- Diagonal(x = 1/(J.mat$sqrt.w))
  #R.full is the error covariance matrix, which is formed as the product
  #diag(sqrt.w)%*%R%*%diag(sqrt.w). The matrix R assumes a variance of 1
  #for all of the binomial responses, while R.full includes the variance from the
  #binomial distribution (in Wolfinger (1993), diag(sqrt.w) is called R_mu).
  R.full <- drop0(symmpart(sqrt.W %*% R %*% sqrt.W))
  R.full.inv <- drop0(symmpart(inv.sqrt.W %*% R_inv %*% inv.sqrt.W))
  #replace vector of responses with vector of updated pseudoresponses
  J.Y <- J.mat$nu
  for (p in unique(J.mat$pat)) {
    J.Y.p[[p]] <- J.Y[pat[[p]]]
  }
  for (PQL.it in 1:control$max.PQL.it) {
    #begin the outer pseudo-likelihood loop
    #    if (PQL.it <= control$EM.to.NR) {
    EM.iter <- control$max.iter.EM[1]
    #   } else {
    #       EM.iter <- control$max.iter.EM[2]
    #   }
    Y.mat <- Matrix(0, EM.iter + 1, n_ybeta)
    G.mat <- Matrix(0, EM.iter + 1, length(reduce.G(G = G, nyear = nyear.score, 
                                                    nteacher = nteacher)))
    R.mat <- Matrix(0, EM.iter + 1, nyear.pseudo * (nyear.pseudo + 1)/2)
    lgLik <- numeric(EM.iter + 1)
    conv <- FALSE
    if (control$verbose) 
      flush.console()
    #There used to be a conditional test in place of if(TRUE): no longer
    #needed
    if (TRUE) {
      #begin the (inner) EM algorithm loop
      #note control$EM.to.NR is fixed at a large value and not relevant
      cat("Beginning EM algorithm\n")
      for (it in 1:1000) {
        #loop will exit before 1000 reached due to stopping conditions
        ptm <- proc.time()[3]
        suppressWarnings(rm(var.eta.hat, temp_mat, temp_mat_R))
        mresid <- as.numeric(J.Y - J.X %*% ybetas)
        cresid <- as.numeric(mresid - J.Z %*% eta.hat)
        yhat <- as.numeric(J.X %*% ybetas + J.Z %*% eta.hat)
        yhat.m <- as.numeric(J.X %*% ybetas)
        #The E-step. Calculate the new random effects
        new.eta <- update.eta(X = J.X, Y = J.Y, Z = J.Z, R.full.inv = R.full.inv, 
                              ybetas = ybetas, G = G, cons.logLik = cons.logLik, Ny = Ny, nstudent = nstudent, 
                              n_eta = n_eta)
        eta.hat <- attr(new.eta, "eta")
        C12 <- attr(new.eta, "C12")
        betacov <- matrix(attr(new.eta, "betacov"),nrow=length(ybetas))
        var.eta.hat <- new.eta
        #temp_mat is the value inside the summand of Equation (14) of
        #Karl, Yang,  Lohr (2013)
        temp_mat <- var.eta.hat + tcrossprod(eta.hat, eta.hat)
        temp_mat_R <- attr(new.eta, "h.inv") + tcrossprod(eta.hat, eta.hat)
        loglikelihood <- attr(new.eta, "likelihood")
        lgLik[it] <- attr(new.eta, "likelihood")
        Y.mat[it, ] <- ybetas
        R.mat[it, ] <- suppressMessages(ltriangle(R_i))
        G.mat[it, ] <- suppressMessages(reduce.G(G, nyear.score, nteacher))
        rm(new.eta)
        thets1 <- c(Y.mat[it - 1, ], R.mat[it - 1, ], G.mat[it - 1, ])
        thets2 <- c(Y.mat[it, ], R.mat[it, ], G.mat[it, ])
        if (it >= 2) {
          # if (persistence == "CP" | persistence == "ZP") {
          #   thetas <- c(ltriangle(as.matrix(R_i))[-LRI], reduce.G(G = G, nyear.score = nyear.score, 
          #                                                         nteacher = nteacher))
          # } else if (persistence == "VP") {
          #   thetas <- c(ltriangle(as.matrix(R_i))[-LRI], reduce.G(G = G, nyear.score = nyear.score, 
          #                                                         nteacher = nteacher), alpha[!((1:nalpha) %in% alpha.diag)])
          # }
          # gradient<-Score(thetas,
          #                 eta = eta.hat, ybetas = ybetas, X = J.X, Y = J.Y, Z = J.Z, sqrt.W = sqrt.W, 
          #                 inv.sqrt.W = inv.sqrt.W, pattern.sum = pattern.sum, con = control, 
          #                 n_ybeta = n_ybeta, n_eta = n_eta, nstudent = nstudent, nteacher = nteacher, 
          #                 Kg = Kg, cons.logLik = cons.logLik, mis.list = mis.list, pattern.parmlist2 = pattern.parmlist2, 
          #                 pattern.count = pattern.count, pattern.length = pattern.length, pattern.Rtemplate = pattern.Rtemplate, 
          #                 pattern.diag = pattern.diag, pattern.key = pattern.key, Ny = Ny, 
          #                 persistence = persistence, P = P, alpha.diag = alpha.diag, nalpha = nalpha, 
          #                 alpha = alpha)
          #stopping condition for EM (inner) loop
          check.lik <- abs(lgLik[it] - lgLik[it - 1])/abs(lgLik[it] + control$tol1) < 
            control$tol1
          if (check.lik | it == (EM.iter)) {
            conv <- TRUE
            EM.conv <- TRUE
            if (control$verbose) {
              cat("\n\n End EM Algorithm.\n")
              cat("\n\niter:", it, "\n")
              cat("log-likelihood:", sprintf("%.7f", lgLik[it]), "\n")
              cat("change in loglik:", sprintf("%.7f", lgLik[it] - lgLik[it - 
                                                                           1]), "\n")
              cat("fixed effects:", round(ybetas, 4), "\n")
              cat("R_i:\n")
              print(round(as.matrix(R_i), 4))
              cat("\n")
              print(round(cov2cor(as.matrix(R_i)), 4))
              cat("\n")
              for (j in 1:nyear.score) {
                cat("\ngamma_teach_year", j, "\n")
                print(round(ltriangle(reduce.G(G, nyear.score, nteacher)[(3 * 
                                                                            j - 2):(3 * j)]), 4))
                cat("\n")
                print(round(cov2cor(ltriangle(reduce.G(G, nyear.score, nteacher)[(3 * 
                                                                                    j - 2):(3 * j)])), 4))
                cat("\n")
                flush.console()
              }
              if (control$school.effects) {
                j <- j + 1
                cat("\n")
                cat("\ngamma_school_effects\n")
                print(round(ltriangle(reduce.G(G, nyear.score, nteacher)[(3 * 
                                                                            j - 2):(3 * j)]), 4))
                cat("\n")
                print(round(cov2cor(ltriangle(reduce.G(G, nyear.score, nteacher)[(3 * 
                                                                                    j - 2):(3 * j)])), 4))
                cat("\n")
                flush.console()
              }
              cat("\n")
              cat("alphas (persistence parameters):\n")
              print(round(alpha, 4))
              rm(j)
            }
            break
          }
        }
        if ((control$verbose) & (it > 1)) {
          cat("\n\niter:", it, "\n")
          cat("log-likelihood:", sprintf("%.7f", lgLik[it]), "\n")
          cat("change in loglik:", sprintf("%.7f", lgLik[it] - lgLik[it - 
                                                                       1]), "\n")
          cat("fixed effects:", round(ybetas, 4), "\n")
          cat("R_i:\n")
          print(round(as.matrix(R_i), 4))
          cat("\n")
          print(round(cov2cor(as.matrix(R_i)), 4))
          cat("\n")
          cat("G:\n")
          for (j in 1:nyear.score) {
            cat("\ngamma_teach_year", j, "\n")
            print(round(ltriangle(reduce.G(G, nyear.score, nteacher)[(3 * 
                                                                        j - 2):(3 * j)]), 4))
            cat("\n")
            print(round(cov2cor(ltriangle(reduce.G(G, nyear.score, nteacher)[(3 * 
                                                                                j - 2):(3 * j)])), 4))
            cat("\n")
            flush.console()
          }
          if (control$school.effects) {
            j <- j + 1
            cat("\n")
            cat("\ngamma_school_effects\n")
            print(round(ltriangle(reduce.G(G, nyear.score, nteacher)[(3 * 
                                                                        j - 2):(3 * j)]), 4))
            cat("\n")
            print(round(cov2cor(ltriangle(reduce.G(G, nyear.score, nteacher)[(3 * 
                                                                                j - 2):(3 * j)])), 4))
            cat("\n")
            flush.console()
          }
          cat("\n")
          cat("alphas (persistence parameters):\n")
          print(round(alpha, 4))
          rm(j)
        }
        flush.console()
        #M step follows
        #M step update for persistence parameters requies a single
        #newton step
        if (persistence == "VP") {
          alpha.parm <- alpha[!((1:nalpha) %in% alpha.diag)]
          alpha.parm.old <- alpha.parm
          s.prev <- numeric(length(alpha.parm))
          # one step is sufficient because score function is linear in alpha
          s <- alpha.score(alpha.parm, alpha = alpha, temp_mat_R = temp_mat_R[seq(1, 
                                                                                  n_eta, 2), seq(1, n_eta, 2)], nalpha = nalpha, alpha.diag = alpha.diag, 
                           P = P, R_inv = R.full.inv[J.mat$response == "score", J.mat$response == 
                                                       "score"], eta.hat = eta.hat[seq(1, n_eta, 2)], ybetas = ybetas, 
                           X = J.X[J.mat$response == "score", ], Y = J.Y[J.mat$response == 
                                                                           "score"])
          j <- symmpart(jacobian(alpha.score, alpha.parm, method = "simple", alpha = alpha, 
                                 temp_mat_R = temp_mat_R[seq(1, n_eta, 2), seq(1, n_eta, 2)], 
                                 nalpha = nalpha, alpha.diag = alpha.diag, P = P, R_inv = R.full.inv[J.mat$response == 
                                                                                                       "score", J.mat$response == "score"], eta.hat = eta.hat[seq(1, 
                                                                                                                                                                  n_eta, 2)], ybetas = ybetas, X = J.X[J.mat$response == "score", 
                                                                                                                                                                                                       ], Y = J.Y[J.mat$response == "score"]))
          if(control$independent.responses){
            indep.key=ltriangle(1:LRI)
            coord.to.delete=indep.key[nrow(indep.key),-ncol(indep.key)]
            for(CTD.ind in 1:length(coord.to.delete)){
              j[coord.to.delete[CTD.ind],]<-j[,coord.to.delete[CTD.ind]]<-numeric(LRI-1)
              j[coord.to.delete[CTD.ind],coord.to.delete[CTD.ind]]<-1
              s[coord.to.delete[CTD.ind]]<-0
            }
          }
          hesprod <- solve(j, s)
          alpha.parm <- alpha.parm - hesprod
          alphan <- alpha
          alphan[!((1:nalpha) %in% alpha.diag)] <- alpha.parm
        }
        # end update alphas
        #M-step update for G. See Equation (14) of Karl, Yang, Lohr(2013) 
        indx <- 1
        Gn <- c(NULL)
        gt <- list()
        for (i in 1:nyear.score) {
          gt[[i]] <- matrix(0, 2, 2)
          for (j in 1:nteacher[i]) {
            gt[[i]] <- gt[[i]] + temp_mat[(2 * (indx - 1) + 1):(2 * indx), 
                                          (2 * (indx - 1) + 1):(2 * indx)]
            indx <- indx + 1
          }
          if(control$independent.responses) gt[[i]][2,1]<- gt[[i]][1,2]<-0             
        }
        if (control$school.effects) {
          gt.school <- matrix(0, 2, 2)
          for (j in 1:nschool_effects) {
            gt.school <- gt.school + temp_mat[(2 * (indx - 1) + 1):(2 * indx), 
                                              (2 * (indx - 1) + 1):(2 * indx)]
            indx <- indx + 1
          }
        }
        Gn <- Matrix(0, 0, 0)
        for (i in 1:nyear.score) {
          Gn <- bdiag(Gn, suppressMessages(kronecker(Diagonal(nteacher[i]), 
                                                     gt[[i]]/nteacher[i])))
        }
        if (control$school.effects) {
          Gn <- bdiag(Gn, suppressMessages(kronecker(Diagonal(nschool_effects), 
                                                     gt.school/nschool_effects)))
        }
        ################7_16_17
        #The code in this if statement is not used in the calculations.
        #This is older code that has been replaced by C++ code. This code is
        #left in the program for running benchmarks against the C++ code.
        if(control$cpp.benchmark){
          pattern.sum <- list()
          J.Y.p <- list()
          t.J.X.p <- list()
          t.J.Z.p <- list()
          J.Y.p.rownumber <- list()
          J.mat$rownumber <- 1:nrow(J.mat)
          pattern.key <- dec2bin(1:(2^nyear.pseudo - 1))
          for (p in unique(J.mat$pat)) {
            pat[[p]] <- which(J.mat$pat == p)
            t.J.X.p[[p]] <- t(J.X[pat[[p]], , drop = FALSE] )
            t.J.Z.p[[p]] <- t(J.Z[pat[[p]], , drop = FALSE])
            J.Y.p[[p]] <- J.Y[pat[[p]]]
            J.Y.p.rownumber[[p]] <- J.mat$rownumber[pat[[p]]]
            pattern.count[[p]] <- length(J.Y.p[[p]])
            pattern.length[[p]] <- sum(pattern.key[p, ])
          }
          rtime<-proc.time()[3]
          
          if(control$REML){
            for (p in unique(J.mat$pat)) {
              pattern.sum[[p]] <- matrix(0, pattern.length[[p]], pattern.length[[p]])
              for (i in 1:(pattern.count[[p]]/pattern.length[[p]])) {
                X.t <- t(t.J.X.p[[p]][,(1 + (i - 1) * pattern.length[[p]]):(i * pattern.length[[p]]) , drop = FALSE])
                Y.t <- J.Y.p[[p]][(1 + (i - 1) * pattern.length[[p]]):(i * pattern.length[[p]])]
                Z.t <- t(t.J.Z.p[[p]][,(1 + (i - 1) * pattern.length[[p]]):(i * pattern.length[[p]]) , drop = FALSE])
                rownumber.t <- J.Y.p.rownumber[[p]][(1 + (i - 1) * pattern.length[[p]]):(i * pattern.length[[p]])]
                temp.t <- Y.t - X.t %*% ybetas
                pattern.sum[[p]] <- pattern.sum[[p]] + inv.sqrt.W[rownumber.t, rownumber.t] %*% (tcrossprod(temp.t) - tcrossprod(temp.t, Z.t %*% eta.hat) - tcrossprod(Z.t %*% eta.hat, temp.t) + Z.t %*%
                                                                                                   tcrossprod(temp_mat, Z.t)+X.t %*%tcrossprod(betacov, X.t)+X.t %*%tcrossprod(C12, Z.t)+Z.t %*%tcrossprod(t(C12), X.t)) %*% inv.sqrt.W[rownumber.t, rownumber.t]
              }
            }
          }else{
            for (p in unique(J.mat$pat)) {
              pattern.sum[[p]] <- matrix(0, pattern.length[[p]], pattern.length[[p]])
              for (i in 1:(pattern.count[[p]]/pattern.length[[p]])) {
                X.t <- t(t.J.X.p[[p]][,(1 + (i - 1) * pattern.length[[p]]):(i * pattern.length[[p]]) , drop = FALSE])
                Y.t <- J.Y.p[[p]][(1 + (i - 1) * pattern.length[[p]]):(i * pattern.length[[p]])]
                Z.t <- t(t.J.Z.p[[p]][,(1 + (i - 1) * pattern.length[[p]]):(i * pattern.length[[p]]) , drop = FALSE])
                rownumber.t <- J.Y.p.rownumber[[p]][(1 + (i - 1) * pattern.length[[p]]):(i * pattern.length[[p]])]
                temp.t <- Y.t - X.t %*% ybetas
                pattern.sum[[p]] <- pattern.sum[[p]] + inv.sqrt.W[rownumber.t, rownumber.t] %*% (tcrossprod(temp.t) - tcrossprod(temp.t, Z.t %*% eta.hat) - tcrossprod(Z.t %*% eta.hat, temp.t) + Z.t %*%
                                                                                                   tcrossprod(temp_mat_R, Z.t)) %*% inv.sqrt.W[rownumber.t, rownumber.t]
              }
            }
          }
          
          cat("R pattern sum:\n")
          print(pattern.sum)
          r.block.runtime= proc.time()[3] - rtime
          R.rmstep.times<-c(R.rmstep.times,r.block.runtime)
          cat("R code time for R-matrix calculations: ", r.block.runtime, " seconds\n")
          
        }
        ################7_16_17
        #pattern.sum stores the sums shown in Equation (16) of Karl, Yang
        #Lohr (2013). If no observations are missing, then there is only one
        #element in the pattern.sum list, a matrix of size equal to one of the
        #blocks of R. Otherwise, there is a different element (a matrix)
        #in the list for each missing data pattern.
        #This function needs to sum the blocks along the diagonal of a matrix
        #of size equal to that of R. Since J.X and J.Z are stored in sparse
        #(triplet) format in R, we need to extract the columns of those triplets
        #and send them to C++, where they will again be built into sparse 
        #matrices. That is the function of the @i, @p, and @x. See the documentation
        #for the R package Matrix for details.
        pattern.sum <- list()
        
        cpptime<-proc.time()[3]
        for (p in unique(J.mat$pat)) {
          if(control$REML){
            pattern.sum[[p]] <- REML_Rm(invsqrtW_ = as.matrix(diag(inv.sqrt.W)), betacov_ = betacov, C12_ = C12,
                                       JYp_ = as.matrix(J.Y.p[[p]]), loopsize_ = pattern.count[[p]]/pattern.length[[p]], 
                                       patternlength_ = pattern.length[[p]], rownumber_ = as.matrix(J.Y.p.rownumber[[p]]), 
                                       ybetas_ = as.matrix(ybetas), etahat_ = as.matrix(eta.hat), tempmatR_ = as.matrix(temp_mat), 
                                       JXpi_ = as.matrix(J.X.p[[p]]@i), JXpp_ = as.matrix(J.X.p[[p]]@p), 
                                       JXpx_ = as.matrix(J.X.p[[p]]@x), JXpdim_ = as.matrix(J.X.p[[p]]@Dim), 
                                       JZpi_ = as.matrix(J.Z.p[[p]]@i), JZpp_ = as.matrix(J.Z.p[[p]]@p), 
                                       JZpx_ = as.matrix(J.Z.p[[p]]@x), JZpdim_ = as.matrix(J.Z.p[[p]]@Dim))
            
          }else{
            
            pattern.sum[[p]] <- R_mstep2(invsqrtW_ = as.matrix(diag(inv.sqrt.W)), 
                                         JYp_ = as.matrix(J.Y.p[[p]]), loopsize_ = pattern.count[[p]]/pattern.length[[p]], 
                                         patternlength_ = pattern.length[[p]], rownumber_ = as.matrix(J.Y.p.rownumber[[p]]), 
                                         ybetas_ = as.matrix(ybetas), etahat_ = as.matrix(eta.hat), tempmatR_ = as.matrix(temp_mat_R), 
                                         JXpi_ = as.matrix(J.X.p[[p]]@i), JXpp_ = as.matrix(J.X.p[[p]]@p), 
                                         JXpx_ = as.matrix(J.X.p[[p]]@x), JXpdim_ = as.matrix(J.X.p[[p]]@Dim), 
                                         JZpi_ = as.matrix(J.Z.p[[p]]@i), JZpp_ = as.matrix(J.Z.p[[p]]@p), 
                                         JZpx_ = as.matrix(J.Z.p[[p]]@x), JZpdim_ = as.matrix(J.Z.p[[p]]@Dim))
            
          }
          
        }
        if(control$cpp.benchmark){
          cat("C++ pattern sum:\n")
          print(pattern.sum)
          cpp.block.runtime= proc.time()[3] - cpptime
          cat("C++ code time for R-matrix calculations: ", cpp.block.runtime, " seconds\n")
          
          
          cpp.rmstep.times<-c(cpp.rmstep.times,cpp.block.runtime)
          
        }
        #extract the unique parameters in the (full) block of R, R_i.
        #M step update for R
        #a NR algorithm (with gradient ascent at the beginning for 
        #stability) is used to solve the score equations from
        #Equation (16) of Karl, Yang, Lohr (2013)
        R_i.parm <- ltriangle(as.matrix(R_i))
        LRI <- length(R_i.parm)
        R_i.parm.constrained <- R_i.parm[-LRI]
        R.cc <- 1
        hes.count <- 1
        s.prev <- numeric(length(R_i.parm))
        while (R.cc > 1e-04) {
          s <- pattern.f.score.constrained(R_i.parm.constrained = R_i.parm.constrained, 
                                           nyear = nyear.pseudo, pattern.parmlist2 = pattern.parmlist2, 
                                           pattern.count = pattern.count, pattern.length = pattern.length, 
                                           pattern.Rtemplate = pattern.Rtemplate, pattern.diag = pattern.diag, 
                                           pattern.key = pattern.key, pattern.sum = pattern.sum)
          j <- symmpart(jacobian(pattern.f.score.constrained, c(R_i.parm.constrained), 
                                 method = "simple", nyear = nyear.pseudo, pattern.parmlist2 = pattern.parmlist2, 
                                 pattern.count = pattern.count, pattern.length = pattern.length, 
                                 pattern.Rtemplate = pattern.Rtemplate, pattern.diag = pattern.diag, 
                                 pattern.key = pattern.key, pattern.sum = pattern.sum))
          if(control$independent.responses){
            indep.key=ltriangle(1:LRI)
            coord.to.delete=indep.key[nrow(indep.key),-ncol(indep.key)]
            for(CTD.ind in 1:length(coord.to.delete)){
              j[coord.to.delete[CTD.ind],]<-j[,coord.to.delete[CTD.ind]]<-numeric(LRI-1)
              j[coord.to.delete[CTD.ind],coord.to.delete[CTD.ind]]<-1
              s[coord.to.delete[CTD.ind]]<-0
            }
          }
          if (it == 1 & (hes.count < 10) & PQL.it <= 3) {
            hesprod <- solve(j + max(c(diag(j), 5)) * diag(LRI - 1), s)
            R.cc <- s %*% s
            if (hes.count == 9) 
              R.cc <- 0
          } else if ((it <= 4) & (hes.count < 30) & PQL.it == 1) {
            hesprod <- solve(j + max(c(diag(j), 5) * ((1 - hes.count/31)^2)) * 
                               diag(LRI - 1), s)
            R.cc <- s %*% s
            if (hes.count == 29) 
              R.cc <- 0
          } else {
            hesprod <- solve(j, s)
            R.cc <- s %*% s
          }
          R_i.parm.constrained <- R_i.parm.constrained - hesprod
          hes.count <- hes.count + 1
          s.prev <- s
        }
        R_i.parm[-LRI] <- R_i.parm.constrained
        R_i <- ltriangle(R_i.parm)
        if(control$independent.responses){
          R_i[nrow(R_i),]<-R_i[,nrow(R_i)]<-0
          R_i[nrow(R_i),nrow(R_i)]<-1
        }
        rm(R_i.parm)
        dimnames(R_i) <- list(NULL, NULL)
        if (length(mis.list) > 0) {
          R <- symmpart(suppressMessages(kronecker(Diagonal(nstudent), R_i))[-mis.list, 
                                                                             -mis.list])
          rinv <- Matrix(0, 0, 0)
          R_i.inv <- chol2inv(chol(R_i))
          for (i in 1:nstudent) {
            if (!any(mis.list %in% seq((i - 1) * nyear.pseudo + 1, i * nyear.pseudo))) {
              rinv <- bdiag(rinv, R_i.inv)
            } else {
              inv.indx <- which(!(seq((i - 1) * nyear.pseudo + 1, i * nyear.pseudo) %in% 
                                    mis.list))
              rinv <- bdiag(rinv, chol2inv(chol(R_i[inv.indx, inv.indx])))
            }
          }
          R_inv <- rinv
        } else {
          R <- symmpart(suppressMessages(kronecker(Diagonal(nstudent), R_i)))
          R_inv <- symmpart(suppressMessages(kronecker(Diagonal(nstudent), 
                                                       chol2inv(chol(R_i)))))
        }
        R.full <- drop0(symmpart(sqrt.W %*% R %*% sqrt.W))
        R.full.inv <- drop0(symmpart(inv.sqrt.W %*% R_inv %*% inv.sqrt.W))
        names(ybetas) <- colnames(J.X)
        G <- Gn         
        #if the VP model is used, then Z need to be updated with the 
        #new persistence parameters
        #See the discussion in Section 3.4 of Karl, Yang, Lohr (2013)
        if (persistence == "VP") {
          alpha <- alphan
          if (!control$school.effects) {
            Z <- Matrix(0, nrow(Z_mat), nteach_effects)
          } else {
            Z <- Matrix(0, nrow(Z_mat), nteach_effects + nschool_effects)
          }
          for (i in 1:nalpha) {
            comp <- which(tril(ltriangle(1:nalpha)) == i, arr.ind = TRUE)
            Z <- Z + alpha[i] * P[[comp[1]]][[comp[2]]]
          }
          if (control$school.effects) {
            colnames(Z) <- eta_effects
            Z <- Z + Z.school.only
          }
          if (!control$school.effects) {
            Z.expand <- Matrix(0, nrow(Z_mat), 2 * nteach_effects)
            Z.expand[, seq(1, 2 * nteach_effects, by = 2)] <- Z
          } else {
            Z.expand <- Matrix(0, nrow(Z_mat), 2 * nteach_effects + 2 * nschool_effects)
            Z.expand[, seq(1, 2 * nteach_effects + 2 * nschool_effects, by = 2)] <- Z
          }
          colnames(Z) <- eta_effects
          # B.Z.expand <- Matrix(0, nrow(B.mat), 2 * nteach_effects) B.Z.expand[, seq(2, 2
          # * nteach_effects, by = 2)] <- B.Z
          J.Z <- rbind(Z.expand, B.Z.expand)
          J.Z <- J.Z[order(J.mat.original$student, J.mat.original$year), 
                     ]
          colnames(J.Z) <- interleave(colnames(Z), colnames(B.Z))
          for (p in unique(J.mat$pat)) {
            J.Z.p[[p]] <- J.Z[pat[[p]], , drop = FALSE]
          }
        }
        G.chol <- chol(G)
        G.inv <- chol2inv(G.chol)
        R.inv.Z <- R.full.inv %*% J.Z
        V.1 <- chol2inv(chol(G.inv + t(J.Z) %*% R.inv.Z))
        tX.Rinv.Z <- t(J.X) %*% R.inv.Z
        tX.Rinv.X <- t(J.X) %*% R.full.inv %*% J.X
        #M-step update for fixed effects vector (both those from the
        #score and the binary submodels)
        ybetas <- as.vector(chol2inv(chol(forceSymmetric(symmpart(tX.Rinv.X - 
                                                                    tX.Rinv.Z %*% V.1 %*% t(tX.Rinv.Z))))) %*% (t(J.X) %*% R.full.inv - 
                                                                                                                  tX.Rinv.Z %*% V.1 %*% t(R.inv.Z)) %*% J.Y)
        if (control$verbose) 
          cat("Iteration Time: ", proc.time()[3] - ptm, " seconds\n")
        flush.console()
      }
    }
    # em loop # End of EM
    nalpha.free <- sum(!((1:nalpha) %in% alpha.diag))
    #The code in this block is not used. If uncommented, it shows how
    #a NR approach may be used to estimate the model instead of an EM 
    #algorithm. It is not used because the NR algorithm is unstable in
    #the presence of strong correlations in G. See the discussion in
    #Section 3 of Karl, Yang, Lohr (2013)
    #        if (control$NR == TRUE & PQL.it > control$EM.to.NR) {
    #            # NR loop
    #            cat("Beginning NR algorithm\n")
    #            convh <- abs(loglikelihood)
    #            nr.counter <- 0
    #            while (as.numeric(convh/abs(loglikelihood)) > 1e-09) {
    #                nr.counter <- nr.counter + 1
    #                new.eta <- update.eta(X = J.X, Y = J.Y, Z = J.Z, R.full.inv = R.full.inv, 
    #                  ybetas = ybetas, G = G, cons.logLik = cons.logLik, Ny = Ny, nstudent = nstudent, 
    #                  n_eta = n_eta)
    #                eta.hat <- attr(new.eta, "eta")
    #                C12 <- attr(new.eta, "C12")
    #                var.eta.hat <- new.eta
    #                temp_mat <- var.eta.hat + tcrossprod(eta.hat, eta.hat)
    #                temp_mat_R <- attr(new.eta, "h.inv") + tcrossprod(eta.hat, eta.hat)
    #                cat("\n", attr(new.eta, "likelihood") - loglikelihood, "\n")
    #                loglikelihood <- attr(new.eta, "likelihood")
    #                # if (nr.counter == 1 | persistence == 'VP') {
    #                if (TRUE) {
    #                  pattern.sum <- list()
    #                  for (p in unique(J.mat$pat)) {
    #                    pattern.sum[[p]] <- Matrix(R_mstep2(invsqrtW_ = as.matrix(diag(inv.sqrt.W)), 
    #                      JYp_ = as.matrix(J.Y.p[[p]]), loopsize_ = pattern.count[[p]]/pattern.length[[p]], 
    #                      patternlength_ = pattern.length[[p]], rownumber_ = as.matrix(J.Y.p.rownumber[[p]]), 
    #                      ybetas_ = as.matrix(ybetas), etahat_ = as.matrix(eta.hat), 
    #                      tempmatR_ = as.matrix(temp_mat_R), JXpi_ = as.matrix(J.X.p[[p]]@i), 
    #                      JXpp_ = as.matrix(J.X.p[[p]]@p), JXpx_ = as.matrix(J.X.p[[p]]@x), 
    #                      JXpdim_ = as.matrix(J.X.p[[p]]@Dim), JZpi_ = as.matrix(J.Z.p[[p]]@i), 
    #                      JZpp_ = as.matrix(J.Z.p[[p]]@p), JZpx_ = as.matrix(J.Z.p[[p]]@x), 
    #                      JZpdim_ = as.matrix(J.Z.p[[p]]@Dim)))
    #                  }
    #                  
    #                }
    #                names(ybetas) <- colnames(J.X)
    #                if (persistence == "CP" | persistence == "ZP") {
    #                  thetas <- c(ltriangle(as.matrix(R_i))[-LRI], reduce.G(G = G, nyear.score = nyear.score, 
    #                    nteacher = nteacher))
    #                } else if (persistence == "VP") {
    #                  thetas <- c(ltriangle(as.matrix(R_i))[-LRI], reduce.G(G = G, nyear.score = nyear.score, 
    #                    nteacher = nteacher), alpha[!((1:nalpha) %in% alpha.diag)])
    #                }
    #                thetas.old <- thetas
    #                ybetas.old <- ybetas
    #                Hessian.inv <- hessian.f()
    #                score.thetas <- Score(thetas = thetas, eta = eta.hat, ybetas = ybetas, 
    #                  X = J.X, Y = J.Y, Z = J.Z, sqrt.W = sqrt.W, inv.sqrt.W = inv.sqrt.W, 
    #                  pattern.sum = pattern.sum, con = control, n_ybeta = n_ybeta, n_eta = n_eta, 
    #                  nstudent = nstudent, nteacher = nteacher, Kg = Kg, cons.logLik = cons.logLik, 
    #                  mis.list = mis.list, pattern.parmlist2 = pattern.parmlist2, pattern.count = pattern.count, 
    #                  pattern.length = pattern.length, pattern.Rtemplate = pattern.Rtemplate, 
    #                  pattern.diag = pattern.diag, pattern.key = pattern.key, Ny = Ny, 
    #                  persistence = persistence, P = P, alpha.diag = alpha.diag, nalpha = nalpha, 
    #                  alpha = alpha)
    #                thetas <- thetas - Hessian.inv %*% score.thetas
    #                convh <- t(score.thetas) %*% Hessian.inv %*% score.thetas
    #                # score.ybetas<-t(J.X)%*%R.full.inv%*%(J.Y-J.X%*%ybetas-J.Z%*%eta.hat)
    #                n_Rparm <- nyear.pseudo * (nyear.pseudo + 1)/2 - 1
    #                G <- thetas[(n_Rparm + 1):(n_Rparm + ncol(G.mat))]
    #                G <- reduce.G(G = G, nyear.score = nyear.score, nteacher = nteacher)
    #                if (persistence == "VP") {
    #                  alpha.parm <- thetas[(n_Rparm + ncol(G.mat) + 1):length(thetas)]
    #                  alpha[!((1:nalpha) %in% alpha.diag)] <- alpha.parm
    #                  if (!control$school.effects) {
    #                    Z <- Matrix(0, nrow(Z_mat), nteach_effects)
    #                  } else {
    #                    Z <- Matrix(0, nrow(Z_mat), nteach_effects + nschool_effects)
    #                  }
    #                  for (i in 1:nalpha) {
    #                    comp <- which(tril(ltriangle(1:nalpha)) == i, arr.ind = TRUE)
    #                    Z <- Z + alpha[i] * P[[comp[1]]][[comp[2]]]
    #                  }
    #                  if (control$school.effects) {
    #                    colnames(Z) <- eta_effects
    #                    Z <- Z + Z.school.only
    #                  }
    #                  if (!control$school.effects) {
    #                    Z.expand <- Matrix(0, nrow(Z_mat), 2 * nteach_effects)
    #                    Z.expand[, seq(1, 2 * nteach_effects, by = 2)] <- Z
    #                  } else {
    #                    Z.expand <- Matrix(0, nrow(Z_mat), 2 * nteach_effects + 2 * nschool_effects)
    #                    Z.expand[, seq(1, 2 * nteach_effects + 2 * nschool_effects, by = 2)] <- Z
    #                  }
    #                  colnames(Z) <- eta_effects
    #                  J.Z <- rbind(Z.expand, B.Z.expand)
    #                  J.Z <- J.Z[order(J.mat.original$student, J.mat.original$year), 
    #                    ]
    #                  colnames(J.Z) <- interleave(colnames(Z), colnames(B.Z))
    #                  for (p in unique(J.mat$pat)) {
    #                    J.Z.p[[p]] <- J.Z[pat[[p]], , drop = FALSE]
    #                  }
    #                  # if(!huge.flag){Z.dense <- as.matrix(Z)} for (p in unique(Z_mat$pat)) {
    #                  # if(!huge.flag){ Z.p[[p]] <- Z.dense[pat[[p]], , drop = FALSE]}else{ Z.p[[p]] <-
    #                  # Z[pat[[p]], , drop = FALSE] } } rm(Z.dense)
    #                }
    #                R_i <- ltriangle(c(as.vector(thetas[1:n_Rparm]), 1))
    #                R_i.parm <- c(as.vector(thetas[1:n_Rparm]), 1)
    #                LRI <- length(R_i.parm)
    #                R_i.parm.constrained <- R_i.parm[-LRI]
    #                if (length(mis.list) > 0) {
    #                  R <- symmpart(suppressMessages(kronecker(suppressMessages(Diagonal(nstudent)), 
    #                    R_i)[-mis.list, -mis.list]))
    #                  rinv <- Matrix(0, 0, 0)
    #                  R_i.inv <- chol2inv(chol(R_i))
    #                  for (i in 1:nstudent) {
    #                    if (!any(mis.list %in% seq((i - 1) * nyear.pseudo + 1, i * nyear.pseudo))) {
    #                      rinv <- bdiag(rinv, R_i.inv)
    #                    } else {
    #                      inv.indx <- which(!(seq((i - 1) * nyear.pseudo + 1, i * nyear.pseudo) %in% 
    #                        mis.list))
    #                      rinv <- bdiag(rinv, chol2inv(chol(R_i[inv.indx, inv.indx])))
    #                    }
    #                  }
    #                  R_inv <- rinv
    #                } else {
    #                  R <- symmpart(suppressMessages(kronecker(suppressMessages(Diagonal(nstudent)), 
    #                    R_i)))
    #                  R_inv <- symmpart(suppressMessages(kronecker(suppressMessages(Diagonal(nstudent)), 
    #                    chol2inv(chol(R_i)))))
    #                }
    #                R.full <- drop0(symmpart(sqrt.W %*% R %*% sqrt.W))
    #                R.full.inv <- drop0(symmpart(inv.sqrt.W %*% R_inv %*% inv.sqrt.W))
    #                G.chol <- chol(G)
    #                G.inv <- chol2inv(G.chol)
    #                R.inv.Z <- R.full.inv %*% J.Z
    #                V.1 <- chol2inv(chol(G.inv + t(J.Z) %*% R.inv.Z))
    #                tX.Rinv.Z <- t(J.X) %*% R.inv.Z
    #                tX.Rinv.X <- t(J.X) %*% R.full.inv %*% J.X
    #                ybetas <- as.vector(chol2inv(chol(symmpart(tX.Rinv.X - tX.Rinv.Z %*% 
    #                  V.1 %*% t(tX.Rinv.Z)))) %*% (t(J.X) %*% R.full.inv - tX.Rinv.Z %*% 
    #                  V.1 %*% t(R.inv.Z)) %*% J.Y)
    #            }
    #            new.eta <- update.eta(X = J.X, Y = J.Y, Z = J.Z, R.full.inv = R.full.inv, 
    #                ybetas = ybetas, G = G, cons.logLik = cons.logLik, Ny = Ny, nstudent = nstudent, 
    #                n_eta = n_eta)
    #            eta.hat <- attr(new.eta, "eta")
    #              C12 <- attr(new.eta, "C12")
    #            var.eta.hat <- new.eta
    #            temp_mat <- var.eta.hat + tcrossprod(eta.hat, eta.hat)
    #            temp_mat_R <- attr(new.eta, "h.inv") + tcrossprod(eta.hat, eta.hat)
    #            cat("\n", attr(new.eta, "likelihood") - loglikelihood, "\n")
    #        }
    # end NR
    #extract the unique model parameters and  compare to those from the
    #previous pseudolikihood iteration
    #The convergence criterion is the same as one of the available
    #options in SAS PROC GLIMMIX
    if (persistence == "CP" | persistence == "ZP") {
      thetas.pql <- c(ltriangle(as.matrix(R_i))[-LRI], reduce.G(G = G, nyear.score = nyear.score, 
                                                                nteacher = nteacher))
    } else if (persistence == "VP") {
      thetas.pql <- c(ltriangle(as.matrix(R_i))[-LRI], reduce.G(G = G, nyear.score = nyear.score, 
                                                                nteacher = nteacher), alpha[!((1:nalpha) %in% alpha.diag)])
    }
    parms.new <- c(ybetas, thetas.pql)
    if (PQL.it > 1) {
      max.change <- 2 * max(abs(parms.new - parms.old)/(abs(parms.new) + abs(parms.old) + 
                                                          1e-06))
      if (control$verbose) {
        cat("\nMaximum relative parameter change between PQL iterations: ", 
            max.change, "\n")
      }
      if (max.change < control$pconv) 
        break
    }
    #Update the pseudo-responses after each pseudo-likelihood iteration
    #This is direct out of Wolfinger and O'Connell (1993).
    #Note: even though the variables are called "PQL", this is actually
    #a pseudo-likelihood approach. They are identical under some circumstances
    thetas.pql.old <- thetas.pql
    ybetas.pql.old <- ybetas
    parms.old <- c(ybetas.pql.old, thetas.pql.old)
    # we use inverse of Wolfinger's W this part of the code transforms the responses
    # (pseudo-likelihood approach) most of this code was taken from the R function
    # glmmPQL
    J.mat[J.mat$response == "outcome", ]$fit <- as.vector(J.X[J.mat$response == 
                                                                "outcome", ] %*% ybetas + J.Z[J.mat$response == "outcome", ] %*% eta.hat)
    J.mat[J.mat$response == "outcome", ]$mu <- as.vector(fam.binom$linkinv(J.mat[J.mat$response == 
                                                                                   "outcome", ]$fit))
    J.mat[J.mat$response == "outcome", ]$mu.eta.val <- fam.binom$mu.eta(J.mat[J.mat$response == 
                                                                                "outcome", ]$fit)
    J.mat[J.mat$response == "outcome", ]$nu <- as.vector(J.mat[J.mat$response == 
                                                                 "outcome", ]$fit) + (J.mat[J.mat$response == "outcome", ]$y - J.mat[J.mat$response == 
                                                                                                                                       "outcome", ]$mu)/J.mat[J.mat$response == "outcome", ]$mu.eta.val
    J.mat[J.mat$response == "outcome", ]$sqrt.w <- sqrt(fam.binom$variance(J.mat[J.mat$response == 
                                                                                   "outcome", ]$mu)/J.mat[J.mat$response == "outcome", ]$mu.eta.val^2)
    J.mat[J.mat$response == "score", ]$sqrt.w <- 1
    J.mat[J.mat$response == "score", ]$nu <- J.mat[J.mat$response == "score", 
                                                   ]$y
    sqrt.W <- Diagonal(x = J.mat$sqrt.w)
    inv.sqrt.W <- Diagonal(x = 1/(J.mat$sqrt.w))
    R.full <- drop0(symmpart(sqrt.W %*% R %*% sqrt.W))
    R.full.inv <- drop0(symmpart(inv.sqrt.W %*% R_inv %*% inv.sqrt.W))
    J.Y <- J.mat$nu
    for (p in unique(J.mat$pat)) {
      J.Y.p[[p]] <- J.Y[pat[[p]]]
    }
  }
  # pql loop
  #end of the outer pseduo-likelihood loop
  #parameter estimates do not changes past this point
  names(ybetas) <- colnames(J.X)
  if (persistence == "CP" | persistence == "ZP") {
    thetas <- c(ltriangle(as.matrix(R_i))[-LRI], reduce.G(G = G, nyear.score = nyear.score, 
                                                          nteacher = nteacher))
  } else if (persistence == "VP") {
    thetas <- c(ltriangle(as.matrix(R_i))[-LRI], reduce.G(G = G, nyear.score = nyear.score, 
                                                          nteacher = nteacher), alpha[!((1:nalpha) %in% alpha.diag)])
  }
  lgLik.hist <- lgLik
  lgLik <- lgLik[it]
  Hessian <- NA
  std_errors <- c(rep(NA, length(thetas)))
  #hessian of the variance components (R and G) calculated with a central
  #difference approximation using the Score() function
  if (control$hessian == TRUE) {
    if (control$verbose) 
      cat("Calculating Hessian of the variance components...")
    flush.console()
    pattern.sum <- list()
    for (p in unique(J.mat$pat)) {
      if(control$REML){
        pattern.sum[[p]] <- Matrix(REML_Rm(invsqrtW_ = as.matrix(diag(inv.sqrt.W)), betacov_ = betacov, C12_ = C12,
                                          JYp_ = as.matrix(J.Y.p[[p]]), loopsize_ = pattern.count[[p]]/pattern.length[[p]], 
                                          patternlength_ = pattern.length[[p]], rownumber_ = as.matrix(J.Y.p.rownumber[[p]]), 
                                          ybetas_ = as.matrix(ybetas), etahat_ = as.matrix(eta.hat), tempmatR_ = as.matrix(temp_mat), 
                                          JXpi_ = as.matrix(J.X.p[[p]]@i), JXpp_ = as.matrix(J.X.p[[p]]@p), 
                                          JXpx_ = as.matrix(J.X.p[[p]]@x), JXpdim_ = as.matrix(J.X.p[[p]]@Dim), 
                                          JZpi_ = as.matrix(J.Z.p[[p]]@i), JZpp_ = as.matrix(J.Z.p[[p]]@p), 
                                          JZpx_ = as.matrix(J.Z.p[[p]]@x), JZpdim_ = as.matrix(J.Z.p[[p]]@Dim)))
        
      }else{
        
        pattern.sum[[p]] <- Matrix(R_mstep2(invsqrtW_ = as.matrix(diag(inv.sqrt.W)), 
                                            JYp_ = as.matrix(J.Y.p[[p]]), loopsize_ = pattern.count[[p]]/pattern.length[[p]], 
                                            patternlength_ = pattern.length[[p]], rownumber_ = as.matrix(J.Y.p.rownumber[[p]]), 
                                            ybetas_ = as.matrix(ybetas), etahat_ = as.matrix(eta.hat), tempmatR_ = as.matrix(temp_mat_R), 
                                            JXpi_ = as.matrix(J.X.p[[p]]@i), JXpp_ = as.matrix(J.X.p[[p]]@p), 
                                            JXpx_ = as.matrix(J.X.p[[p]]@x), JXpdim_ = as.matrix(J.X.p[[p]]@Dim), 
                                            JZpi_ = as.matrix(J.Z.p[[p]]@i), JZpp_ = as.matrix(J.Z.p[[p]]@p), 
                                            JZpx_ = as.matrix(J.Z.p[[p]]@x), JZpdim_ = as.matrix(J.Z.p[[p]]@Dim)))
        
      }
    }
    
    if (control$hes.method == "richardson") {
      Hessian <- ltriangle(ltriangle(jacobian(Score, thetas, eta = eta.hat, 
                                              ybetas = ybetas, X = J.X, Y = J.Y, Z = J.Z, sqrt.W = sqrt.W, inv.sqrt.W = inv.sqrt.W, 
                                              pattern.sum = pattern.sum, con = control, n_ybeta = n_ybeta, n_eta = n_eta, 
                                              nstudent = nstudent, nteacher = nteacher, Kg = Kg, cons.logLik = cons.logLik, 
                                              mis.list = mis.list, pattern.parmlist2 = pattern.parmlist2, pattern.count = pattern.count, 
                                              pattern.length = pattern.length, pattern.Rtemplate = pattern.Rtemplate, 
                                              pattern.diag = pattern.diag, pattern.key = pattern.key, Ny = Ny, 
                                              persistence = persistence, P = P, alpha.diag = alpha.diag, nalpha = nalpha, 
                                              alpha = alpha)))
    } else {
      Hessian <- ltriangle(ltriangle(jacobian(Score, thetas, method = "simple", 
                                              eta = eta.hat, ybetas = ybetas, X = J.X, Y = J.Y, Z = J.Z, sqrt.W = sqrt.W, 
                                              inv.sqrt.W = inv.sqrt.W, pattern.sum = pattern.sum, con = control, 
                                              n_ybeta = n_ybeta, n_eta = n_eta, nstudent = nstudent, nteacher = nteacher, 
                                              Kg = Kg, cons.logLik = cons.logLik, mis.list = mis.list, pattern.parmlist2 = pattern.parmlist2, 
                                              pattern.count = pattern.count, pattern.length = pattern.length, pattern.Rtemplate = pattern.Rtemplate, 
                                              pattern.diag = pattern.diag, pattern.key = pattern.key, Ny = Ny, 
                                              persistence = persistence, P = P, alpha.diag = alpha.diag, nalpha = nalpha, 
                                              alpha = alpha)))
    }
    std_errors <- try(c(sqrt(diag(solve(Hessian)))), silent = TRUE)
    hes.warn <- FALSE
    if (any(eigen(Hessian)$values <= 0)) {
      if (control$verbose) 
        cat("Warning: Hessian not PD", "\n")
      std_errors <- c(rep(NA, length(thetas)))
      hes.warn <- TRUE
    }
  }
  #Henderson mixed model equations produce standard errors for
  #fixed effects and EBLUPs
  c.temp <- crossprod(J.X, R.full.inv) %*% J.Z
  c.1 <- rbind(crossprod(J.X, R.full.inv) %*% J.X, t(c.temp))
  G.inv <- chol2inv(chol(G))
  c.2 <- rbind(c.temp, H.eta(G.inv, J.Z, R.full.inv))
  C_inv <- cbind(c.1, c.2)
  C <- solve(C_inv)
  eblup_stderror <- sqrt(diag(C)[-c(1:n_ybeta)])
  ybetas_stderror <- sqrt(diag(C)[1:n_ybeta])
  C<-as.matrix(C)
  rownames(C)<-colnames(C)<-c(names(ybetas),colnames(J.Z))
  rm( C_inv, c.2, c.1, c.temp)
  eblup <- as.matrix(cbind(eta.hat, eblup_stderror))
  eblup <- as.data.frame(eblup)
  eblup <- as.data.frame(cbind(colnames(J.Z), eblup))
  #remainder of the code is just formatting to return to user
  colnames(eblup) <- c("effect", "EBLUP", "std_error")
  t_lab <- as.vector(NULL)
  r_lab <- as.vector(NULL)
  for (j in 1:nyear.score) {
    # ne <- (Kg[j] * (Kg[j] + 1))/2 y <- c(NULL) x <- c(NULL) for (k in 1:Kg[j]) { x
    # <- c(x, k:Kg[j]) y <- c(y, rep(k, (Kg[j] - k + 1))) }
    t_lab <- c(t_lab, paste("teacher effect from year", j, " [1,1] (score variance)", 
                            sep = ""))
    t_lab <- c(t_lab, paste("teacher effect from year", j, " [2,1] (score-outcome covariance)", 
                            sep = ""))
    t_lab <- c(t_lab, paste("teacher effect from year", j, " [2,2] (outcome variance)", 
                            sep = ""))
  }
  if (control$school.effects) {
    t_lab <- c(t_lab, rep("School Effect", 3))
  }
  y <- c(NULL)
  x <- c(NULL)
  for (k in 1:nyear.pseudo) {
    x <- c(x, k:nyear.pseudo)
    y <- c(y, rep(k, (nyear.pseudo - k + 1)))
  }
  r_lab <- paste("error covariance", ":[", x, ",", y, "]", sep = "")
  # rm(j, ne)
  alpha.label <- c(NULL)
  for (i in 1:(nyear.score - 1)) {
    for (j in (i + 1):(nyear.score)) {
      alpha.label <- c(alpha.label, paste("alpha_", j, i, sep = ""))
    }
  }
  if (persistence == "CP" | persistence == "ZP") {
    effect_la <- c(names(ybetas), r_lab, t_lab)
  } else if (persistence == "VP") {
    effect_la <- c(names(ybetas), r_lab, t_lab, alpha.label)
  }
  constant.element <- which(effect_la == paste("error covariance", ":[", nyear.pseudo, 
                                               ",", nyear.pseudo, "]", sep = ""))
  if (control$hessian == TRUE) {
    parameters <- round(cbind(append(c(ybetas, thetas), 1, after = constant.element - 
                                       1), append(c(ybetas_stderror, std_errors), NA, after = constant.element - 
                                                    1)), 4)
    colnames(parameters) <- c("Estimate", "Standard Error")
    rownames(parameters) <- as.character(effect_la)
  }
  if (control$hessian == FALSE) {
    parameters <- round(cbind(append(c(ybetas, thetas), 1, after = constant.element - 
                                       1), c(ybetas_stderror, rep(NA, length(thetas) + 1))), 4)
    colnames(parameters) <- c("Estimate", "Standard Error")
    rownames(parameters) <- as.character(effect_la)
  }
  if (control$verbose) 
    cat("done.\n")
  mresid <- as.numeric(J.Y - J.X %*% ybetas)
  cresid <- as.numeric(mresid - J.Z %*% eta.hat)
  yhat <- as.numeric(J.X %*% ybetas + J.Z %*% eta.hat)
  yhat.m <- as.numeric(J.X %*% ybetas)
  # rchol <- try(chol(R.full.inv)) yhat.s <- try(as.vector(rchol %*% (yhat)))
  # sresid <- try(as.vector(rchol %*% J.Y - yhat.s)) print(system.time(
  # Vchol<-chol(J.Z%*%G%*%t(J.Z)+R.full))) print(system.time(
  # sresid<-solve(Vchol,J.Y-res$yhat) ))
  gam_t <- list()
  for (i in 1:nyear.score) {
    gam_t[[i]] <- as.matrix(ltriangle(reduce.G(G, nyear.score = nyear.score, 
                                               nteacher = nteacher)[(1 + 3 * (i - 1)):(3 * i)]))
    colnames(gam_t[[i]]) <- rownames(gam_t[[i]]) <- c("score", "outcome")
  }
  persistence_parameters = ltriangle(alpha)
  persistence_parameters[upper.tri(persistence_parameters)] <- NA
  if (!control$school.effects) {
    school.subset <- NULL
  } else {
    school.subset <- eblup[(2 * nteach_effects + 1):n_eta, ]
  }
  teach.cov <- lapply(gam_t, function(x) round(x, 4))
  names(J.mat)[names(J.mat) == "y"] <- "y.combined.original"
  R_i <- as.matrix(R_i)
  rownames(R_i) <- colnames(R_i) <- c(paste("year", 1:(ncol(R_i) - 1), "score", 
                                            sep = "_"), "outcome")
  parameters <- cbind(parameters, z = round(parameters[, 1]/parameters[, 2], 2))
  parameters[, 3][grepl("variance", rownames(parameters)) & !grepl("covariance", 
                                                                   rownames(parameters))] <- NA
  parameters <- cbind(parameters, pvalue = round(2 * (1 - pnorm(abs(parameters[, 
                                                                               3]))), 3))
  parameters <- cbind(parameters, lower.95.CI = round(parameters[, 1] - parameters[, 
                                                                                   2] * qnorm(0.975), 3))
  parameters <- cbind(parameters, upper.95.CI = round(parameters[, 1] + parameters[, 
                                                                                   2] * qnorm(0.975), 3))
  parameters[, 5:6][grepl("variance", rownames(parameters)) & !grepl("covariance", 
                                                                     rownames(parameters))] <- NA
  res <- list(loglik = lgLik, teach.effects = eblup[1:(2 * nteach_effects), ], 
              school.effects = school.subset, parameters = parameters, Hessian = Hessian, 
              R_i = R_i, teach.cov = gam_t, mresid = mresid, cresid = cresid, y.combined = J.Y, 
              y.combined.hat = yhat, y.response.type = J.mat$response, y.year = J.mat$year, 
              num.obs = Ny, num.student = nstudent, num.year = nyear.score, num.teach = nteacher, 
              persistence = control$persistence, persistence_parameters = persistence_parameters, 
              X = J.X, Z = J.Z, G = G, R = R,C=C, R.full = R.full, sqrt.W = J.mat$sqrt.w, eblup = eblup, 
              fixed.effects = ybetas, joined.table = J.mat)
  try(res$teach.effects$teacher_year <- as.numeric(substr(gsub(".*year", "", res$teach.effects[, 
                                                                                               1]), 1, 1)))
  # res$teach.effects$teacher_year <- key[match(res$teach.effects$teacher_year,
  # key[, 2]), 1] res$teach.effects$effect_year <-
  # key[match(res$teach.effects$effect_year, key[, 2]), 1]
  class(res) <- "RealVAMS"
  cat("Total Time: ", proc.time()[3] - ptm.total, " seconds\n")
  if(control$cpp.benchmark){
    cat("Run times for identical R and C++ code during each iteration\n")
    df.bench=data.frame(Cpp.times=cpp.rmstep.times,R.times=R.rmstep.times)
    print(df.bench)
    print(summary(df.bench))
  }
  return(res)
}