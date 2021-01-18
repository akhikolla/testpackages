##################################################################*
##################################################################*
create_multiple_minimax <-function(model, FIM,  multipars, is.only.w, only_w_varlist = NULL){
  #sqrt(.Machine$double.xmin)
  ## model: a character denotes the model
  ## we should define this because for very small numbers near zero we have NaN, instead we return -Inf

  log2 <- function(x){
    out <- suppressWarnings(log(x))
    if (is.na(out))
      out <- -1e24
    return(out)
  }

  if (model == "4pl"){
    lambda <- multipars$lambda
    delta <- multipars$delta
    s <- 4
    crfunc <- function(param, q, npred){

      if (!is.only.w){
        lq <- length(q)
        n_seg <- lq/(npred + 1)
        x_ind <- 1:(npred * n_seg)
        w_ind <- (x_ind[length(x_ind)]+1):lq
        x <- q[x_ind]
        w <- q[w_ind]
      } else{
        w <- q
        x <- only_w_varlist$x
      }
      ################*
      ## compute the information matrix and determine whether it is exactly singular
      FIM_val <- FIM(x=x, w = w, param = param)
      logdet_FIM <- det2(FIM_val, logarithm=TRUE)
      ### we calculate the logarithm to be sure that the information matrix is not singular
      if (length(w) < length(param) || logdet_FIM == 1e+24)
        singular_FIM <- TRUE else
          singular_FIM <- FALSE
      ###############*

      ###############*
      ## calculating the first part in equation of page 262
      if (lambda[1] != 0){
        if (!singular_FIM)
          part1 <- lambda[1]/s * logdet_FIM else{
            return(1e+24)
          }
      }else
        part1 <- 0
      ###############*

      ###############*
      ### computing the generalized inverse of Fisher information matrix
      if (lambda[2] != 0 || lambda[3] != 0){
        #if (singular_FIM)
        FIM_inv <- mpginv(FIM_val, tol = multipars$tol) #else
        #sqrt(.Machine$double.xmin)
        #FIM_inv <- solve(FIM_val, tol =sqrt(.Machine$double.eps))
        #det( FIM_inv)
        #MASS::ginv(FIM_val, tol = sqrt(.Machine$double.eps))
      }
      ###############*

      ###############*
      ## part2 ED50
      if (lambda[2] == 0)
        part2 <- 0 else{
          ED50_prime <- matrix(c(0, param[3]/param[2]^2, -1/param[2], 0), 1, 4)
          var_ED50 <- ED50_prime %*% FIM_inv %*% t(ED50_prime)
          #see equation 4 of Hyun and Wong (2015)
          ## we dont have negative var!!
          if (var_ED50 > 0){
            # if (lambda[2] == 1)
            #   part2 <- var_ED50 else
            part2 <- lambda[2] * log(var_ED50)
          } else
            return(1e+24) # because we are finding the minimum

        }
      ###############*

      ###############*
      ## part3 MED50
      if (lambda[3] == 0)
        part3 <- 0 else {
          ## computing the MED prime
          ###############* MED prime and var_MED
          if (param[2] > 0)
            MED_prime <- matrix(c(-1/((param[1]+delta)*param[2]),
                                  (param[3] - log(-delta/(param[1]+delta)))/param[2]^2,
                                  -1/param[2], 0),1, 4)
          if (param[2] < 0)
            MED_prime <- matrix(c(1/((param[1]-delta)*param[2]),
                                  (param[3] - log((param[1]-delta)/delta))/param[2]^2,
                                  -1/param[2], 0), 1, 4)
          var_MED <- MED_prime %*% FIM_inv %*% t(MED_prime)
          ###########################################*
          ## Eq. 5 Hyun and Wong (2015)
          ## we dont have negative var!!

          if (var_MED > 0){
            # if (lambda[3] == 1)
            #   part3 <- var_MED else
            part3 <- lambda[3] * log(var_MED)
          } else
            return(1e+24) # because we are finding the minimum

        }
      ###############*
      # we find the minimum here
      locally_crfunc <- -part1 + part2 + part3
      return(locally_crfunc)
    }

    # only for 4PL!
    # equation 6 of Multiple-Objective Optimal Designs for Studying the Dose Response Function and Interesting Dose Levels
    # Hyun and Wong (2015)
    # directional derivative considerations and show that for a fixed lambda,
    # the sensitivity function for the locally multiple-objective
    d_multi_x <- function(x1, FIM, x, w, param, lambda, delta){
      npar <- 4 #  four parameter logistic
      ###############* find the generalized inverse
      FIM_val = FIM(x = x, w = w, param = param)
      #FIM_inv <- solve(FIM_val, tol = .Machine$double.xmin)
      # FIM_inv <- mpginv(FIM_val, tol =  sqrt(.Machine$double.xmin))
      FIM_inv <- mpginv(FIM_val, tol =  multipars$tol)
      #sqrt(.Machine$double.eps)
      ###############*

      ####################### writign g
      constant1 <- exp(param[2]*x1 + param[3])
      constant2 <- 1/(1+constant1)
      g_x <- matrix(c(constant2, - param[1] * x1 * constant1*constant2^2, -param[1] * constant1 * constant2^2, 1), 1, 4)
      ###############*

      ##ED50 prime!
      ED50_prime <- matrix(c(0, param[3]/param[2]^2, -1/param[2], 0), 1, 4)
      ##MED prime!
      if(param[2]>0)
        MED_prime <- matrix(c(-1/((param[1] + delta) * param[2]),
                              (param[3] - log(-delta/(param[1] + delta)))/param[2]^2,
                              -1/param[2], 0), 1, 4)
      if(param[2]<0)
        MED_prime <- matrix(c(1/((param[1]-delta) * param[2]),
                              (param[3] - log((param[1] - delta)/delta))/param[2]^2,
                              -1/param[2], 0), 1, 4)

      ## finally we can write d(x, xi)
      if (lambda[2] == 1){
        d_x_xi <-  g_x %*% FIM_inv %*% t(ED50_prime)  %*% ED50_prime %*% FIM_inv %*% t(g_x) - (ED50_prime %*% FIM_inv %*% t(ED50_prime))
      }else
        # d_x_xi <-  ((g_x %*% FIM_inv %*% t(ED50_prime))^2) - (ED50_prime %*% FIM_inv %*% t(ED50_prime)) else
        if (lambda[3] == 1){
          d_x_xi<-  g_x %*% FIM_inv %*% t(MED_prime)  %*% MED_prime %*% FIM_inv %*% t(g_x) - (MED_prime %*% FIM_inv %*% t(MED_prime))
          #d_x_xi <-  ((g_x %*% FIM_inv %*% t(MED_prime))^2) - (MED_prime %*% FIM_inv %*% t(MED_prime))
        }else{
          d_x_xi <- lambda[1]/npar * g_x %*% FIM_inv %*% t(g_x) +
            lambda[2] * ((g_x %*% FIM_inv %*% t(ED50_prime))^2)/(ED50_prime %*% FIM_inv %*% t(ED50_prime)) +
            lambda[3] * ((g_x %*% FIM_inv %*% t(MED_prime))^2)/(MED_prime %*% FIM_inv %*% t(MED_prime)) - 1
        }
      return(d_x_xi)
    }

    ## lambda and delta are from the global enironment
    PsiMulti_x <- function(x1, mu, FIM, x, w, answering){
      if (length(mu) != dim(answering)[1])
        stop("The number of measures is not equal to the number of elements of answering set.")

      if (typeof(FIM) != "closure")
        stop("'FIM' must be of type 'closure'.")


      CardOfRegion <- dim(answering)[2] ## cardinality of region of uncertainty
      n_mu <- dim(answering)[1]

      ## each row is  tr(M^{-1}(\xi, mu_j) %*% I(x, \mu_j))
      # so sum of each row of 'Psi_Point_answering' is \int tr(M^{-1}(\xi, mu) %*% I(x_i, \mu))
      Psi_Point_answering <- matrix(NA, 1,  n_mu)

      ##the sum of eah row minus number of parameter is psi
      i <- 1
      for(j in 1:n_mu){
        Psi_Point_answering[i, j]  <-
          mu[j] *
          d_multi_x(x1 = x1, FIM = FIM, x = x, w = w, param = answering[j, ], lambda = lambda, delta = delta)
      }

      PsiAtEachPoint <- rowSums(Psi_Point_answering) #- CardOfRegion
      #x can be any support points, p783 King and Wong(2004). So we choose the first one
      PsiFunction <- PsiAtEachPoint[1]
      return(PsiFunction)

    }

    ## lambda and delta are from the global enironmen
    PsiMulti_Mu <- function(mu, FIM,  x, w, answering, PenaltyCoeff){
      if(length(mu) != dim(answering)[1])
        stop("The number of measures is not equal to the number of elements of answering set.")
      if(typeof(FIM) != "closure")
        stop("'FIM' must be of type 'closure'.")


      CardOfRegion <- dim(answering)[2] ##cardinal of region of uncertainty
      n_mu <- dim(answering)[1] ## number of mu, measures


      npred <- length(x)/length(w) # number of independent variables.
      one_point_mat <- matrix(x, length(w), npred)

      ##The value of Psi at each design x1 = points[i] and each element of answering
      Psi_Point_answering <- matrix(NA, length(w),  n_mu)

      ##the sum of eah row minus number of parameter is psi
      for(i in 1:length(w)){
        for(j in 1:n_mu){
          Psi_Point_answering[i, j]  <-
            mu[j] *
            d_multi_x(x1 = as.vector(one_point_mat[i,]), FIM = FIM, x = x, w = w,
                      param = answering[j, ], lambda = lambda,  delta = delta)
        }
      }
      #Psi at each points x = poi
      #Psi at each point
      PsiAtEachPoint <- round(rowSums(Psi_Point_answering) - CardOfRegion, 5)
      ##the value of Psi as each design point of design m(x,nu)

      # x can be any support points, p783 King and Wong(2004). So we choose the first one
      PsiEqualityPenalty <- sum(PenaltyCoeff*pmax(PsiAtEachPoint, 0)^2 + PenaltyCoeff*pmax(-PsiAtEachPoint, 0)^2)
      PsiFunction <-  PsiEqualityPenalty + PenaltyCoeff *  (sum(mu) -1)^2

      return(PsiFunction)
    }
  }#4pl

  # we dont need psi_Mu for locally version but we kepp it
  return(list(crfunc = crfunc, PsiMulti_Mu = PsiMulti_Mu, PsiMulti_x = PsiMulti_x))
}
##################################################################*
##################################################################*
mpginv <- function (a, tol = sqrt(.Machine$double.eps)){
  #### Moore-Penrose Matrix Inverse
  if (length(dim(a)) > 2 || !(is.numeric(a) || is.complex(a)))
    stop("a must be a numeric or complex matrix")
  if (!is.matrix(a))
    a <- as.matrix(a)
  asvd <- La.svd(a)
  if (is.complex(a)) {
    asvd$u <- Conj(asvd$u)
    asvd$v <- t(Conj(asvd$vt))
  }
  else {
    asvd$v <- t(asvd$vt)
  }
  Positive <- asvd$d > max(tol * asvd$d[1], 0)
  if (!any(Positive))
    array(0, dim(a)[2:1])
  else asvd$v[, Positive] %*% ((1/asvd$d[Positive]) * t(asvd$u[ , Positive]))
}
##################################################################*
##################################################################*
# create_multiple_bayes <-function(model, fimfunc2, multiple_arg){
#   ## model: a character strings denotes the model. Currently, only "FIM_logistic_4par"
#   ## multiple_arg is a list of necessary parameters for the model
#
#
#
#   ## we should define the following because for very small numbers near zero we have NaN, instead we return -Inf
#   log2 <- function(x){
#     out <- suppressWarnings(log(x))
#     if (is.na(out))
#       out <- 0
#     return(out)
#   }
#
#   if (model == "FIM_logistic_4par"){
#     if (all(multiple_arg$lambda != 0)){
#       crfunc <- function(param, x, w, multiple){
#         #param <- c(1.563, param, .137)
#         s <- 4
#         lambda <- multiple$lambda
#         delta <- multiple$delta
#         ###############*
#         ## we compute the information matrix and determine whether it is exactly singular
#         FIM_val <- fimfunc2(x=x, w = w, param = param)
#         logdet_FIM <- det2(FIM_val, logarithm=TRUE)
#         part1 <- lambda[1]/s * logdet_FIM
#         FIM_inv <- mpginv(FIM_val, tol = sqrt(.Machine$double.xmin))
#         #part2
#         ED50_prime <- matrix(c(0, param[3]/param[2]^2, -1/param[2], 0), 1, 4)
#         var_ED50 <- ED50_prime %*% FIM_inv %*% t(ED50_prime)
#
#         if(var_ED50>0)
#           part2 <-  lambda[2] * log(var_ED50) else
#             part2 <- 1e+24
#         #part3
#         if (param[2] > 0)
#           MED_prime <- matrix(c(-1/((param[1]+delta)*param[2]),
#                                 (param[3] - log(-delta/(param[1]+delta)))/param[2]^2,
#                                 -1/param[2], 0),1, 4)
#         if (param[2] < 0)
#           MED_prime <- matrix(c(1/((param[1]-delta)*param[2]),
#                                 (param[3] - log((param[1]-delta)/delta))/param[2]^2,
#                                 -1/param[2], 0), 1, 4)
#         var_MED <- MED_prime %*% FIM_inv %*% t(MED_prime)
#         if(var_MED>0)
#           part3 <- lambda[3] * log(var_MED) else
#             part3 <- 1e+24
#         locally_crfunc <- -part1 + part2 + part3
#         return(locally_crfunc)
#       }
#       ### sensitivity function
#       d_multi_x <- function(x1, FIM, x, w, param, multiple){
#         npar <- 4 #  four parameter logistic
#         lambda <- multiple$lambda
#         delta <- multiple$delta
#         #################### find the generalized inverse
#         FIM_val = FIM(x = x, w = w, param = param)
#         FIM_inv <- mpginv(FIM_val, tol =  sqrt(.Machine$double.eps))
#         ##################################################*
#
#         ####################### writign g
#         constant1 <- exp(param[2]*x1 + param[3])
#         constant2 <- 1/(1+constant1)
#         g_x <- matrix(c(constant2, - param[1] * x1 * constant1*constant2^2, -param[1] * constant1 * constant2^2, 1), 1, 4)
#         ##################################################*
#         ##ED50 prime!
#         ED50_prime <- matrix(c(0, param[3]/param[2]^2, -1/param[2], 0), 1, 4)
#         ##MED prime!
#         if(param[2]>0)
#           MED_prime <- matrix(c(-1/((param[1] + delta) * param[2]),
#                                 (param[3] - log(-delta/(param[1] + delta)))/param[2]^2,
#                                 -1/param[2], 0), 1, 4)
#         if(param[2]<0)
#           MED_prime <- matrix(c(1/((param[1]-delta) * param[2]),
#                                 (param[3] - log((param[1] - delta)/delta))/param[2]^2,
#                                 -1/param[2], 0), 1, 4)
#         d_x_xi <- lambda[1]/npar * g_x %*% FIM_inv %*% t(g_x) +
#           lambda[2] * ((g_x %*% FIM_inv %*% t(ED50_prime))^2)/(ED50_prime %*% FIM_inv %*% t(ED50_prime)) +
#           lambda[3] * ((g_x %*% FIM_inv %*% t(MED_prime))^2)/(MED_prime %*% FIM_inv %*% t(MED_prime)) - 1
#         return(d_x_xi)
#       }
#     }
#     if (round(multiple_arg$lambda[1], 5) == 1){
#       crfunc <- function(param, x, w, multiple){
#         s <- 4
#         lambda <- multiple$lambda
#         delta <- multiple$delta
#         ##########################################################*
#         ## compute the information matrix and determine whether it is exactly singular
#         FIM_val <- fimfunc2(x=x, w = w, param = param)
#         logdet_FIM <- det2(FIM_val, logarithm=TRUE)
#         part1 <- lambda[1]/s * logdet_FIM
#         locally_crfunc <- -part1
#         return(locally_crfunc)
#       }
#       # only for D-optimal design
#       d_multi_x <- function(x1, FIM, x, w, param, multiple){
#
#         npar <- 4 #  four parameter logistic
#         lambda <- multiple$lambda
#         delta <- multiple$delta
#         #################### find the generalized inverse
#         FIM_val = FIM(x = x, w = w, param = param)
#         FIM_inv <- mpginv(FIM_val, tol =  sqrt(.Machine$double.eps))
#         ##################################################*
#         ####################### writign g
#         constant1 <- exp(param[2]*x1 + param[3])
#         constant2 <- 1/(1+constant1)
#         g_x <- matrix(c(constant2, - param[1] * x1 * constant1*constant2^2, -param[1] * constant1 * constant2^2, 1), 1, 4)
#         ##################################################*
#         d_x_xi <- lambda[1]/npar * g_x %*% FIM_inv %*% t(g_x)
#         return(d_x_xi)
#       }
#     }
#     if (round(multiple_arg$lambda[2], 5) == 1){
#       crfunc <- function(param, x, w, multiple){
#         FIM_val <- fimfunc2(x=x, w = w, param = param)
#         lambda <- multiple$lambda
#         delta <- multiple$delta
#         FIM_inv <- mpginv(FIM_val, tol = sqrt(.Machine$double.xmin))
#         #part2
#         ED50_prime <- matrix(c(0, param[3]/param[2]^2, -1/param[2], 0), 1, 4)
#         var_ED50 <- ED50_prime %*% FIM_inv %*% t(ED50_prime)
#         part2 <-  lambda[2] * log(var_ED50)
#         locally_crfunc <- part2
#         return(locally_crfunc)
#       }
#       d_multi_x <- function(x1, FIM, x, w, param, multiple){
#         npar <- 4 #  four parameter logistic
#         lambda <- multiple$lambda
#         delta <- multiple$delta
#         #################### find the generalized inverse
#         FIM_val = FIM(x = x, w = w, param = param)
#         #FIM_inv <- solve(FIM_val, tol = .Machine$double.xmin)
#         # FIM_inv <- mpginv(FIM_val, tol =  sqrt(.Machine$double.xmin))
#         FIM_inv <- mpginv(FIM_val, tol =  sqrt(.Machine$double.eps))
#         ##################################################*
#         ####################### writign g
#         constant1 <- exp(param[2]*x1 + param[3])
#         constant2 <- 1/(1+constant1)
#         g_x <- matrix(c(constant2, - param[1] * x1 * constant1*constant2^2, -param[1] * constant1 * constant2^2, 1), 1, 4)
#         ##################################################*
#         ##ED50 prime!
#         ED50_prime <- matrix(c(0, param[3]/param[2]^2, -1/param[2], 0), 1, 4)
#         d_x_xi <-  g_x %*% FIM_inv %*% t(ED50_prime)  %*% ED50_prime %*% FIM_inv %*% t(g_x) - (ED50_prime %*% FIM_inv %*% t(ED50_prime))
#         return(d_x_xi)
#       }
#     }
#
#     if (round(multiple_arg$lambda[3], 5) == 1){
#       crfunc <- function(param, x, w, multiple){
#         FIM_val <- fimfunc2(x=x, w = w, param = param)
#         lambda <- multiple$lambda
#         delta <- multiple$delta
#         FIM_inv <- mpginv(FIM_val, tol = sqrt(.Machine$double.xmin))
#         #part3
#         if (param[2] > 0)
#           MED_prime <- matrix(c(-1/((param[1]+delta)*param[2]),
#                                 (param[3] - log(-delta/(param[1]+delta)))/param[2]^2,
#                                 -1/param[2], 0),1, 4)
#         if (param[2] < 0)
#           MED_prime <- matrix(c(1/((param[1]-delta)*param[2]),
#                                 (param[3] - log((param[1]-delta)/delta))/param[2]^2,
#                                 -1/param[2], 0), 1, 4)
#         var_MED <- MED_prime %*% FIM_inv %*% t(MED_prime)
#         part3 <- lambda[3] * log(var_MED)
#
#         locally_crfunc <- part3
#
#         return(locally_crfunc)
#       }
#       d_multi_x <- function(x1, FIM, x, w, param, multiple){
#
#         npar <- 4 #  four parameter logistic
#         lambda <- multiple$lambda
#         delta <- multiple$delta
#         #################### find the generalized inverse
#         FIM_val = FIM(x = x, w = w, param = param)
#         #FIM_inv <- solve(FIM_val, tol = .Machine$double.xmin)
#         # FIM_inv <- mpginv(FIM_val, tol =  sqrt(.Machine$double.xmin))
#         FIM_inv <- mpginv(FIM_val, tol =  sqrt(.Machine$double.eps))
#         ##################################################*
#         ####################### writign g
#         constant1 <- exp(param[2]*x1 + param[3])
#         constant2 <- 1/(1+constant1)
#         g_x <- matrix(c(constant2, - param[1] * x1 * constant1*constant2^2, -param[1] * constant1 * constant2^2, 1), 1, 4)
#         ##################################################*
#         ##MED prime!
#         if(param[2]>0)
#           MED_prime <- matrix(c(-1/((param[1] + delta) * param[2]),
#                                 (param[3] - log(-delta/(param[1] + delta)))/param[2]^2,
#                                 -1/param[2], 0), 1, 4)
#         if(param[2]<0)
#           MED_prime <- matrix(c(1/((param[1]-delta) * param[2]),
#                                 (param[3] - log((param[1] - delta)/delta))/param[2]^2,
#                                 -1/param[2], 0), 1, 4)
#         d_x_xi<-  g_x %*% FIM_inv %*% t(MED_prime)  %*% MED_prime %*% FIM_inv %*% t(g_x) - (MED_prime %*% FIM_inv %*% t(MED_prime))
#         return(d_x_xi)
#       }
#     }
#   }
#
#
#   Psi_x_integrand_multiple <- function(x1,  FIM,  x, w, param, prior_func, multiple){
#
#     deriv_integrand <- apply(param, 2, FUN = function(col_par)d_multi_x(x1 = x1, FIM = FIM, x = x,  w=w, param = col_par, multiple = multiple)) * prior_func(t(param))
#
#     dim(deriv_integrand) <- c(1, length(deriv_integrand))
#     return(deriv_integrand)
#   }
#   Psi_x_bayes_multiple <- function(x1, prior_func, FIM,  x, w, lp, up, npar, truncated_standard, multiple){
#     out <- hcubature(f = Psi_x_integrand_multiple, lowerLimit = lp, upperLimit = up, vectorInterface = TRUE,
#                      x = x, w = w, prior_func = prior_func, FIM = FIM, x1 = x1,
#                      multiple = multiple,
#                      tol = 1e-5, maxEval = 50000)
#
#     return(out$integral/truncated_standard)
#   }
#   # we dont need psi_Mu for locally version but we kepp it
#   return(list(crfunc = crfunc,  sensitivity = Psi_x_bayes_multiple))
# }
