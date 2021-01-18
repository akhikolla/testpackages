######################################################################################################*
######################################################################################################*
parvars_rightside <-function(fixedpars){
  ## return the character strings after modifiyin the
  NA_id <- which(is.na(fixedpars))
  if (length(NA_id) > 0){
    # no vectorization
    fixedpars[NA_id] <- paste("param[", 1:length(NA_id), "]", sep = "")
    char1 <- paste(" <-", fixedpars)
    # with vectorization
    fixedpars[NA_id] <- paste("param[, ", 1:length(NA_id), "]", sep = "")
    char2 <- paste(" <-", fixedpars)
  } else{
    # no vectorization
    fixedpars <- paste("param[", 1:length(fixedpars), "]", sep = "")
    char1 <- paste(" <-", fixedpars)
    # with vectorization
    fixedpars <- paste("param[, ", 1:length(fixedpars), "]", sep = "")
    char2 <- paste(" <-", fixedpars)
  }
  return(list(right_vec = char1, right_mat = char2))
}
######################################################################################################*
######################################################################################################*
#' @importFrom stats deriv
#' @importFrom stats as.formula
create_grad <- function (formula, predvars, parvars, paramvectorized = FALSE, fixedpars){
  #  if paramvectorized = TRUE, then we vectorized with respect to param not x, otherwise x
  ## number of parameters
  npar <- length(parvars)
  npred <- length(predvars)
  rightside <- parvars_rightside(fixedpars = fixedpars)
  # number of design points
  ## derviative with respect to parameters
  #DD <- deriv(expr = formula, namevec = parvars, function.arg = c(predvars, parvars))
  formula <- as.formula(formula)
  DD <- deriv(expr = formula, namevec = parvars)
  aDD <- as.character(DD)
  if (!paramvectorized){
    to <- paste("(", paste("npoints/", npred, "*" , 1:npred, sep  = ""), ")", sep  =  "")
    from <- c(1, paste("(", to[-length(to)], " +  1)", sep  = ""))
    temp1 <- "    "
    for(i in 1:length(predvars))
      temp1 <- paste(temp1, paste(predvars[i], " <- points[", from[i], ":", to[i], "]\n    ", sep = ""), sep = "")
    # from <- from + k
    # to <- to + k
    grad <- NA ## to define the variable in the global environment and avoid  check Note
    gradtext <- "grad <- function(points, param){\n npoints <- length(points) \n"
    gradtext <- paste(gradtext, temp1, sep = "")
    ### for parameters
    #gradtext <- paste(gradtext,  parvars[1], " <- param[1]", " \n    ", sep = "")
    gradtext <- paste(gradtext,  parvars[1], rightside$right_vec[1], " \n    ", sep = "")
    if (npar > 1) {
      for (j in 2:npar) {
        #gradtext <- paste(gradtext, parvars[j], " <- param[", j, "]", " \n  ", sep = "")
        gradtext <- paste(gradtext, parvars[j], rightside$right_vec[j], " \n  ", sep = "")
      }
    }
  }else{
    temp1 <- "    "
    for(i in 1:length(predvars))
      temp1 <- paste(temp1, paste(predvars[i], " <- points[", i, "]\n    ", sep = ""), sep = "")
    gradtext <- "grad <- function(points, param){ \n"
    gradtext <- paste(gradtext, temp1, sep = "")
    #gradtext <- paste(gradtext,  parvars[1], " <- param[, 1]", " \n    ", sep = "")
    gradtext <- paste(gradtext,  parvars[1], rightside$right_mat[1], " \n    ", sep = "")
    if (npar > 1) {
      for (j in 2:npar) {
        #gradtext <- paste(gradtext, parvars[j], " <- param[, ", j, "]", " \n  ", sep = "")
        gradtext <- paste(gradtext, parvars[j], rightside$right_mat[j], " \n  ", sep = "")
      }
    }
  }## end of paramvectorized = TRUE
  gradtext <- paste(gradtext, substr(x = aDD, start = 2, stop = nchar(aDD)), sep = "")
  eval(parse(text = gradtext))
  return(grad)
}



## fro exampe three points
# (1):(npoint/ndes * 1)
# (npoint/ndes +  1):(npoint/ndes * (ndes-1))
# (npoint/ndes * 2 +  1):(npoint/ndes * ndes)
######################################################################################################*
######################################################################################################*
create_mu <- function (formula, predvars, parvars, paramvectorized = FALSE, fixedpars){
  ## number of parameters
  formula <- as.formula(formula)
  npar <- length(parvars)
  rightside <- parvars_rightside(fixedpars = fixedpars)
  # number of design points
  npred <- length(predvars)
  if (!paramvectorized){
    to <- paste("(", paste("npoints/", npred, "*" , 1:npred, sep  = ""), ")", sep  =  "")
    from <- c(1, paste("(", to[-length(to)], " +  1)", sep  = ""))
    temp1 <- " "

    for(i in 1:npred)
      temp1 <- paste(temp1, paste(predvars[i], " <- points[", from[i], ":", to[i], "]\n ", sep = ""), sep = "")

    ## derviative with respect to parameters
    formula_char <- as.character(formula)
    formula_char <- formula_char[length(formula_char)]
    mu <- NA ## to define the variable in the global environment and avoid  check Note
    mutext <- "mu <- function(points, param){\n npoints <- length(points) \n"
    mutext <- paste(mutext, temp1, sep = "")
    ### for parameters
    #mutext <- paste(mutext, parvars[1], " <- param[1]", " \n ", sep = "")
    mutext <- paste(mutext, parvars[1], rightside$right_vec[1], " \n ", sep = "")

    if (npar > 1) {
      for (j in 2:npar) {
        #mutext <- paste(mutext, parvars[j], " <- param[", j, "]", " \n ", sep = "")
        mutext <- paste(mutext, parvars[j], rightside$right_vec[j], " \n ", sep = "")
      }
    }
  }else{
    temp1 <- " "
    for(i in 1:length(predvars))
      temp1 <- paste(temp1, paste(predvars[i], " <- points[", i, "]\n ", sep = ""), sep = "")
    ## derviative with respect to parameters
    formula_char <- as.character(formula)
    formula_char <- formula_char[length(formula_char)]
    mutext <- "mu <- function(points, param){ \n"
    mutext <- paste(mutext, temp1, sep = "")
    ### for parameters
    #mutext <- paste(mutext, parvars[1], " <- param[, 1]", " \n ", sep = "")
    mutext <- paste(mutext, parvars[1], rightside$right_mat[1], " \n ", sep = "")
    if (npar > 1) {
      for (j in 2:npar) {
        #mutext <- paste(mutext, parvars[j], " <- param[, ",j, "]", " \n ", sep = "")
        mutext <- paste(mutext, parvars[j], rightside$right_mat[j], " \n ", sep = "")
      }
    }
  }
  mutext <- paste(mutext,  formula_char, "\n}", sep = "")
  eval(parse(text = mutext))
  return(mu)
}



######################################################################################################*
######################################################################################################*
fim <- function(x, w, param, grad, mu, family, paramvectorized = FALSE){
  # ## create the functions that returnthe value of gradient and mu at each point
  # grad<- create_grad(formula = formula, predvars = predvars, parvars = parvars)
  # mu <- create_mu(formula = formula, predvars = predvars, parvars = parvars)
  if (paramvectorized && is.vector(param))
    stop("Bug: 'param' should be matrix not vector")
  if (!paramvectorized && is.matrix(param))
    stop("Bug: 'param' should be vector not matrix")
  if (is.null(family))
    stop("Bug: 'family' is NULL")

  kk <- length(w)
  if (!paramvectorized){
    grad_des <- attributes(grad(points = x, param = param))$gradient
    mu_des <- mu(points = x, param = param)
    var_des <- family$variance(mu_des)

    #grad_des_w <- grad_des %*% diag(x=sqrt(w/var_des), nrow = length(w))
    #FIM_mat <- grad_des_w %*% t(grad_des_w)
    #FIM_mat <- crossprod(.Internal(diag(x=sqrt(w/var_des), nrow = kk, ncol = kk)) %*% grad_des)
    FIM_mat <- crossprod(diag(x=sqrt(w/var_des), nrow = kk, ncol = kk) %*% grad_des)

    #FIM_mixed_inhibition(S = x[1:4], I = x[5:8], w = w, param = param)
    #(f%*%x) %*% (f%*%x)

  }else{
    npoint <- length(w)
    x_mat <- matrix(x, nrow = npoint)
    ## a list as length of length(w)
    # each row of [[i]] is the gradient with respect to each value of the set of parameters
    # the number of rows is each value of the parameters
    grad_des <- lapply(1:npoint,  FUN = function(j)attributes(grad(point = x_mat[j, ], param = param))$gradient)
    grad_des <- do.call(cbind, grad_des)
    #every row is for one set of parameters
    mu_des <- lapply(1:npoint,  FUN = function(j)mu(point = x_mat[j, ], param = param))
    mu_des <- do.call(cbind, mu_des)
    ## every row is the evaulated mu at design points (columns are th e design points)

    ## mu_des is of length dim(param)[1] * npoint
    if(!is.null(family))
      var_des   <- family$variance(mu_des)

    if (is.vector(var_des))
      var_des <- matrix(var_des, ncol = kk)
    #FIM_mat <- lapply(1:dim(param)[1], function(j)crossprod(.Internal(diag(x=sqrt(w/var_des[j, ]), nrow = kk, ncol = kk)) %*% matrix(grad_des[j, ], nrow = npoint, byrow = TRUE)))
    FIM_mat <- lapply(1:dim(param)[1], function(j)crossprod(diag(x=sqrt(w/var_des[j, ]), nrow = kk, ncol = kk) %*% matrix(grad_des[j, ], nrow = npoint, byrow = TRUE)))
    ## the output is a list of informataion matrix with respect to each values of set of parameters
    # the length is equal to the length of param_set
  }


  ## it is a list if design FIM matrices when  the param is a matrix (vectorized with respect to the parameters)
  # or a design matrix when param is vector

  # for (i in 1:dim(param)[1]){
  #   if(any(abs(FIM_logistic(x = x, w=w,param[i, ]) - FIM_mat[[i]]) > 1e-12))
  #     browser()
  # }
  return(FIM_mat)
}
######################################################################################################*
######################################################################################################*
create_prob <- function(prob, parvars, predvars){
  predvars <- gsub(" ", "", predvars, fixed = TRUE)
  allvars <- all.vars(prob)
  parvars_list <- strsplit(parvars, "=")
  fixedpars <- sapply(1:length(parvars_list), FUN = function(j)parvars_list[[j]][2])
  fixedpars <- gsub(" ", "", fixedpars, fixed = TRUE)
  parvars <- sapply(parvars_list, '[[', 1)
  parvars <- gsub(" ", "", parvars, fixed = TRUE)
  #num_unknown_param <- sum(is.na(fixedpars))

  rightside <- parvars_rightside(fixedpars = fixedpars)



  npar <- length(parvars)
  npred <- length(predvars)
  ######### pedictors
  to <- paste("(", paste("npoints/", npred, "*" , 1:npred, sep  = ""), ")", sep  =  "")
  from <- c(1, paste("(", to[-length(to)], " +  1)", sep  = ""))
  temp1 <- "    "
  for(i in 1:length(predvars))
    temp1 <- paste(temp1, paste(predvars[i], " <- x[", from[i], ":", to[i], "]\n    ", sep = ""), sep = "")
  ## create the function as text
  prob2 <- NA ## to define the variable in the global environment and avoid R CMD check Note
  probtext <- "prob2 <- function(x, param){\n npoints <- length(x) \n"
  probtext <- paste(probtext, temp1, sep = "")
  ### for parameters
  probtext <- paste(probtext,  parvars[1], rightside$right_vec[1], " \n    ", sep = "")
  if (npar > 1) {
    for (j in 2:npar) {
      #gradtext <- paste(gradtext, parvars[j], " <- param[", j, "]", " \n  ", sep = "")
      probtext <- paste(probtext, parvars[j], rightside$right_vec[j], " \n  ", sep = "")
    }
  }

  prob <- as.formula(prob)
  probtext <- paste(probtext, "\n out <- ", prob, "\n", "return(out)\n}")
  eval(parse(text = probtext))
  return(prob2)
}




# prob <- function(x, param){
#   npoint <- length(x)/2
#   x1 <- x[1:npoint]
#   x2 <- x[(npoint+1):(npoint*2)]
#   out <- 1-1/(1+exp(param[1] + param[2] * x1 + param[3] * x2 + param[4] * x1 * x2))
#   return(out)
# }


# prob_formula <- function(x, param){
#   npoints <- length(x)
#   x1 <- x[1:(npoints/2*1)]
#   x2 <- x[((npoints/2*1) +  1):(npoints/2*2)]
#   b0 <- param[1]
#   b1 <- param[2]
#   b2 <- param[3]
#   b3 <- param[4]
#
#   out <-  1 - 1/(1 + exp(b0 + b1 * x1 + b2 * x2 + b3 * x1 * x2))
#   return(out)
# }

######################################################################################################*
######################################################################################################*



######################################################################################################*
######################################################################################################*

