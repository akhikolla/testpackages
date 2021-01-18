# ###############################################################################################################*
# ###############################################################################################################*
print_optima_char <- function(mat, type, merge_tol){
  # if (npred == 1)
  #   return(x)
  if (type == "minimax"){
    optimachar <- "minima"
    optimumchar <- "minimum"
  }
  if (type == "standardized"){
    optimachar <- "maxima"
    ptimumchar <- "maximum"
  }

  npar <- ncol(mat) - 2
  answering <- mat[, ncol(mat)]
  answering[which(answering == 1)] <- "YES"
  answering[which(answering == 0)] <- "NO"
  answering <- c("In Answering Set?*", answering)
  answering  <- format(answering , width = max(nchar(answering)), justify = "centre")
  crt <- mat[, ncol(mat)-1]

  global_opt <- which.max(crt)[1]
  par_mat <- mat[, -c(ncol(mat), ncol(mat)-1), drop = FALSE]
  par_mat <- sprintf("%.5f", round(par_mat,5))
  crt <- sprintf("%.5f", round(crt,5))
  # warning: do not change the posistion of two following lines
  crt[global_opt] <- paste(crt[global_opt], "**", sep = "")
  crt <- c("Criterion Value", crt)



  parchar <- paste("Par", 1:npar, sep = " ")
  charwidth <- max(nchar(c(mat, parchar, crt)))


  crt <- format(crt, width = charwidth, justify = "centre")


  par_mat <- matrix(par_mat, ncol = npar)
  par_mat <- rbind(parchar, par_mat)
  par_mat <- sapply(X = 1:nrow(par_mat), function(j)paste("(", paste(par_mat[j, ], collapse = ","), ")", sep = ""))
  par_mat <- format(par_mat, width = charwidth, justify = "centre")



  par_mat <- paste("| ", par_mat, " | ", crt, " | ", answering, " |", "\n", sep = "", collapse = "")
  par_mat <- paste("\nAll the local ",  optimachar,
                   " over the parameter space found at the input design:\n\n", par_mat,
                   "* The Asnwering set is built by picking the local ", optimachar,  " nearest to the global ", optimumchar, " subject to a tolerance of ", merge_tol,
                   "\n** Global ",  optimumchar, sep = "", collapse = "")


  return(par_mat)
}

# ###############################################################################################################*
# ###############################################################################################################*
# control.answering <- function(n_seg = 6, merge_tol = .005){
#   if (!is.numeric(n_seg) || n_seg <= 0)
#     stop("value of 'n_seg' must be > 0")
#   if (!is.numeric(merge_tol) || merge_tol <= 0)
#     stop("value of 'merge_tol' must be > 0")
#   return(list(n_seg = n_seg, merge_tol = merge_tol))
# }
# ###############################################################################################################*
# ###############################################################################################################*
# ###############################################################################################################*
Calculate_Cost_minimax <- function(mat, fixed_arg){
  # mat: is the matrix of positions
  # fixed_arg: passed to calculate
  if(!is.matrix(mat))
    stop("'mat' must be matrix")
  n_cost <- dim(mat)[1]
  if (fixed_arg$type != "robust")
    inner_optima_store <- matrix(NA, ncol = length(fixed_arg$lp), nrow=  n_cost) else
      inner_optima_store <- matrix(NA, ncol = dim(fixed_arg$parset)[2], nrow=  n_cost)
  store <- vector("list", n_cost) #temporarily
  cost <- vector("double", n_cost)
  if(fixed_arg$equal_weight)
    w_equal <- rep(1/fixed_arg$k, fixed_arg$k)
  #if (!fixed_arg$parallel){
  for(i in 1:dim(mat)[1]){
    if (!fixed_arg$is.only.w)
      x <- mat[i, fixed_arg$x_id] else
        x <- NULL
      if(!fixed_arg$equal_weight)
        w <- mat[i, fixed_arg$w_id] else
          w <- w_equal
        if(fixed_arg$sym){
          x_w <- ICA_extract_x_w(x = x, w = w,
                                 sym_point = fixed_arg$sym_point)
          x <- x_w$x
          w <- x_w$w
        }
        optim_func <- fixed_arg$optim_func
        store[[i]] <- optim_func(fn = fixed_arg$crfunc,
                                 x = x,
                                 w = w,
                                 lower = fixed_arg$lp,
                                 upper = fixed_arg$up,
                                 fixedpar = fixed_arg$fixedpar,
                                 fixedpar_id = fixed_arg$fixedpar_id,
                                 npred= fixed_arg$npred)
        MinRowId <- which.min(store[[i]]$minima[, fixed_arg$CostColumnId])
        cost[i] <- store[[i]]$minima[MinRowId , fixed_arg$CostColumnId]
        inner_optima_store[i, ] <- store[[i]]$minima[MinRowId , -fixed_arg$CostColumnId]
        x <- w <-MinRowId <- NA
  }
  #}else{
  #################### parallel
  # store <- foreach(x = 1:dim(mat)[1])%dopar%{
  #     i <- x
  #     if (!fixed_arg$is.only.w)
  #       x <- mat[i, fixed_arg$x_id] else
  #         x <- NULL
  #       if(!fixed_arg$equal_weight)
  #         w <- mat[i, fixed_arg$w_id] else
  #           w <- w_equal
  #         if(fixed_arg$sym){
  #           x_w <- ICA_extract_x_w(x = x, w = w,
  #                                  sym_point = fixed_arg$sym_point)
  #           x <- x_w$x
  #           w <- x_w$w
  #         }
  #         temp <- fixed_arg$optim_func(fn = fixed_arg$crfunc,
  #                                      x = x,
  #                                      w = w,
  #                                      lower = fixed_arg$lp,
  #                                      upper = fixed_arg$up,
  #                                      fixedpar = fixed_arg$fixedpar,
  #                                      fixedpar_id = fixed_arg$fixedpar_id,
  #                                      npred= fixed_arg$npred)
  #         MinRowId <- which.min(temp$minima[, fixed_arg$CostColumnId])
  #         temp$cost <- temp$minima[MinRowId , fixed_arg$CostColumnId]
  #         temp$inner_optima_store <- temp$minima[MinRowId , -fixed_arg$CostColumnId]
  #         return(temp)
  #   }
  # }
  #################### parallel

  ##we multiply the cost function by negative for maximin and minimax optimal design!
  # if (fixed_arg$parallel){
  #   cost <-  sapply(store, "[[", "cost")
  #   inner_optima_store <-  t(sapply(store, "[[", "inner_optima_store"))
  # }
  if(fixed_arg$type != "locally" && fixed_arg$type != "robust"){
    cost <- -cost
    #inner_optima_store <- matrix(NA, nrow = length(cost))
  }

  nfeval <-  sum(sapply(store, "[[", "counts"))
  return(list(cost = cost, nfeval = nfeval, inner_optima = inner_optima_store))
}
# ###############################################################################################################*
# ###############################################################################################################*
find_measure <- function(npar, x, w, answering, FIM, Psi_Mu){
  ## retrun the measure mu for given x, w and A(\xi)
  # FIM is function
  if (!is.matrix(answering))
    stop("'answering' must be a matrix")
  if (dim(answering)[1] == 1){
    mu <- 1
    System <- "One equation"
  }else{
    n_mu <-  dim(answering)[1] ## number of  measures
    if (n_mu <= 50)
      mu <- optim(par = rep(1/dim(answering)[1], dim(answering)[1]),
                  fn=Psi_Mu, lower=rep(0, n_mu), upper=rep(1, n_mu),
                  control = list(factr =  1e+7), method="L-BFGS-B",
                  answering = answering, x = x, w = w, FIM = FIM,
                  PenaltyCoeff = 50)$par

    # mu <- optimx(par = rep(1/dim(answering)[1], dim(answering)[1]), fn = Psi_Mu,
    #        lower=rep(0, n_mu), upper=rep(1, n_mu),  answering = answering, x = x, w = w, FIM = FIM,
    #        method="bobyqa",
    #        PenaltyCoeff = 50, itnmax = 1000)

    if (n_mu > 50)
      mu <- nloptr::nloptr(x0= rep(1/dim(answering)[1], dim(answering)[1]),
                           eval_f = Psi_Mu, lb = rep(0, n_mu), ub = rep(1, n_mu),
                           opts = list(algorithm = "NLOPT_LN_BOBYQA",
                                       xtol_rel = 1e-4,
                                       ftol_rel = 1e-4,
                                       maxeval = 100),
                           answering = answering, x = x, w = w, FIM = FIM,
                           PenaltyCoeff = 50)$solution

  }
  if (round(sum(mu), 3) != 0)
    mu <- mu/sum(mu)
  return(list(mu = mu))
}
# ###############################################################################################################*
# ###############################################################################################################*
find_nearest <- function(values, tol, compare_with_minimum){
  # returns the indices of "values" that are nearer tothe  minimum or the maximim with 'tol'
  if(compare_with_minimum)
    global_index <- which.min(values) else
      global_index <- which.max(values)
    global <- values[global_index]
    keep <- which(values - global < tol)
    return(keep)
}
# ###############################################################################################################*
# ###############################################################################################################*
# create_criterion_minimax <- function(FIM, type, localdes = NULL){
#   if (type == "minimax") {
#     crfunc_minimax <- function(param, q, npred) {
#       lq <- length(q) # q is the design points and weights
#       pieces <- lq / (npred + 1)
#       x_ind <- 1:(npred * pieces)
#       w_ind <- (x_ind[length(x_ind)] + 1):lq
#       x <- q[x_ind]
#       w <- q[w_ind]
#       if (is.vector(param)){
#         # here the fim only accepts param as vector
#         minimax_crfunc <-  det2(FIM(x = x, w = w, param = param), logarithm = TRUE) - 5000 * (sum(w) - 1) ^ 2   ## -(-det+pen) = det-pen
#       }
#
#       if (is.matrix(param)){
#         FIM_list <- FIM(x = x, w = w, param = param)
#         minimax_crfunc <- sapply(1:dim(param)[1], function(j) det2(FIM_list[[j]], logarithm = TRUE))
#         minimax_crfunc <- minimax_crfunc - 5000 * (sum(w) - 1) ^ 2
#       }
#       return(minimax_crfunc)
#       localdes = NULL
#     }
#     # we dont need minus for minimax because we want to maximize the -log(det) or minimze log(det)
#   }
#   if (type == "standardized" && is.null(localdes))
#     stop("'localdes' must be given for standardized optimal designs")
#
#   if (type == "standardized") {
#     crfunc_standardized <- function(q, param, npred) {
#       lq <- length(q)
#       pieces <- lq/(npred + 1)
#       x_ind <- 1:(npred * pieces)
#       w_ind <- (x_ind[length(x_ind)] + 1):lq
#       x <- q[x_ind]
#       w <- q[w_ind]
#       denominator <- localdes(param = param)
#       numerator <- det2(FIM(x = x, w = w, param = param), logarithm = FALSE)
#       fraction <- (numerator/denominator)
#       if (round(fraction, 7) > 1)
#         stop("locally D-efficiency of the non-optimal design is higher than the locally optimal design for",
#              paste("iniparam = c(", paste(round(param, 5), collapse = ","), ").",sep = ""),
#              "\nProbably the generated design by 'localdes' is not the true locally D-optimal design.\n",
#              "Check 'localdes' function and be sure that it returns the true locally optimal design")
#       if (npar %% 2 != 0) {
#         maximin_crfunc <- (fraction)^(1/npar)
#       }else{
#         maximin_crfunc <-  ifelse(fraction < 0, 0,(fraction)^(1/npar))
#       }
#       return(maximin_crfunc)
#     }
#   }
#   if (type[1] == "locally") {
#     crfunc_locally <- function(param, q, npred) {
#       lq <- length(q)
#       pieces <- lq / (npred + 1)
#       x_ind <- 1:(npred * pieces)
#       w_ind <- (x_ind[length(x_ind)] + 1):lq
#       x <- q[x_ind]
#       w <- q[w_ind]
#       locally_det <-
#         -det2(FIM(x = x, w = w, param = param), logarithm = TRUE) + 5000 * (sum(w) - 1) ^ 2
#       if (locally_det == -1e24)
#         ## becuase for locally the 'locally_det' will be -1e24 and spoil the algorithm!
#         locally_det <- 1e24
#       return(locally_det)
#     }
#   }
#   crfunc <- switch(type, "minimax" = crfunc_minimax, "standardized" = crfunc_standardized, "locally" = crfunc_locally)
#   return(crfunc)
# }
#############################################################################################################*
#############################################################################################################*
paste_mat <- function(mat){
  # cat matrix
  temp <- ""
  for(i in 1:dim(mat)[1])
    temp <- paste(temp, "\n", paste(mat[i, ], sep = "", collapse = "  "),  sep = "")
  return(temp)
}

###############################################################################################################*
#############################################################################################################*


#############################################################################################################*
#############################################################################################################*
# check_maxeval <- function(fn, lower, upper, maxeval, ...){
#   ## this function check whether a higher value of maxeval will give a better value or not
#   if (missing(fn))
#     stop("'fn' is missing")
#   if (missing(fn))
#     stop("'lower' is missing")
#   if (missing(upper))
#     stop("'upper' is missing")
#   if (missing(maxeval))
#     stop("'maxeval' is missing")
#   ## WARNINGS: fn should have argument 'param'
#   fn1 <- function(param) fn(param, ...)
#
#   # fn1(param = c(1, 1))
#   opt <-list()
#
#   maxeval_vec <- maxeval + c(100, 200, 300, 500, 600, 800, 900, 1000, 1500, 2000, 2500, 3000, 4000, 5000, 6000, 10000, 20000, 30000)
#   for(i in 1:length(maxeval_vec))
#     opt[[i]] <- directL(fn = fn1, lower = lower, upper = upper,
#                         nl.info = FALSE,
#                         control=list(xtol_rel=.Machine$double.eps,
#                                      maxeval = maxeval_vec[i]))
#   val <- sapply(opt, "[[", 'value')
#
#   ## all equal?
#   all.eq <- (abs(max(val) - min(val)) < .Machine$double.eps ^ 0.5)
#
#   ############################################*
#   ## we produce a warnings if the points that directL missed is NOT on the vertices, because if it was on any of the vertices
#   ## cause no problem as our optimization function cosider this case as well.
#   vertices <- make_vertices(lower = lower, upper = upper)
#   vertices_val <- find_on_points(fn = fn1, points = vertices)
#   vertices_val <- min(vertices_val$minima[, dim(vertices_val$minima)[2]])
#
#
#
#   is.min_on_vertices <- ((vertices_val - min(val)) <  .Machine$double.eps ^ 0.5)
#   ##################################################*
#
#   if (!all.eq && !is.min_on_vertices){
#     temp_id <- which(abs(val - val[length(val)]) < .Machine$double.eps ^ 0.5)
#     maxeval_recom <- maxeval_vec[temp_id[1]]
#     msg <- paste("'inner_maxeval' recommended to be larger than ",   maxeval_recom, "\nPlease also consider increasing 'n_seg' if the value of 'inner_maxeval' is large.")
#     warning(msg, call. = FALSE)
#   }else{
#     maxeval_recom <- maxeval
#     msg <- NULL
#   }
#   return(list(msg = msg, maxeval = maxeval_recom, opt = opt))
# }

#############################################################################################################*
#############################################################################################################*
##maybe removed later
#############################################################################################################*
#############################################################################################################*
# check_inner <- function(obj, maxeval){
#   ### obj is an object of class "ICA"
#   ## only when type = "minimax" or "standardized"
#   it <- length(obj$evol)
#   obj$evol[[it]]$x
#   out <- nloptr::direct(fn = obj$arg$crfunc, lower = obj$arg$lp, upper = obj$arg$up,
#                         q = c(obj$evol[[it]]$x, obj$evol[[it]]$w),
#                         npred = length(obj$arg$lx),
#                         control=list(xtol_rel = .Machine$double.eps,
#                                      maxeval = maxeval) )
#   return(out)
# }
#############################################################################################################*
#############################################################################################################*
check_minimax_args <- function(lp, up, type, localdes, localdes_par, parvars, fimfunc, crt.minimax.control, crtfunc){


  #### minimax version
  if (missing(lp))
    stop("\"lp\" is  missing")
  if (missing(up))
    stop("\"up\" is missing")
  if (length(lp) != length(up))
    stop("length of \"lp\" is not equal to length of \"up\"")
  if (!type[1] %in% c("minimax", "standardized", "locally"))
    stop("\"type\" must be \"minimax\", \"standardized\" or  \"locally\"")
  if (type[1] == "standardized" && is.null(localdes))
    stop("Please specify the 'localdes'")
  if (type[1] == "standardized"){

    ## check what localdes takes as arguments
    if (!is.null(fimfunc)){
      if (!any(formalArgs(localdes) == "param"))
        stop("'localdes' must have an argument named 'param'")
    }else{
      if (!all(formalArgs(localdes) %in% c(parvars, "fimfunc")))
        stop("arguments of 'localdes' must be fimfunc, ", paste(parvars, collapse = ", "))
    }

    ## check what local des returns
    if (crt.minimax.control$inner_space == "discrete")
      localdes_check <- localdes_par(param = crt.minimax.control$param_set[1, ]) else
        localdes_check <- localdes_par(param = (lp + up)/2)
      if (!is.list(localdes_check))
        stop("'localdes' must be a list")
      if (!all(names(localdes_check) %in% c("x", "w")))
        stop("'localdes' must return a list with components 'x' and 'w'")
  }
  # if (type[1] == "locally")
  #   iniparas <- lp
  #   return(iniparas)
}
#############################################################################################################*
#############################################################################################################*
# localdes <- mixed_inhibition_formula
# parvars <- c("V", "Km", "Kic", "Kiu")

create_localdes <- function(parvars, localdes){
  npar <- length(parvars)
  localdes_formula <- NA ## to define the variable in the global environment and avoid R CMD check Note
  localdes_char <- "localdes_formula <- function(param)\n{\n  "
  ### for parameters
  localdes_char <- paste(localdes_char,  parvars[1], " <- param[1]", " \n  ", sep = "")
  if (npar > 1) {
    for (j in 2:npar) {
      localdes_char <- paste(localdes_char, parvars[j], " <- param[", j, "]", " \n  ", sep = "")
    }

  }
  localdes_char_func <- paste("out <- localdes(", paste(parvars, parvars, sep = " = ", collapse = ", "), ")", sep = "")
  localdes_char <- paste(localdes_char, localdes_char_func, "\n  return(out)\n}",  sep = "")
  #cat(localdes_char)
  eval(parse(text = localdes_char))
  return(localdes_formula)
}
#############################################################################################################*
#############################################################################################################*
create_optim_func <- function(type, lp_nofixed, up_nofixed, crt.minimax.control, discrete_set, robpars, inipars, is.only.w){
  # defining the nloptr directL
  if (type != "robust")
    vertices_inner <- make_vertices(lower = lp_nofixed, upper = up_nofixed)
  optim_nloptr <- function(fn, lower, upper, w, x, fixedpar, fixedpar_id,  npred){
    if (is.null(crt.minimax.control$x0))
      x0 <- (lower + upper)/2 else
        x0 <- crt.minimax.control$x0

      out_nloptr <- nloptr::nloptr(x0= x0, eval_f = fn, lb = lower, ub = upper,
                                   opts = crt.minimax.control$optslist,
                                   q = c(x, w), fixedpar = fixedpar,
                                   fixedpar_id = fixedpar_id, npred = npred)


      out <- find_on_points(fn = fn,
                            points = vertices_inner,
                            q = c(x, w),
                            fixedpar = fixedpar,
                            fixedpar_id = fixedpar_id,
                            npred = npred)

      counts <- out$count + out_nloptr$iterations
      #counts <- out$count + out_nloptr$iter
      minima <- out$minima
      minima_nloptr <- c(out_nloptr$solution, out_nloptr$objective)
      #minima_nloptr <- c(out_nloptr$par, out_nloptr$value)
      minima <- rbind(minima, minima_nloptr)
      return(list(minima =minima, counts = counts))
  }
  ########################*

  # for locally optimal design we dont need an inner optimization
  if (type == "locally")
    optim_locally <- function(fn, lower, upper, w, x, fixedpar, fixedpar_id,  npred){
      log_det <- fn(param = inipars, q=c(x, w),  fixedpar = fixedpar,
                    fixedpar_id = fixedpar_id, npred = npred)
      counts <- 1
      minima <- matrix(c(inipars, log_det), nrow = 1)
      return(list(minima =minima, counts = counts))
    }
  ########################*
  # for locally optimal design we dont need an inner optimization
  if (type == "robust"){
    if (is.null(robpars))
      stop("BUG: 'robpars is NULL when in 'create_optim_func'")
    robust_value <- function(fn, lower, upper, w, x, fixedpar, fixedpar_id,  npred){

      log_det <- fn(param = robpars$parset, q=c(x, w),  fixedpar = fixedpar,
                    fixedpar_id = fixedpar_id, npred = npred)
      counts <- 1
      minima <- cbind(robpars$parset, log_det)
      return(list(minima =minima, counts = counts))
    }
  }
  ########################*

  ## set the optim_func that will be used for calculating the cost values
  if (type != "locally" && type != "robust"){
    if (crt.minimax.control$inner_space == "continuous"){
      optim_func <- optim_nloptr
    }
    if (crt.minimax.control$inner_space == "vertices"){
      vertices <- make_vertices(lower = lp_nofixed, upper = up_nofixed)
      optim_func <- function(fn, lower, upper, w, x,
                             fixedpar, fixedpar_id, npred){
        out <- find_on_points(fn = fn,
                              points = vertices,
                              q = c(x, w),
                              fixedpar = fixedpar,
                              fixedpar_id = fixedpar_id,
                              npred = npred)
        minima <- out$minima
        counts <- out$counts
        return(list(minima =minima, counts = counts))
      }
    }

    if (crt.minimax.control$inner_space == "discrete"){
      if (is.null(discrete_set))
        stop("Bug: when inner_space = 'discrete', 'discrete_set must be given'")

      optim_func <- function(fn, lower, upper, w, x, fixedpar, fixedpar_id, npred){
        cost  <- fn(param = discrete_set, q = c(x, w), fixedpar = fixedpar, fixedpar_id = fixedpar_id,npred = npred)
        minima <- cbind(discrete_set, cost, deparse.level = 0)
        counts <- length(discrete_set)[1]
        return(list(minima =minima, counts = counts))
      }
    }

  }
  if (type == "locally")
    optim_func <- optim_locally
  if(type == "robust")
    optim_func <- robust_value
  return(optim_func)
}
#############################################################################################################*
#############################################################################################################*
create_criterion_minimax <- function(FIM, type, lp, up, localdes = NULL, npar, robpars = NULL, crt_type, multipars = NULL, compound = NULL, is.only.w, only_w_varlist = NULL, user_crtfunc2){
  if (crt_type == "D"){
    if (type == "minimax") {
      crfunc_minimax <- function(param, q, npred) {
        if (!is.only.w){
          lq <- length(q) # q is the design points and weights
          pieces <- lq / (npred + 1)
          x_ind <- 1:(npred * pieces)
          w_ind <- (x_ind[length(x_ind)] + 1):lq
          x <- q[x_ind]
          w <- q[w_ind]
        }else{
          w <- q
          x <- only_w_varlist$x
        }
        if (is.vector(param)){
          # here the fim only accepts param as vector
          minimax_crfunc <-  det2(FIM(x = x, w = w, param = param), logarithm = TRUE) - 5000 * (sum(w) - 1) ^ 2   ## -(-det+pen) = det-pen
        }
        if (is.matrix(param)){
          FIM_list <- FIM(x = x, w = w, param = param)
          if (!is.list(FIM_list))
            stop("Bug: FIM must return list of matrices when param is a matrix (discrete space)")
          minimax_crfunc <- sapply(1:length(FIM_list), function(j) det2(FIM_list[[j]], logarithm = TRUE))
          #minimax_crfunc <- sapply(1:dim(param)[1], function(j) det2(FIM_list[[j]], logarithm = TRUE))
          minimax_crfunc <- minimax_crfunc - 5000 * (sum(w) - 1) ^ 2
        }
        return(minimax_crfunc)
        localdes = NULL
      }



      # we dont need minus for minimax because we want to maximize the -log(det) or minimze log(det)
    }
    if (type == "standardized" && is.null(localdes))
      stop("'localdes' must be given for standardized optimal designs")

    if (type == "standardized") {
      crfunc_standardized <- function(q, param, npred) {
        if (!is.only.w){
          lq <- length(q)
          pieces <- lq/(npred + 1)
          x_ind <- 1:(npred * pieces)
          w_ind <- (x_ind[length(x_ind)] + 1):lq
          x <- q[x_ind]
          w <- q[w_ind]
        }else{
          w <- q
          x <- only_w_varlist$x
        }
        if (is.matrix(param)){
          # localdes_check<- localdes(param = param[1, ])
          # if (!is.list(localdes_check))
          #   stop("'localdes' must be a list")
          # if (!all(names(localdes_check) %in% c("x", "w")))
          #    stop("'localdes' must return a list with components 'x' and 'w'")
          #denominator <- apply(param, 1, localdes)
          denominator <- sapply(1:dim(param)[1], FUN = function(j){localdes_res <- localdes(param = param[j, ]); det2(FIM(x = localdes_res$x, w = localdes_res$w, param = param[j, ])[[1]], logarithm = FALSE)})
          if (length(denominator) != dim(param)[1])
            stop("'localdes' for each vector of param must return one value")
          numerator <- sapply(1:dim(param)[1], FUN = function(j)det2(FIM(x = x, w = w, param = param[j, ])[[1]], logarithm = FALSE))
          if (length(denominator) != length(numerator))
            stop("Bug: 'length of nominator is not equal o the denominator in standardized criterion'")
          eff <- (numerator/denominator)
          if (any(is.nan(eff)))
            stop("The criterion (D-efficiency) value is 'NaN' for", paste("iniparam = c(", paste(round(param[which(is.nan(eff))[1],], 5), collapse = ","), ").",sep = ""),
                 "\nCheck the Fisher information matrix, number of design points, lx, ux and .... for possible unmatched design parameters.")
          if (any(round(eff, 7) > 1))
            stop("Efficiency value ", eff[which(round(eff , 7) > 1)[1]],
                 " is greater than one for for ",
                 paste("iniparam = c(", paste(round(param[which(round(eff , 7) > 1)[1],], 5), collapse = ","), ").",sep = ""),
                 "\n  Probably 'localdes' does not return the true locally optimal design for te given iniparam.")

        }
        if (is.vector(param)){
          localdes_res <- localdes(param = param)
          # if (!is.list(localdes_res ))
          #   stop("'localdes' must be a list")
          # if (!all(names(localdes_res) %in% c("x", "w")))
          #   stop("'localdes' must return a list with components 'x' and 'w'")
          denominator <- det2(FIM(x = localdes_res$x, w = localdes_res$w, param = param), logarithm = FALSE)
          numerator <- det2(FIM(x = x, w = w, param = param), logarithm = FALSE)
          eff <- (numerator/denominator)
          if (is.nan(eff))
            stop("The criterion (D-efficiency) value is 'NaN' for ", paste("iniparam = c(", paste(round(param, 5), collapse = ","), ").",sep = ""),
                 "\nCheck the Fisher information matrix, number of design points, lx, ux and .... for possible unmatched design parameters.")
          if (round(eff, 7) > 1)
            stop("Efficiency value ", eff,
                 " is greater than one. This results in a wrong conclusion that the non-optimal design is better than the true optimal design when ",
                 paste("iniparam = c(", paste(round(param, 5), collapse = ","), ").",sep = ""),
                 "\n  Probably 'localdes' does not return the true locally optimal designfor iniparam.")
        }
        if (npar %% 2 != 0) {
          eff <- (eff)^(1/npar)
        }else{
          eff <-  ifelse(eff < 0, 0,(eff)^(1/npar))
        }
        return(eff)
      }
    }
    if (type[1] == "locally") {
      crfunc_locally <- function(param, q, npred) {
        if (!is.only.w){
          lq <- length(q)
          pieces <- lq / (npred + 1)
          x_ind <- 1:(npred * pieces)
          w_ind <- (x_ind[length(x_ind)] + 1):lq
          x <- q[x_ind]
          w <- q[w_ind]
        } else{
          w <- q
          x <- only_w_varlist$x
        }
        locally_det <- -det2(FIM(x = x, w = w, param = param), logarithm = TRUE) + 5000 * (sum(w) - 1) ^ 2
        # if (!is.finite(locally_det))
        #   locally_det <- 1e+24
        return(locally_det)
      }
    }

    if (type == "robust") {
      if (is.null(robpars$prob) || is.null(robpars$parset))
        stop("'prob' and 'parset' must be given in the list 'robpars'")
      # FIM_list <- FIM(x = x, w = w, param = param)
      # if (!is.list(FIM_list))
      #   stop("Bug: FIM must return list of matrices when param is a matrix (ave design)")
      # ave_crfunc <- sapply(1:length(FIM_list), function(j) robpars$prob[j] * det2(FIM_list[[j]], logarithm = TRUE))
      # minimax_crfunc <- minimax_crfunc - 5000 * (sum(w) - 1) ^ 2

      crfunc_ave <- function(param, q, npred) {
        ## param here is a matrix: parset: a set of parameters
        if (!is.only.w){
          lq <- length(q) # q is the design points and weights
          pieces <- lq / (npred + 1)
          x_ind <- 1:(npred * pieces)
          w_ind <- (x_ind[length(x_ind)] + 1):lq
          x <- q[x_ind]
          w <- q[w_ind]
        } else{
          w <- q
          x <- only_w_varlist$x
        }
        on_average_crfunc <- sum(
          sapply(1:length(robpars$prob), FUN= function(j)
            robpars$prob[j] * -det2(FIM(x = x, w = w, param = param[j, ]), logarithm = TRUE))
        ) + 5000 * (sum(w) - 1) ^ 2   ## -(-det+pen) = det-pen
        return(on_average_crfunc)
      }
    }


    crfunc <- switch(type, "minimax" = crfunc_minimax, "standardized" = crfunc_standardized, "locally" = crfunc_locally, "robust" = crfunc_ave)
  }
  if (crt_type == "user"){
    if (type == "standarized" & crt_type == "user")
      stop("'standardized' criterion is not defined for the user given criterion.")
    if (type == "minimax") {
      crfunc_minimax <- function(param, q, npred) {
        if (!is.only.w){
          lq <- length(q) # q is the design points and weights
          pieces <- lq / (npred + 1)
          x_ind <- 1:(npred * pieces)
          w_ind <- (x_ind[length(x_ind)] + 1):lq
          x <- q[x_ind]
          w <- q[w_ind]
        }else{
          w <- q
          x <- only_w_varlist$x
        }
        if (is.vector(param)){
          # here the fim only accepts param as vector
          minimax_crfunc <-  -user_crtfunc2(x = x, w = w, param = param, fimfunc = FIM)
        }
        if (is.matrix(param)){
          minimax_crfunc <- sapply(1:dim(param)[1], function(j) -user_crtfunc2(x = x, w = w, param = param[j, ], fimfunc = FIM))

          # #stop("BUG:there is no parallelization is designed for user matrix")
          # FIM_list <- user_FIM_paramvectorized(x = x, w = w, param = param)
          # if (!is.list(FIM_list))
          #   stop("Bug: FIM must return list of matrices when param is a matrix (discrete space)")
          # minimax_crfunc <- sapply(1:length(FIM_list), function(j) -user_crtfunc2(x = x, w = w, param = param, fimfunc = FIM_list[[j]]))
        }
        return(minimax_crfunc)
      }
    }
    if (type[1] == "locally") {
      crfunc_locally <- function(param, q, npred) {
        if (!is.only.w){
          lq <- length(q)
          pieces <- lq / (npred + 1)
          x_ind <- 1:(npred * pieces)
          w_ind <- (x_ind[length(x_ind)] + 1):lq
          x <- q[x_ind]
          w <- q[w_ind]
        } else{
          w <- q
          x <- only_w_varlist$x
        }
        locally_out <- user_crtfunc2(x = x, w = w, param = param, fimfunc = FIM)
        return(locally_out)
      }
    }

    if (type == "robust") {
      if (is.null(robpars$prob) || is.null(robpars$parset))
        stop("'prob' and 'parset' must be given in the list 'robpars'")

      crfunc_ave <- function(param, q, npred) {
        ## param here is a matrix: parset: a set of parameters
        if (!is.only.w){
          lq <- length(q) # q is the design points and weights
          pieces <- lq / (npred + 1)
          x_ind <- 1:(npred * pieces)
          w_ind <- (x_ind[length(x_ind)] + 1):lq
          x <- q[x_ind]
          w <- q[w_ind]
        } else{
          w <- q
          x <- only_w_varlist$x
        }
        on_average_crfunc <- sum(
          sapply(1:length(robpars$prob), FUN= function(j)
            robpars$prob[j] * user_crtfunc2(x = x, w = w, param = param[j, ], fimfunc = FIM))
        )
        return(on_average_crfunc)
      }
    }
    crfunc <- switch(type, "minimax" = crfunc_minimax,  "locally" = crfunc_locally, "robust" = crfunc_ave)
  }
  if (crt_type == "multiple"){
    temp1<- create_multiple_minimax(model = "4pl", FIM = FIM, multipars = multipars, is.only.w = is.only.w, only_w_varlist = only_w_varlist)
    crfunc <- temp1$crfunc
  }

  if (crt_type == "DPA"){
    crfunc <- function(param, q, npred){
      if (!is.only.w){
        lq <- length(q) # q is the design points and weights
        pieces <- lq / (npred + 1)
        x_ind <- 1:(npred * pieces)
        w_ind <- (x_ind[length(x_ind)] + 1):lq
        x <- q[x_ind]
        w <- q[w_ind]
      } else{
        w <- q
        x <- only_w_varlist$x
      }
      if (compound$alpha != 0)
        bcrfunc1 <- -(compound$alpha/4  * det2(FIM(x = x, w = w, param = param),logarithm= TRUE) +
                        (1 -compound$alpha) * log(sum(w * compound$prob(x = x, param = param)))) else
                          bcrfunc1 <-  -log(sum(w * compound$prob(x = x, param = param)))

                        return(bcrfunc1)
    }
    if (crt_type == "DPM")

      crfunc <- function(param, q, npred){
        if (!is.only.w){
          lq <- length(q) # q is the design points and weights
          pieces <- lq / (npred + 1)
          x_ind <- 1:(npred * pieces)
          w_ind <- (x_ind[length(x_ind)] + 1):lq
          x <- q[x_ind]
          w <- q[w_ind]
        } else{
          w <- q
          x <- only_w_varlist$x
        }
        bcrfunc1 <- -(compound$alpha/4  * det2(FIM(x = x, w = w, param = param),logarithm= TRUE) +
                        (1 -compound$alpha) * log(min(compound$prob(x = x, param = param))))
        return(bcrfunc1)
      }

  }

  if (type != "locally" &  type != "robust" ){
    # here we search if one of the parameters are fixed. then we pass it to the optimization function in the inner problem because otherwise it may casue an error.
    any_fixed <- sapply(1:length(lp), function(i) lp [i] == up[i]) # is a vector
    if (any(any_fixed)){
      is_fixed <- TRUE
      fixedpar_id <- which(any_fixed)
      fixedpar <- lp[fixedpar_id]
      lp_nofixed <- lp[-fixedpar_id]
      up_nofixed <- up[-fixedpar_id]
      ## warning: lp and up are channged here if 'any_fix == TRUE'
    }else{
      fixedpar <- NA
      fixedpar_id <- NA
      is_fixed <- FALSE
      lp_nofixed <- lp
      up_nofixed <- up
    }
  }else{
    fixedpar <- NA
    fixedpar_id <- NA
    is_fixed <- FALSE
    lp_nofixed <- lp
    up_nofixed <- up
  }

  ##################################################################*
  ### re-defimimg crfunc to handle fixed parameters.
  if (is_fixed){
    crfunc2 <- function(param, q, fixedpar = NA, fixedpar_id = NA, npred){
      param1 <- insert_fixed(param = param, fixedpar = fixedpar, fixedpar_id = fixedpar_id)
      out <- crfunc(param = param1, q = q, npred = npred)
      return(out)
    }
  }else{
    crfunc2 <- function(param, q, fixedpar = NA, fixedpar_id = NA, npred){
      out <- crfunc(param = param, q = q, npred = npred)
      return(out)
    }
  }
  #####################################################################*
  return(list(crfunc = crfunc2, fixedpar = fixedpar, fixedpar_id = fixedpar_id, is_fixed = is_fixed, lp_nofixed = lp_nofixed, up_nofixed = up_nofixed))
}
#############################################################################################################*
#############################################################################################################*
make_grid <- function(lp, up, n.grid, on_perimeter = FALSE){
  seq2 <- Vectorize(seq.default, vectorize.args = c("from", "to"))
  seq_res <- seq2(lp, up, length.out =  n.grid)
  seq_res <- lapply(seq_len(ncol(seq_res)), function(i) seq_res[,i])
  all_grid <- expand.grid(seq_res)
  all_grid <- as.matrix(all_grid)
  if (!on_perimeter)
    param_set <- all_grid else{
      on_surface <- apply(all_grid, 1, function(col) any(col == lp) || any(col == up))
      surface_grid <- all_grid[on_surface,, drop = FALSE]
      param_set <- as.matrix(surface_grid)
      param_set <- unique(param_set)
    }
  if (ncol(param_set) != length(lp))
    stop("Bug: each row of the param_set should be an initial estimates of the prameters")
  return(param_set)
}
#############################################################################################################*
#############################################################################################################*
insert_fixed <- function(fixedpar, fixedpar_id, param){
  if (!is.matrix(param) & !is.vector(param))
    stop("Bug: 'param' must be a matrix or vector when inserting a column or vector!")
  if (is.vector(param)){
    param_out <- rep(NA, length(param) + length(fixedpar))
    param_out[fixedpar_id] <- fixedpar
    param_out[-fixedpar_id] <- param
  }
  if (is.matrix(param)){
    param_out  <- matrix(NA, ncol =  dim(param)[2] + length(fixedpar), nrow = dim(param)[1])
    param_out[, fixedpar_id] <- fixedpar
    param_out[, -fixedpar_id] <- param
  }
  return(param_out)
}

create_psy_minimax <- function(crt_type, multipars, compound, user_sensfunc){
  if (crt_type == "D"){
    Psi_Mu <- function(mu, FIM,  x, w,  answering, PenaltyCoeff){
      # mu is a vector of measure that usually Psi_mu is optimized with respect to.
      # FIM: is the Fisher information matrix as function
      # x: vector of design points
      # w: vector of design weights
      # answering: the matrix of elements of answering set. Each row is an element.
      # answering set are '\mu'in (3) in ICA paper.
      # return the value of left hand side of (4) in ICA paper.
      # if mu = 1  and answering has only one row, then we computing the equivalence theorem left hand side for locally D_optimal design.
      # also can be used for multidimensional models like enzyme kinetic
      if(length(mu) != dim(answering)[1])
        stop("The number of measures is not equal to the number of elements of answering set.")
      if(typeof(FIM) != "closure")
        stop("'FIM' must be of type 'closure'.")
      CardOfRegion <- dim(answering)[2] ##cardinal of region of uncertainty
      n_mu <- dim(answering)[1] ## number of mu, measures

      npred <- length(x)/length(w) # number of independent variables.
      one_point_mat <- matrix(x, length(w), npred)

      Psi_Point_answering <- matrix(NA, length(w),  n_mu)
      ## each row is  tr(M^{-1}(\xi, mu_j) %*% I(x_i, \mu_j))
      # so sum of each row of 'Psi_Point_answering' is \int tr(M^{-1}(\xi, mu) %*% I(x_i, \mu))
      for(i in 1:length(w)){
        for(j in 1:n_mu){
          Psi_Point_answering[i, j]  <-
            mu[j] *
            sum(diag(solve(FIM(x = x, w = w, param = answering[j, ]), tol = .Machine$double.xmin) %*%
                       FIM(x = as.vector(one_point_mat[i,]), w = 1, param = answering[j, ])
            ))
        }
      }
      #Psi at each points x = poi
      #Psi at each Point
      PsiAtEachPoint <- round(rowSums(Psi_Point_answering) - CardOfRegion, 5)

      ##now we should creat the penalty function. Psi must be zero at each supoport point
      PsiEqualityPenalty <- sum(PenaltyCoeff*pmax(PsiAtEachPoint, 0)^2 + PenaltyCoeff*pmax(-PsiAtEachPoint, 0)^2)
      PsiFunction <-  PsiEqualityPenalty + PenaltyCoeff *  (sum(mu) -1)^2
      return(PsiFunction)
    }
    Psi_x <- function(x1, mu, FIM,  x, w,  answering){

      ## here x is a degenerate design that putt all its mass on x.
      # x1 is one point
      # This function is required to check the equivalence theorem by ploting and also find the Efficiency lower bound
      # mu is a vector of measure. For locally optimal design mu = 1
      # FIM: is the Fisher information matrix
      # x: vector of design points
      # w: vector of design weights
      # answering: the matrix of elements of answering set. Each row is an element.
      # answering set are '\mu'in (3) in ICA paper.
      # return the value of left hand side of (4) in ICA paper as a function of x.
      # if mu = 1  and answering has only one row, then we computing the equivalence theorem left hand side for locally D_optimal design.
      # also can be used for multidimensional models like enzyme kinetic
      if(length(mu) != dim(answering)[1])
        stop("The number of measures is not equal to the number of elements of answering set.")
      if(typeof(FIM) != "closure")
        stop("'FIM' must be of type 'closure'.")
      CardOfRegion <- dim(answering)[2] ## cardinality of region of uncertainty
      n_mu <- dim(answering)[1]
      ## each row is  tr(M^{-1}(\xi, mu_j) %*% I(x, \mu_j))
      # so sum of each row of 'Psi_Point_answering' is \int tr(M^{-1}(\xi, mu) %*% I(x_i, \mu))
      Psi_Point_answering <- matrix(NA, 1,  n_mu)
      i <- 1 ## we need Psi at only point x, so we have only one row
      for(j in 1:n_mu){
        Psi_Point_answering[i,j]  <- mu[j] *
          sum(diag(solve(FIM(x = x, w = w, param = answering[j, ]), tol = .Machine$double.xmin) %*%
                     FIM(x = x1, w = 1, param = answering[j, ])
          ))
      }
      PsiAtEachPoint <- rowSums(Psi_Point_answering) - CardOfRegion
      PsiFunction <- PsiAtEachPoint[1]
      return(PsiFunction)
    }
    Psi_xy <- function(x1, y1, mu, FIM,  x, w,  answering){
      ## WARNINGS: do not change names of 'x1' and 'y1' here unless you check the vectorize in 'PlotPsi_x'
      ## there we have 'Vectorize(FUN = Psi_x, vectorize.args=c("x1", "y1"))'
      ##here x is a degenerate matrix that putt all its mass on x.
      ## this function is used for plotting the equivalence theorem equaltion for model with two independent variables.
      # the function is exactly as psy_x, only with two arument.
      # 'Point' will be handeled by the FIM function of the models itself, see 'common_mulit_dimensional_design.R'

      if(length(mu) != dim(answering)[1])
        stop("The number of measures is not equal to the number of elements of answering set.")
      if(typeof(FIM) != "closure")
        stop("'FIM' must be of type 'closure'.")

      CardOfRegion <- dim(answering)[2] # dimension of the region of uncertainty
      n_mu <- dim(answering)[1]

      Psi_Point_answering <- matrix(NA, 1,  n_mu)
      i <- 1
      for(j in 1:n_mu){
        Psi_Point_answering[i,j]  <-
          mu[j] *
          sum(diag(solve(FIM(x = x, w = w, param = answering[j, ]), tol = .Machine$double.xmin) %*%
                     FIM(x = c(x1, y1), w = 1, param = answering[j, ])))
      }

      PsiAtEachPoint <- round(rowSums(Psi_Point_answering) - CardOfRegion, 5)
      PsiFunction <- PsiAtEachPoint[1]

      return(-PsiFunction)
    }
  }
  if (crt_type == "multiple"){
    temp_multi <- create_multiple_minimax(model = "4pl", FIM = NA, multipars = multipars)
    Psi_x <- temp_multi$PsiMulti_x
    Psi_Mu <- temp_multi$PsiMulti_Mu
    Psi_xy <- NULL
  }

  if (crt_type == "user"){
    Psi_Mu <- function(mu, FIM,  x, w,  answering, PenaltyCoeff){
      # mu is a vector of measure that usually Psi_mu is optimized with respect to.
      # FIM: is the Fisher information matrix as function
      # x: vector of design points
      # w: vector of design weights
      # answering: the matrix of elements of answering set. Each row is an element.
      # answering set are '\mu'in (3) in ICA paper.
      # return the value of left hand side of (4) in ICA paper.
      # if mu = 1  and answering has only one row, then we computing the equivalence theorem left hand side for locally D_optimal design.
      # also can be used for multidimensional models like enzyme kinetic

      if(length(mu) != dim(answering)[1])
        stop("The number of measures is not equal to the number of elements of answering set.")
      if(typeof(FIM) != "closure")
        stop("'FIM' must be of type 'closure'.")
      CardOfRegion <- dim(answering)[2] ##cardinal of region of uncertainty
      n_mu <- dim(answering)[1] ## number of mu, measures

      npred <- length(x)/length(w) # number of independent variables.
      one_point_mat <- matrix(x, length(w), npred)

      Psi_Point_answering <- matrix(NA, length(w),  n_mu)
      ## each row is  tr(M^{-1}(\xi, mu_j) %*% I(x_i, \mu_j))
      # so sum of each row of 'Psi_Point_answering' is \int tr(M^{-1}(\xi, mu) %*% I(x_i, \mu))
      for(i in 1:length(w)){
        for(j in 1:n_mu){
          Psi_Point_answering[i, j]  <-
            mu[j] * user_sensfunc(xi_x = as.vector(one_point_mat[i,]), x = x, w = w, param = answering[j, ], fimfunc = FIM)
        }
      }
      #Psi at each points x = poi
      #Psi at each Point
      PsiAtEachPoint <- rowSums(Psi_Point_answering)

      ##now we should creat the penalty function. Psi must be zero at each supoport point
      PsiEqualityPenalty <- sum(PenaltyCoeff*pmax(PsiAtEachPoint, 0)^2 + PenaltyCoeff*pmax(-PsiAtEachPoint, 0)^2)
      PsiFunction <-  PsiEqualityPenalty + PenaltyCoeff *  (sum(mu) -1)^2
      return(PsiFunction)
    }
    Psi_x <- function(x1, mu, FIM,  x, w,  answering){
      if(length(mu) != dim(answering)[1])
        stop("The number of measures is not equal to the number of elements of answering set.")
      if(typeof(FIM) != "closure")
        stop("'FIM' must be of type 'closure'.")

      n_mu <- dim(answering)[1]
      ## each row is  tr(M^{-1}(\xi, mu_j) %*% I(x, \mu_j))
      # so sum of each row of 'Psi_Point_answering' is \int tr(M^{-1}(\xi, mu) %*% I(x_i, \mu))
      Psi_Point_answering <- matrix(NA, 1,  n_mu)
      i <- 1 ## we need Psi at only point x, so we have only one row
      for(j in 1:n_mu){
        Psi_Point_answering[i,j]  <- mu[j] *
          user_sensfunc(xi_x = x1, x = x, w = w, param = answering[j, ], fimfunc = FIM)
      }

      PsiAtEachPoint <- rowSums(Psi_Point_answering)
      PsiFunction <- PsiAtEachPoint[1]
      return(PsiFunction)
    }
    Psi_xy <- function(x1, y1, mu, FIM,  x, w,  answering){
      ## WARNINGS: do not change names of 'x1' and 'y1' here unless you check the vectorize in 'PlotPsi_x'
      ## there we have 'Vectorize(FUN = Psi_x, vectorize.args=c("x1", "y1"))'
      ##here x is a degenerate matrix that putt all its mass on x.
      ## this function is used for plotting the equivalence theorem equaltion for model with two independent variables.
      # the function is exactly as psy_x, only with two arument.
      # 'Point' will be handeled by the FIM function of the models itself, see 'common_mulit_dimensional_design.R'
      if(length(mu) != dim(answering)[1])
        stop("The number of measures is not equal to the number of elements of answering set.")
      if(typeof(FIM) != "closure")
        stop("'FIM' must be of type 'closure'.")

      n_mu <- dim(answering)[1]

      Psi_Point_answering <- matrix(NA, 1,  n_mu)
      i <- 1
      for(j in 1:n_mu){
        Psi_Point_answering[i,j]  <-
          mu[j] *
          user_sensfunc(xi_x = c(x1, y1), x = x, w = w, param = answering[j, ], fimfunc = FIM)
      }

      PsiAtEachPoint <- rowSums(Psi_Point_answering)
      PsiFunction <- PsiAtEachPoint[1]

      return(-PsiFunction)
    }
  }
  if (crt_type == "DPA" || crt_type == "DPAM"){

    Psi_x <- function(x1, mu, FIM,  x, w,  answering){

      if(length(mu) != dim(answering)[1])
        stop("The number of measures is not equal to the number of elements of answering set.")
      if(typeof(FIM) != "closure")
        stop("'FIM' must be of type 'closure'.")

      n_mu <- dim(answering)[1]

      Psi_Point_answering <- matrix(NA, 1,  n_mu)
      i <- 1
      for(j in 1:n_mu){
        FIM_x <- FIM(x = x, w = w, par = answering[j, ])
        FIM_x1 <- FIM(x = x1, w = 1, par = answering[j, ])
        if (crt_type == "DPA")
          if (compound$alpha != 0)
            Psi_Point_answering[i,j] <- mu[j] * compound$alpha/compound$npar * sum(diag(solve(FIM_x) %*% FIM_x1)) +
          (1-compound$alpha) * (compound$prob(x1, answering[j, ])- sum(w * compound$prob(x, answering[j, ])))/sum(w * compound$prob(x, answering[j, ])) else
            Psi_Point_answering[i,j] <- mu[j] * (compound$prob(x1, answering[j, ])- sum(w * compound$prob(x, answering[j, ])))/sum(w * compound$prob(x, answering[j, ]))
        if (crt_type == "DPM")
          Psi_Point_answering[i,j] <- mu[j] * compound$alpha/compound$npar * sum(diag(solve(FIM_x) %*% FIM_x1)) +
          (1-compound$alpha) * (compound$prob(x1, answering[, j])- min(compound$prob(x, answering[, j])))/min(compound$prob(x, answering[j, ]))
      }
      PsiAtEachPoint <- round(rowSums(Psi_Point_answering) - compound$alpha, 5)
      PsiFunction <- PsiAtEachPoint[1]
      return(PsiFunction)
    }

    Psi_xy <- function(x1,y1,  mu, FIM,  x, w,  answering){

      if(length(mu) != dim(answering)[1])
        stop("The number of measures is not equal to the number of elements of answering set.")
      if(typeof(FIM) != "closure")
        stop("'FIM' must be of type 'closure'.")

      n_mu <- dim(answering)[1]

      Psi_Point_answering <- matrix(NA, 1,  n_mu)
      i <- 1
      for(j in 1:n_mu){
        FIM_x <- FIM(x = x, w = w, par = answering[j, ])
        FIM_x1 <- FIM(x = c(x1, y1), w = 1, par = answering[j, ])
        if (crt_type == "DPA")
          if (compound$alpha != 0)
            Psi_Point_answering[i,j] <- mu[j] * compound$alpha/compound$npar * sum(diag(solve(FIM_x) %*% FIM_x1)) + (1-compound$alpha) * (compound$prob(c(x1, y1), answering[j, ])- sum(w * compound$prob(x, answering[j, ])))/sum(w * compound$prob(x, answering[j, ])) else
              Psi_Point_answering[i,j] <- mu[j] * (compound$prob(c(x1, y1), answering[j, ])- sum(w * compound$prob(x, answering[j, ])))/sum(w * compound$prob(x, answering[j, ]))

        if (crt_type == "DPM")
          Psi_Point_answering[i,j] <- mu[j] * compound$alpha/compound$npar * sum(diag(solve(FIM_x) %*% FIM_x1)) +
          (1-compound$alpha) * (compound$prob(c(x1, y1), answering[, j])- min(compound$prob(x, answering[, j])))/min(compound$prob(x, answering[j, ]))
      }
      PsiAtEachPoint <- round(rowSums(Psi_Point_answering) - compound$alpha, 5)
      PsiFunction <- PsiAtEachPoint[1]
      return(-PsiFunction)
    }
    ## add later if you decide to find minimax compound optimal design
    Psi_Mu <- NULL
  }
  Psi_x_minus_minimax <- function(x1, mu, FIM,  x, w,  answering){
    -Psi_x(x1 = x1, mu = mu, FIM = FIM,  x = x, w = w,  answering = answering)
  }
  return(list(Psi_x = Psi_x, Psi_Mu = Psi_Mu, Psi_xy = Psi_xy, Psi_x_minus_minimax = Psi_x_minus_minimax))
}
#############################################################################################################*
#############################################################################################################*
#' @importFrom grDevices  rainbow
PlotPsi_x <- function(x, w, lower, upper, Psi_x, FIM, answering, mu, plot_3d  = "lattice"){
  # plot the equivalanece theorem, for both one and two dimensional.
  # x: vector of the point of the design
  # w: vector of the weights of the design
  # lower: lower bound
  # psi_x psi_x as a function of x. in multiobjective optimal design ('multi_ICA'), 'PsiMulti_x' is 'given as Psi_x'
  # FIM: fisher information matrix function
  # AsnweringSet: answering set matrix, each row is one element of the A(\xi) or N(\xi)
  # vector of the measure. For locally optimal design mu = 1.
  # plot_3d: which package should be used to plot the 3d plots, 'rgl' or 'lattice'


  if(length(lower) == 1){
    xPlot <- seq(lower,upper, length.out = 1000)

    PsiPlot<- sapply(1:length(xPlot), FUN = function(j) Psi_x(x1 = xPlot[j], mu = mu, FIM = FIM, x = x,
                                                              w = w, answering = answering))
    plot(xPlot, PsiPlot, type = "l",
         col = "blue",   xlab = "Design Interval",
         ylab = "c", xaxt = "n")
    abline(h = 0, v = c(x) ,col = "grey", lty = 3)

    Point_y<- sapply(1:length(x), FUN = function(j) Psi_x(x1 = x[j], mu = mu, FIM = FIM, x = x,
                                                          w = w, answering = answering))

    points(x = x,  y = Point_y, col = "red" ,pch = 16, cex = 1)

    axis(side = 1, at = c(lower, x, upper, 0),
         labels = round(c(lower, x, upper, 0), 4))
    # varlist$Psi_x_minus_minimax(x1 = 2, mu = c(.5006285, .4993715),
    #                             answering = matrix(c(-2, 2, 1.5, 1.5), nrow = 2),
    #                             FIM = varlist$fimfunc_sens,
    #                             x = c(-2, .0033, 2), w = c(.274, .452, .274))
  }

  if(length(lower) == 2){
    xlab_3d <- "x1"
    ylab_3d <- "x2"
    len <- 40
    xPlot <- seq(lower[1],upper[1], length.out = len)
    yPlot <- seq(lower[2],upper[2], length.out = len)
    ncol = 100
    color1 <- rev(rainbow(ncol, start = 0/6, end = 4/6))
    ## Psi_x here is Psi_xy
    Psi_xy1 <- Vectorize(FUN = Psi_x, vectorize.args=c("x1", "y1"))
    ###does not work, we can njot vectorize
    z <- -outer(X = xPlot,
                Y =yPlot,
                FUN = Psi_xy1,
                answering = answering,
                mu=mu,
                FIM = FIM,
                w = w,
                x = x)

    if (plot_3d[1] == "lattice"){
      if (requireNamespace("lattice", quietly = TRUE)){
        wireframe_dat <-expand.grid(x = xPlot, y = yPlot)
        wireframe_dat$z <- as.vector(z)
        p1 <- lattice::wireframe( z ~ x * y, data = wireframe_dat,
                                  col.regions=  color1,
                                  drap = TRUE,
                                  scales = list(arrows = FALSE),
                                  shade = FALSE,
                                  xlab = xlab_3d, ylab = ylab_3d,
                                  pretty = TRUE,
                                  zlab = "c",
                                  #zlab = expression(c(x, xi, mu)),
                                  #main = "Optimality Verification Plot",
                                  colorkey = list(col = color1, tick.number = 16))
        #lattice::print.trellis(p1)
        print(p1)
      }else {
        warning("Package 'lattice' is not installed in your system. The 3D derivation plot can not be plotted unless this packages is installed.")
      }
    }
    if(plot_3d[1] == "rgl"){

      if (requireNamespace("rgl", quietly = TRUE)){
        rgl::.check3d()



        zcol  = cut(z, ncol)

        rgl::persp3d(x = xPlot, y = yPlot, z = z, col = color1[zcol],
                     smooth = FALSE,
                     alpha = .8,
                     shininess    = 128,
                     xlab = xlab_3d, ylab = ylab_3d, zlab = "c(x, \\xi, \\mu)")


        Point_mat <- matrix(round(x, 3), ncol = length(lower), nrow = length(x)/length(lower))
        ##adding poin
        Point_y <- sapply(1:dim(Point_mat)[1], FUN = function(j) Psi_x(x1 = Point_mat[j, 1], y1 = Point_mat[j, 2], mu = mu, FIM = FIM, x=x,
                                                                       w = w, answering = answering))
        rgl::points3d(x = Point_mat[, 1],
                      y = Point_mat[, 2],
                      z = Point_y,
                      col = "darkred",
                      size = 11,
                      point_antialias = TRUE)

        text_point <- paste("(", xlab_3d, "=", Point_mat[, 1], ", ", ylab_3d, "=", Point_mat[, 2], "; c=",  Point_y, ")", sep ="")


        rgl::text3d(x = Point_mat[, 1],
                    y = Point_mat[, 2],
                    z = Point_y + .1,
                    texts = text_point,
                    col = "darkred",
                    font = 4)

        rgl::grid3d("z")
      }else{
        warning("Package 'rgl' is not installed in your system. The 3D derivation plot can not be plotted unless this packages is installed.")
      }
    }

    ###contour plot
    xyz <- expand.grid(x = xPlot, y = yPlot)
    xyz$z <- as.vector(z)

    if (requireNamespace("lattice", quietly = TRUE)){
      p2 <- lattice::contourplot( z ~ x * y, data = xyz,
                                  xlab = xlab_3d,
                                  ylab = ylab_3d,
                                  #main = "Optimality Verification Contour Plot",
                                  region = TRUE,
                                  col.regions = color1,
                                  colorkey = list(col = color1, tick.number = 16),
                                  cuts = 13)
      #lattice::print.trellis(p2)
      print(p2)
    }else{
      warning("Package 'lattice' is not installed in your system. The countor plot can not be plotted.")
    }
  }
}




#############################################################################################################*
#############################################################################################################*
PlotELB <- function(Iter ,
                    ELB,
                    plot_main = TRUE,
                    ...){
  # Iter: the iterations of ELB
  # ELB: a vector of Efficiency lower bounds

  cex.main <- .9
  prec <- 8

  if(plot_main)
    main1 <- paste("D-Efficiency Lower Bound (ELB): ", round(ELB[length(ELB)], prec), paste = "") else
      main1 = NULL

  plot(x = Iter, y = ELB ,
       xlim = c(Iter[1], Iter[length(Iter)]), xlab = "Iteration", ylab = "D-Efficiency Lower Bound (ELB)", type = "s",
       main = main1,
       cex.main = cex.main,
       xaxt = "n",...)


  if(Iter[length(Iter)] < 5)
    axis(side = 1, at = c(Iter),
         labels = c(Iter)) else{
           pretty_plot <- pretty(Iter)
           axis(side = 1, at = pretty_plot,
                labels = pretty_plot)
         }
}
###############################################################################################################*
###############################################################################################################*
#' @importFrom nloptr nloptr
minimax_inner <- function(formula,
                          predvars, parvars,
                          family = "gaussian",
                          lx,
                          ux,
                          lp,
                          up,
                          iter,
                          k,
                          fimfunc = NULL,
                          ICA.control =  list(),
                          sens.control = list(),
                          sens.minimax.control = list(),
                          crt.minimax.control = list(),
                          #locally.control = list(...),
                          type = c("minimax", "locally", "standardized", "robust"),
                          initial = NULL,
                          localdes = NULL,
                          npar,
                          robpars = list(),
                          crt_type = c("D", "multiple", "DPA"),
                          multipars = list(),
                          plot_3d = c("lattice", "rgl"),
                          compound = list(prob = NULL, alpha = NULL, npar = NULL),
                          only_w_varlist = list(x = NULL),
                          user_crtfunc = NULL,
                          user_sensfunc = NULL,
                          ...) {
  time1 <- proc.time()
  #   param_set: a matrix that denotes the fixed values for the parameters and is required when inner_space = "discrete". Each row of the matrix is the values of the components of the parameters,
  #    The number of columns should be equal to length(lp).
  if (!is.null(only_w_varlist$x))
    is.only.w <- TRUE else
      is.only.w <- FALSE
    # ################################################*
    # ### do not change the seed
    # if (exists(".Random.seed")) {
    #   OldSeed <- get(".Random.seed", envir = .GlobalEnv)
    #   on.exit(assign(".Random.seed", OldSeed, envir = .GlobalEnv))
    # }
    # ################################################*
    if (is.null(crt.minimax.control$inner_space))
      stop("BUG: set the 'crt.minimax.control$inner_space' in the outer function")
    if (!crt.minimax.control$inner_space %in% c("discrete", "continuous", "robust_set", "locally"))
      stop("BUG: 'inner_space must be either 'discrete' or 'continuous' or 'robust_set' or 'locally'")
    ########################################################*
    #### dealing with FIM and formula besides checking some common argument
    if (type == "minimax" & crt.minimax.control$inner_space == "discrete" & crt_type != "user")
      paramvectorized <- TRUE else
        paramvectorized <- FALSE
      funcs_formula <- check_common_args(fimfunc = fimfunc, formula = formula,
                                         predvars = predvars, parvars = parvars,
                                         family = family, lx =lx, ux = ux, iter = iter, k = k,
                                         paramvectorized = paramvectorized,
                                         prior = NULL, x = only_w_varlist$x,
                                         user_crtfunc = user_crtfunc, user_sensfunc = user_sensfunc)
      if(missing(formula)){
        # to handle ...
        #crt_type == "user"
        if (crt.minimax.control$inner_space == "continuous" || crt.minimax.control$inner_space == "robust_set" || crt.minimax.control$inner_space == "locally" || crt_type == "user")
          fimfunc2 <- function(x, w, param)
            fimfunc(x = x, w = w, param = param,...)
        if (crt.minimax.control$inner_space == "discrete" & crt_type != "user")
          fimfunc2 <- function(x, w, param)
            list(fimfunc(x = x, w = w, param = param,...))

        fimfunc_sens <- function(x, w, param)
          fimfunc(x = x, w = w, param = param,...) ## it can not be equal to fimfunc2 when inner_space is discrete
      } else{
        if (type != "robust")
          if (length(lp) != length(parvars))
            stop("length of 'lp' is not equal to the length of 'parvars'")
        # fim_localdes <- funcs_formula$fimfunc_formula
        fimfunc2 <- funcs_formula$fimfunc_formula ## can be vectorized with respect to parameters!
        fimfunc_sens <- funcs_formula$fimfunc_sens_formula
      }

      if (type[1] == "standardized"){
        if (!missing(formula))
          localdes_par <- create_localdes(parvars = parvars, localdes = localdes) else
            localdes_par <- localdes
          #localdes_par <- function(param)
          #  localdes(param = param, fimfunc = fim_localdes)

      }else
        localdes_par <- NULL

      if (type == "minimax" || type == "standardized" || type == "locally" )
        checkminimax <- check_minimax_args(lp = lp, up = up, type = type, localdes = localdes, localdes_par = localdes_par, parvars = parvars, fimfunc = fimfunc, crt.minimax.control = crt.minimax.control)

      #######################################################*
      ICA.control <- do.call("ICA.control", ICA.control)
      ICA.control <- add_fixed_ICA.control(ICA.control.list = ICA.control)

      sens.minimax.control <- do.call("sens.minimax.control", sens.minimax.control)
      sens.control <- do.call("sens.control", sens.control)
      #############################################################################*
      ## if only one point design was requested, then the weight can only be one
      if (k == 1)
        ICA.control$equal_weight <- TRUE
      if (length(lx) != 1 && ICA.control$sym)
        stop("currently symetric property only can be applied to models with one variable")
      #############################################################################*


      # global variables needed for tin the creation of crfunc
      npred <- length(lx)
      if (is.null(npar))
        stop("BUG: please provide 'npar'")
      #   npar <- length(lp)
      #if (type == "minimax" || type == "standardized" || type == "locally" ){

      # if (crt_type == "user" & type == "minimax" & crt.minimax.control$inner_space == "discrete")
      #   user_FIM_paramvectorized <- create_FIM_paramvectorizes(parvars = parvars, FIM = fimfunc2)

      temp1 <- create_criterion_minimax(FIM = fimfunc2, type = type[1],
                                        localdes = localdes_par, lp = lp,
                                        up = up, npar = npar, robpars = robpars,
                                        crt_type = crt_type[1], multipars = multipars,
                                        compound = compound, is.only.w = is.only.w,
                                        only_w_varlist = only_w_varlist,
                                        user_crtfunc2 = funcs_formula$user_crtfunc)

      temp2 <- create_criterion_minimax(FIM = fimfunc_sens, type = type[1],
                                        localdes = localdes_par, lp = lp,
                                        up = up, npar = npar, robpars = robpars,
                                        crt_type = crt_type[1], multipars = multipars,
                                        compound = compound, is.only.w = is.only.w,
                                        only_w_varlist = only_w_varlist, user_crtfunc2 = funcs_formula$user_crtfunc)

      #}
      #if (type == "robust")
      #  create_criterion_ave(FIM, prob, parset)


      #############################################################################*
      ## return ld and ud, the lower and upper bound of the design weighs and points
      temp3 <- return_ld_ud (sym = ICA.control$sym, equal_weight = ICA.control$equal_weight, k = k, npred = npred, lx = lx, ux = ux, is.only.w  =  is.only.w)
      initial <- check_initial(initial = initial, ld = temp3$ld, ud = temp3$ud)
      #############################################################################*
      # Psi_Mu, Psi_x, Psi_xy
      Psi_funcs <- create_psy_minimax(crt_type = crt_type, multipars = multipars,
                                      compound = compound,
                                      user_sensfunc = funcs_formula$user_sensfunc)
      #############################################################################*

      #############################################################################*
      ## making the arg list
      ## the variables that will be added to control not by user, but by mica
      arg <- list(lx = lx, ux = ux, lp = lp, up = up, k = k, npar = npar,
                  ld = temp3$ld, ud = temp3$ud, type = type[1], #localdes = localdes_par,
                  initial = initial, ICA.control = ICA.control,
                  sens.minimax.control = sens.minimax.control,
                  crt.minimax.control = crt.minimax.control,
                  sens.control = sens.control,
                  FIM = fimfunc2,  FIM_sens = fimfunc_sens,
                  crfunc_sens = temp2$crfunc,
                  crfunc = temp1$crfunc,
                  fixedpar = temp1$fixedpar, fixedpar_id = temp1$fixedpar_id,
                  is_fixed = temp1$is_fixed,
                  lp_nofixed = temp1$lp_nofixed, up_nofixed = temp1$up_nofixed,
                  robpars = robpars, Psi_funcs = Psi_funcs,
                  plot_3d = plot_3d[1],
                  localdes = localdes,
                  compound = compound,
                  is.only.w =  is.only.w,
                  only_w_varlist = only_w_varlist,
                  time_start = time1
                  # SK@03052020
                  ,family=family
                  ,grad=funcs_formula$grad
                  ,mu=funcs_formula$mu
                  ,paramvectorized=paramvectorized
                  )
      #lx_sens = lx_sens,
      #ux_sens = ux_sens)
      if (type == "locally")
        arg$inipars <- lp

      ## updating will be added to arg in iterate
      ## because when the inner_space is discrete crfunc_sen is not vectorized with respect to the parameters

      ### sensitivity function required for cheking the equivalence theorem

      #   arg$Psi_x <- temp_der$Psi_x
      #   arg$Psi_Mu <- temp_der$Psi_Mu
      #   if (length(lx) == 2)
      #     arg$Psi_xy <- temp_der$Psi_xy
      # }
      ##  Psi_x works for both one, two and three dimensional
      ## but Psi_xy is needed for plotting becasue the function should have two arguments
      #############################################################################*

      ICA_object <- list(arg = arg, evol = NULL)
      #class(ICA_object) <- c("list", "minimax") # 06202020@seongho
      class(ICA_object) <- c("minimax")
      out <- update.minimax(object = ICA_object, iter = iter)
      return(out)

      #list(out=out,object=ICA_object,arg=arg)
}

######################################################################################################*
######################################################################################################*
sensminimax_inner <- function (formula,
                               predvars, parvars,
                               family = gaussian(),
                               x, w,
                               lx, ux,
                               lp, up,
                               fimfunc = NULL,
                               sens.control = list(),
                               sens.minimax.control = list(),
                               type = c("minimax", "locally", "standardized", "robust"),
                               localdes = NULL,
                               plot_3d = c("lattice", "rgl"),
                               plot_sens = TRUE,
                               varlist = list(),
                               calledfrom = c("sensfuncs", "iter", "plot"),
                               npar = NULL,
                               crt.minimax.control = list(),
                               calculate_criterion = TRUE,
                               robpars = list(),
                               crt_type = c("D", "multiple", "DPA"),
                               multipars = list(),
                               silent = FALSE,
                               calculate_sens = TRUE,
                               compound = list(prob = NULL, alpha = NULL, npar = NULL),
                               user_crtfunc = NULL,
                               user_sensfunc = NULL,
                               #only_w_varlist = list(x = NULL),
                               ...){
  #calculate_sens is for when you call the function from plot and only
  #want to calculate the criterion!
  time1 <- proc.time()
  if (calledfrom[1]  == "sensfuncs"){
    # if (!is.null( only_w_varlist$x))
    #   is.only.w <- TRUE else
    #     is.only.w <- FALSE

    ## you should create the varlist!! as in iter function
    # ################################################*
    #### minimax version
    if (length(lx) > 2)
      plot_sens <- FALSE
    if (!is.logical(calculate_criterion))
      stop("'calculate_criterion' must be logical")
    if (!is.null(npar))
      if(!is.numeric(npar) || npar <= 0)
        stop("'npar' must be positive numeric")
    if (!is.character(plot_3d[1]))
      stop("'plot_3d' must be a character string")
    if (!(plot_3d[1] %in% c("lattice", "rgl")))
      stop("'plot_3d' must be either 'lattice' or 'rgl'")


    ########################################################*
    #### dealing with FIM and formula besides checking some common argument
    funcs_formula <- check_common_args(fimfunc = fimfunc, formula = formula,
                                       predvars = predvars, parvars = parvars,
                                       family = family, lx =lx, ux = ux,
                                       iter = 1, k = length(w),
                                       paramvectorized = FALSE, prior = NULL,
                                       x =NULL, user_crtfunc = user_crtfunc, user_sensfunc = user_sensfunc)
    if(missing(formula)){
      # to handle ...
      fimfunc2 <- function(x, w, param)
        fimfunc(x = x, w = w, param = param,...)

      fimfunc_sens <- fimfunc2
      fim_localdes <- fimfunc2
    }else{
      if (type != "robust")
        if (length(lp) != length(parvars))
          stop("length of 'lp' is not equal to the length of 'parvars'")
      fimfunc2 <- funcs_formula$fimfunc_formula
      fimfunc_sens <- funcs_formula$fimfunc_sens_formula
      fim_localdes <- funcs_formula$fimfunc_sens_formula
    }

    if (type[1] == "standardized"){
      if (!missing(formula))
        localdes_par <- create_localdes(parvars = parvars, localdes = localdes) else
          localdes_par <- localdes
        # localdes_par <- function(param)
        #   localdes(param = param, fimfunc = fim_localdes)
    }else
      localdes_par <- NULL

    if (is.null(npar))
      npar <- length(lp)

    sens.minimax.control <- do.call("sens.minimax.control", sens.minimax.control)
    sens.control <- do.call("sens.control", sens.control)
    crt.minimax.control <- do.call("crt.minimax.control", crt.minimax.control)
    crt.minimax.control$inner_space <- "continuous"
    if (type == "minimax" || type == "standardized" || type == "locally" )
      checkminimax <- check_minimax_args(lp = lp, up = up, type = type[1], localdes = localdes, localdes_par = localdes_par, parvars = parvars, fimfunc = fimfunc, crt.minimax.control = crt.minimax.control)
    temp2 <- create_criterion_minimax(FIM = fimfunc_sens, type = type[1], localdes = localdes_par,
                                      lp = lp, up = up, npar = npar, robpars = robpars,
                                      crt_type = crt_type[1], multipars = multipars,
                                      compound = compound, is.only.w = FALSE,
                                      only_w_varlist = NULL, user_crtfunc2 = funcs_formula$user_crtfunc)

    ###############################################################################*
    # required for finding the answering set for verification
    # copied from iterate
    #if (length(lp) <= 2)
    optim_starting <- function(fn, lower, upper, w, x, fixedpar, fixedpar_id,  npred){
      # if (!is.only.w)
      #   q1 <- c(x, w) else
      #     q1 <- w
      out <- optim2(fn = fn, lower = lower, upper = upper,
                    n_seg = sens.minimax.control$n_seg,
                    q = c(x, w),
                    fixedpar = fixedpar, fixedpar_id = fixedpar_id,
                    npred= npred)
      minima <- out$minima
      counts <- out$counts
      return(list(minima =minima, counts = counts))
    }
    ######################################################################*
    # Psi_Mu, Psi_x, Psi_xy
    Psi_funcs <- create_psy_minimax(crt_type = crt_type,
                                    multipars = multipars,
                                    compound = compound,
                                    user_sensfunc = funcs_formula$user_sensfunc)
    varlist <-list(fixedpar = temp2$fixedpar, fixedpar_id = temp2$fixedpar_id,
                   npred =  length(lx),
                   crfunc_sens = temp2$crfunc,
                   #Psi_x_minus_minimax = Psi_x_minus_minimax,
                   lp_nofixed = temp2$lp_nofixed, up_nofixed = temp2$up_nofixed,
                   plot_3d = plot_3d[1],
                   optim_starting = optim_starting,
                   fimfunc_sens =  fimfunc_sens,
                   npar = npar,
                   Psi_x_minus_minimax = Psi_funcs$Psi_x_minus_minimax, Psi_x = Psi_funcs$Psi_x,
                   Psi_xy = Psi_funcs$Psi_xy, Psi_Mu = Psi_funcs$Psi_Mu)

  }
  if (calledfrom[1] ==  "iter")
    if (is.null(varlist) || is.null(npar))
      stop("Bug: 'varlist' and 'npar' must be given when you call the sensitivity function from '")

  # if (calledfrom[1] ==  "iter")
  #   if (is.null(varlist) || is.null(npar) || calculate_criterion)
  #     stop("Bug: 'varlist' and 'npar' must be given when you call the sensitivity function from 'iter'\n'calculate_criterion' can not be TRUE because it is overdoing!")

  # finding the answering set, measure and ELB
  # coppied from iterate
  if (calculate_sens){
    if (type[1] != "locally" & type[1] != "robust"){
      if(!silent){
        # if (calledfrom == "sensfuncs")
        #   cat("\n********************************************************************\n")
        cat("Please be patient! It may take long for minimax designs..... \nDecrease the value of 'n_seg' in 'sens.minimax.control' for higher speed. \nCalculating ELB..................................................\n")
      }
      ## finding the answering set, measure
      ######################################################################*
      # find all local minima
      output <- varlist$optim_starting(fn = varlist$crfunc_sens, lower = varlist$lp_nofixed, upper = varlist$up_nofixed, w = w, x = x,
                                       fixedpar = varlist$fixedpar, fixedpar_id = varlist$fixedpar_id, npred = varlist$npred)
      # total_nfeval <- total_nfeval + output$counts
      all_optima <- output$minima
      ########################################################################*

      ########################################################################*
      ### find the nearest elements in all_optima to construct the answering set larter
      near_ind <- find_nearest(values = all_optima[, dim(all_optima)[2]],
                               tol = sens.minimax.control$merge_tol,
                               compare_with_minimum = TRUE)
      ########################################################################*


      ########################################################################*
      ## add the columns of fixedpar to all local minima
      if (any(!is.na(varlist$fixedpar))){
        ## warnings: all_optima also contain the cost values in last column!
        ## create the columns of fixedpar
        fixedpar_col <- matrix(varlist$fixedpar, dim(all_optima)[1], length(varlist$fixedpar), byrow = TRUE)
        ## how many columns the all_optima with the fixed param will have
        all_dim <- 1:(dim(fixedpar_col)[2] + dim(all_optima)[2])
        all_optima <- cbind(all_optima, fixedpar_col)[, order(c(setdiff(all_dim, varlist$fixedpar_id), varlist$fixedpar_id)), drop = FALSE]
        #Modifying the all local_minima
      }
      ########################################################################*

      ########################################################################*
      if (type[1] == "minimax")
        all_optima[, ncol(all_optima)] <- -all_optima[, ncol(all_optima)]
      # construct the answering set and alos answering set cost
      answering <- all_optima[near_ind, , drop = FALSE]
      answering_cost <- answering[, dim(answering)[2], drop = FALSE]
      # if (type[1] == "minimax")
      #   answering_cost <- -answering_cost
      answering <- answering[, -dim(answering)[2], drop = FALSE]
      ########################################################################*

      ########################################################################*
      # find the measure
      mu <- find_measure(npar = varlist$npar, x = x, w = w,
                         answering = answering,
                         FIM = varlist$fimfunc_sens,
                         Psi_Mu = varlist$Psi_Mu)$mu
      ########################################################################*

      ########################################################################*
      ## we save all minima as well
      parchar <- paste("Par", 1:(dim(all_optima)[2]-1), sep = "")
      all_optima <- cbind(all_optima, rep(FALSE, nrow(all_optima)))
      colnames(all_optima) <- c(parchar, "Crtiterion_Value", "Answering_Set")
      if (type[1] == "minimax")
        rownames(all_optima) <- paste("Maximum", 1:dim(all_optima)[1], sep = "")
      if (type[1] == "standardized")
        rownames(all_optima) <- paste("Minimum", 1:dim(all_optima)[1], sep = "")
      #all_optima_cost <- all_optima[, dim(all_optima)[2], drop = FALSE]
      #all_optima <- all_optima[, -dim(all_optima)[2], drop = FALSE]
      all_optima[near_ind, ncol(all_optima)] <- TRUE
      ########################################################################*


    }else{
      if (type[1] == "locally"){
        answering <- matrix(up, nrow = 1) ## we need it for find measure and check. they use answering set
        mu <- 1 # we need it for check and plot
      }
      if (type[1] == "robust"){
        answering <- robpars$parset
        mu <- robpars$prob
      }
      all_optima <- NA
    }
    ##########################################################################*
    # find the maximum of derivative function
    if (is.null(sens.control$x0))
      x0 <- (lx + ux)/2 else
        x0 <- sens.control$x0
      OptimalityCheck <- nloptr::nloptr(x0 = x0,
                                        eval_f = varlist$Psi_x_minus_minimax,
                                        lb = lx,
                                        ub = ux,
                                        opts = sens.control$optslist,
                                        mu = mu,
                                        answering = answering,
                                        x = x,
                                        w = w,
                                        FIM = varlist$fimfunc_sens)
      ##sometimes the optimization can not detect maximum  in the bound, so here we add the cost values on the bound
      vertices_outer <- make_vertices(lower = lx, upper = ux)
      check_vertices <- find_on_points(fn = varlist$Psi_x,
                                       points = vertices_outer,
                                       mu = mu,
                                       answering = answering,
                                       x = x,
                                       w = w,
                                       FIM = varlist$fimfunc_sens)
      check_vertices <- check_vertices$minima[, dim(check_vertices$minima)[2]]
      if (any(is.nan(check_vertices))){
        stop("'NaN produced when computing sensitivity function at ",
             paste(vertices_outer[which(is.nan(check_vertices)), ], collapse = ", "),
             ". Check your design space.")
      }

      ## minus because for optimality check we minimize the minus derivative function
      max_deriv <- c(-OptimalityCheck$objective, check_vertices)
      max_deriv <- max(max_deriv)
      ##########################################################################*
      # D-efficiency lower bound
      #exp(-max_deriv/varlist$npar)
      ELB <- varlist$npar/(varlist$npar + max_deriv)
      if (ELB < 0)
        warning("'ELB' can not be negative.\n1) Provide the number of model parameters by argument 'npar' when any of them are fixed. \n2) Increase the value of 'maxeval' in 'optslist'")
      if ((!silent))
        cat("    Maximum of the sensitivity function is ", max_deriv, "\n    Efficiency lower bound (ELB) is ", ELB, "\n")

      if (length(lx) == 2)
        Psi_plot <- varlist$Psi_xy else
          Psi_plot <- varlist$Psi_x
      if (plot_sens){
        if(!silent)
          cat("Plotting the sensitivity function................................\n")
        PlotPsi_x(lower = lx, upper = ux, Psi_x = Psi_plot, FIM  = varlist$fimfunc_sens,
                  mu = mu, x = x, w = w, plot_3d = varlist$plot_3d, answering = answering)
      }
      object <- list(type = type[1],
                     optima = all_optima,
                     #all_optima_cost = all_optima_cost,
                     #answering = answering,
                     #answering_cost = answering_cost,
                     mu = mu,
                     max_deriv = max_deriv,
                     ELB = ELB,
                     merge_tol = sens.minimax.control$merge_tol)
  }
  if (calculate_criterion){
    if(!silent)
      cat("Evaluating the criterion.........................................\n")
    optim_func <- create_optim_func(type = type[1], lp_nofixed = varlist$lp_nofixed, up_nofixed = varlist$up_nofixed, crt.minimax.control = crt.minimax.control, discrete_set = NULL, robpars = robpars, inipars = lp)
    temp3 <- optim_func(fn = varlist$crfunc,
                        x = x, w = w,
                        lower = varlist$lp_nofixed,
                        upper = varlist$up_nofixed,
                        fixedpar = varlist$fixedpar,
                        fixedpar_id = varlist$fixedpar_id,
                        npred= length(lx))


    MinRowId <- which.min(temp3$minima[, dim(temp3$minima)[2]])
    crtval  <- temp3$minima[MinRowId , dim(temp3$minima)[2]]
    if (type[1] == "minimax")
      crtval <- -crtval
    if (!calculate_sens)
      object <- list()

    object$crtval  <- crtval


  }
  #if (calledfrom[1] == "sensfuncs")

  # if (type == "locally")
  #   object$crtval <- crfunc(param = lp, q = c(x, w), npred = 1)
  # if (type == "minimax")
  #   object$crtval <- max(answering_cost)
  # if(type == "standardized")
  #   object$crtval <- min(answering_cost)
  time2 <- proc.time() - time1
  if(!silent){
    if (calculate_criterion & calculate_sens)
      cat("    Criterion value is",  object$crtval, "\nVerification is done in",time2[3], "seconds!", "\nAdjust the control parameters for higher speed.\n")
    if (calculate_criterion & !calculate_sens)
      cat("    Criterion value is",  object$crtval)
    if (!calculate_criterion & calculate_sens)
      cat("\nVerification is done in",time2[3], "seconds!", "\nAdjust the control parameters for higher speed.\n")
  }
  #, "\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n")
  object$time <- time2[3]
  #class(object) <- c("list", "sensminimax") # 06202020@seongho
  class(object) <- c("sensminimax")
  if (calculate_criterion || calculate_sens) # to avoid error
    return(object)
}

##########################################################################*

