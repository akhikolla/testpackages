## --------------- Degree error handling ---------------
handler1 <- function(){
  e <- parent.frame()

  #variable names in parent frame

  vnames <- objects(, envir = e)
  # "sourcing" the parent.frame
  for(n in vnames) assign(n, get(n, e))

  if(!(degree %% 1 == 0))
    stop('Polynomial degree must be integer.')

  aux <- match(c("formula", "data"),
               names(Call), nomatch = 0)

  if (aux[1] == 0) stop("A formula argument is required")
  if (aux[2] == 0) stop("A dataset argument is required")
  e$aux <- aux
}

## --------------- Frailty handling ---------------
handler2 <- function(){
  e <- parent.frame()
  #variable names in parent frame

  vnames <- objects(, envir = e)
  # "sourcing" the parent.frame
  for(n in vnames) assign(n, get(n, e))

  id <- NULL
  if (!is.null(attr(temp$formula, "specials")$frailty)) {
    frailty_idx <- attr(temp$formula, "specials")$frailty
    id <- model.matrix(formula)[, attr(temp$formula, "specials")$frailty]
    dist <- 1 ## gamma
  }
  else if (!is.null(attr(temp$formula, "specials")$frailty.gamma)) {
    frailty_idx <- attr(temp$formula, "specials")$frailty.gamma
    id <- model.matrix(formula)[, attr(temp$formula, "specials")$frailty.gamma]
    dist <- 1 ## gamma
  }
  else if (!is.null(attr(temp$formula, "specials")$frailty.gauss)) {
    frailty_idx <- attr(temp$formula, "specials")$frailty.gauss
    id <- model.matrix(formula)[, attr(temp$formula, "specials")$frailty.gauss]
    dist <- 2 ## gauss
  }
  else if (!is.null(attr(temp$formula, "specials")$frailty.t)) {
    frailty_idx <- attr(temp$formula, "specials")$frailty.t
    id <- model.matrix(formula)[, attr(temp$formula, "specials")$frailty.t]
    dist <- 3 ## t-student
  }
  else{
    dist <- 0
    frailty_idx <- NULL
  }
  e$dist <- dist
  e$id <- id
  e$frailty_idx <- frailty_idx
}

## --------------- Priors handling ---------------
handler3 <- function(){
  e <- parent.frame()
  #variable names in parent frame

  vnames <- objects(, envir = e)
  # "sourcing" the parent.frame
  for(n in vnames) assign(n, get(n, e))

  betap <- try(lapply(priors$beta, read_prior), silent = T)
  # locationp <- list()
  # scalep <- list()

  # priordist = priordist,
  #             priorpars = priorpars,
  #             priordistbeta = priordistbeta,
  #             priordistlocation = priordistlocation,
  #             priordistscale = priordistscale,
  #             priorparsbeta = priorparsbeta,

  # for(i in 1:length(betap)){
  #   if(is.na(betap[[i]][2])){
  #     locationp[[i]] <- read_prior(priors$location_beta[i])
  #     scalep[[i]] <- read_prior(priors$scale_beta[i])
  #   }
  #   else{
  #     locationp[[i]] <- c("normal", betap[[i]][2], .000001)
  #     scalep[[i]] <- c("normal", betap[[i]][3], .000001)
  #   }
  # }

  e$priordist_beta <- sapply(betap, `[[`, 1)
  e$location_beta <- sapply(betap, `[[`,2)
  e$scale_beta <- sapply(betap, `[[`,3)

  gammap <- try(read_prior(priors$gamma), silent = T)
  # if(is.na(gammap[2])){
  #   hyperp <- read_prior(priors$hyper_gamma)
  # }
  # else{
  #   hyperp <- c("normal", gammap[3], .000001)
  # }

  frailtyp <- try(read_prior(priors$frailty), silent = T)
  if(is.na(frailtyp[2])){
    frailtyp <- c("gamma", ".1", ".1")
  }

  e$priordist <- c(gamma = gammap[1],
                   frailty = frailtyp[1]
                   )

  e$priorpars <- as.numeric(c(hyper1 = gammap[2],
                   hyper2 = gammap[3],
                   frailty1 = frailtyp[2],
                   frailty2 = frailtyp[3]
                   ))
}

## --------------- Extra args error handling ---------------
##  ... arguments directly passed to `rstan::stan`, handles typos
## like "chans=4".
handler4 <- function(){
  ## stan arguments
  e <- parent.frame()
  #variable names in parent frame

  vnames <- objects(, envir = e)
  # "sourcing" the parent.frame
  for(n in vnames) assign(n, get(n, e))

  if (length(stanArgs)) {
    ifelse(approach == 0,
           stanformals <- c(names(formals(rstan::stan)),
                            "seed", "check_data", "sample_file",
                            ~~"algorithm", "verbose", "hessian", "as_vector",
                            "draws", "constrained", "save_iterations",
                            "refresh", "init_alpha", "tol_obj",
                            "tol_rel_obj", "tol_grad", "tol_rel_grad",
                            "tol_param", "history_size"),
           stanformals <- names(formals(rstan::stan))) #legal arg names
    aux <- pmatch(names(stanArgs), stanformals, nomatch = 0)

    if (any(aux == 0))
      stop(gettextf("Argument %s not matched", names(stanArgs)[aux==0]))
  }
}

  ## --------------- Model Frame error handling ---------------
handler5 <- function(){
  e <- parent.frame()
  #variable names in parent frame

  vnames <- objects(, envir = e)
  # "sourcing" the parent.frame
  for(n in vnames) assign(n, get(n, e))

  if (nrow(mf) == 0) stop("Only missing observations")
  if (!inherits(Y, "Surv")) stop("Response must be a survival object")
  if (type!='right' && type!='counting')
    stop(paste("Proportional hazards model doesn't support \"", type,
               "\" survival data", sep=''))
  if (length(attr(Terms, '))variables')) > 2) { # a ~1 formula has length 2
    ytemp <- terms.inner(formula)[1:2]
    xtemp <- terms.inner(formula)[-c(1,2)]
    if (any(!is.na(match(xtemp, ytemp))))
      warning("a variable appears on both the left and right sides of
                the formula")
  }
}


