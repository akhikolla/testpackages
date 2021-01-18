# Train.R

################################################################
# Since FieldTranscriptome3 (a predecessor of FIT),
# the recipe for constructing models is a part of the library interface.
#
# Models are constructed from params and coefs.
# Params and coefs can be computed thru the following 3 phases:
#
# - init:  chooses initial values for the model params
# - optim: tunes the params
# - fit:   computes the model coefs using the tuned params
#
# init -> optim (-> optim ..) -> fit => (params,coefs) => model
#
################################################################

################################
# init

# returns a list (length genes.n) of lists (lengths envs.n) of models
# - both lists should be named

Train$init <- function(exprs, weights, attribute.data, weather.data,
                       envs, method, init.data, data.step, time.step) {
  genes <- colnames(exprs)
  if (!is.character(envs)) stop('init.gridsearch(): composite env not supported yet.')

  names(envs) <- envs
  names(genes) <- genes
  models <- if (method == 'manual') {
    Train$init.manual(exprs, weights, attribute.data, weather.data, envs, init.data, genes,
                      data.step, time.step)
  } else if (method == 'gridsearch') {
    Train$init.gridsearch(exprs, weights, attribute.data, weather.data, envs, init.data, genes,
                          data.step, time.step)
  }
}

Train$init.manual <- function(exprs, weights, attribute.data, weather.data, envs, init.data, genes, data.step, time.step) {
  manual.params <- init.data$params
  
  if(is.null(init.data$response.type)){
    manual.reponse.type <- rep(1, length(genes))
    names(manual.response.type) <- genes
  }else if(length(init.data$response.type)==1){
    manual.response.type <- rep(init.data$response.type, length(genes))
    names(manual.response.type) <- genes
  }else{
    manual.response.type <- init.data$response.type
  }
  
  if(is.null(init.data$input.mean)){
    manual.input.mean <- c(age=mean(attribute.data$age, na.rm=T), sapply(envs, function(e) mean(weather.data[[e]], na.rm=T)))
  }else{
    manual.input.mean <- init.data$input.mean
  }
  if(is.null(init.data$input.sd)){
    manual.input.sd <- c(age=sd(attribute.data$age, na.rm=T), sapply(envs, function(e) sd(weather.data[[e]], na.rm=T)))
  }else{
    manual.input.sd <- init.data$input.sd
  }
  
  if (!all(vapply(genes, function(g) g %in% names(manual.params), TRUE)))
    stop('init.manual(): init params missing for some genes')

  cat('# manually setting init params\n')
  models <- lapply(genes, function(g) {
    lapply(envs, function(e) {
      log <- c('manual')
      params <- manual.params[[g]][Model$param.names(e)]
      
      attribute.data.norm <- Model$normalize.attribute(attribute.data, manual.input.mean['age'], manual.input.sd['age'])
      weather.data.norm <- Model$normalize.weather(weather.data, e, manual.response.type[g], 
                                                   manual.input.mean[e], manual.input.sd[e])
      
      fit <- Train$fitLasso(params, e, exprs[, g], weights[, g], attribute.data.norm, weather.data.norm,
                            data.step, time.step)
      coefs <- Train$coefsLasso(fit, e)
      dev <- Train$devLm(params, e, exprs[, g], weights[, g], attribute.data.norm, weather.data.norm,
                         data.step, time.step)
      response.type <- manual.response.type[g]
      input.mean <- manual.input.mean[c('age', e)]
      input.sd <- manual.input.sd[c('age', e)]
      
      Model$new(g, e, time.step, log, params, dev, coefs, response.type, input.mean, input.sd)
    })
  })
}

Train$init.gridsearch <- function(exprs, weights, attribute.data, weather.data, envs, init.data, genes, data.step, time.step) {
  grid.coords <- init.data

  age.mean <- mean(attribute.data$age, na.rm=T)
  age.sd <- sd(attribute.data$age, na.rm=T)
  if(age.sd == 0) age.sd <- 1
  weather.mean <- sapply(envs, function(e) mean(weather.data[[e]], na.rm=T))
  weather.sd <- sapply(envs, function(e) sd(weather.data[[e]], na.rm=T))
  weather.sd[weather.sd==0] <- 1
  input.mean <- c(age=age.mean, weather.mean)
  input.sd <- c(age=age.sd, weather.sd)
  
  grid.coords.posi <- grid.coords
  grid.coords.nega <- grid.coords
  for(env in envs){
    threshold.name <- sprintf('env.%s.threshold', env)
    grid.coords.posi[[threshold.name]] <- (grid.coords.posi[[threshold.name]] - weather.mean[env]) / weather.sd[env]
    grid.coords.nega[[threshold.name]] <- -(grid.coords.nega[[threshold.name]] - weather.mean[env]) / weather.sd[env]
  }

  attribute.data <- Model$normalize.attribute(attribute.data, age.mean, age.sd)
  weather.data.posi <- Model$normalize.weather(weather.data, envs, 1, weather.mean, weather.sd)
  weather.data.nega <- Model$normalize.weather(weather.data, envs, -1, weather.mean, weather.sd)
  
  init.params.and.devs.posi <- Train$initParamsAndDevs(exprs, weights, attribute.data,
                                                      weather.data.posi, envs, grid.coords.posi,
                                                      data.step, time.step)
  init.params.and.devs.nega <- Train$initParamsAndDevs(exprs, weights, attribute.data,
                                                      weather.data.nega, envs, grid.coords.nega,
                                                      data.step, time.step)
  
  models <- lapply(genes, function(g) {
    lapply(envs, function(e) {
      log       <- c('gridsearch')
      dev.posi  <- init.params.and.devs.posi[[e]][['init.devs']][[g]]
      dev.nega  <- init.params.and.devs.nega[[e]][['init.devs']][[g]]
      if(is.na(dev.nega) || dev.posi < dev.nega){
        params <- init.params.and.devs.posi[[e]][['init.params']][, g]
        dev    <- dev.posi
        coefs  <- Train$coefsLm(params, e, exprs[, g], weights[, g], attribute.data, weather.data.posi, data.step, time.step)
        Model$new(g, e, time.step, log, params, dev, coefs, 1, input.mean, input.sd)
      }else{
        params <- init.params.and.devs.nega[[e]][['init.params']][, g]
        dev    <- dev.nega
        coefs  <- Train$coefsLm(params, e, exprs[, g], weights[, g], attribute.data, weather.data.nega, data.step, time.step)        
        Model$new(g, e, time.step, log, params, dev, coefs, -1, input.mean, input.sd)
      }
    })
  })
}

################################
# optim
# - refine models by performing a step in recipe$optim

# list of supported options for each method
Train$optim.default.opts <- list(
  none  = list(),
  lm    = list(maxit=1500, nfolds=-1), # nfolds for lm is simply ignored
  lasso = list(maxit=1000, nfolds=10)
  )

Train$optim <- function(exprs, weights, attribute.data, weather.data, models, method,
                        data.step, time.step, maxit = NULL, nfolds = NULL,
                        min.expressed.rate = 0.01, gate.open.min=0) {
  if      (tolower(method) == 'none')   { cat('# *** (none)\n'); return(models) }
  else if (tolower(method) == 'lm')     { cat('# *** Lm:\n');   method <- 'lm' }
  else if (tolower(method) == 'lasso')  { cat('# *** Lasso:\n'); method <- 'lasso' }
  else stop('# !! optimization method ', method, ' not supported. Aborting..')

  opts <- Train$optim.default.opts[[method]]
  if (!is.null(maxit) && is.numeric(maxit)) opts$maxit <- maxit
  if (!is.null(nfolds) && is.numeric(nfolds)) opts$nfolds <- nfolds

  fn <- if (method == 'lasso') {
          function(params, env, s.expr, s.weight, s.attrs, weather.data, data.step, time.step) {
            fit <- Train$fitLasso(params, env, s.expr, s.weight, s.attrs, weather.data, data.step, time.step, opts$nfolds)
            Train$devLasso(fit)
          }
        } else devLm # call C wrapper directly
  genes        <- colnames(exprs)
  genes.n      <- length(genes)
  best.devs    <- rep(-1, genes.n) # deviances are never negative (note that NA < 1 is NA)
  best.envs    <- character(genes.n)
  best.params  <- vector(mode = 'list', length = genes.n)

  # exprs are logs of (offset + expresseion count)s where offset is a const
  # across a run (normally 1, or 0.1).
  # We here assume that a gene with var(expr) < thres.expr is an unexpressed gene,
  # and set its model as: expr = log(offset) + 0*inputs
  var.exprs <- rep(-1, genes.n)
  is.unexpressed <- logical(genes.n)

  names(genes)       <- genes
  names(best.devs)   <- genes
  names(best.envs)   <- genes
  names(best.params) <- genes
  names(var.exprs) <- genes
  names(is.unexpressed) <- genes

  # upper bound of gate.threshold
  gate.threshold.max <- cos(pi * min(max(gate.open.min, 0), 1440) / 1440)
  
  for (g in genes) {
    cat('# optimizing', g, '\n# | ')
    var.exprs[[g]]  <- var(exprs[, g])
    ## 2015-08-12
    ## v0.06: "unexpressed" genes were defined like this:
    ## if (var.exprs[[g]] < thres.expr) is.unexpressed[[g]] <- TRUE
    ## v0.07:
    expr.valBg <- median(exprs[, g])
    expressed.n.thres <- floor(min.expressed.rate * length(exprs[, g]))
    if (sum (expr.valBg != exprs[, g]) < expressed.n.thres)
      is.unexpressed[[g]] <- TRUE

    if (is.unexpressed[[g]]) {
      cat('unexpressed | => (', 'unexpressed', ',', var.exprs[[g]], ')\n')
    } else {
      envs <- names(models[[g]])
      for (e in envs) {
        cat(e)
        # normalize
        attribute.data.norm <- models[[g]][[e]]$normalize.attribute(attribute.data)
        weather.data.norm <- models[[g]][[e]]$normalize.weather(weather.data)
        
        res <- Train$doOptim(models[[g]][[e]][['params']], fn,
                             e, exprs[, g], weights[, g], attribute.data.norm, weather.data.norm,
                             data.step, time.step, maxit = opts$maxit, gate.threshold.max=gate.threshold.max)
        if (res$value < best.devs[[g]] || best.devs[[g]] < 0) {
          best.envs[[g]]   <- e
          best.devs[[g]]   <- res$value
          best.params[[g]] <- res$par
        }
        # use nested if, as R switch function is awkward..
        s <- res$convergence
        if (s == 0)       cat (' o | ') # normal
        else if (s == 1)  cat (' * | ') # maxit
        else if (s == 10) cat (' # | ') # degenerate
        else cat (' ? | ')
      }
      cat('=> (', best.envs[[g]], ',', best.devs[[g]], ')\n')
    }
  }
  models <- lapply(genes, function(g) {
    if (is.unexpressed[[g]]) {
      e <- 'unexpressed'
      log <- 'unexpressed'
      names.params <- Model$param.names(e)
      params <- rep(NA_real_, length(names.params))
      names(params) <- names.params
      names.coefs <- Model$coef.names(e)
      coefs <- numeric(length(names.coefs))
      coefs[[1]] <- expr.valBg # intercept, others = 0
      names(coefs) <- Model$coef.names(e)
      
      response.type <- 1
      input.mean <- numeric(length(names.coefs))
      input.sd <- rep(1, length(names.coefs))
      
      m <- Model$new(g, e, time.step, log, params, var.exprs[[g]], coefs, response.type, input.mean, input.sd)
    } else {
      e      <- best.envs[[g]]
      cur    <- models[[g]][[e]]
      log    <- c(cur$log, method)
      params <- best.params[[g]]
      dev    <- best.devs[[g]]
      
      fit    <- Train$fitLasso(params, e, exprs[, g], weights[, g], 
                               cur$normalize.attribute(attribute.data), 
                               cur$normalize.weather(weather.data),
                               data.step, time.step)
      coefs  <- Train$coefsLasso(fit, e)

      m <- Model$new(g, e, time.step, log, params, dev, coefs, cur$response.type, cur$input.mean, cur$input.sd)
    }
    # keep models a list (of length genes) of lists (of lengths envs)
    # although envs will be length 1 from the first iteration of the optim phase
    m <- list(m)
    names(m) <- e
    m
  })
}

################################
# fit
Train$fit <- function(exprs, weights, attribute.data, weather.data, models, method,
                      data.step, time.step) {
  if (method != 'fit.lm' && method != 'fit.lasso')
    stop('unknown fitting method: ', method)

  genes <- colnames(exprs)
  names(genes) <- genes
  models <- lapply(genes, function(g) {
    envs <- names(models[[g]])
    names(envs) <- envs
    lapply(envs, function(e) {
      cur    <- models[[g]][[e]]
      log    <- c(cur$log, method)
      params <- cur$params
      if (cur$env == 'unexpressed') {
        # models with e == 'unexpressed' has coefs/dev set
        coefs <- cur$coefs
        dev <- cur$dev
      } else {
        attribute.data.norm <- cur$normalize.attribute(attribute.data)
        weather.data.norm <- cur$normalize.weather(weather.data)
        
        if (method == 'fit.lm') {
          coefs <- Train$coefsLm(params, e, exprs[, g], weights[, g], attribute.data.norm, weather.data.norm, data.step, time.step)
          dev   <- Train$devLm(params, e, exprs[, g], weights[, g], attribute.data.norm, weather.data.norm, data.step, time.step)
        } else if (method == 'fit.lasso') {
          fit    <- Train$fitLasso(params, e, exprs[, g], weights[, g], attribute.data.norm, weather.data.norm, data.step, time.step)
          coefs  <- Train$coefsLasso(fit, e)
          dev    <- Train$devLm(params, e, exprs[, g], weights[, g], attribute.data.norm, weather.data.norm, data.step, time.step)
        } else stop('cannot happen')
      }
      Model$new(g, e, time.step, log, params, dev, coefs, cur$response.type, cur$input.mean, cur$input.sd)
    })
  })
}


################################################################
# some lower level stuff

Train$initParamsAndDevs <- function(exprs, weights, attribute.data, weather.data,
                                    env.factors, grid.coords, data.step, time.step) {
  initParamsAndDevs(exprs, weights, attribute.data, weather.data, env.factors, grid.coords,
                    data.step, time.step)
}

Train$doOptim <- function(params, fn, env, expr, weight, attribute.data, weather.data,
                          data.step, time.step, maxit, gate.threshold.max=1) {
  param.bounds <- Model$param.bounds(env, attribute.data, weather.data, data.step, gate.threshold.max)
  bounds.name <- colnames(param.bounds)
  params[bounds.name] <- pmin(pmax(params[bounds.name], param.bounds[1, bounds.name]), param.bounds[2, bounds.name])
  params.period.log <- Model$period.log(params, env)
  
  fn.bounded <- function(params, env, expr, weight, attribute.data, weather.data, data.step, time.step, param.bounds){
    params.inv <- Model$period.log.inv(params, env)
    if(
      prod(vapply(colnames(param.bounds), function(bound.name){
        param.bounds[1,bound.name] <= params.inv[bound.name] & params.inv[bound.name] < param.bounds[2,bound.name]
      }, TRUE))
    ){
      fn(params.inv, env, expr, weight, attribute.data, weather.data, data.step, time.step)
    }else{
      Inf
    }
  }
  
  res <- stats::optim(params.period.log, fn.bounded, gr = NULL,
        env, expr, weight, attribute.data, weather.data, data.step, time.step, param.bounds, 
        control = list(trace=0, maxit=maxit))
  res$par <- Model$period.log.inv(res$par, env)

  res
}

# coefsLm() and devLm() are light enough so that we do not provide Train$fitLm()
Train$coefsLm <- function(params, env, expr, weight, attribute.data, weather.data, data.step, time.step) {
  coefs <- coefsLm(params, env, expr, weight, attribute.data, weather.data, data.step, time.step)
  names(coefs) <- Model$coef.names(env)
  coefs
}
Train$devLm <- function(params, env, expr, weight, attribute.data, weather.data, data.step, time.step) {
  devLm(params, env, expr, weight, attribute.data, weather.data, data.step, time.step)
}

Train$fitLasso <- function(params, env, expr, weight, attribute.data, weather.data,
                           data.step, time.step) {
  # Lasso (use gglasso)
  # Group Lasso
  
  # weight
  weight.sq <- sqrt(weight)
  expr.weighted <- weight.sq * expr
  inputs <- inputVars(params, env, attribute.data, weather.data, data.step, time.step)[, -1]
  inputs.weighted <- diag(weight.sq) %*% inputs
  
  # prepare 
  group.index <- c(1, 2, 3, 3, 4, 4, 5:(ncol(inputs)-2))

  deg.freedom <- sapply(unique(group.index), 
                        function(gi) sqrt(sum(group.index==gi)) + sum(Model$coef.deg.freedom(env)[-1][group.index==gi]-1))
  
#  inputs.weighted.sd <- apply(inputs.weighted, 2, sd)
#  inputs.weighted.sd[inputs.weighted.sd==0] <- 1
#  inputs.weighted.norm <- apply(inputs.weighted, 2, function(a) a - mean(a)) / matrix(inputs.weighted.sd, nrow(inputs.weighted), ncol(inputs.weighted), byrow=T)
#  beta.ols <- lm(expr.weighted~inputs.weighted.norm)$coefficients[-1]
  beta.ols <- lm(expr~inputs, weights = weight)$coefficients[-1]
  beta.ols[is.na(beta.ols) | beta.ols==0] <- 1e-5
  multiplier <- deg.freedom / sapply(unique(group.index), function(gi) sqrt(sum(beta.ols[group.index==gi]^2))^2)
  
#  fit <- gglasso::gglasso(inputs.weighted.norm, expr.weighted, group=group.index, pf=multiplier)
  fit <- gglasso::gglasso(inputs.weighted, expr.weighted, group=group.index, pf=multiplier)
  coefficients <- fit$beta
  intercepts <- mean(expr) - colMeans(inputs)%*%coefficients
  
  lambda <- fit$lambda

  lambda.n <- length(lambda)
  sample.n <- nrow(inputs)
  errors <- matrix(0, sample.n, lambda.n)
  for(i in 1:lambda.n){
    active.set <- which(coefficients[, i] != 0)
    if(length(active.set)>0){
      inputs.inv <- MASS::ginv(t(inputs.weighted[, active.set])%*%inputs.weighted[, active.set])
      for(j in 1:sample.n){
        chi <- 
          inputs.inv + inputs.inv %*% inputs.weighted[j, active.set] %*% t(inputs.weighted[j, active.set]) %*% inputs.inv / 
          as.numeric(1 - t(inputs.weighted[j, active.set])%*%inputs.inv%*%inputs.weighted[j, active.set])
        errors[j,i] <- (1 + sum(inputs.weighted[j, active.set]%*%t(inputs.weighted[j, active.set]) * chi))^2 * 
          (expr.weighted[j] - inputs.weighted[j,] %*% coefficients[, i] - intercepts[i])^2
      }
    }else{
      errors[, i] <- (expr.weighted - inputs.weighted %*% coefficients[, i] - intercepts[i])^2
    }
  }
  cvm <- colMeans(errors)
  
  coefs <- rep(0, ncol(inputs)+1)
  lambda.idx <- which(lambda==max(lambda[cvm <= min(cvm)+mad(errors[, which.min(cvm)])/sqrt(sample.n)]))
  nzero.idx <- which(coefficients[, lambda.idx]!=0)
  
  if(length(nzero.idx) > 0){
    coefs[c(1,nzero.idx+1)] <- lm(expr~inputs[, nzero.idx], weights = weight)$coefficients
    coefs[is.na(coefs)] <- 0
  }else{
    coefs[1] <- sum(expr.weighted) / sum(weight.sq)
  }
  list(coefs=coefs, cvm=cvm[lambda.idx])
}

# Note: does not name the coefs
Train$coefsLasso <- function(fit, env) {
  coefs <- fit$coefs
  names(coefs) <- Model$coef.names(env)
  coefs
}

Train$devLasso <- function(fit) {
  fit$cvm
}
