# Model.R

######################################################################
# Statistical model of transcriptomic dynamics                       #
######################################################################

################################################################
### Model definition
#
# For each gene type i:
#
#   log2(expression) = Norm(mu, sigma^2)
#
# where mu = a + b1*D + c*type + b2*C + b3*D*C + b4*E + b5*D*E
# and (a,b1,..,b5,c) are some regression parameters.
#
# (C,D,E,..) represent effects of environmental factors
# (ages of samples, time of the day etc.) and type takes into account
# the effect of gene types. See FIT.pdf for the details.
#
# Various parameters carry the following indices:
#   mu_{i,s}: expression of gene i in the sample s (training data)
#   a_{i}           (s independent intercept)
#   c_{i} * type_{i}
#   b1_{i} * D_{s}
#   b2_{i} * C_{i,s}
#   b3_{i} * D_{s} * C_{i,s}
#   b4_{i,e} * E_{ies}
#   b5_{i,k,e} * D_{s} * E_{i,s,e}
# where
#   i:1..ngenes      (~3x10^4),
#   e:1..length(env) (up to ~6; just 1 in FIT.R)
#   s:1..nsamples    (~500(2008), ~100(2009), size of microarrays)
#

################################################################
# An object representation of the way we construct and refine models
Model$Recipe <- setRefClass(
  'Model$Recipe',
  fields = list(
    # limitation: at the moment, an env must be a single weather factor
    envs      = 'character', # the keyword 'all' or a vec of 'wind', 'temperature' etc.
    init      = 'character', # (gridsearch | manual)
    optim     = 'character', # (none | lm | lasso)+
    fit       = 'character', # (fit.lm | fit.lasso)
    init.data = 'list',      # init=gridsearch => grid.params | init=manual => init.params
    time.step = 'integer', 
    gate.open.min = 'numeric', 
    opts      = 'list'
    )
  )

################################################################
# The linear model for transcriptome in a field environment
#
# Note:
# - In FIT2, params were a 7*2 matrix (2: init, optim),
#   but perhaps that was a bad idea..
# - It appears to be simpler to construct instances of models at the end of
#   each phase (init, optim1, optim2, .., fit), and restrict params to be a 7-vector
# - If needed, user shall save the entire model at the end of each phase.
Model$LogNormal <- setRefClass(
  Class = 'LogNormal',

  fields = list(
    gene    = 'character',      # gene name
    env     = 'character',      # environmental factor (restricted to only one factor)
    time.step = 'integer',      # sampling interval of meteorological data (integral multiple of data.step in Weather object)
    log     = 'character',      # recipe
    params  = 'numeric',        # parameters (named numeric vector)
    dev     = 'numeric',        # deviance
    coefs   = 'numeric',        # regression coefficients (named numeric vector)
    response.type = 'numeric',  # response to environmental stimuli over the threshold (1) or under the threshold (-1)
    input.mean = 'numeric',     # mean values of inputs
    input.sd = 'numeric'        # standard deviations of inputs
    ),

  methods = list(
    # In the R5 class system, initialize() has to work with zero args
    # to make a class instance to be copy()able.
    # (copy() first calls initialize() and then copies field values.)
    initialize = function(gene = character(0), env = character(0), time.step = 1,
                          log = character(0), params = NULL, dev = -1, coefs = NULL,
                          response.type = 1, input.mean = NULL, input.sd = NULL) {
      if (is.null(names(params)) || !all(names(params) == Model$param.names(env))
          || (dev >=0
              && (is.null(names(coefs)) || !all(names(coefs) == Model$coef.names(env)))))
        stop("Protocol bug: invalid names/ordering for params.")
      gene      <<- gene
      env       <<- env
      time.step <<- time.step
      log       <<- log
      params    <<- params
      dev       <<- dev
      coefs     <<- coefs
      response.type <<- sign(response.type)
      if(is.null(input.mean)){
        input.mean <<- rep(0, 1+length(env))
        names(input.mean) <<- c('age', env)
      }else{
        input.mean <<- input.mean
      }
      if(is.null(input.sd)){
        input.sd <<- rep(1, 1+length(env))
        names(input.sd) <<- c('age', env)
      }else{
        input.sd <<- input.sd
        input.sd[input.sd==0] <<- 1
      }
    },

    # normalize inputs
    normalize.attribute = function(attribute.data){
      Model$normalize.attribute(attribute.data, input.mean['age'], input.sd['age'])
    },
    normalize.weather = function(weather.data){
      Model$normalize.weather(weather.data, env, response.type, input.mean[env], input.sd[env])
    },

    # predict gene expression using attributes and meteorological data
    predict = function(attribute.data, weather.data, data.step) {
      # normally, mu = (1,D,type,C,E,..)*(optimized regression coefs)
      # but inputs are NA when env == 'unexpressed' and
      # we return the value of the first coef (intercept) as mu
      mu <- if (env == 'unexpressed')
        coefs[['intercept']]
      else {
        input <- inputVars(params, env,
                           normalize.attribute(attribute.data),
                           normalize.weather(weather.data),
                           data.step, time.step)
        as.vector(input %*% coefs)
      }
    },

    deviance = function(exprs, attribute.data, weather.data, data.step) {
      diff <- exprs[, gene] - predict(attribute.data, weather.data, data.step)
      sum(diff ^ 2)
    }
    )
  )

# Wrapper for Model$LogNormal$new()
# - We provide this in case that we want optional args for the ctor.
#   (It is difficult to have non-optional args in R5 ctors.)
# args
# - see the comments to Model$LogNormal$fields
Model$new <- function(gene, env, time.step, log, params, dev, coefs, response.type, input.mean, input.sd) {
  Model$LogNormal(gene, env, time.step, log, params, dev, coefs, response.type, input.mean, input.sd)
}

# DOC:
# - caller of Model$new() must order params in this order:
Model$param.names <- function(env) {
  vapply(env, function(e) {
    c(sprintf('env.%s.period',     e),
      sprintf('env.%s.amplitude',  e),
      sprintf('env.%s.threshold',  e),
      sprintf('gate.%s.phase',     e),
      sprintf('gate.%s.amplitude', e),
      sprintf('gate.%s.threshold', e))
  }, rep("", 6))
}

# DOC:
# - caller of Model$new() must order coefs in this order:
Model$coef.names <- function(env) {
  c('intercept',                      # 1
    'coef.age',                       # D
    'coef.genotype',                  # N
    'coef.clock.cos',                 # C
    'coef.clock.sin',                 # C
    'coef.ageClock.cos',              # D*C
    'coef.ageClock.sin',              # D*C
    vapply(env, function(e) {         # E(e) and D*E(e)
      c(sprintf('coef.env.%s', e),
        sprintf('coef.ageEnv.%s', e))
    }, rep("", 2)))
}

# Normalizer of attribute:
Model$normalize.attribute <- function(attribute.data, age.mean, age.sd){
  attribute.data[, 'age'] <- (attribute.data[, 'age'] - age.mean) / age.sd
  attribute.data
}

# Normalizer of weather:
Model$normalize.weather <- function(weather.data, env, response.type, weather.mean, weather.sd){
  for(e in env){
    weather.data[, e] <- response.type * (weather.data[, e] - weather.mean[e]) / weather.sd[e]
  }
  weather.data
}

# Upper and lower bounds of parameters:
Model$param.bounds <- function(env, attribute.data, weather.data, data.step, gate.threshold.max=1){
  period.max <- min(attribute.data$times.pickup) -1
  env.n <- length(env)

  threshold.bound <- sapply(env, function(e)
    c(min(weather.data[, env]), max(weather.data[, env]))
  )
  colnames(threshold.bound) <- sprintf('env.%s.threshold', env)

  cbind(
    matrix(c(1, period.max), 2, env.n, dimnames = list(NULL, sprintf('env.%s.period', env))),
    threshold.bound,
    matrix(c(0, 1440), 2, env.n, dimnames = list(NULL, sprintf('gate.%s.phase', env))),
    matrix(c(-1, gate.threshold.max), 2, env.n, dimnames = list(NULL, sprintf('gate.%s.threshold', env)))
  )
}

# convert period to log-scale
Model$period.log <- function(params, env){
  params[sprintf('env.%s.period', env)] <- log(params[sprintf('env.%s.period', env)])

  params
}

Model$period.log.inv <- function(params, env){
  params[sprintf('env.%s.period', env)] <- exp(params[sprintf('env.%s.period', env)])

  params
}

# degree of freedom
Model$coef.deg.freedom <- function(env){
  c(rep(1, 7), rep(7, 2*length(env)))
}
