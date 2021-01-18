Posterior_HaldDP = R6::R6Class(
  "Posterior_HaldDP",
  public = list(
    q = NULL,
    s = NULL,
    alpha = NULL,
    r = NULL,
    lambda_i = NULL,
    xi = NULL,
    xi_prop = NULL,
    iters = 0,
    maxIters = 0,
    chunkSize = 2000,
    initialize = function(n_iter,
                          namesSources,
                          namesTimes,
                          namesLocations,
                          namesTypes) {
      self$maxIters = n_iter
      self$alpha = to.tensor(
        NA,
        dim = list(
          'source' = namesSources,
          'time' = namesTimes,
          'location' = namesLocations,
          'iter' = 1:n_iter
        )
      )
      class(self$alpha) = append(class(self$alpha), 'array')
      self$q = to.tensor(NA,
                         dim = list('type' = namesTypes,
                                    'iter' = 1:n_iter))
      class(self$q) = append(class(self$q), 'array')
      self$s = to.tensor(NA,
                         dim = list('type' = namesTypes,
                                    'iter' = 1:n_iter))
      class(self$s) = append(class(self$s), 'array')
      self$r = to.tensor(
        NA,
        dim = list(
          'type' = namesTypes,
          'source' = namesSources,
          'time' = namesTimes,
          'iter' = 1:n_iter
        )
      )
      class(self$r) = append(class(self$r), 'array')
    },
    sample = function(model) {
      # Extend the posterior if we need to
      self$iters = self$iters + 1
      if (self$iters > self$maxIters)
        self$extend(self$chunkSize)

      # Save chain state
      for (time in 1:dim(self$alpha)[2]) {
        for (location in 1:dim(self$alpha)[3]) {
          self$alpha[, time, location, self$iters] <-
            model$alphaNodes[[time]][[location]]$getData()
        }
        for (source in 1:dim(self$alpha)[1]) {
          # TODO: change to an lapply?
          self$r[, source, time, self$iters] <-
            model$rNodes[[time]][[source]]$getData()
        }
      }
      self$q[, self$iters] <-
        model$qNodes$getData()
      self$s[, self$iters] <- model$qNodes$s
    },
    extend = function(iters) {
      # Extends the posterior file by iters iterations
      rnames = list(1:(self$maxIters + iters))
      self$alpha = arrayextend(self$alpha,
                               along = 4,
                               size = iters,
                               newdimnames = rnames)
      self$q     = arrayextend(self$q,
                               along = 2,
                               size = iters,
                               newdimnames = rnames)
      self$s     = arrayextend(self$s,
                               along = 2,
                               size = iters,
                               newdimnames = rnames)
      self$r     = arrayextend(self$r,
                               along = 4,
                               size = iters,
                               newdimnames = rnames)
      self$maxIters = self$maxIters + iters
    },
    calc_lambda_i = function(n_iter,
                             nTimes,
                             nLocations,
                             nTypes,
                             namesTimes,
                             namesLocations,
                             namesTypes,
                             namesIters,
                             k) {
      k = as.tensor(k)
      names(k) = c('source', 'time')
      self$lambda_i = self$q * mul.tensor(
        self$r,
        i = 'source',
        self$alpha * k,
        j = 'source',
        by = c('time', 'iter')
      )
      self$lambda_i = reorder.tensor(self$lambda_i, c('type', 'time', 'location', 'iter'))
    },
    calc_xi = function(n_iter,
                             nSources,
                             nTimes,
                             nLocations,
                             namesSources,
                             namesTimes,
                             namesLocations,
                             namesIters,
                             k) {
      # tensorA implementation
      k = as.tensor(k)
      names(k) = c('source', 'time')
      self$xi = self$alpha * mul.tensor(self$r, 'type', self$q, 'type', by =
                                                'iter') * k
    },
    calc_xi_prop = function(n_iter,
                                  nSources,
                                  nTimes,
                                  nLocations,
                                  namesSources,
                                  namesTimes,
                                  namesLocations,
                                  namesIters,
                                  k) {
      if (is.null(self$xi))
        self$calc_xi(
          n_iter,
          nSources,
          nTimes,
          nLocations,
          namesSources,
          namesTimes,
          namesLocations,
          namesIters,
          k
        )
      self$xi_prop = to.tensor(apply(self$xi, c('time', 'location', 'iter'), function(x)
        x / sum(x)))
      names(self$xi_prop)[1] = 'source'
    }
  )
)
