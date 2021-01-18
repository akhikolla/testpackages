#' #' Calculate Information Criteria
#' #'
#' #' @description Calculate information criteria for \code{D*M} models. This function can be called with
#' #' the resObserved argument or the resDecision and, in case of D*M analyses, also the nondecision analysis.
#' #' If resObserved is not supplied a call to \code{\link{estObserved}} will be made.
#' #'
#' #' @param resObserved an object of class \code{D*M}, specifically output of the function \code{estObserved}.
#' #' @param resDecision an object of class \code{D*M}, specifically output of the function \code{estDstarM}.
#' #' @param resND an object of class \code{D*M}, specifically output of the function \code{estND}.
#' #' @param data The data used to estimate \code{resDecision}.
#' #' @param npar The number of parameters estimated.
#' #'
#' #' @details Calculates several information criteria.
#' #'
#' #' @return a list (S3 object of class 'D*M') that contains:
#' #'
#' #' \item{AIC}{Akaike Information Criterion.}
#' #' \item{AICc}{AIC with a correction for finite sample sizes.}
#' #' \item{BIC}{Bayesian Information Criterion, also known as Schwarz criterion.}
#' #' \item{npar}{Number of parameters.}
#' #' \item{logLik}{The natural log of the likelihood.}
#' #'
#' #' @details Calculate information criteria for D*M models. The results should be interpreted
#' #' carefully! D*M models are \strong{not} estimated by maximizing a likelihood, and
#' #' therefore more complex models do not always necessarily have a higher likelihood than
#' #' less complex models. Also, adding parameters to the decision model might improve the fit
#' #' of the decision model, but can worsen the fit of the nondecision model. As a result the total
#' #' fit decreases when adding unneeded parameters. This is shown in the example. Caution is
#' #' adviced when interpreting likelihood based information criteria.
#' #'
#' #' @examples
#' #'rm(list = ls())
#' #'set.seed(42)
#' #'# Simulate data, fit 3 models, compare models.
#' #'# simulate data with three stimuli of different difficulty.
#' #'# this implies different drift rates across conditions.
#' #'# define a time grid. A more reasonable stepsize is .01; this is just for speed.
#' #'tt = seq(0, 5, .1)
#' #'pars = c(.8, 2, .5, .5, .5, # condition 1
#' #'         .8, 3, .5, .5, .5, # condition 2
#' #'         .8, 4, .5, .5, .5) # condition 3
#' #'pdfND = dbeta(tt, 10, 30)
#' #'# simulate data
#' #'dat = simData(n = 3e3, pars = pars, tt = tt, pdfND = pdfND)
#' #'# define restriction matrices
#' #'restr1 = restr2 = restr3 = matrix(1:5, 5, 3)
#' #'## restr1 allows nothing to differ over conditions
#' #'restr2[2, 2:3] = 6:7 # allow drift rates to differ
#' #'restr3[1:3, 2:3] = 6:11 # allow all decision model parameters to differ
#' #'# fix variance parameters
#' #'fixed = matrix(c('sz1', .5, 'sv1', .5), 2, 2)
#' #'\dontrun{
#' #'# Run D*M analysis
#' #'resD1 = estDstarM(dat = dat, tt = tt, restr = restr1, fixed = fixed, Optim = list(parallelType = 1))
#' #'resD2 = estDstarM(dat = dat, tt = tt, restr = restr2, fixed = fixed, Optim = list(parallelType = 1))
#' #'resD3 = estDstarM(dat = dat, tt = tt, restr = restr3, fixed = fixed, Optim = list(parallelType = 1))
#' #'# Estimate nondecision density - The warnings can be ignored.
#' #'# Quantiles for nondecision densities cannot be calculated due to the big time grid (by .1).
#' #'resND1 = estND(resD1, Optim = list(parallelType = 1))
#' #'resND2 = estND(resD2, Optim = list(parallelType = 1))
#' #'resND3 = estND(resD3, Optim = list(parallelType = 1))
#' #'# Estimate observed densities
#' #'resObs1 = estObserved(resD1, resND1, data = dat)
#' #'resObs2 = estObserved(resD2, resND2, data = dat)
#' #'resObs3 = estObserved(resD3, resND3, data = dat)
#' #'# Compare optimizer fitness
#' #'lstObs = list(resObs1, resObs2, resObs3)
#' #'sapply(lstObs, function(x) x$fit)          # model 2 performs best!
#' #'sapply(lstObs, function(x) x$fit$chisq$sum)# model 2 performs best!
#' #'# model 3 has lower Decision fit due to overfitting the decision model
#' #'sapply(lstObs, function(x) x$fit$Decision)
#' #'
#' #'# However, model 3 has worse nondecision fit compared to model 2 due to overfitting!
#' #'sapply(lstObs, function(x) x$fit$ND)
#' #'# Calculate information criteria
#' #'IC1 = calcIC(resObs1, resD1, data = dat)
#' #'IC2 = calcIC(resObs2, resD2, data = dat)
#' #'IC3 = calcIC(resObs3, resD3, data = dat)
#' #'sapply(list(IC1, IC2, IC3), function(x) x$logLik)
#' #'# Do likelihood ratio tests
#' #'anova(IC1, IC2, IC3) # unorthodox output, model 2 is best?
#' #'which.min(c(IC1$AIC, IC2$AIC, IC3$AIC)) # model 2 is best
#' #'which.min(c(IC1$BIC, IC2$BIC, IC3$BIC)) # model 2 is best
#' #'# We can do the same for traditional analyses
#' #'# likelihood based methods do work for traditional analyses
#' #'# define restriction matrices
#' #'restrC1 = restrC2 = restrC3 = matrix(1:7, 7, 3)
#' #'# restr1 allows nothing to differ over conditions
#' #'restrC2[2, 2:3] = 8:9 # allow drift rates to differ
#' #'restrC3[c(1:2, 4), 2:3] = 8:13 # allow all decision model parameters to differ
#' #'# Do traditional analyses
#' #'resC1 = estDstarM(dat = dat, tt = tt, restr = restrC1, fixed = fixed,
#' #'DstarM = FALSE, Optim = list(parallelType = 1))
#' #'resC2 = estDstarM(dat = dat, tt = tt, restr = restrC2, fixed = fixed,
#' #'DstarM = FALSE, Optim = list(parallelType = 1))
#' #'resC3 = estDstarM(dat = dat, tt = tt, restr = restrC3, fixed = fixed,
#' #'DstarM = FALSE, Optim = list(parallelType = 1))
#' #'# resObs does not need to be calculated now since resC* contains the full model
#' #'# Estimate observed densities
#' #'resObsC1 = estObserved(resC1, data = dat)
#' #'resObsC2 = estObserved(resC2, data = dat)
#' #'resObsC3 = estObserved(resC3, data = dat)
#' #'# Compare optimizer fitness
#' #'lstObsC = list(resObsC1, resObsC2, resObsC3)
#' #'sapply(lstObsC, function(x) x$fit)
#' #'sapply(lstObsC, function(x) x$fit$chisq$sum)
#' #'# Calculate information criteria
#' #'ICC1 = calcIC(resObserved = resObsC1, data = dat)
#' #'ICC2 = calcIC(resDecision = resC2, data = dat) # both input methods are possible
#' #'ICC3 = calcIC(resObserved = resObsC3, data = dat)
#' #'sapply(list(ICC1, ICC2, ICC3), function(x) x$logLik)
#' #'# Do likelihood ratio tests
#' #'# correct model is retrieved by all methods
#' #'anova(ICC1, ICC2, ICC3)
#' #'which.min(c(ICC1$AIC, ICC2$AIC, ICC3$AIC))
#' #'which.min(c(ICC1$BIC, ICC2$BIC, ICC3$BIC))
#' #'}
#'
#'
#'
#' #' @export
#' calcIC = function(resObserved, resDecision, resND, data, npar) {
#'   # error control
#'   if (missing(resObserved)) {
#'     if (!missing(resDecision)) {
#'       if (!is.DstarM(resDecision)) {
#'         stop('Argument resDecision must be of class 'DstarM'', call. = FALSE)
#'       }
#'       if (!resDecision$DstarM) { # traditional analyses
#'         resObserved = estObserved(resDecision = resDecision)
#'       } else {
#'         if (missing(resND)) {
#'           stop('Input was a decision model of a 'D*M' analysis, please also supply a nondecision analysis.', call. = FALSE)
#'         } else {
#'           if (!is.DstarM(resND)) {
#'             stop('Argument resND must be of class 'DstarM'', call. = FALSE)
#'           }
#'           resObserved = estObserved(resDecision = resDecision, resND = resND)
#'         }
#'       }
#'     } else {
#'       stop('Please supply at least argument resObserved or argument resDecision.', call. = FALSE)
#'     }
#'   } else {
#'     if (!is.DstarM(resObserved)) {
#'       stop('Argument resObserved must be of class 'DstarM'', call. = FALSE)
#'     }
#'   }
#'   if (any(!(c('rt', 'condition') %in% names(data)))) {
#'     stop('Data must be a dataframe with a column rt containing observed reaction times and a condition column (named condition)',
#'          call. = FALSE)
#'   }
#'   if (missing(npar)) {
#'     # npar = estimated parameters - fixed parameters
#'     npar = resObserved$npar
#'   }
#'   ncondition = dim(resObserved$obs)[2L] / 2
#'   if (ncondition != length(unique(data$condition))) {
#'     stop(sprintf('Number of conditions in resObserved (%d) does not match number of conditions in data (%d)',
#'                  ncondition, length(unique(data$condition))))
#'   }
#'   # helper matrix
#'   mm = matrix(0, ncondition * 2, ncondition)
#'   mm[1:dim(mm)[1L] + dim(mm)[1L] * rep(1:dim(mm)[2L] - 1, each = 2)] = 1
#'   # round observed rt to grid of tt
#'   tt = resObserved$tt
#'   by = tt[2] - tt[1]
#'   rt = mround(data$rt, by)
#'   tbRt = table(rt)
#'   # condition index
#'   uniqC = unique(data$condition)
#'   idxC = lapply(uniqC, function(find, where) which(where == find), where = data$condition)
#'   tbLst = vector('list', length = length(idxC))
#'   for (i in seq_along(idxC)) {
#'     rt = rt[idxC[[i]]]
#'     tb = table(rt)
#'     tbLst[[i]] = cbind(as.numeric(names(tb)), tb)
#'   }
#'
#'   # sum densities of two condition response pairs to obtain
#'   # density of one condition, normalize it and take the log of it
#'   totDensity = resObserved$obs %*% mm
#'   # cor = apply(totDensity, 2, simpson, x = tt)
#'   # totDensity = totDensity %*% (diag(dim(totDensity)[2L]) / cor)
#'   logDensity = log(totDensity)
#'   logDensity[is.infinite(logDensity)] = NA
#'
#'   # calc log likelihood
#'   logLik = 0
#'   #browser() # check if works properly
#'   for (i in seq_along(idxC)) {
#'     # rounding to deal with small floats (see print(tt, 20))
#'     idxTt = which(round(tt, 10) %in% round(tbLst[[i]][, 1], 10))
#'     logLik = logLik + sum(logDensity[idxTt, i] * tbLst[[i]][, 2], na.rm = TRUE)
#'   }
#'   n = dim(data)[1L]
#'   if (!missing(resDecision)) {
#'     if (resDecision$n != n) {
#'       warning(sprintf('Decision model was estimated on a dataset with a different number of observations (%d) than the number of observations in argument data (%d)?',
#'                       resDecision$n, n))
#'     }
#'   }
#'   # AIC, AICc, BIC, HQC all from wiki
#'   AIC = 2 * npar - 2 * logLik
#'   AICc = AIC + 2 * npar * (npar + 1) / (n - npar - 1)
#'   BIC = log(n) * npar - 2 * logLik
#'   HQC = log(log(n)) * npar - 2 * logLik
#'
#'   res = list(AIC = AIC, AICc = AICc, BIC = BIC,
#'              HQC = HQC, npar = npar, logLik = logLik)
#'   class(res) = 'DstarM'
#'   warning('D*M models are not estimated by maximizing the likelihood! Results of these functions are experimental, be careful when interpreting the results. See ?calcIC for more information.')
#'   return(res)
#' }
#'
#' #' @export
#' anova.DstarM = function(...) {
#'   dots = list(...)
#'   nModel = length(dots)
#'   # helper function for error control
#'   .f2 = function(x) names(x)[[1]] == 'AIC'
#'   errIdx1 = unlist(lapply(dots, is.DstarM))
#'   errIdx2 = unlist(lapply(dots, .f2))
#'   errIdxCombined = which(errIdx1 & errIdx2)
#'   # indices to remove
#'   rmIdx = which(!(seq_along(dots) %in% errIdxCombined))
#'   # gets input argument names
#'   argNames = sapply(substitute(list(...))[-1L], deparse)
#'   rmNames = argNames[rmIdx]
#'   if (length(rmNames)) {
#'     if (length(rmNames) == 1) {
#'       msg = c('Model', 'is', 'does')
#'     } else {
#'       msg = c('Models', 'are', 'do')
#'       if(length(rmIdx) == 2) {
#'         rmNames = paste(rmNames, collapse = ' and ')
#'       } else {
#'         len = length(rmNames)
#'         rmNames = paste0(paste(rmNames[-len], collapse = ', '), ', and ', rmNames[len])
#'       }
#'     }
#'     warning(sprintf('%s %s %s not of class D*M or %s not have AIC as first element and will be omitted.',
#'                     msg[1], rmNames, msg[2], msg[3]), call. = FALSE, immediate. = TRUE)
#'     # remove models
#'     dots = dots[-rmIdx]
#'   }
#'   output = matrix(NA, nrow = length(dots), ncol = 6)
#'   colnames(output) = c('Model', 'logLik', 'nPar', '-2 logLik', 'df', 'p')
#'   # abbreviate names if they are too long
#'   if (any(sapply(argNames, nchar) > 10)) {
#'     argNames = abbreviate(argNames, minlength = 10)
#'   }
#'   # sort list from lowest number of parameters to highest
#'   npar = sapply(dots, `[[`, 'npar')
#'   idxNpars = sort.int(npar, index.return = TRUE)$ix
#'   dots = dots[idxNpars]
#'
#'   logLik = sapply(dots, `[[`, 'logLik')
#'   npar = sapply(dots, `[[`, 'npar')
#'   L.ratio = c(NA, 2*(logLik[2:length(dots)] - logLik[1:(length(dots) - 1)]))
#'   df = c(NA, npar[2:length(dots)] - npar[1:(length(dots) - 1)])
#'   pvals = ifelse(df == 0 | L.ratio < 0, NA, stats::pchisq(L.ratio, df, lower.tail = FALSE))
#'
#'   # collect output
#'   output = cbind(logLik, npar, L.ratio, df, pvals)
#'   dimnames(output) = list(argNames,
#'                           c('logLik', 'npar', 'L.ratio', 'df', 'p.values'))
#'   class(output) = 'DstarM'
#'   return(output)
#' }
