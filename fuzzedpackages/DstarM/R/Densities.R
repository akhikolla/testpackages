#' Calculate model density for a given set of parameters
#'
#' @param t Time grid for density to be calculated on.
#' @param pars Parameter vector where (if \code{DstarM == TRUE}) the first index contains the boundary parameter, the second contains the drift speed, the third contains the relative starting point, the fourth contains a proportion of the maximum size of the variance on the relative starting point, the fifth contains the standard deviation of the drift speed.
#' if \code{DstarM == FALSE} then third index of \code{pars} contains the Ter, the fifth the drift speed, the the sixth contains a proportion of the maximum size of the variance on the relative starting point, the fifth contains the standard deviation of the drift speed, and the seventh contains a proportion of the maximum variance of the Ter.
#' @param boundary For which response option will the density be calculated? Either 'upper' or 'lower'.
#' @param DstarM Logical, see \code{pars}.
#' @param prec Precision with which the density is calculated. Corresponds roughly to the number of decimals accurately calculated.
#' @param ... Other arguments, see \code{\link{dLBA}}
#'
#'
#' @return A numeric vector of length \code{length(t)} containing a density.
#' @details
#' These functions are examples of what \code{fun.density} should look like.
#' \code{Voss.density} is an adaptation of \code{\link{ddiffusion}},
#' \code{LBA.density} is an adaptation of \code{\link{dLBA}}, and
#' \code{wiener.density} is an adaptation of \code{\link{dwiener}}.
#' To improve speed one can remove error handling.
#' Normally error handling is useful, however
#' because differential evolution can result in an incredible number of
#' function evaluations (more than 10.000) it is recommended to omit error handling in custom
#' density functions. \code{estDstarM} will apply some internal error checks
#' (see \code{\link{testFun}}) on the density functions before starting differential
#' evolution. A version of \code{ddifusion} without error handling can be found in
#' the source code (commented out to pass R check). Note that for in \code{Voss.density}
#' if DstarM == FALSE nondecision parameters are implemented manually and might differ
#' from from how they are implemented in other packages. The parameter \code{t0}
#' specifies the mean of a uniform distribution and \code{st0} specifies the relative
#' size of this uniform distribution. To obtain the lower and upper range of the
#' uniform distribution calculate a = t0 - t0*st0, and b = t0 + t0*st0.
#'
#' @examples
#' t = seq(0, .75, .01)
#' V.pars = c(1, 2, .5, .5, .5)
#' L.pars = c(1, .5, 2, 1, 1, 1)
#' W.pars = V.pars[1:3]
#' V1 = Voss.density(t = t, pars = V.pars, boundary = 'upper', DstarM = TRUE)
#' V2 = Voss.density(t = t, pars = V.pars, boundary = 'lower', DstarM = TRUE)
#' L1 = LBA.density(t = t, pars = L.pars, boundary = 'upper', DstarM = TRUE)
#' L2 = LBA.density(t = t, pars = L.pars, boundary = 'lower', DstarM = TRUE)
#' W1 = Wiener.density(t = t, pars = W.pars, boundary = 'upper', DstarM = TRUE)
#' W2 = Wiener.density(t = t, pars = W.pars, boundary = 'lower', DstarM = TRUE)
#' densities = cbind(V1, V2, L1, L2, W1, W2)
#' matplot(t, densities, type = 'b', ylab = 'Density', lty = 1, las = 1, bty = 'n',
#'         col = rep(1:3, each = 2), pch = c(0, 15, 1, 16, 2, 17), cex = .8,
#'         main = 'Model densities')
#' legend('topright', legend = c('Voss', 'LBA', 'RWiener'), lty = 1,
#'        pch = 15:17, col = 1:3, bty = 'n')
#'
#'
#'

#' @importFrom rtdists ddiffusion

# get Voss densities; returns vector of length t
#' @export
Voss.density = function(t, pars, boundary, DstarM = TRUE, prec = 3) {
  if (DstarM) {
    sz = pars[4L] * 2 * pars[1L] * min(c(pars[3L], 1 - pars[3L])) # rescale sz
    pars[3L] = pars[3L] * pars[1L] # rescale z
    dist = abs(ddiffusion(rt = t, response = boundary,
                                   a = pars[1L], v = pars[2L], t0 = 0, z = pars[3],
                                   d = 0, sz = sz, sv = pars[5L], st0 = 0,
                                   precision = prec))
  } else {
    sz = pars[5L] * 2 * pars[1L] * min(c(pars[4L], 1 - pars[4L])) # rescale sz
    pars[4L] = pars[4L] * pars[1L] # rescale z
    st0 = pars[7L]# * 2 * pars[3L] # rescale st0 to % * max size
    t0 <- pars[3L]
    dist = abs(ddiffusion(rt = t, response = boundary,
                                   a = pars[1L], v = pars[2L], t0 = t0, z = pars[4L], # pars[3]
                                   d = 0, sz = sz, sv = pars[6L], st0 = st0, # st0
                                   precision = prec))
    # dist = abs(rtdists::ddiffusion(rt = t, response = boundary,
    #                                a = pars[1L], v = pars[2L], t0 = 0, z = pars[4L], # pars[3]
    #                                d = 0, sz = sz, sv = pars[6L], st0 = 0, # st0
    #                                precision = prec))
    # by = t[2L] - t[1L]
    # a <- t0
    # b <- t0 + st0
    # if (!(abs(a - b) <= by)) {
    #   ND = rev(stats::dunif(t, t0, t0 + st0))
    #   if (any(ND != 0)) {
    #     # dist = zapsmall(customConvolveO(dist, by * ND)[seq_along(t)], 13)
    #   }
    # }
  }
  return(dist)
}

# This function uess customDdiffusion rather than rtdists::ddiffusion
# Voss.density = function(t, pars, boundary, DstarM = TRUE, prec = 3) {
#   if (DstarM) {
#     sz = pars[4L] * 2 * min(c(pars[3L], 1L - pars[3L])) # rescale sz
#     dist = abs(customDdiffusion(t = t, boundary = boundary,
#                                 a = pars[1L], v = pars[2L], t0 = 0, z = pars[3],
#                                 d = 0L, sz = sz, sv = pars[5L], st0 = 0,
#                                 precision = prec))
#   } else {
#     sz = pars[5L] * 2 * min(c(pars[4L], 1 - pars[4L])) # rescale sz
#     st0 = pars[7L] * 2 * pars[3L] # rescale st0 to % * max size
#     dist = abs(customDdiffusion(t = t, boundary = boundary,
#                                 a = pars[1L], v = pars[2L], t0 = 0, z = pars[4L], # pars[3]
#                                 d = 0, sz = sz, sv = pars[6L], st0 = 0, # st0
#                                 precision = prec))
#     ND = rev(stats::dunif(t, pars[3L] - .5*st0, pars[3L] + .5*st0 + ((t[2L] - t[1L]) / 100) ))
#     if (any(ND != 0)) {
#       dist = zapsmall(customConvolveO(abs(dist), ND)[seq_along(t)], 13)
#     }
#   }
#   return(dist)
# }
# custom wrapper around dfastdm_b ~3x faster than rtdists::ddiffusion
# customDdiffusion = function(t, boundary = c("upper", "lower"), a, v, t0, z = 0.5,
#                      d = 0, sz = 0, sv = 0, st0 = 0, precision = 3) {
#   pl = c(a, v, t0, d, sz, sv, st0, z)
#   i = ifelse(boundary == 'upper', 2L, 1L)
#   output = .C(rtdists:::dfastdm_b, as.integer(length(t)), as.vector(pl),
#               t, as.double(precision), i, rep.int(0, length(t)))
#   return(output[[6L]])
# }

#' @rdname Voss.density
#' @export
# LBA density
LBA.density = function(t, pars, boundary, DstarM = TRUE, ...) {
  mean_v = pars[3:4 + !DstarM]
  sd_v = pars[5:6 + !DstarM]
  A = pars[2] * pars[1]
  if (DstarM) {
    pdf = rtdists::dLBA(t, ifelse(boundary == 'upper', 1, 2), A = A, b = pars[2L],
                        t0 = 0, mean_v = mean_v, sd_v = sd_v, ..., silent = TRUE)
  } else {
    pdf = rtdists::dLBA(t, ifelse(boundary == 'upper', 1, 2), A = A, b = pars[2L],
                        t0 = pars[3L], mean_v = mean_v, sd_v = sd_v, st0 = pars[8L], silent = TRUE, ...)
  }
  return(pdf)
}

#' @rdname Voss.density
#' @export
Wiener.density = function(t, pars, boundary, DstarM) {
  if (DstarM) {
    dist = RWiener::dwiener(q = t, alpha = pars[1], tau = 1e-100,
                            beta = pars[3], delta = pars[2], resp = rep(boundary, length(t)))
  } else {
    dist = RWiener::dwiener(q = t, alpha = pars[1], tau = pars[4],
                            beta = pars[3], delta = pars[2], resp = rep(boundary, length(t)))
  }
  return(abs(dist))
}

# custom wrapper around dwiener_c - DEPRECATED
# customDwiener = function (q, alpha, tau, beta, delta, resp = "upper", give_log = FALSE) {
#   d = vector("double", length = length(q))
#   if (resp == 'upper') {
#     for (i in seq_along(q)) d[i] = .Call(RWiener:::dwiener_c, q[i], alpha, tau, beta, delta, give_log)
#   } else {
#     q = -q
#     for (i in seq_along(q)) d[i] = .Call(RWiener:::dwiener_c, q[i], alpha, tau, beta, delta, give_log)
#   }
#   return(d)
# }





