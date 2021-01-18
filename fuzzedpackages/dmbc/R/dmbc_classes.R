#' An S4 class to represent the data to use in a DMBC model.
#'
#' @slot diss A list whose elements are the dissimilarity matrices corresponding
#'   to the judgments expressed by the \emph{S} subjects/raters. These matrices
#'   must be defined as a \code{dist} object.
#' @slot n A length-one character vector representing the number of objects
#'   compared by each subject.
#' @slot S A length-one numeric vector representing the number of subjects.
#' @slot family A length-one character vector representing the type of data to
#'   analyze. Currently, it accepts only the 'binomial' value, but future
#'   developments will include the possibility to analyze continuous,
#'   multinomial and count data.
#'
#' @name dmbc_data-class
#' @rdname dmbc_data-class
#' @aliases dmbc_data
#'
#' @author Sergio Venturini \email{sergio.venturini@unito.it}
#'
#' @references
#'   Venturini, S., Piccarreta, R. (2019), "A Bayesian Approach for Model-Based
#'   Clustering of Several Binary Dissimilarity Matrices: the \pkg{dmbc}
#'   Package in \code{R}", Technical report.
#'
#' @examples
#' showClass("dmbc_data")
#'
#' @exportClass dmbc_data
setClass(Class = "dmbc_data",
  slots = c(
    diss = "list",
    n = "numeric",
    S = "numeric",
    family = "character"
  )
)

#' Create an instance of the \code{dmbc_data} class using new/initialize.
#'
#' @param .Object Prototype object from the class \code{\link{dmbc_data}}.
#' @param diss A list whose elements are the dissimilarity matrices corresponding
#'   to the judgments expressed by the \emph{S} subjects/raters. These matrices
#'   must be defined as a \code{dist} object.
#' @param n A length-one character vector representing the number of objects
#'   compared by each subject.
#' @param S A length-one numeric vector representing the number of subjects.
#' @param family A length-one character vector representing the type of data to
#'   analyze. Currently, it accepts only the 'binomial' value, but future
#'   developments will include the possibility to analyze continuous,
#'   multinomial and count data.
#'
#' @author Sergio Venturini \email{sergio.venturini@unito.it}
#'
#' @aliases initialize,dmbc_data-method
#' @aliases dmbc_data-initialize
#' 
#' @importFrom methods initialize
#' @exportMethod initialize
setMethod("initialize", "dmbc_data",
  function(
    .Object,
    diss = list(),
    n = numeric(),
    S = numeric(),
    family = character()
  )
  {
    .Object@diss <- diss
    .Object@n <- n
    .Object@S <- S
    # .Object@family <- family
    .Object@family <- "binomial"
    .Object
  }
)

#' Show an instance of the \code{dmbc_data} class.
#'
#' @author Sergio Venturini \email{sergio.venturini@unito.it}
#'
#' @param object An object of class \code{\link{dmbc_data}}.
#'
#' @aliases show,dmbc_data-method
#' @aliases dmbc_data-show
#' 
#' @importFrom methods show
#' @exportMethod show
setMethod("show",
  "dmbc_data",
  function(object) {
    cat("Observed dissimilarity matrices to use in a DMBC analysis\n")
    cat("Number of objects (n):", object@n, "\n")
    cat("Number of subjects (S):", object@S, "\n")
    cat("Family:", object@family, "\n")
  }
)

#' Provide a summary of a \code{dmbc_data} class instance.
#'
#' @param object An object of class \code{\link{dmbc_data}}.
#'
#' @author Sergio Venturini \email{sergio.venturini@unito.it}
#'
#' @aliases summary,dmbc_data-method
#' @aliases dmbc_data-summary
#' 
#' @exportMethod summary
setMethod("summary",
  "dmbc_data",
    function(object) {
      show(object)
      cat("Observed dissimilarities:\n")

      wd <- as.integer(options("width"))
      for (s in 1:object@S) {
        dash_string_L <- strrep("-", floor((wd - nchar(paste0("Subject ", object@S)) - 2)/2))
        dash_string_R <- strrep("-", floor((wd - nchar(paste0("Subject ", object@S)) - 2)/2) - nchar(s) + 1)
        cat(dash_string_L, " ", "Subject ", s, " ", dash_string_R, "\n", sep = "")
        diss_nm <- attr(object@diss[[s]], "Labels")
        print_matrix(as.matrix(object@diss[[s]]), rownm = diss_nm, colnm = diss_nm,
          ndigits = ifelse(object@family != "normal", 0, 2),
          colwidth = ifelse(is.null(diss_nm), 5, max(nchar(diss_nm))),
          isint = (object@family != "normal"), between_cols = 1)
      }
    }
)

#' Provide a graphical summary of a \code{dmbc_data} class instance.
#'
#' @param x An object of class \code{\link{dmbc_data}}.
#' @param colors A character vector providing the colors to use in the plot.
#' @param font A length-one numeric vector for the font to use for text.
#'   Can be a vector. \code{NA} values (the default) mean use \code{par("font")}.
#' @param cex.font A length-one numeric vector for the character expansion
#'   factor. \code{NULL} and \code{NA} are equivalent to \code{1.0}. This is an
#'   absolute measure, not scaled by \code{par("cex")} or by setting
#''   \code{par("mfrow")} or \code{par("mfcol")}. Can be a vector.
#' @param ... Further arguments to pass on (currently ignored).
#'
#' @author Sergio Venturini \email{sergio.venturini@unito.it}
#'
#' @aliases plot,dmbc_data-method
#' @aliases dmbc_data-plot
#' 
#' @exportMethod plot
#' 
#' @examples
#' data(simdiss)
#' library(bayesplot)
#' cols <- color_scheme_set("brightblue")
#' plot(simdiss, colors = unlist(cols)[c(1, 6)], font = 1, cex.font = 0.75)
setMethod("plot",
  signature(x = "dmbc_data"),
  function(x, colors = c("white", "black"), font = NA, cex.font = NA, ...) {
    if (x@family == "binomial") {
      plot_dmbc_data <- function(ncat, colors, crit, INPUT, MD_SOL, MD_DIM,
        labxaxis = "", labyaxis = "") {
        nsat <- length((1:nrow(INPUT))[crit])
        if (nsat < 2) {
          a <- matrix(0, ncol = ncol(INPUT), nrow = 1)
        }
        if (nsat > 1) {
          a <- as.matrix(INPUT[crit, ][order(MD_SOL[crit, MD_DIM]), ])
        }
        graphics::image(x = 1:nrow(a), z = a, y = 1:ncol(a), axes = FALSE,
          zlim = c(1, ncat), xlim = c(1, nsat), xlab = NA, ylab = NA, col = colors)
         
        graphics::box()
        graphics::axis(side = 1, tck = -.03, labels = NA)
        graphics::axis(side = 2, tck = -.03, labels = NA)
        graphics::axis(side = 1, lwd = 0, line = -.9, cex.axis = 0.6)
        graphics::axis(side = 2, lwd = 0, line = -.6, las = 1, cex.axis = 0.6)
        graphics::mtext(side = 1, labxaxis, line = 2, cex = 0.6)
        graphics::mtext(side = 2, labyaxis, line = 2.5, cex = 0.6)
      }

      D <- x@diss
      n <- x@n
      S <- x@S
      col_bn <- colors[1:2]
      nr <- floor(sqrt(S))
      nc <- S %/% nr
      nr <- ifelse(S %% nr, nr + 1, nr)
      opar <- graphics::par(c("mfrow", "mar", "oma"))
      on.exit(par(opar))
      graphics::par(mfrow = c(nr, nc), mar = c(1, .75, .75, .5) + 0.1, oma = c(1, 0, 0, 0))
      for(s in 1:S) {
        use <- as.matrix(D[[s]]) + 1
        plot_dmbc_data(2, col_bn, use[, 1] > 0, use, as.matrix((1:nrow(use))), 1)
        mtext(paste0("Subject ", s), side = 3, font = font, line = 0.05, cex = cex.font)
      }
    } else {
      stop("no plot methods implemented yet for non-binary data.")
    }
  }
)

#' An S4 class to represent a DMBC model.
#'
#' @slot p A length-one character vector representing the number of dimensions
#'   of the latent space to use in the MDS analysis.
#' @slot G A length-one numeric vector representing the number of clusters to
#'   partition the subjects into.
#' @slot family A length-one character vector representing the type of data to
#'   analyze. Currently, it accepts only the 'binomial' value, but future
#'   developments will include the possibility to analyze continuous,
#'   multinomial and count data.
#'
#' @name dmbc_model-class
#' @rdname dmbc_model-class
#' @aliases dmbc_model
#'
#' @references
#'   Venturini, S., Piccarreta, R. (2019), "A Bayesian Approach for Model-Based
#'   Clustering of Several Binary Dissimilarity Matrices: the \pkg{dmbc}
#'   Package in \code{R}", Technical report.
#'
#' @examples
#' showClass("dmbc_model")
#'
#' @exportClass dmbc_model
setClass(Class = "dmbc_model",
  slots = c(
    p = "numeric",
    G = "numeric",
    family = "character"
  )
)

#' Create an instance of the \code{dmbc_model} class using new/initialize.
#'
#' @param .Object Prototype object from the class \code{\link{dmbc_model}}.
#' @param p A length-one character vector representing the number of dimensions
#'   of the latent space to use in the MDS analysis.
#' @param G A length-one numeric vector representing the number of clusters to
#'   partition the subjects into.
#' @param family A length-one character vector representing the type of data to
#'   analyze. Currently, it accepts only the 'binomial' value, but future
#'   developments will include the possibility to analyze continuous,
#'   multinomial and count data.
#'
#' @author Sergio Venturini \email{sergio.venturini@unito.it}
#'
#' @aliases initialize,dmbc_model-method
#' @aliases dmbc_model-initialize
#' 
#' @importFrom methods initialize
#' @exportMethod initialize
setMethod("initialize", "dmbc_model",
  function(
    .Object,
    p = numeric(),
    G = numeric(),
    family = character()
  )
  {
    .Object@p <- p
    .Object@G <- G
    .Object@family <- family
    .Object
  }
)

#' Show an instance of the \code{dmbc_model} class.
#'
#' @param object An object of class \code{\link{dmbc_model}}.
#'
#' @author Sergio Venturini \email{sergio.venturini@unito.it}
#'
#' @aliases show,dmbc_model-method
#' @aliases dmbc_model-show
#' 
#' @importFrom methods show
#' @exportMethod show
setMethod("show",
  "dmbc_model",
  function(object) {
    cat("Dissimilarity Model Based Clustering definition\n")
    cat("Number of latent dimensions (p):", object@p, "\n")
    cat("Number of clusters (G):", object@G, "\n")
    cat("Family:", object@family, "\n")
  }
)

#' An S4 class to represent the results of fitting DMBC model.
#'
#' @description
#'   An S4 class to represent the results of fitting DMBC model using a single
#'   Markov Chain Monte Carlo chain.
#'
#' @slot z.chain An object of class \code{array}; posterior draws from
#'   the MCMC algorithm for the (untransformed) latent configuration \eqn{Z}.
#' @slot z.chain.p An object of class \code{array}; posterior draws from
#'   the MCMC algorithm for the (Procrustes-transformed) latent configuration
#'   \eqn{Z}.
#' @slot alpha.chain An object of class \code{matrix}; posterior draws
#'   from the MCMC algorithm for the \eqn{\alpha} parameters.
#' @slot eta.chain An object of class \code{matrix}; posterior draws
#'   from the MCMC algorithm for the \eqn{\eta} parameters.
#' @slot sigma2.chain An object of class \code{matrix}; posterior draws
#'   from the MCMC algorithm for the \eqn{\sigma^2} parameters.
#' @slot lambda.chain An object of class \code{matrix}; posterior draws
#'   from the MCMC algorithm for the \eqn{\lambda} parameters.
#' @slot prob.chain An object of class \code{array}; posterior draws
#'   from the MCMC algorithm for the cluster membership probabilities.
#' @slot x.ind.chain An object of class \code{array}; posterior draws
#'   from the MCMC algorithm for the cluster membership indicators.
#' @slot x.chain An object of class \code{matrix}; posterior draws from
#'   the MCMC algorithm for the cluster membership labels.
#' @slot accept An object of class \code{matrix}; final acceptance rates
#'   for the MCMC algorithm.
#' @slot diss An object of class \code{list}; list of observed
#'   dissimilarity matrices.
#' @slot dens An object of class \code{list}; list of log-likelihood,
#'   log-prior and log-posterior values at each iteration of the MCMC simulation.
#' @slot control An object of class \code{list}; list of the control
#'   parameters (number of burnin and sample iterations, number of MCMC chains,
#'   etc.). See \code{\link{dmbc_control}()} for more information.
#' @slot prior An object of class \code{list}; list of the prior
#'   hyperparameters. See \code{\link{dmbc_prior}()} for more information.
#' @slot dim An object of class \code{list}; list of dimensions for
#'   the estimated model, i.e. number of objects (\emph{n}), number of latent
#'   dimensions (\emph{p}), number of clusters (\emph{G}), and number of
#'   subjects (\emph{S}).
#' @slot model An object of class \code{\link{dmbc_model}}.
#'
#' @name dmbc_fit-class
#' @rdname dmbc_fit-class
#'
#' @references
#'   Venturini, S., Piccarreta, R. (2019), "A Bayesian Approach for Model-Based
#'   Clustering of Several Binary Dissimilarity Matrices: the \pkg{dmbc}
#'   Package in \code{R}", Technical report.
#'
#' @examples
#' showClass("dmbc_fit")
#'
#' @exportClass dmbc_fit
setClass(Class = "dmbc_fit",
	slots = c(
		z.chain = "array",
		z.chain.p = "array",
		alpha.chain = "matrix",
		eta.chain = "matrix",
		sigma2.chain = "matrix",
		lambda.chain = "matrix",
		prob.chain = "array",
		x.ind.chain = "array",
		x.chain = "matrix",
		accept = "matrix",
		diss = "list",
		dens = "list",
		control = "list",
    prior = "list",
		dim = "list",
    model = "dmbc_model"
	)
)

#' Create an instance of the \code{dmbc_fit} class using new/initialize.
#'
#' @param .Object Prototype object from the class \code{\link{dmbc_fit}}.
#' @param z.chain An object of class \code{array}; posterior draws from
#'   the MCMC algorithm for the (untransformed) latent configuration \eqn{Z}.
#' @param z.chain.p An object of class \code{array}; posterior draws from
#'   the MCMC algorithm for the (Procrustes-transformed) latent configuration
#'   \eqn{Z}.
#' @param alpha.chain An object of class \code{matrix}; posterior draws
#'   from the MCMC algorithm for the \eqn{\alpha} parameters.
#' @param eta.chain An object of class \code{matrix}; posterior draws
#'   from the MCMC algorithm for the \eqn{\eta} parameters.
#' @param sigma2.chain An object of class \code{matrix}; posterior draws
#'   from the MCMC algorithm for the \eqn{\sigma^2} parameters.
#' @param lambda.chain An object of class \code{matrix}; posterior draws
#'   from the MCMC algorithm for the \eqn{\lambda} parameters.
#' @param prob.chain An object of class \code{array}; posterior draws
#'   from the MCMC algorithm for the cluster membership probabilities.
#' @param x.ind.chain An object of class \code{array}; posterior draws
#'   from the MCMC algorithm for the cluster membership indicators.
#' @param x.chain An object of class \code{matrix}; posterior draws from
#'   the MCMC algorithm for the cluster membership labels.
#' @param accept An object of class \code{matrix}; final acceptance rates
#'   for the MCMC algorithm.
#' @param diss An object of class \code{list}; list of observed
#'   dissimilarity matrices.
#' @param dens An object of class \code{list}; list of log-likelihood,
#'   log-prior and log-posterior values at each iteration of the MCMC simulation.
#' @param control An object of class \code{list}; list of the control
#'   parameters (number of burnin and sample iterations, number of MCMC chains,
#'   etc.). See \code{\link{dmbc_control}()} for more information.
#' @param prior An object of class \code{list}; list of the prior
#'   hyperparameters. See \code{\link{dmbc_prior}()} for more information.
#' @param dim An object of class \code{list}; list of dimensions for
#'   the estimated model, i.e. number of objects (\emph{n}), number of latent
#'   dimensions (\emph{p}), number of clusters (\emph{G}), and number of
#'   subjects (\emph{S}).
#' @param model An object of class \code{\link{dmbc_model}}.
#'
#' @author Sergio Venturini \email{sergio.venturini@unito.it}
#'
#' @aliases initialize,dmbc_fit-method
#' @aliases dmbc_fit-initialize
#' 
#' @importFrom methods initialize
#' @exportMethod initialize
setMethod("initialize",
  "dmbc_fit",
		function(
			.Object,
			z.chain = array(),
			z.chain.p = array(),
			alpha.chain = matrix(),
			eta.chain = matrix(),
			sigma2.chain = matrix(),
			lambda.chain = matrix(),
			prob.chain = array(),
			x.ind.chain = array(),
			x.chain = matrix(),
			accept = matrix(),
			diss = list(),
			dens = list(),
			control = list(),
      prior = list(),
			dim = list(),
      model = NA
		)
		{
			.Object@z.chain <- z.chain
			.Object@z.chain.p <- z.chain.p
			.Object@alpha.chain <- alpha.chain
			.Object@eta.chain <- eta.chain
			.Object@sigma2.chain <- sigma2.chain
			.Object@lambda.chain <- lambda.chain
			.Object@prob.chain <- prob.chain
			.Object@x.ind.chain <- x.ind.chain
			.Object@x.chain <- x.chain
			.Object@accept <- accept
			.Object@diss <- diss
			.Object@dens <- dens
      .Object@control <- control
      .Object@prior <- prior
			.Object@dim <- dim
      .Object@model <- model
			.Object
		}
)

#' Show an instance of the \code{dmbc_fit} class.
#'
#' @param object An object of class \code{\link{dmbc_fit}}.
#'
#' @author Sergio Venturini \email{sergio.venturini@unito.it}
#'
#' @aliases show,dmbc_fit-method
#' @aliases dmbc_fit-show
#' 
#' @importFrom methods show
#' @exportMethod show
setMethod("show",
  "dmbc_fit",
  function(object) {
    cat("Dissimilarity Model Based Clustering simulated chain\n")
    cat("Number of latent dimensions (p):", object@model@p, "\n")
    cat("Number of clusters (G):", object@model@G, "\n")
    cat("Family:", object@model@family, "\n")
    cat("\n")
    cat("To get a summary of the object, use the 'summary()' function.")
  }
)

#' Provide a summary of a \code{dmbc_fit} class instance.
#'
#' @param object An object of class \code{\link{dmbc_fit}}.
#' @param include.burnin A length-one logical vector. If \code{TRUE} the
#'   burnin iterations (if available) are included in the summary.
#' @param summary.Z A length-one logical vector. If \code{TRUE} the summary
#'   also includes the latent configuration coordinates.
#' @param ... Further arguments to pass on (currently ignored).
#'
#' @author Sergio Venturini \email{sergio.venturini@unito.it}
#'
#' @aliases summary,dmbc_fit-method
#' @aliases dmbc_fit-summary
#' 
#' @exportMethod summary
setMethod("summary",
	"dmbc_fit",
    function(object, include.burnin = FALSE, summary.Z = FALSE, ...) {
      control <- object@control

      n <- object@dim[["n"]]
      p <- object@dim[["p"]]
      G <- object@dim[["G"]]
      S <- object@dim[["S"]]

      res.coda <- dmbc_fit_to_mcmc(object, include.burnin = include.burnin, verbose = FALSE)

      if (summary.Z) {
        res.coda <- res.coda[, (n*p*G + 1):(2*n*p*G + 4*G)]
      } else {
        res.coda <- res.coda[, (2*n*p*G + 1):(2*n*p*G + 4*G)]
      }

      out <- summary(res.coda)

      return(out)
    }
)

#' @export
setGeneric("subset", function(x) standardGeneric("subset"))

#' Subsetting a \code{dmbc_fit} object.
#'
#' @param x An object of class \code{\link{dmbc_fit}}.
#' @param pars An optional character vector of parameter names. If neither 
#'   \code{pars} nor \code{regex_pars} is specified, the default is to use all
#'   parameters.
#' @param regex_pars An optional \code{\link[=grep]{regular expression}} to use
#'   for parameter selection. Can be specified instead of \code{pars} or in addition
#'   to \code{pars}.
#' @param ... Further arguments to pass on (currently ignored).
#'
#' @author Sergio Venturini \email{sergio.venturini@unito.it}
#'
#' @aliases subset,dmbc_fit-method
#' @aliases dmbc_fit-subset
#' 
#' @export
setMethod("subset",
  "dmbc_fit",
  function(x, pars = character(), regex_pars = character(), ...) {
    x_mcmc <- dmbc_fit_to_mcmc(x, include.burnin = TRUE, verbose = FALSE)

    parnames <- colnames(x_mcmc)
    pars <- select_pars(explicit = pars, patterns = regex_pars, complete = parnames)
    
    out <- x_mcmc[, pars]
    
    return(out)
  }
)

#' Provide a graphical summary of a \code{dmbc_fit} class instance.
#'
#' @param x An object of class \code{\link{dmbc_fit}}.
#' @param what A length-one character vector providing the plot type to produce.
#'   Admissible values are those provided by the \pkg{\link{bayesplot}} package,
#'   that is: \code{acf}, \code{areas}, \code{dens}, \code{hex}, \code{hist},
#'   \code{intervals}, \code{neff}, \code{pairs}, \code{parcoord}, \code{recover},
#'   \code{rhat}, \code{scatter}, \code{trace}, \code{violin} or \code{combo}.
#'   In particular, \code{combo} allows to mix different plot types. For more
#'   details see the documentation of the \pkg{\link{bayesplot}} package,
#'   starting from \code{\link[=MCMC-overview]{this overview page}}.
#' @param pars An optional character vector of parameter names. If neither 
#'   \code{pars} nor \code{regex_pars} is specified, the default is to use all parameters.
#' @param regex_pars An optional \code{\link[=grep]{regular expression}} to use for
#'   parameter selection. Can be specified instead of \code{pars} or in addition to
#'   \code{pars}.
#' @param include.burnin A length-one logical vector. If \code{TRUE} the
#'   burnin iterations (if available) are included in the summary.
#' @param combo A character vector providing the plot types to combine (see
#'   \code{\link[bayesplot]{mcmc_combo}}).
#' @param ... Further arguments to pass on.
#'
#' @author Sergio Venturini \email{sergio.venturini@unito.it}
#'
#' @aliases plot,dmbc_fit-method
#' @aliases dmbc_fit-plot
#' 
#' @exportMethod plot
setMethod("plot",
  signature(x = "dmbc_fit"),
  function(x, what = "trace", pars = character(), regex_pars = "lambda", include.burnin = FALSE,
    combo = NULL, ...) {
    stopifnot(is.character(pars),
              is.character(regex_pars),
              is.character(what))
    
    if (!(what %in% unlist(all_plots_list, use.names = FALSE)))
      stop("the plot type specified is not available.")

    x_mcmc <- dmbc_fit_to_mcmc(x, include.burnin = include.burnin, verbose = FALSE)

    control <- x@control

    if (what %in% acf_plot_list) {
      if (what == "acf") {
        p <- bayesplot::mcmc_acf(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
      } else if (what == "acf_bar") {
        p <- bayesplot::mcmc_acf_bar(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
      }
    }

    if (what %in% areas_plot_list) {
      if (what == "areas") {
        p <- bayesplot::mcmc_areas(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
      } else if (what == "areas_ridges") {
        p <- bayesplot::mcmc_areas_ridges(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
      }
    }

    if (what %in% dens_plot_list) {
      if (what == "dens") {
        p <- bayesplot::mcmc_dens(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
      } else if (what == "dens_overlay") {
        p <- bayesplot::mcmc_dens_overlay(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
      } else if (what == "dens_chains") {
        p <- bayesplot::mcmc_dens_chains(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
      }
    }

    if (what %in% hex_plot_list) {
      if (what == "hex") {
        p <- bayesplot::mcmc_hex(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
      }
    }

    if (what %in% hist_plot_list) {
      if (what == "hist") {
        p <- bayesplot::mcmc_hist(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
      } else if (what == "hist_by_chain") {
        p <- bayesplot::mcmc_hist_by_chain(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
      }
    }

    if (what %in% intervals_plot_list) {
      if (what == "intervals") {
        p <- bayesplot::mcmc_intervals(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
      }
    }

    if (what %in% neff_plot_list) {
      x_sub <- subset(x, pars = pars, regex_pars = regex_pars)
      nsample <- floor((control[["burnin"]] + control[["nsim"]])/control[["thin"]])
      neff <- coda::effectiveSize(x_sub)
      ratio <- neff/nsample
      if (what == "neff") {
        p <- bayesplot::mcmc_neff(ratio = ratio, ...)
      } else if (what == "neff_hist") {
        p <- bayesplot::mcmc_neff_hist(ratio = ratio, ...)
      }
    }

    if (what %in% pairs_plot_list) {
      if (what == "pairs") {
        p <- bayesplot::mcmc_pairs(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
      }
    }

    if (what %in% parcoord_plot_list) {
      if (what == "parcoord") {
        p <- bayesplot::mcmc_parcoord(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
      }
    }

    if (what %in% recover_plot_list) {
      x_sub <- subset(x, pars = pars, regex_pars = regex_pars)
      if (what == "recover_hist") {
        p <- bayesplot::mcmc_recover_hist(x = x_sub, ...)
      } else if (what == "recover_intervals") {
        p <- bayesplot::mcmc_recover_intervals(x = x_sub, ...)
      } else if (what == "recover_scatter") {
        p <- bayesplot::mcmc_recover_scatter(x = x_sub, ...)
      }
    }

    if (what %in% rhat_plot_list) {
      x_sub <- subset(x, pars = pars, regex_pars = regex_pars)
      rhat <- coda::gelman.diag(x_sub, multivariate = FALSE)$psrf[, 1]
      if (what == "rhat") {
        p <- bayesplot::mcmc_rhat(rhat = rhat, ...)
      } else if (what == "rhat_hist") {
        p <- bayesplot::mcmc_rhat_hist(rhat = rhat, ...)
      }
    }

    if (what %in% scatter_plot_list) {
      if (what == "scatter") {
        p <- bayesplot::mcmc_scatter(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
      }
    }

    if (what %in% trace_plot_list) {
      if (what == "trace") {
        p <- bayesplot::mcmc_trace(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
      } else if (what == "trace_highlight") {
        p <- bayesplot::mcmc_trace_highlight(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
      }
    }

    if (what %in% violin_plot_list) {
      if (what == "violin") {
        p <- bayesplot::mcmc_violin(x = x_mcmc, pars = pars, regex_pars = regex_pars, ...)
      }
    }

    if (what == "combo") {
      if (!is.null(combo)) {
        p <- bayesplot::mcmc_combo(x = x_mcmc, pars = pars, regex_pars = regex_pars, combo = combo, ...)
      } else {
        stop("to produce an 'mcmc_combo' plot, the 'combo' option must be specified.")
      }
    }

    p
  }
)

#' An S4 class to represent the results of fitting DMBC model.
#'
#' @description
#'   An S4 class to represent the results of fitting DMBC model using multiple
#'   Markov Chain Monte Carlo chains.
#'
#' @slot results An object of class \code{list}; list of \code{dmbc_fit}
#'   objects corresponding to the parallel MCMC chains simulated during the
#'   estimation.
#'
#' @name dmbc_fit_list-class
#' @rdname dmbc_fit_list-class
#' @aliases dmbc_fit_list
#'
#' @seealso
#' \code{\link{dmbc_fit}} for more details on the components of each element of
#'   the list.
#'
#' @references
#'   Venturini, S., Piccarreta, R. (2019), "A Bayesian Approach for Model-Based
#'   Clustering of Several Binary Dissimilarity Matrices: the \pkg{dmbc}
#'   Package in \code{R}", Technical report.
#'
#' @examples
#' showClass("dmbc_fit_list")
#'
#' @exportClass dmbc_fit_list
setClass(Class = "dmbc_fit_list",
  slots = c(
    results = "list"
  )
)

#' Create an instance of the \code{dmbc_fit_list} class using new/initialize.
#'
#' @param .Object Prototype object from the class \code{\link{dmbc_fit_list}}.
#' @param results A list whose elements are the \code{dmbc_fit} objects for
#'   each simulated chain.
#'
#' @author Sergio Venturini \email{sergio.venturini@unito.it}
#'
#' @aliases initialize,dmbc_fit_list-method
#' @aliases dmbc_fit_list-initialize
#' 
#' @importFrom methods initialize
#' @exportMethod initialize
setMethod("initialize", "dmbc_fit_list",
  function(
    .Object,
    results = list()
  )
  {
    .Object@results <- results
    .Object
  }
)

#' Show an instance of the \code{dmbc_fit_list} class.
#'
#' @param object An object of class \code{\link{dmbc_fit_list}}.
#'
#' @author Sergio Venturini \email{sergio.venturini@unito.it}
#'
#' @aliases show,dmbc_fit_list-method
#' @aliases dmbc_fit_list-show
#' 
#' @importFrom methods show
#' @exportMethod show
setMethod("show",
  "dmbc_fit_list",
  function(object) {
    cat("List of Dissimilarity Model Based Clustering simulated chains\n")
    cat("Number of simulated chains:", length(object@results), "\n")
    cat("Number of latent dimensions (p):", object@results[[1]]@model@p, "\n")
    cat("Number of clusters (G):", object@results[[1]]@model@G, "\n")
    cat("Family:", object@results[[1]]@model@family, "\n")
    cat("\n")
    cat("To get a summary of the object, use the 'summary()' function.")
  }
)

#' Provide a summary of a \code{dmbc_fit_list} class instance.
#'
#' @param object An object of class \code{\link{dmbc_fit_list}}.
#' @param include.burnin A length-one logical vector. If \code{TRUE} the
#'   burnin iterations (if available) are included in the summary.
#' @param summary.Z A length-one logical vector. If \code{TRUE} the summary
#'   also includes the latent configuration coordinates.
#' @param ... Further arguments to pass on (currently ignored).
#'
#' @author Sergio Venturini \email{sergio.venturini@unito.it}
#'
#' @aliases summary,dmbc_fit_list-method
#' @aliases dmbc_fit_list-summary
#' 
#' @exportMethod summary
setMethod("summary",
  "dmbc_fit_list",
    function(object, include.burnin = FALSE, summary.Z = FALSE, ...) {
      control <- object@results[[1]]@control
      nchains <- control[["nchains"]]

      n <- object@results[[1]]@dim[["n"]]
      p <- object@results[[1]]@dim[["p"]]
      G <- object@results[[1]]@dim[["G"]]
      S <- object@results[[1]]@dim[["S"]]

      res.coda <- dmbc_fit_list_to_mcmc.list(object, include.burnin = include.burnin, verbose = FALSE)

      for (c in 1:nchains) {
        if (summary.Z) {
          res.coda[[c]] <- res.coda[[c]][, (n*p*G + 1):(2*n*p*G + 4*G)]
        } else {
          res.coda[[c]] <- res.coda[[c]][, (2*n*p*G + 1):(2*n*p*G + 4*G)]
        }
      }

      out <- summary(res.coda)

      return(out)
    }
)

#' Subsetting a \code{dmbc_fit_list} object.
#'
#' @param x An object of class \code{\link{dmbc_fit_list}}.
#' @param pars An optional character vector of parameter names. If neither 
#'   \code{pars} nor \code{regex_pars} is specified, the default is to use all
#'   parameters.
#' @param regex_pars An optional \code{\link[=grep]{regular expression}} to use
#'   for parameter selection. Can be specified instead of \code{pars} or in addition
#'   to \code{pars}.
#' @param ... Further arguments to pass on (currently ignored).
#'
#' @author Sergio Venturini \email{sergio.venturini@unito.it}
#'
#' @aliases subset,dmbc_fit_list-method
#' @aliases dmbc_fit_list-subset
#' 
#' @export
setMethod("subset",
  "dmbc_fit_list",
  function(x, pars = character(), regex_pars = character(), ...) {
    x_mcmc.list <- dmbc_fit_list_to_mcmc.list(x, include.burnin = TRUE, verbose = FALSE)

    parnames <- colnames(x_mcmc.list[[1]])
    pars <- select_pars(explicit = pars, patterns = regex_pars, complete = parnames)
    
    out <- coda::mcmc.list(lapply(x_mcmc.list, function(chain) chain[, pars]))
    
    return(out)
  }
)

#' Provide a graphical summary of a \code{dmbc_fit_list} class instance.
#'
#' @param x An object of class \code{\link{dmbc_fit_list}}.
#' @param what A length-one character vector providing the plot type to produce.
#'   Admissible values are those provided by the \pkg{\link{bayesplot}} package,
#'   that is: \code{acf}, \code{areas}, \code{dens}, \code{hex}, \code{hist},
#'   \code{intervals}, \code{neff}, \code{pairs}, \code{parcoord}, \code{recover},
#'   \code{rhat}, \code{scatter}, \code{trace}, \code{violin} or \code{combo}.
#'   In particular, \code{combo} allows to mix different plot types. For more
#'   details see the documentation of the \pkg{\link{bayesplot}} package,
#'   starting from \code{\link[=MCMC-overview]{this overview page}}.
#' @param pars An optional character vector of parameter names. If neither 
#'   \code{pars} nor \code{regex_pars} is specified, the default is to use all parameters.
#' @param regex_pars An optional \code{\link[=grep]{regular expression}} to use for
#'   parameter selection. Can be specified instead of \code{pars} or in addition to
#'   \code{pars}.
#' @param include.burnin A length-one logical vector. If \code{TRUE} the
#'   burnin iterations (if available) are included in the summary.
#' @param combo A character vector providing the plot types to combine (see
#'   \code{\link[bayesplot]{mcmc_combo}}).
#' @param ... Further arguments to pass on.
#'
#' @author Sergio Venturini \email{sergio.venturini@unito.it}
#'
#' @aliases plot,dmbc_fit_list-method
#' @aliases dmbc_fit_list-plot
#' 
#' @exportMethod plot
setMethod("plot",
	signature(x = "dmbc_fit_list"),
	function(x, what = "trace", pars = character(), regex_pars = "lambda", include.burnin = FALSE,
    combo = NULL, ...) {
    stopifnot(is.character(pars),
              is.character(regex_pars),
              is.character(what))
		
    if (!(what %in% unlist(all_plots_list, use.names = FALSE)))
      stop("the plot type specified is not available.")

    x_mcmc.list <- dmbc_fit_list_to_mcmc.list(x, include.burnin = include.burnin, verbose = FALSE)

    control <- x@results[[1]]@control

    if (what %in% acf_plot_list) {
      if (what == "acf") {
        p <- bayesplot::mcmc_acf(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, ...)
      } else if (what == "acf_bar") {
        p <- bayesplot::mcmc_acf_bar(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, ...)
      }
    }

    if (what %in% areas_plot_list) {
      if (what == "areas") {
        p <- bayesplot::mcmc_areas(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, ...)
      } else if (what == "areas_ridges") {
        p <- bayesplot::mcmc_areas_ridges(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, ...)
      }
    }

    if (what %in% dens_plot_list) {
      if (what == "dens") {
        p <- bayesplot::mcmc_dens(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, ...)
      } else if (what == "dens_overlay") {
        p <- bayesplot::mcmc_dens_overlay(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, ...)
      } else if (what == "dens_chains") {
        p <- bayesplot::mcmc_dens_chains(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, ...)
      }
    }

    if (what %in% hex_plot_list) {
      if (what == "hex") {
        p <- bayesplot::mcmc_hex(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, ...)
      }
    }

    if (what %in% hist_plot_list) {
      if (what == "hist") {
        p <- bayesplot::mcmc_hist(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, ...)
      } else if (what == "hist_by_chain") {
        p <- bayesplot::mcmc_hist_by_chain(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, ...)
      }
    }

    if (what %in% intervals_plot_list) {
      if (what == "intervals") {
        p <- bayesplot::mcmc_intervals(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, ...)
      }
    }

    if (what %in% neff_plot_list) {
      x_sub <- subset(x, pars = pars, regex_pars = regex_pars)
      nsample <- control[["nchains"]]*floor((control[["burnin"]] + control[["nsim"]])/control[["thin"]])
      neff <- coda::effectiveSize(x_sub)
      ratio <- neff/nsample
      if (what == "neff") {
        p <- bayesplot::mcmc_neff(ratio = ratio, ...)
      } else if (what == "neff_hist") {
        p <- bayesplot::mcmc_neff_hist(ratio = ratio, ...)
      }
    }

    if (what %in% pairs_plot_list) {
      if (what == "pairs") {
        p <- bayesplot::mcmc_pairs(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, ...)
      }
    }

    if (what %in% parcoord_plot_list) {
      if (what == "parcoord") {
        p <- bayesplot::mcmc_parcoord(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, ...)
      }
    }

    if (what %in% recover_plot_list) {
      x_sub <- subset(x, pars = pars, regex_pars = regex_pars)
      if (what == "recover_hist") {
        p <- bayesplot::mcmc_recover_hist(x = x_sub, ...)
      } else if (what == "recover_intervals") {
        p <- bayesplot::mcmc_recover_intervals(x = x_sub, ...)
      } else if (what == "recover_scatter") {
        p <- bayesplot::mcmc_recover_scatter(x = x_sub, ...)
      }
    }

    if (what %in% rhat_plot_list) {
      x_sub <- subset(x, pars = pars, regex_pars = regex_pars)
      rhat <- coda::gelman.diag(x_sub, multivariate = FALSE)$psrf[, 1]
      if (what == "rhat") {
        p <- bayesplot::mcmc_rhat(rhat = rhat, ...)
      } else if (what == "rhat_hist") {
        p <- bayesplot::mcmc_rhat_hist(rhat = rhat, ...)
      }
    }

    if (what %in% scatter_plot_list) {
      if (what == "scatter") {
        p <- bayesplot::mcmc_scatter(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, ...)
      }
    }

    if (what %in% trace_plot_list) {
      if (what == "trace") {
        p <- bayesplot::mcmc_trace(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, ...)
      } else if (what == "trace_highlight") {
        p <- bayesplot::mcmc_trace_highlight(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, ...)
      }
    }

    if (what %in% violin_plot_list) {
      if (what == "violin") {
        p <- bayesplot::mcmc_violin(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, ...)
      }
    }

    if (what == "combo") {
      if (!is.null(combo)) {
        p <- bayesplot::mcmc_combo(x = x_mcmc.list, pars = pars, regex_pars = regex_pars, combo = combo, ...)
      } else {
        stop("to produce an 'mcmc_combo' plot, the 'combo' option must be specified.")
      }
    }

    p
	}
)

#' An S4 class to represent the comparison of a set of DMBC models.
#'
#' @description
#'   An S4 class to represent the comparison of a set of DMBC models through
#'   the dissimilarity model-based clustering information criterion (DCIC).
#'
#' @slot logprior An object of class \code{matrix} providing the
#'   log-prior values corresponding to different values of \emph{p} and
#'   \emph{G}.
#' @slot logmlik An object of class \code{matrix} providing the
#'   marginal log-likelihood values corresponding to different values of
#'   \emph{p} and \emph{G}.
#' @slot logcorrfact An object of class \code{matrix} providing the
#'   logarithm of the correction factors corresponding to different values of
#'   \emph{p} and \emph{G}.
#' @slot DCIC An object of class \code{matrix} providing the values
#'   of the DCIC index corresponding to different values of \emph{p} and
#'   \emph{G}.
#' @slot post.est An object of class \code{list}; named list with
#'   elements representing the parameter estimates corresponding to different
#'   values of \emph{p} and \emph{G}.
#' @slot est A length-one character vector representing the estimate
#'   type used in computing the DCIC index. Possible values are \code{mean},
#'   \code{median}, \code{ml} and \code{map}. See \code{\link{dmbc_ic}()} for
#'   more details about these values.
#' @slot res_last_p An object of class \code{list}; list of
#'   \code{dmbc_fit_list} objects with the results of fitting the DMBC
#'   models corresponding to the last value of \emph{p}. This is needed in case
#'   of an update of the DCIC calculations using additional \emph{p} and/or
#'   \emph{G} values.
#'
#' @name dmbc_ic-class
#' @rdname dmbc_ic-class
#' @aliases dmbc_ic
#'
#' @references
#'   Venturini, S., Piccarreta, R. (2019), "A Bayesian Approach for Model-Based
#'   Clustering of Several Binary Dissimilarity Matrices: the \pkg{dmbc}
#'   Package in \code{R}", Technical report.
#'
#' @examples
#' showClass("dmbc_ic")
#'
#' @exportClass dmbc_ic
setClass(Class = "dmbc_ic",
	slots = c(
		logprior = "matrix",
		logmlik = "matrix",
		logcorrfact = "matrix",
		DCIC = "matrix",
		post.est = "list",
		est = "character",
		res_last_p = "list"
	)
)

#' Create an instance of the \code{dmbc_ic} class using new/initialize.
#'
#' @param .Object Prototype object from the class \code{\link{dmbc_ic}}.
#' @param logprior An object of class \code{matrix} providing the
#'   log-prior values corresponding to different values of \emph{p} and
#'   \emph{G}.
#' @param logmlik An object of class \code{matrix} providing the
#'   marginal log-likelihood values corresponding to different values of
#'   \emph{p} and \emph{G}.
#' @param logcorrfact An object of class \code{matrix} providing the
#'   logarithm of the correction factors corresponding to different values of
#'   \emph{p} and \emph{G}.
#' @param DCIC An object of class \code{matrix} providing the values
#'   of the DCIC index corresponding to different values of \emph{p} and
#'   \emph{G}.
#' @param post.est An object of class \code{list}; named list with
#'   elements representing the parameter estimates corresponding to different
#'   values of \emph{p} and \emph{G}.
#' @param est A length-one character vector representing the estimate
#'   type used in computing the DCIC index. Possible values are \code{mean},
#'   \code{median}, \code{ml} and \code{map}. See \code{\link{dmbc_ic}()} for
#'   more details about these values.
#' @param res_last_p An object of class \code{list}; list of
#'   \code{dmbc_fit_list} objects with the results of fitting the DMBC
#'   models corresponding to the last value of \emph{p}. This is needed in case
#'   of an update of the DCIC calculations using additional \emph{p} and/or
#'   \emph{G} values.
#'
#' @author Sergio Venturini \email{sergio.venturini@unito.it}
#'
#' @aliases initialize,dmbc_ic-method
#' @aliases dmbc_ic-initialize
#' 
#' @importFrom methods initialize
#' @exportMethod initialize
setMethod("initialize", "dmbc_ic",
	function(
		.Object,
		logprior = matrix(),
		logmlik = matrix(),
		logcorrfact = matrix(),
		DCIC = matrix(),
		post.est = list(),
		est = character(),
		res_last_p = list()
	)
	{
		.Object@logprior <- logprior
		.Object@logmlik <- logmlik
		.Object@logcorrfact <- logcorrfact
		.Object@DCIC <- DCIC
		.Object@post.est <- post.est
		.Object@est <- est
    .Object@res_last_p <- res_last_p
		.Object
	}
)

#' Show an instance of the \code{dmbc_ic} class.
#'
#' @param object An object of class \code{\link{dmbc_ic}}.
#'
#' @author Sergio Venturini \email{sergio.venturini@unito.it}
#'
#' @aliases show,dmbc_ic-method
#' @aliases dmbc_ic-show
#' 
#' @importFrom methods show
#' @exportMethod show
setMethod("show",
  "dmbc_ic",
    function(object) {
      pmax <- nrow(object@DCIC)
      Gmax <- ncol(object@DCIC)
      last <- length(object@res_last_p)
      family <- object@res_last_p[[1]]@results[[1]]@model@family

      cat("Dissimilarity Model Based Clustering -- Model choice\n")
      cat("Latent dimensions (p) requested: from 1 to", pmax, "\n")
      cat("Number of clusters (G) requested: from 1 to", Gmax, "\n")
      cat("Family:", family, "\n")
      cat("Dissimilarity Model Based Clustering Information Criterion (DCIC):\n")

      dcic_rownm <- paste0("p = ", 1:pmax)
      dcic_colnm <- paste0("G = ", 1:Gmax)
      print_matrix(object@DCIC, dcic_rownm, dcic_colnm, colwidth = 8)
    }
)

#' Provide a summary of a \code{dmbc_ic} class instance.
#'
#' @param object An object of class \code{\link{dmbc_ic}}.
#' @param p An optional length-one numeric vector providing the number of
#'   latent space dimension to use in the summary.
#' @param G An optional length-one numeric vector providing the number of
#'   clusters to use in the summary.
#'
#' @author Sergio Venturini \email{sergio.venturini@unito.it}
#'
#' @aliases summary,dmbc_ic-method
#' @aliases dmbc_ic-summary
#' 
#' @exportMethod summary
setMethod("summary",
  "dmbc_ic",
    function(object, p = NULL, G = NULL) {
      pmax <- nrow(object@DCIC)
      Gmax <- ncol(object@DCIC)
      dcic_rownm <- paste0("p = ", 1:pmax)
      dcic_colnm <- paste0("G = ", 1:Gmax)

      if (!is.null(p))
        if (p > pmax) stop("the number of dimensions (p) specified is not available.")
      if (!is.null(G))
        if (G > Gmax) stop("the number of clusters (G) specified is not available.")

      show(object)
      cat("\n")
      cat("Dissimilarity Model Based Clustering -- Logarithm of marginal likelihood values:\n")
      print_matrix(object@logmlik, dcic_rownm, dcic_colnm)
      cat("\n")
      cat("Dissimilarity Model Based Clustering -- Logarithm of prior values:\n")
      print_matrix(object@logprior, dcic_rownm, dcic_colnm)
      cat("\n")
      cat("Dissimilarity Model Based Clustering -- Logarithm of correction factors:\n")
      print_matrix(object@logcorrfact, dcic_rownm, dcic_colnm)

      if (is.null(p) || is.null(G)) {
        G_best <- ifelse((which.min(t(object@DCIC)) %% Gmax) == 0, Gmax, (which.min(t(object@DCIC)) %% Gmax))
        p_best <- which.min(object@DCIC[, G_best])
      } else {
        G_best <- G
        p_best <- p
      }
      best_lab <- paste("p = ", p_best, " -- G = ", G_best, sep = "")
      best_idx <- which(names(object@post.est) == best_lab)
      best_res <- object@post.est[[best_idx]]
      n <- dim(best_res$z.m)[1]
      cat("\n")
      est <- switch(object@est,
        mean = "posterior mean",
        median = "posterior median",
        ml = "maximum likelihood",
        map = "maximum a posteriori"
      )
      cat("Dissimilarity Model Based Clustering -- Parameter estimates (", est, ") -- p = ", p_best,
        " - G = ", G_best, "\n", sep = "")
      wd <- as.integer(options("width"))
      for (g in 1:G_best) {
        dash_string_L <- strrep("-", floor((wd - nchar(paste0("Cluster ", Gmax)) - 2)/2))
        dash_string_R <- strrep("-", floor((wd - nchar(paste0("Cluster ", Gmax)) - 2)/2) - nchar(g) + 1)
        cat(dash_string_L, " ", "Cluster ", g, " ", dash_string_R, "\n", sep = "")

        cat("(1) Z:\n")
        print_matrix(as.matrix(best_res$z.m[, , g]), paste0("n = ", 1:n), paste0("Dimension ", 1:p_best),
          colwidth = 11, ndigits = 4, shift = 7)
        cat("(2) alpha:", round(best_res$alpha.m[g], digits = 4), "\n")
        cat("(3) eta:", round(best_res$eta.m[g], digits = 4), "\n")
        cat("(4) sigma2:", round(best_res$sigma2.m[g], digits = 4), "\n")
        cat("(5) lambda:", round(best_res$lambda.m[g], digits = 4), "\n")
      }
    }
)

#' Provide a graphical summary of a \code{dmbc_ic} class instance.
#'
#' @param x An object of class \code{\link{dmbc_ic}}.
#' @param size A length-two numeric vector providing the optional sizes of
#'   points and lines in the plot.
#' @param ... Further arguments to pass on (currently ignored).
#'
#' @author Sergio Venturini \email{sergio.venturini@unito.it}
#'
#' @aliases plot,dmbc_ic-method
#' @aliases dmbc_ic-plot
#' 
#' @exportMethod plot
setMethod("plot",
	signature(x = "dmbc_ic"),
	function(x, size = NULL, ...) {
    data <- reshape2::melt(data = x@DCIC,
                           varnames = c("p", "G"),
                           value.name = "DCIC")
    data$p <- factor(data$p)
    n_p <- nlevels(data$p)
    geom_args <- list()
    if (is.null(size)) {
      geom_args$size <- c(3, 1)
    } else {
      if (length(size) != 2)
        stop("the 'size' argument must be a numeric vector with two elements.")
      geom_args$size <- size
    }
    # mapping <- ggplot2::aes_(x = ~ G, y = ~ DCIC, color = ~ p, group = ~ p, shape = ~ p)
    mapping <- ggplot2::aes(x = G, y = DCIC, color = p, group = p, shape = p)
    graph <- ggplot2::ggplot(data, mapping) + scale_shape_manual(values = 1:n_p) +
      bayesplot::bayesplot_theme_get()

    graph <- graph + ggplot2::geom_point(size = geom_args$size[1])
    if (n_p > 1) {
      graph <- graph + ggplot2::geom_line(size = geom_args$size[2])
    }

    graph <- graph +
      ggplot2::scale_color_manual("p", values = choose_colors(n_p))

    graph <- graph +
      ggplot2::scale_x_continuous(breaks = 1:ncol(x@DCIC)) +
      bayesplot::legend_move(ifelse(n_p > 1, "right", "none")) +
      bayesplot::xaxis_title(on = TRUE) +
      bayesplot::yaxis_title(on = TRUE)

    graph
	}
)

#' Provide an update of a \code{dmbc_ic} class instance.
#'
#' @param object An object of class \code{\link{dmbc_ic}}.
#' @param pmax A length-one numeric vector indicating the maximum number of
#'   dimensions of the latent space to consider.
#' @param Gmax A length-one numeric vector indicating the maximum number of
#'   cluster to consider.
#' @param ... Further arguments to pass on (currently ignored).
#'
#' @author Sergio Venturini \email{sergio.venturini@unito.it}
#'
#' @aliases update,dmbc_ic-method
#' @aliases dmbc_ic-update
#' 
#' @return A \code{dmbc_ic} object.
#'
#' @seealso \code{\link{dmbc}()} for fitting a DMBC model.
#' @seealso \code{\link{dmbc_ic}} for a description of the elements included
#'   in the returned object.
#'
#' @references
#'   Venturini, S., Piccarreta, R. (2019), "A Bayesian Approach for Model-Based
#'   Clustering of Several Binary Dissimilarity Matrices: the \pkg{dmbc}
#'   Package in \code{R}", Technical report.
#'
#' @examples
#' \dontrun{
#' data(simdiss, package = "dmbc")
#'
#' pmax <- 2
#' Gmax <- 2
#' prm.prop <- list(z = 1.5, alpha = .75)
#' burnin <- 2000
#' nsim <- 1000
#' seed <- 1809
#' 
#' set.seed(seed)
#' 
#' control <- list(burnin = burnin, nsim = nsim, z.prop = prm.prop[["z"]],
#'   alpha.prop = prm.prop[["alpha"]], random.start = TRUE, verbose = TRUE,
#'   thin = 10, store.burnin = TRUE)
#' sim.ic <- dmbc_IC(data = simdiss, pmax = pmax, Gmax = Gmax, control = control,
#'   est = "mean")
#' 
#' pmax <- pmax + 1
#' Gmax <- Gmax + 2
#' new.ic <- update(sim.ic, pmax = pmax, Gmax = Gmax)
#' new.ic
#' 
#' # plot the results
#' library(bayesplot)
#' library(ggplot2)
#' color_scheme_set("mix-yellow-blue")
#' p <- plot(new.ic, size = c(4, 1.5))
#' p + panel_bg(fill = "gray90", color = NA)
#' }
#'
#' @importFrom methods new
#' @importFrom stats4 update
#' @exportMethod update
setMethod("update", "dmbc_ic",
  function(object, pmax = NULL, Gmax = NULL, ...) {
    control <- object@res_last_p[[1]]@results[[1]]@control
    verbose <- control[["verbose"]]
    burnin <- control[["burnin"]]
    nsim <- control[["nsim"]]
    prior <- object@res_last_p[[1]]@results[[1]]@prior

    pmax.old <- nrow(object@DCIC)
    Gmax.old <- ncol(object@DCIC)
    if (is.null(pmax) & is.null(Gmax))
      stop("pmax and Gmax cannot be both null.")
    if (is.null(pmax)) {
      warning("you did not provide a pmax value, so it is set to its previous value.",
        call. = FALSE, immediate. = TRUE)
      pmax <- pmax.old
    }
    if (is.null(Gmax)) {
      warning("you did not provide a Gmax value, so it is set to its previous value.",
        call. = FALSE, immediate. = TRUE)
      Gmax <- Gmax.old
    }
    if ((pmax <= pmax.old) & (Gmax <= Gmax.old))
      stop("at least one of the new values for pmax and Gmax must be larger than the previous ones.")
    if (pmax < pmax.old) {
      warning("the pmax value (", pmax, ") is smaller than the previous one (", pmax.old, ") and hence it is set at that value.", call. = FALSE, immediate. = TRUE)
      pmax <- pmax.old
    }
    if (Gmax < Gmax.old) {
      warning("the Gmax value (", Gmax, ") is smaller than the previous one (", Gmax.old, ") and hence it is set at that value.", call. = FALSE, immediate. = TRUE)
      Gmax <- Gmax.old
    }

    D <- object@res_last_p[[1]]@results[[1]]@diss
    data <- new("dmbc_data", diss = D, n = object@res_last_p[[1]]@results[[1]]@dim[["n"]],
      S = object@res_last_p[[1]]@results[[1]]@dim[["S"]])
    est <- object@est
    logprior <- logmlik <- logcorrfact <- dcic <- matrix(NA, nrow = pmax, ncol = Gmax)
    logprior[1:pmax.old, 1:Gmax.old] <- object@logprior
    logmlik[1:pmax.old, 1:Gmax.old] <- object@logmlik
    logcorrfact[1:pmax.old, 1:Gmax.old] <- object@logcorrfact
    dcic[1:pmax.old, 1:Gmax.old] <- object@DCIC
    res_list <- list()
    res_save <- list()
    res_all <- object@post.est
    res_last_p <- object@res_last_p
    res.i <- pmax.old*Gmax.old + 1
    p.seq <- if (pmax == pmax.old) pmax.old else seq(from = (pmax.old + 1), to = pmax, by = 1)
    G.seq <- if (Gmax == Gmax.old) Gmax.old else seq(from = (Gmax.old + 1), to = Gmax, by = 1)
    
    if (Gmax > Gmax.old) {
      for (G.i in G.seq) {
        for (p.i in 1:pmax.old) {
          if (verbose)
            message("--- p = ", p.i, " -- G = ", G.i, " ---")
          
          prior.i <- update_prior(prior, p.i, G.i)
          res <- dmbc(data = data, p = p.i, G = G.i, control = control, prior = prior.i)
          res_list[[p.i]] <- res
          if (p.i == pmax.old)
            res_last_p[[G.i]] <- res
          
          if (est == "mean") {
            est.tmp <- dmbc_get_postmean(res)
          } else if (est == "median") {
            est.tmp <- dmbc_get_postmedian(res)
          } else if (est == "ml") {
            est.tmp <- dmbc_get_ml(res)
          } else if (est == "map") {
            est.tmp <- dmbc_get_map(res)
          }
          z.m <- est.tmp$z
          alpha.m <- est.tmp$alpha
          eta.m <- est.tmp$eta
          sigma2.m <- est.tmp$sigma2
          lambda.m <- est.tmp$lambda
          class.m <- dmbc_get_postmean(res)$cluster
          res_all[[res.i]] <- list(z.m = z.m, alpha.m = alpha.m, eta.m = eta.m, sigma2.m = sigma2.m,
            lambda.m = lambda.m, class.m = class.m)
          names(res_all)[res.i] <- paste("p = ", p.i, " -- G = ", G.i, sep = "")

          res.i <- res.i + 1
          
          logprior[p.i, G.i] <- log_marg_prior(res, z.m)
          logmlik[p.i, G.i] <- log_marg_lik(res, z.m)
          if (p.i > 1)
            logcorrfact[p.i, G.i] <- log_corr_fact(res_list[[p.i - 1]], z.m)
          dcic[p.i, G.i] <-  -2*(logprior[p.i, G.i] + logmlik[p.i, G.i] +
            ifelse(p.i > 1, sum(logcorrfact[2:p.i, G.i], na.rm = TRUE), 0))
          
          if (verbose) {
            message(" ")
          }
        }
        res_list <- list()
      }
    }

    if (pmax > pmax.old) {
      for (G.i in 1:Gmax) {
        for (p.i in p.seq) {
          if (verbose)
            message("--- p = ", p.i, " -- G = ", G.i, " ---")
          
          prior.i <- update_prior(prior, p.i, G.i)
          res <- dmbc(data = data, p = p.i, G = G.i, control = control, prior = prior.i)
          res_list[[p.i]] <- res
          if (p.i == p.seq[1])
            res_list[[p.i - 1]] <- res_last_p[[G.i]]
          if (p.i == pmax)
            res_save[[G.i]] <- res
          
          if (est == "mean") {
            est.tmp <- dmbc_get_postmean(res)
          } else if (est == "median") {
            est.tmp <- dmbc_get_postmedian(res)
          } else if (est == "ml") {
            est.tmp <- dmbc_get_ml(res)
          } else if (est == "map") {
            est.tmp <- dmbc_get_map(res)
          }
          z.m <- est.tmp$z
          alpha.m <- est.tmp$alpha
          eta.m <- est.tmp$eta
          sigma2.m <- est.tmp$sigma2
          lambda.m <- est.tmp$lambda
          class.m <- dmbc_get_postmean(res)$cluster
          res_all[[res.i]] <- list(z.m = z.m, alpha.m = alpha.m, eta.m = eta.m, sigma2.m = sigma2.m,
            lambda.m = lambda.m, class.m = class.m)
          names(res_all)[res.i] <- paste("p = ", p.i, " -- G = ", G.i, sep = "")

          res.i <- res.i + 1
          
          logprior[p.i, G.i] <- log_marg_prior(res, z.m)
          logmlik[p.i, G.i] <- log_marg_lik(res, z.m)
          logcorrfact[p.i, G.i] <- log_corr_fact(res_list[[p.i - 1]], z.m)
          dcic[p.i, G.i] <-  -2*(logprior[p.i, G.i] + logmlik[p.i, G.i] +
            ifelse(p.i > 1, sum(logcorrfact[2:p.i, G.i], na.rm = TRUE), 0))
          
          if (verbose) {
            message(" ")
          }
        }
        res_list <- list()
      }
    }
    
    out <- new("dmbc_ic",
      logprior = logprior,
      logmlik = logmlik,
      logcorrfact = logcorrfact,
      DCIC = dcic,
      post.est = res_all,
      est = est,
      res_last_p = res_save
    )

    return(out)
  }
)

#' An S4 class to represent the latent configuration estimate for a DMBC model.
#'
#' @description
#'   An S4 class to represent the the latent configuration estimate for a DMBC
#'   model.
#'
#' @slot Z.est An array containing the estimate of the latent
#'   configuration for a DMBC model.
#' @slot Z.sd An array containing the standard deviation of the latent
#'   configuration for a DMBC model.
#' @slot cluster A numeric vector providing the estimated group
#'   membership for the \emph{S} subjects in the data.
#' @slot est A length-one character vector providing the estimate type
#'   returned in \code{Z.est}. Possible values are \code{mean} (posterior
#'   mean), \code{median} (posterior median), \code{ml} (maximum likelihood)
#'   and \code{map} (maximum-a-posteriori).
#' @slot n A length-one numeric vector providing the number of objects.
#' @slot p A length-one numeric vector providing the number of latent
#'   dimensions.
#' @slot S A length-one numeric vector providing the number of subjects.
#' @slot G A length-one numeric vector providing the number of clusters.
#' @slot family An object of class \code{list}; named list with
#'   elements representing the parameter estimates corresponding to different
#'   values of \emph{p} and \emph{G}.
#' @slot chain A length-one numeric vector representing the ID of
#'   the MCMC chain used to compute the estimates.
#' @slot labels A character vector for the (optional) strings to use
#'   in the plots for labeling the objects.
#'
#' @name dmbc_config-class
#' @rdname dmbc_config-class
#' @aliases dmbc_config
#'
#' @references
#'   Venturini, S., Piccarreta, R. (2019), "A Bayesian Approach for Model-Based
#'   Clustering of Several Binary Dissimilarity Matrices: the \pkg{dmbc}
#'   Package in \code{R}", Technical report.
#'
#' @examples
#' showClass("dmbc_config")
#'
#' @exportClass dmbc_config
setClass(Class = "dmbc_config",
  slots = c(
    Z.est = "array",
    Z.sd = "array",
    cluster = "numeric",
    est = "character",
    n = "numeric",
    p = "numeric",
    S = "numeric",
    G = "numeric",
    family = "character",
    chain = "numeric",
    labels = "character"
  )
)

#' Create an instance of the \code{dmbc_config} class using new/initialize.
#'
#' @param .Object Prototype object from the class \code{\link{dmbc_config}}.
#' @param Z.est An array containing the estimate of the latent
#'   configuration for a DMBC model.
#' @param Z.sd An array containing the standard deviation of the latent
#'   configuration for a DMBC model.
#' @param cluster A numeric vector providing the estimated group
#'   membership for the \emph{S} subjects in the data.
#' @param est A length-one character vector providing the estimate type
#'   returned in \code{Z.est}. Possible values are \code{mean} (posterior
#'   mean), \code{median} (posterior median), \code{ml} (maximum likelihood)
#'   and \code{map} (maximum-a-posteriori).
#' @param n A length-one numeric vector providing the number of objects.
#' @param p A length-one numeric vector providing the number of latent
#'   dimensions.
#' @param S A length-one numeric vector providing the number of subjects.
#' @param G A length-one numeric vector providing the number of clusters.
#' @param family An object of class \code{list}; named list with
#'   elements representing the parameter estimates corresponding to different
#'   values of \emph{p} and \emph{G}.
#' @param chain A length-one numeric vector representing the ID of
#'   the MCMC chain used to compute the estimates.
#' @param labels A character vector for the (optional) strings to use
#'   in the plots for labeling the objects.
#'
#' @author Sergio Venturini \email{sergio.venturini@unito.it}
#'
#' @aliases initialize,dmbc_config-method
#' @aliases dmbc_config-initialize
#' 
#' @importFrom methods initialize
#' @exportMethod initialize
setMethod("initialize", "dmbc_config",
  function(
    .Object,
    Z.est = array(),
    Z.sd = array(),
    cluster = numeric(),
    est = character(),
    n = numeric(),
    S = numeric(),
    p = numeric(),
    G = numeric(),
    family = character(),
    chain = numeric(),
    labels = character()
  )
  {
    .Object@Z.est <- Z.est
    .Object@Z.sd <- Z.sd
    .Object@cluster <- cluster
    .Object@est <- est
    .Object@n <- n
    .Object@S <- S
    .Object@p <- p
    .Object@G <- G
    .Object@family <- family
    .Object@chain <- chain
    .Object@labels <- labels
    .Object
  }
)

#' Show an instance of the \code{dmbc_config} class.
#'
#' @param object An object of class \code{\link{dmbc_config}}.
#'
#' @author Sergio Venturini \email{sergio.venturini@unito.it}
#'
#' @aliases show,dmbc_config-method
#' @aliases dmbc_config-show
#' 
#' @importFrom methods show
#' @exportMethod show
setMethod("show",
  "dmbc_config",
    function(object) {
      est <- switch(object@est,
        mean = "posterior mean",
        median = "posterior median",
        ml = "maximum likelihood",
        map = "maximum a posteriori"
      )
      cat("Dissimilarity Model Based Clustering -- Latent configuration estimates\n")
      cat("Number of objects (n):", object@n, "\n")
      cat("Number of latent dimensions (p):", object@p, "\n")
      cat("Number of subjects (S):", object@S, "\n")
      cat("Number of clusters (G):", object@G, "\n")
      cat("Family:", object@family, "\n")
      cat("Chain ID:", object@chain, "\n")
      cat("Estimate type:", est, "\n")
      cat("Cluster sizes:\n")
      cl <- factor(object@cluster, levels = 1:object@G)
      cl_tbl <- table(cl)
      cl_sizes <- as.matrix(trunc(cl_tbl))
      rownames(cl_sizes) <- paste0("Cluster ", levels(cl))
      colnames(cl_sizes) <- "Size"
      print_matrix(cl_sizes, rownames(cl_sizes), colnames(cl_sizes), ndigits = 0, colwidth = 6, isint = TRUE)
      cat("Subjects assigned to each cluster:\n")
      for (g in 1:object@G) {
        subj_clg <- which(object@cluster == g)
        subj_str <- subj_clg[1]
        if (length(subj_clg) > 1) {
          for (h in 2:length(subj_clg)) {
            subj_str <- paste0(subj_str, ", ", subj_clg[h])
          }
        }
        cat("Cluster ", g, ": ", subj_str, "\n", sep = "")
      }
    }
)

#' Provide a summary of a \code{dmbc_config} class instance.
#'
#' @param object An object of class \code{\link{dmbc_config}}.
#'
#' @author Sergio Venturini \email{sergio.venturini@unito.it}
#'
#' @aliases summary,dmbc_config-method
#' @aliases dmbc_config-summary
#' 
#' @exportMethod summary
setMethod("summary",
  "dmbc_config",
    function(object) {
      show(object)
      cat("Latent configuration estimates (Z):\n")

      wd <- as.integer(options("width"))
      for (g in 1:object@G) {
        dash_string_L <- strrep("-", floor((wd - nchar(paste0("Cluster ", object@G)) - 2)/2))
        dash_string_R <- strrep("-", floor((wd - nchar(paste0("Cluster ", object@G)) - 2)/2) - nchar(g) + 1)
        cat(dash_string_L, " ", "Cluster ", g, " ", dash_string_R, "\n", sep = "")
        print_matrix(as.matrix(object@Z.est[, , g]),
          if (length(object@labels)) object@labels else paste0("i = ", 1:object@n),
          paste0("Dimension ", 1:object@p), colwidth = 11, ndigits = 4)
      }
    }
)

#' Provide a graphical summary of a \code{dmbc_config} class instance.
#'
#' @param x An object of class \code{\link{dmbc_config}}.
#' @param size A length-two numeric vector providing the optional sizes of
#'   points and lines in the plot.
#' @param size_lbl A length-one numeric vector providing the size of labels.
#' @param nudge_x A length-one numeric vector providing the optional horizontal
#'   adjustment to nudge labels by.
#' @param nudge_y A length-one numeric vector providing the optional vertical
#'   adjustment to nudge labels by.
#' @param label_objects A length-one logical vector. If \code{TRUE}, labels are
#'   added to the plot.
#' @param ... Further arguments to pass on (currently ignored).
#'
#' @author Sergio Venturini \email{sergio.venturini@unito.it}
#'
#' @aliases plot,dmbc_config-method
#' @aliases dmbc_config-plot
#' 
#' @exportMethod plot
setMethod("plot",
  signature(x = "dmbc_config"),
  function(x, size = NULL, size_lbl = NULL, nudge_x = 0, nudge_y = 0, label_objects = TRUE, ...) {
    n <- x@n
    p <- x@p
    G <- x@G
    data <- prepare_data_to_plot(x)
    n_G <- nlevels(data$G)
    geom_args <- list()
    if (is.null(size)) {
      geom_args$size <- 2
    } else {
      geom_args$size <- size
    }
    if (is.null(size_lbl)) {
      geom_args$size_lbl <- 3
    } else {
      geom_args$size_lbl <- size_lbl
    }

    mapping <- ggplot2::aes_(x = ~ p_1, y = ~ p_2, color = ~ G)
    graph <- ggplot2::ggplot(data, mapping) + bayesplot::bayesplot_theme_get()
    graph <- graph + ggplot2::geom_point(size = geom_args$size)
 
    if (p <= 2) {
      if (p == 1) {
        min_val <- tapply(data[, 1], data$G, min, na.rm = TRUE)*1.15
        max_val <- tapply(data[, 1], data$G, max, na.rm = TRUE)*1.15
        lbl_cl <- unique(data$cl)
        lbl_cl_rep <- margin.table(table(data$G, data$cl), margin = 2)/n
        if (any(lbl_cl_rep > 1)) {
          lbl_cl <- rep(lbl_cl, lbl_cl_rep[match(lbl_cl, names(lbl_cl_rep))])
        }
        lbl_data <- data.frame(min_val = min_val, max_val = max_val, G = factor(1:G, levels = 1:G),
          cl = lbl_cl)
        graph <- graph +
          ggplot2::geom_segment(data = lbl_data,
            mapping = ggplot2::aes_(x = ~ min_val, y = ~ min_val, xend = ~ max_val, yend = ~ max_val),
            arrow = ggplot2::arrow(angle = 45, length = ggplot2::unit(0.05, "npc"), type = "open"), size = .5)
      }
      graph <- graph + ggplot2::facet_wrap(~ G + cl,
        labeller = ggplot2::label_bquote(cols = paste("Cluster ", .(G), ", ", italic(S)[.(G)], " = ", .(cl))),
        scales = "free")
    } else {
      graph <- graph + ggplot2::facet_grid(p_vs + p_i + p_j ~ G + cl,
        labeller = ggplot2::label_bquote(
                    cols = paste("Cluster ", .(G), ", ", italic(S)[.(G)], " = ", .(cl)),
                    rows = paste("Dimension ", .(p_j), " vs. ", "Dimension ", .(p_i))),
        scales = "free")
    }

    if (label_objects) {
      graph <- graph +
        ggrepel::geom_text_repel(ggplot2::aes(label = lbl), size = geom_args$size_lbl,
          segment.size = .25, nudge_x = nudge_x, nudge_y = nudge_y)
    }

    graph <- graph +
      ggplot2::scale_color_manual(G, values = choose_colors(n_G))
    graph <- graph +
      bayesplot::legend_move("none") +
      bayesplot::xaxis_title(on = (p == 2)) +
      bayesplot::yaxis_title(on = (p == 2)) +
      force_axes_in_facets()

    if (p == 2) {
      graph <- graph + ggplot2::xlab(expression(italic(p)[1]))
      graph <- graph + ggplot2::ylab(expression(italic(p)[2]))
    }

    graph
  }
)

#' Extract the final cluster memberships from a \code{dmbc_config} class instance.
#'
#' @param object An object of class \code{\link{dmbc_config}}.
#' @param newdata An object of no explicit specification (currently ignored).
#' @param ... Further arguments to pass on (currently ignored).
#'
#' @author Sergio Venturini \email{sergio.venturini@unito.it}
#'
#' @aliases clusters,dmbc_config-method
#' @aliases dmbc_config-clusters
#' 
#' @importFrom modeltools clusters
#' @exportMethod clusters
setMethod("clusters",
  signature(object = "dmbc_config", newdata = "ANY"),
  function(object, newdata = NULL, ...) {
    object@cluster
  }
)
