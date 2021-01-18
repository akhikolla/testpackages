#' Model-Based Clustering of Several Dissimilarity Matrices.
#'
#' @docType package
#'
#' @name dmbc-package
#' @aliases dmbc-pkg
#' @aliases dmbc-p
#'
#' @description
#' The \pkg{dmbc} package implements a Bayesian algorithm for clustering a set
#' of dissimilarity matrices within a model-based framework. In particular,
#' we consider the case where \emph{S} matrices are available, each
#' describing the dissimilarities among \emph{n} objects, possibly expressed by
#' \emph{S} subjects (judges), or measured under different experimental conditions,
#' or with reference to different characteristics of the objects them- selves.
#' Specifically, we focus on binary dissimilarities, taking values 0 or 1
#' depending on whether or not two objects are deemed as similar, with the goal
#' of analyzing such data using multidimensional scaling (MDS). Differently
#' from the standard MDS algorithms, we are interested in partitioning the
#' dissimilarity matrices into clusters and, simultaneously, to extract a
#' specific MDS configuration for each cluster. The parameter estimates
#' are derived using a hybrid Metropolis-Gibbs Markov Chain Monte Carlo
#' algorithm. We also include a BIC-like criterion for jointly selecting the
#' optimal number of clusters and latent space dimensions.
#'
#' For efficiency reasons, the core computations in the package are implemented
#' using the \code{C} programming language and the \pkg{RcppArmadillo} package.
#'
#' The \pkg{dmbc} package also supports the simulation of multiple chains
#' through the support of the \pkg{parallel} package.
#'
#' Plotting functionalities are imported from the nice \pkg{bayesplot} package.
#' Currently, the package includes methods for binary data only. In future
#' releases routines will be added specifically for continuous (i.e. normal),
#' multinomial and count data.
#'
#' @section \pkg{dmbc} classes:
#' The \pkg{dmbc} package defines the following new classes:
#' \itemize{
#'   \item{\code{\link{dmbc_data}}: }{defines the data to use in a DMBC model.}
#'   \item{\code{\link{dmbc_model}}: }{defines a DMBC model.}
#'   \item{\code{\link{dmbc_fit}}: }{defines the results of a DMBC analysis
#'     for a single MCMC chain.}
#'   \item{\code{\link{dmbc_fit_list}}: }{defines the results of a DMBC analysis
#'     for multiple MCMC chains.}
#'   \item{\code{\link{dmbc_ic}}: }{defines the results of the computation of
#'     the information criterion for a DMBC analysis.}
#'   \item{\code{\link{dmbc_config}}: }{defines the estimate of the latent
#'     configuration for a DMBC analysis.}
#' }
#' The package includes \code{print}, \code{summary} and \code{plot} methods
#'   for each one of these classes.
#'
#' @section Resources:
#' \itemize{
#'  \item{\strong{Bug reports}:}{
#'  If you have noticed a bug that needs to be fixed, please let us know at the
#'  \pkg{dmbc} issue tracker on GitHub:
#'
#'  \url{https://github.com/sergioventurini/dmbc/issues/}.
#'  }
#'  \item{\strong{General questions and help}:}{
#'  To ask a question about \pkg{dmbc} send and email to:
#'
#'  \email{sergio.venturini@unito.it}.
#' }
#' }
#'
#' @seealso \code{\link[bayesplot]{theme_default}} for the default ggplot theme
#'  used by \pkg{bayesplot}.
#' @seealso \code{\link[bayesplot]{bayesplot-colors}} to set or view the color
#'  scheme used for plotting with \pkg{bayesplot}.
#' @seealso \code{\link[ggplot2]{ggsave}} in \pkg{ggplot2} for saving plots.
#'
#' @references
#'   Venturini, S., Piccarreta, R. (2019), "A Bayesian Approach for Model-Based
#'   Clustering of Several Binary Dissimilarity Matrices: the \pkg{dmbc}
#'   Package in \code{R}", Technical report.
#'
#' @useDynLib dmbc
#' @import bayesplot
#' @import coda
#' @import ggplot2
#' @import graphics
#' @import stats
NULL

#' Simulated binary dissimilarity matrices.
#'
#' A dataset containing a list of simulated binary dissimilarity matrices.
#'
#' @usage data(simdiss)
#'
#' @format A \code{\link{dmbc_data}} object whose \code{diss} element is a list
#'   of 10 binary dissimilarity matrices. Each matrix is defined as a \code{dist}
#'   object measuring the agreement among 16 different units.
#'
#' @examples
#' data(simdiss)
#' library(bayesplot)
#' cols <- color_scheme_set("brightblue")
#' plot(simdiss, colors = unlist(cols)[c(1, 6)], font = 1, cex.font = 0.75)
"simdiss"

#' List of binary dissimilarity matrices among 18 animals.
#'
#' To illustrate the MDS analysis of sorting data, Takane et al. (2009) refer
#'   to judgments on the similarity between \emph{n} = 18 animals expressed by
#'   \emph{S} = 20 subjects. Each subject was asked to divide the animals into
#'   as many groups as needed, based on their similarity. We converted these
#'   values to 0 or 1 depending on whether a pair of animals is placed or not
#'   in the same group by a subject.
#'
#' @usage data(animals)
#'
#' @format{
#'   A \code{\link{dmbc_data}} object whose \code{diss} element is a list of 20
#'   binary dissimilarity matrices. Each matrix is defined as a \code{dist}
#'   object measuring whether each pair of the 18 animals has is placed in the
#'   same group (1) or not (0).
#'
#'   The \code{dist} objects have rows and columns that are named as follows:
#'   \describe{
#'     \item{be}{bear}
#'     \item{cm}{camel}
#'     \item{ct}{cat}
#'     \item{cw}{cow}
#'     \item{dg}{dog}
#'     \item{el}{elephant}
#'     \item{gf}{giraffe}
#'     \item{fx}{fox}
#'     \item{hs}{horse}
#'     \item{li}{lion}
#'     \item{mk}{monkey}
#'     \item{ms}{mouse}
#'     \item{pg}{pig}
#'     \item{rb}{rabbit}
#'     \item{sh}{sheep}
#'     \item{sq}{squirrel}
#'     \item{tg}{tiger}
#'     \item{wf}{wolf}
#'   }
#' }
#'
#' @references
#' Takane, Y., Jung, S., Takane, Y. O. (2009). "Multidimensional Scaling". In
#'   Millsap, R. E., Maydeu-Olivares, A. (eds.), The SAGE Handbook of
#'   Quantitative Methods in Psychology, chapter 10, pp. 217–242,.
#'
#' @examples
#' data(animals)
#' library(bayesplot)
#' cols <- color_scheme_set("teal")
#' plot(animals, colors = unlist(cols)[c(1, 6)], font = 1, cex.font = 0.75)
"animals"

#' List of binary dissimilarity matrices among 15 kinship terms.
#'
#' @description{
#' Rosenberg and Kim (1975) designed an experiment to analyze the perceived
#'   similarities of 15 kinship terms.
#'
#'   Here, we consider the data relative to 85 females made available in
#'   Rosenberg (1982). Each subject was asked to group the kinship terms
#'   according to the perceived similarity. Thus, \emph{S} = 85 binary
#'   dissimilarity matrices are available whose elements (0 or 1) indicate
#'   whether or not two kinship terms were grouped together by each individual.
#' }
#'
#' @usage data(kinship)
#'
#' @format{
#'   A \code{\link{dmbc_data}} object whose \code{diss} element is a list of 85
#'   binary dissimilarity matrices. Each matrix is defined as a \code{dist}
#'   object measuring whether each pair of the 15 kinship terms is judged as
#'   similar (1) or not (0).
#'
#'   The \code{dist} objects have rows and columns that are named as follows:
#'   \describe{
#'     \item{GrF}{grandfather}
#'     \item{GrM}{grandmother}
#'     \item{GrD}{granddaughter}
#'     \item{GrS}{grandson}
#'     \item{Bro}{brother}
#'     \item{Sis}{sister}
#'     \item{Fat}{father}
#'     \item{Mot}{mother}
#'     \item{Dau}{daughter}
#'     \item{Son}{son}
#'     \item{Nep}{nephew}
#'     \item{Nie}{niece}
#'     \item{Cou}{cousin}
#'     \item{Aun}{aunt}
#'     \item{Unc}{uncle}
#'   }
#' }
#'
#' @references{
#'   Rosenberg, S. (1982). The method of sorting in multivariate research with
#'     applications selected from cognitive psychology and person perception. In
#'     N Hirschberg, LG Humphreys (eds.), Multivariate Applications in the Social
#'     Sciences, pp. 117–142. Erlbaum., Hillsdale, NJ.
#'
#'   Rosenberg, S., Kim, M. P. (1975). The method of sorting as a data-gathering
#'     procedure in multivariate research. Multivariate Behavioral Research, 10.
#' }
#'
#' @examples
#' data(kinship)
#' library(bayesplot)
#' cols <- color_scheme_set("mix-red-blue")
#' plot(kinship, colors = unlist(cols)[c(1, 6)], font = 1, cex.font = 0.75)
"kinship"
