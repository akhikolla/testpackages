#' @useDynLib roben, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#' @docType package
#' @keywords overview
#' @name roben-package
#' @title Robust Bayesian Variable Selection for Gene-Environment Interactions
#' @aliases roben-package
#' @description In this package, we provide a set of robust Bayesian variable selection methods tailored for interaction analysis. A Bayesian formulation of the least absolute deviation (LAD) regression has been adopted to accommodate data contamination and long-tailed distributions in the response/ phenotype. The default method (the proposed method) conducts variable selection by accounting for structural sparsity. In particular, the spike--and--slab priors are imposed on both individual and group levels to identify important main and interaction effects (bi-level/ sparse-group selection).
#'
#' In addition to the default method, users can also choose different selection structures (group-level-only or individual-level-only), methods without spike--and--slab priors and non-robust methods. In total, \emph{roben} provides 12 different methods (6 robust and 6 non-robust). Among them, robust methods with spike--and--slab priors and the robust method for bi-level selection have been developed for the first time. Please read the Details below for how to configure the method used.
#'
#' @details The user friendly, integrated interface \strong{roben()} allows users to flexibly choose the fitting methods they prefer. There are three arguments in roben() that control the fitting method:
#' \tabular{rl}{
#' robust: \tab whether to use robust methods. \cr\cr
#' sparse: \tab whether to use the spike-and-slab priors to create sparsity. \cr\cr
#' structure: \tab structural identification. Three choices are available: \cr \tab "sparsegroup", "group" and “individual”.
#' }
#'
#' The function roben() returns a roben object that contains the posterior estimates of each coefficients.
#' S3 generic functions GxESelection(), predict() and print() are implemented for roben objects.
#' GxESelection() takes a roben object and returns the variable selection results.
#' predict() takes a roben object and returns the predicted values for new observations.
#'
#' @references
#' Ren, J., Zhou, F., Li, X., Ma, S., Jiang, Y. and Wu, C. (2020). Robust Bayesian variable selection for gene-environment interactions.
#'
#' Wu, C., and Ma, S. (2015). A selective review of robust variable selection with applications in bioinformatics.
#' {\emph{Briefings in Bioinformatics}, 16(5), 873–883} \url{https://doi.org/10.1093/bib/bbu046}
#'
#' Zhou, F., Ren, J., Lu, X., Ma, S. and Wu, C. (2020). Gene–Environment Interaction: a Variable Selection Perspective. Epistasis. Methods in Molecular Biology.
#' {\emph{Humana Press} (Accepted)} \url{https://arxiv.org/abs/2003.02930}
#'
#' Ren, J., Zhou, F., Li, X., Chen, Q., Zhang, H., Ma, S., Jiang, Y. and Wu, C. (2020) Semi-parametric Bayesian variable selection for gene-environment interactions.
#' {\emph{Statistics in Medicine}, 39: 617– 638} \url{https://doi.org/10.1002/sim.8434}
#'
#' Ren, J., Zhou, F., Li, X., Wu, C. and Jiang, Y. (2019) spinBayes: Semi-Parametric Gene-Environment Interaction via Bayesian Variable Selection.
#' R package version 0.1.0. \url{https://CRAN.R-project.org/package=spinBayes}
#'
#' Wu, C., Jiang, Y., Ren, J., Cui, Y. and Ma, S. (2018). Dissecting gene-environment interactions: A penalized robust approach accounting for hierarchical structures.
#' {\emph{Statistics in Medicine}, 37:437–456} \url{https://doi.org/10.1002/sim.7518}
#'
#' Wu, C., Shi, X., Cui, Y. and Ma, S. (2015). A penalized robust semiparametric approach for gene-environment interactions.
#' {\emph{Statistics in Medicine}, 34 (30): 4016–4030} \url{https://doi.org/10.1002/sim.6609}
#'
#' Wu, C., Cui, Y., and Ma, S. (2014). Integrative analysis of gene–environment interactions under a multi–response partially linear varying coefficient model.
#' {\emph{Statistics in Medicine}, 33(28), 4988–4998} \url{https://doi.org/10.1002/sim.6287}
#'
#' Wu, C., Zhong, P.S. and Cui, Y. (2018). Additive varying–coefficient model for nonlinear gene–environment interactions.
#' {\emph{Statistical Applications in Genetics and Molecular Biology}, 17(2)} \url{https://doi.org/10.1515/sagmb-2017-0008}
#'
#' Wu, C., Zhong, P.S. and Cui, Y. (2013). High dimensional variable selection for gene-environment interactions.
#' {\emph{Technical Report. Michigan State University.}}
#'
#' @seealso \code{\link{roben}}
NULL
