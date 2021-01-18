#' @useDynLib springer, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#' @docType package
#' @keywords overview
#' @name springer-package
#' @title Sparse Group Variable Selection for Gene-Environment Interactions in the Longitudinal Study
#' @aliases springer-package
#' @description
#' In this package, we provide a set of regularized variable selection methods tailored for longitudinal studies of gene- environment interactions. The proposed method conducts sparse group variable selection by accounting for bi-level sparsity. Specifically, the individual and group level penalties have been simultaneously imposed to identify important main and interaction effects under three working correlation structures (exchangeable , AR-1 and independence), based on either the quadratic inference function (QIF) or generalized estimating equation (GEE). In addition, only the individual or group level selection in the longitudinal setting can also be conducted using springer. In total, springer provides 18 (=3\eqn{\times}3\eqn{\times}2) methods. Among them, sparse group variable selection for longitudinal studies have been developed for the first time. Please read the Details below for how to configure the method used.
#'
#' @details Users can flexibly choose the methods to fit the model by specifying the three arguments in the user interface \strong{springer()}:
#' \tabular{rl}{
#' func: \tab the framework to obtain the score equation.  Two choices are available: \cr \tab "GEE" and "QIF". \cr\cr
#' corr: \tab working correlation.  Three choices are available: \cr \tab "exchangeable", "AR-1" and "independence". \cr\cr
#' structure: \tab structural identification. Three choices are available: \cr \tab "bilevel", "group" and “individual”.
#' }
#'
#' The function springer() returns a springer object that contains the estimated coefficients.
#'
#' @references
#'
#' Zhou, F., Lu, X., Ren, J., Ma, S. and Wu, C. (2020+). Sparse Group Variable Selection for Gene-Environment Interactions in the Longitudinal Study.
#'
#' Zhou, F., Ren, J., Lu, X., Ma, S. and Wu, C. (2020). Gene–Environment Interaction: a Variable Selection Perspective. Epistasis. Methods in Molecular Biology.
#' {\emph{Humana Press} (Accepted)} \url{https://arxiv.org/abs/2003.02930}
#'
#' Zhou, F., Ren, J.,  Li, G., Jiang, Y., Li, X., Wang, W. and Wu, C. (2019). Penalized Variable Selection for Lipid–Environment Interactions in a Longitudinal Lipidomics Study.
#' {\emph{Genes}, 10(12), 1002} \url{https://doi.org/10.3390/genes10121002}
#'
#' Zhou, F., Ren, J., Li, X., Wu, C. and Jiang, Y. (2019) interep: Interaction Analysis of Repeated Measure Data.
#' R package version 0.3.1. \url{https://CRAN.R-project.org/package=interep}
#'
#' Ren, J., Du, Y., Li, S., Ma, S., Jiang, Y. and Wu, C. (2019). Robust network-based regularization and variable selection for high-dimensional genomic data in cancer prognosis.
#' {\emph{Genetic epidemiology}, 43(3), 276-291} \url{https://doi.org/10.1002/gepi.22194}
#'
#' Wu, C., Zhang, Q., Jiang, Y. and Ma, S. (2018). Robust network-based analysis of the associations between (epi) genetic measurements.
#' {\emph{Journal of multivariate analysis}, 168, 119-130} \url{https://doi.org/10.1016/j.jmva.2018.06.009}
#'
#' Wu, C., Jiang, Y., Ren, J., Cui, Y. and Ma, S. (2018). Dissecting gene-environment interactions: A penalized robust approach accounting for hierarchical structures.
#' {\emph{Statistics in Medicine}, 37:437–456} \url{https://doi.org/10.1002/sim.7518}
#'
#' Wu, C., Zhong, P.S. and Cui, Y. (2018). Additive varying-coefficient model for nonlinear gene-environment interactions.
#' {\emph{Statistical Applications in Genetics and Molecular Biology}, 17(2)} \url{https://doi.org/10.1515/sagmb-2017-0008}
#'
#' Ren, J., He, T., Li, Y., Liu, S., Du, Y., Jiang, Y. and Wu, C. (2017). Network-based regularization for high dimensional SNP data in the case–control study of Type 2 diabetes.
#' \href{https://doi.org/10.1186/s12863-017-0495-5}{\emph{BMC genetics}, 18(1), 44}
#'
#' Wu, C., Shi, X., Cui, Y. and Ma, S. (2015). A penalized robust semiparametric approach for gene-environment interactions.
#' {\emph{Statistics in Medicine}, 34 (30): 4016–4030} \url{https://doi.org/10.1002/sim.6609}
#'
#' Wu, C., Cui, Y., and Ma, S. (2014). Integrative analysis of gene–environment interactions under a multi–response partially linear varying coefficient model.
#' {\emph{Statistics in Medicine}, 33(28), 4988–4998} \url{https://doi.org/10.1002/sim.6287}
#'
#' Wu, C. and Cui, Y. (2014). Boosting signals in gene-based association studies via efficient SNP selection.
#' {\emph{Briefings in bioinformatics}, 15(2), 279-291} \url{https://doi.org/10.1093/bib/bbs087}
#'
#' Wu, C., Zhong, P.S. and Cui, Y. (2013). High dimensional variable selection for gene-environment interactions.
#' {\emph{Technical Report. Michigan State University.}}
#'
#' @seealso \code{\link{springer}}
NULL
