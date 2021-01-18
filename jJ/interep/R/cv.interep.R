#' @useDynLib interep, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#' k-folds cross-validation for interep
#'
#' This function does k-fold cross-validation for interep and returns the optimal value of lambda.
#' @param e matrix of environment factors.
#' @param g matrix of omics factors. In the case study, the omics measurements are lipidomics data.
#' @param y the longitudinal response.
#' @param beta0 the intial value for the coefficient vector.
#' @param lambda1 a user-supplied sequence of \eqn{\lambda_{1}} values, which serves as a tuning parameter for individual predictors.
#' @param lambda2 a user-supplied sequence of \eqn{\lambda_{2}} values, which serves as a tuning parameter for interactions.
#' @param nfolds the number of folds for cross-validation.
#' @param corre the working correlation structure that is used in the estimation algorithm. interep provides three choices for the
#' working correlation structure: "a" as AR-1", "i" as "independence" and "e" as "exchangeable".
#' @param pmethod the penalization method. "mixed" refers to MCP penalty to individual main effects and group MCP penalty to interactions; "individual" means MCP penalty to all effects.
#' @param maxits the maximum number of iterations that is used in the estimation algorithm.
#' @details
#' When dealing with predictors with both main effects and interactions, this function returns two optimal tuning parameters,
#' \eqn{\lambda_{1}} and \eqn{\lambda_{2}}; when there are only main effects in the predictors, this function returns \eqn{\lambda_{1}},
#' which is the optimal tuning parameter for individual predictors containing main effects.
#' @return an object of class "cv.interep" is returned, which is a list with components:
#' \item{lam1}{the optimal \eqn{\lambda_{1}}.}
#' \item{lam2}{the optimal \eqn{\lambda_{2}}.}
#' @references
#' Zhou, F., Ren, J., Li, G., Jiang, Y., Li, X., Wang, W.and Wu, C. (2019). Penalized variable selection for Lipid--environment interactions in a longitudinal lipidomics study.
#' \href{https://www.mdpi.com/2073-4425/10/12/1002/htm}{\emph{Genes}, 10(12), 1002}
#'
#' Zhou, F., Ren, J., Lu, X., Ma, S. and Wu, C. (2020) Gene–Environment Interaction: a Variable Selection Perspective.
#' \href{https://arxiv.org/abs/2003.02930}{\emph{Epistasis}, Methods in Molecular Biology. Humana Press. (Accepted)}
#'
#' Ren, J., Zhou, F., Li, X., Chen, Q., Zhang, H., Ma, S., Jiang,Y. and Wu, C. (2020). Semi-parametric Bayesian variable selection for Gene-Environment interactions.
#' \href{https://doi.org/10.1002/sim.8434}{\emph{Statistics in Medicine}, 39(5): 617–638}
#'
#' Wu, C., Zhou, F., Ren, J., Li, X., Jiang, Y., Ma, S. (2019). A Selective Review of Multi-Level Omics Data Integration Using Variable Selection.
#' \href{https://doi.org/10.3390/ht8010004}{\emph{High-Throughput}, 8(1)}
#'
#' Ren, J., Du, Y., Li, S., Ma, S., Jiang, Y. and Wu, C. (2019). Robust network-based regularization and variable selection for high-dimensional genomic data in cancer prognosis.
#' \href{https://doi.org/10.1002/gepi.22194}{\emph{Genetic epidemiology}, 43(3), 276-291}
#'
#' Ren, J., Jung, L., Du, Y., Wu, C., Jiang, Y. and Liu, J. (2019). regnet: Network-Based Regularization for Generalized Linear Models.
#' \href{https://cran.r-project.org/package=regnet}{\emph{R package}, version 0.4.0}
#'
#' Wu, C., Zhang, Q., Jiang, Y. and Ma, S. (2018). Robust network-based analysis of the associations between (epi) genetic measurements.
#' \href{https://doi.org/10.1016/j.jmva.2018.06.009}{\emph{Journal of multivariate analysis}, 168, 119-130}
#'
#' Wu, C., Zhong, P.-S., and Cui, Y. (2018). Additive varying-coefficient model for nonlinear gene-environment interactions.
#' \href{https://doi.org/10.1515/sagmb-2017-0008}{\emph{ Statistical Applications in Genetics and Molecular Biology}, 17(2)}
#'
#' Wu, C., Jiang, Y., Ren, J., Cui, Y., Ma, S. (2018). Dissecting gene-environment interactions: A penalized robust approach accounting for hierarchical structures.
#' \href{https://doi.org/10.1002/sim.7518}{\emph{Statistics in Medicine}, 37:437–456}
#'
#' Ren, J., He, T., Li, Y., Liu, S., Du, Y., Jiang, Y. and Wu, C. (2017). Network-based regularization for high dimensional SNP data in the case–control study of Type 2 diabetes.
#' \href{https://doi.org/10.1186/s12863-017-0495-5}{\emph{BMC genetics}, 18(1), 44}
#'
#' Jiang, Y., Huang, Y., Du, Y., Zhao, Y., Ren, J., Ma, S., & Wu, C. (2017). Identification of prognostic genes and pathways in lung adenocarcinoma using a Bayesian approach.
#' \href{https://www.researchgate.net/profile/Cen_Wu/publication/310828910_Identification_of_Prognostic_Genes_and_Pathways_in_Lung_Adenocarcinoma_Using_a_Bayesian_Approach/links/5cfbdac692851c874c5947f6/Identification-of-Prognostic-Genes-and-Pathways-in-Lung-Adenocarcinoma-Using-a-Bayesian-Approach.pdf}{\emph{Cancer Inform}, 1(7)}
#'
#' Wu, C., and Ma, S. (2015). A selective review of robust variable selection with applications in bioinformatics.
#' \href{https://doi.org/10.1093/bib/bbu046}{\emph{Briefings in Bioinformatics}, 16(5), 873–883}
#'
#' Wu, C., Shi, X., Cui, Y. and Ma, S. (2015). A penalized robust semiparametric approach for gene-environment interactions.
#' \href{https://doi.org/10.1002/sim.6609}{\emph{Statistics in Medicine}, 34 (30): 4016–4030}
#'
#' Wu, C., Cui, Y., and Ma, S. (2014). Integrative analysis of gene–environment interactions under a multi–response partially linear varying coefficient model.
#' \href{https://doi.org/10.1002/sim.6287}{\emph{Statistics in Medicine}, 33(28), 4988–4998}
#'
#' Wu, C. and Cui, Y. (2013). A novel method for identifying nonlinear gene–environment interactions in case–control association studies.
#' \href{https://doi.org/10.1007/s00439-013-1350-z}{\emph{Human Genetics}, 132(12):1413–1425}
#'
#' Wu, C. and Cui, Y. (2013). Boosting signals in gene–based association studies via efficient SNP selection.
#' \href{https://doi.org/10.1093/bib/bbs087}{\emph{Briefings in Bioinformatics}, 15(2):279–291}
#'
#' Wu, C., Zhong, P.S. and Cui, Y. (2013). High dimensional variable selection for gene-environment interactions.
#' \href{https://www.stt.msu.edu/Links/Research_Memoranda/RM/RM_701.pdf}{\emph{Technical Report}, Michigan State University.}
#'
#' Wu, C., Li, S., and Cui, Y. (2012). Genetic Association Studies: An Information Content Perspective.
#' \href{https://doi.org/10.2174/138920212803251382}{\emph{Current Genomics}, 13(7),  566–573}
#'
#' @export



cv.interep <- function(e, g, y, beta0, lambda1, lambda2, nfolds, corre, pmethod, maxits){
  n=dim(y)[1]
  q=dim(e)[2]
  p1=dim(g)[2]
  len1=length(lambda1)
  len2=length(lambda2)

  folds=rep(1,n/nfolds)
  for (i in 2:nfolds) {
    folds=c(folds,rep(i,n/nfolds))
  }
  folds=sample(folds)
  pred=matrix(0,len1,len2)
  for (i in 1:len1) {
    lam1=lambda1[i]
    for (j in 1:len2){
      lam2=lambda2[j]
      mse=0
      for (cv in 1:nfolds) {
        #Segement your data by fold using the which() function
        testIndexes <- which(folds==cv,arr.ind=TRUE)
        e.test=e[testIndexes,]
        g.test=g[testIndexes,]
        y.test=y[testIndexes,]
        e.train=e[-testIndexes,]
        g.train=g[-testIndexes,]
        y.train=y[-testIndexes,]
        e.train=as.matrix(e.train)
        g.train=as.matrix(g.train)
        y.train=as.matrix(y.train)
        beta=interep(e.train, g.train, y.train, beta0,corre,pmethod,lam1,lam2,maxits)
        x.test=cbind(e.test,g.test)
        for (i1 in 1:p1) {
          for (j1 in 1:q) {
            x.test=cbind(x.test,e.test[,j1]*g.test[,i1])
          }
        }
        x.test=scale(x.test)
        data.test=reformat(y.test, x.test)
        x.test=data.test$x
        y.test=data.test$y
        mu=x.test%*%beta
        mse=mse+mean((y.test-mu)^2)
      }
      pred[i,j]=mse/nfolds
    }
  }
  lamb1=lambda1[which(pred==min(pred),arr.ind = TRUE)[1]]
  lamb2=lambda2[which(pred==min(pred),arr.ind = TRUE)[2]]
  return(list("lam1"=lamb1,"lam2"=lamb2))
}
