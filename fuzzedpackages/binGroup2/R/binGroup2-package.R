#' @title binGroup2: Identification and Estimation using Group Testing
#' 
#' @description Methods for the group testing identification and 
#' estimation problems.
#' 
#' @details
#' Methods for identification of positive items in group testing designs: 
#' Operating characteristics (e.g., expected number of tests) are calculated 
#' for commonly used hierarchical and array-based algorithms. Optimal testing 
#' configurations for an algorithm can be found as well. Please see Hitt et al. 
#' (2019) for specific details.
#' 
#' Methods for estimation and inference for proportions in group testing 
#' designs: For estimating one proportion or the difference of proportions, 
#' confidence interval methods are included that account for different pool 
#' sizes. Functions for hypothesis testing of proportions, calculation of 
#' power, and calculation of the expected width of confidence intervals are 
#' also included. Furthermore, regression methods and simulation of group 
#' testing data are implemented for simple pooling, halving, and array testing 
#' designs. 
#' 
#' The \code{binGroup2} package is based upon the \code{binGroup} package that 
#' was originally designed for the group testing estimation problem. Over time, 
#' additional functions for estimation and for the group testing identification 
#' problem were included. Due to the diverse styles resulting from these 
#' additions, we have created \code{binGroup2} as a way to unify functions in 
#' a coherent structure and incorporate additional functions for identification. 
#' The name “binGroup” originates from the assumption in basic estimation for 
#' group testing that the number of positive groups has a binomial distribution. 
#' While more advanced estimation methods no longer make this assumption, we 
#' continue with the \code{binGroup} name for consistency.
#' 
#' Bilder (2019a,b) provide introductions to group testing. These papers and 
#' additional details about group testing are available at 
#' \url{http://chrisbilder.com/grouptesting}.
#' 
#' This research was supported by the National Institutes of Health under 
#' grant R01 AI121351.
#' 
#' \subsection{Identification}{
#' The binGroup2 package focuses on the group testing identification problem 
#' using hierarchical and array-based group testing algorithms.
#' 
#' The \code{\link{OTC1}} function implements a number of group testing 
#' algorithms, described in Hitt et al. (2019), which calculate the operating 
#' characteristics and find the optimal testing configuration over a range of 
#' possible initial group sizes and/or testing configurations (sets of 
#' subsequent group sizes). The \code{\link{OTC2}} function does the same with 
#' a multiplex assay that tests for two diseases.
#' 
#' The \code{\link{operatingCharacteristics1}} (\code{\link{opChar1}}) and 
#' \code{\link{operatingCharacteristics2}} (\code{\link{opChar2}}) functions 
#' calculate operating characteristics for a specified testing configuration 
#' with assays that test for one and two diseases, respectively.
#' 
#' These functions allow the sensitivity and specificity to differ across 
#' stages of testing. This means that the accuracy of the diagnostic test can 
#' differ for stages in a hierarchical testing algorithm or between 
#' row/column testing and individual testing in an array testing algorithm.}
#' 
#' \subsection{Estimation}{
#' The binGroup2 package also provides functions for estimation and 
#' inference for proportions in group testing designs.
#' 
#' The \code{\link{propCI}} function calculates the point estimate and 
#' confidence intervals for a single proportion from group testing data. 
#' The \code{\link{propDiffCI}} function does the same for the difference of 
#' proportions. A number of confidence interval methods are available for 
#' groups of equal or different sizes.
#' 
#' The \code{\link{gtWidth}} function calculates the expected width of 
#' confidence intervals in group testing. The \code{\link{gtTest}} function 
#' calculates p-values for hypothesis tests of single proportions. The 
#' \code{\link{gtPower}} function calculates power to reject a hypothesis.
#' 
#' The \code{\link{designPower}} function iterates either the number of groups 
#' or group size in a one-parameter group testing design until a pre-specified 
#' power level is achieved. The \code{\link{designEst}} function finds the 
#' optimal group size corresponding to the minimal mean-squared error of the 
#' point estimator.
#' 
#' The \code{\link{gtReg}} function implements regression methods and the 
#' \code{\link{gtSim}} function simulates group testing data for simple 
#' pooling, halving, and array testing designs.}
#' 
#' @references 
#' \insertRef{Altman1994a}{binGroup2}
#' 
#' \insertRef{Altman1994b}{binGroup2}
#' 
#' \insertRef{Biggerstaff2008}{binGroup2}
#' 
#' \insertRef{Bilder2010a}{binGroup2}
#' 
#' \insertRef{Bilder2019}{binGroup2}
#' 
#' \insertRef{Bilder2019est}{binGroup2}
#' 
#' \insertRef{Bilder2019id}{binGroup2}
#' 
#' \insertRef{Black2012}{binGroup2}
#' 
#' \insertRef{Black2015}{binGroup2}
#' 
#' \insertRef{Graff1972}{binGroup2}
#' 
#' \insertRef{Hepworth1996}{binGroup2}
#' 
#' \insertRef{Hepworth2017}{binGroup2}
#' 
#' \insertRef{Hitt2019}{binGroup2}
#' 
#' \insertRef{Hou2019}{binGroup2}
#' 
#' \insertRef{Malinovsky2016}{binGroup2}
#' 
#' \insertRef{McMahan2012a}{binGroup2}
#' 
#' \insertRef{McMahan2012b}{binGroup2}
#' 
#' \insertRef{Schaarschmidt2007}{binGroup2}
#' 
#' \insertRef{Swallow1985}{binGroup2}
#' 
#' \insertRef{Tebbs2004}{binGroup2}
#' 
#' \insertRef{Vansteelandt2000}{binGroup2}
#' 
#' \insertRef{Xie2001}{binGroup2}
#' 
#' @examples 
#' # Estimated running time for all examples was calculated 
#' #   using a computer with 16 GB of RAM and one core of 
#' #   an Intel i7-6500U processor. Please take this into 
#' #   account when interpreting the run times given.
#' 
#' # 1) Identification using hierarchical and array-based group testing 
#' #   algorithms with an assay that tests for one disease.
#' 
#' # 1.1) Find the optimal testing configuration over a range of initial 
#' #   group sizes, using informative three-stage hierarchical testing, where 
#' #   p denotes the overall prevalence of disease;
#' #   Se denotes the sensitivity of the diagnostic test; 
#' #   Sp denotes the specificity of the diagnostic test;
#' #   group.sz denotes the range of initial pool sizes for consideration; and
#' #   obj.fn specifies the objective functions for which to find results.
#' 
#' # This example takes approximately 25 seconds to run.
#' \donttest{
#' set.seed(1002)
#' results1 <- OTC1(algorithm="ID3", p=0.01, Se=0.95, Sp=0.95, 
#'                  group.sz=3:30, obj.fn=c("ET", "MAR"), alpha=2)
#' summary(results1)}
#' 
#' # 1.2) Find the optimal testing configuration using non-informative
#' # array testing without master pooling.
#' # The sensitivity and specificity differ for row/column testing and 
#' #   individual testing.
#' 
#' # This example takes approximately 15 seconds to run.
#' \donttest{
#' results2 <- OTC1(algorithm="A2", p=0.05, Se=c(0.95, 0.99), 
#'                  Sp=c(0.95, 0.98), group.sz=3:20, obj.fn=c("ET", "MAR"))
#' summary(results2)}
#' 
#' # 1.3) Calculate the operating characteristics using informative
#' #   two-stage hierarchical (Dorfman) testing, implemented via the 
#' #   pool-specific optimal Dorfman (PSOD) method described in 
#' #   McMahan et al. (2012a).
#' # Hierarchical testing configurations are specified by a matrix 
#' #   in the hier.config argument. The rows of the matrix correspond 
#' #   to the stages of the hierarchical testing algorithm, the columns 
#' #   correspond to the individuals to be tested, and the cell values 
#' #   correspond to the group number of each individual at each stage.
#' config.mat <- matrix(data=c(rep(1, 5), rep(2, 4), 3, 1:10), 
#'                      nrow=2, ncol=10, byrow=TRUE)
#' set.seed(8791)
#' results3 <- opChar1(algorithm="ID2", p=0.02, Se=0.95, Sp=0.99, 
#'                     hier.config=config.mat, alpha=0.5)
#' summary(results3)
#' 
#' # 1.4) Calculate the operating characteristics using non-informative
#' #   four-stage hierarchical testing. 
#' config.mat <- matrix(data=c(rep(1, 15), rep(c(1, 2, 3), each=5), 
#'                             rep(1, 3), rep(2, 2), rep(3, 3), rep(4, 2), 
#'                             rep(5, 4), 6, 1:15), 
#'                      nrow=4, ncol=15, byrow=TRUE)
#' results4 <- opChar1(algorithm="D4", p=0.008, Se=0.96, Sp=0.98, 
#'                     hier.config=config.mat, a=c(1, 4, 6, 9, 11, 15))
#' summary(results4)
#' 
#' 
#' # 2) Identification using hierarchical and array-based group testing 
#' #   algorithms with a multiplex assay that tests for two diseases.
#' 
#' # 2.1) Find the optimal testing configuration using non-informative 
#' #   two-stage hierarchical testing, given
#' #   p.vec, a vector of overall joint probabilities of disease; 
#' #   Se, a vector of sensitivity values for each disease; and 
#' #   Sp, a vector of specificity values for each disease. 
#' # Se and Sp can also be specified as a matrix, where one value 
#' #   is specified for each disease at each stage of testing.
#' results5 <- OTC2(algorithm="D2", p.vec=c(0.90, 0.04, 0.04, 0.02), 
#'                  Se=c(0.99, 0.99), Sp=c(0.99, 0.99), group.sz=3:50)
#' summary(results5)
#' 
#' # 2.2) Calculate the operating characteristics for informative
#' #   five-stage hierarchical testing, given
#' #   alpha.vec, a vector of shape parameters for the Dirichlet distribution; 
#' #   Se, a matrix of sensitivity values; and 
#' #   Sp, a matrix of specificity values.
#' Se <- matrix(data=rep(0.95, 10), nrow=2, ncol=5, byrow=TRUE)
#' Sp <- matrix(data=rep(0.99, 10), nrow=2, ncol=5, byrow=TRUE)
#' config.mat <- matrix(data=c(rep(1, 24), rep(1, 18), rep(2, 6), 
#'                             rep(1, 9), rep(2, 9), rep(3, 4), 4, 5, 
#'                             rep(1, 6), rep(2, 3), rep(3, 5), rep(4, 4), 
#'                             rep(5, 3), 6, rep(NA, 2), 1:21, rep(NA, 3)), 
#'                      nrow=5, ncol=24, byrow=TRUE)
#' results6 <- opChar2(algorithm="ID5", alpha=c(18.25, 0.75, 0.75, 0.25),
#'                     Se=Se, Sp=Sp, hier.config=config.mat)
#' summary(results6)
#' 
#' # 3) Estimation of the overall disease prevalence and calculation 
#' #   of confidence intervals.
#' 
#' # 3.1) Suppose 3 groups out of 24 test positively. 
#' #   Each group has a size of 7.
#' propCI(x=3, m=7, n=24, ci.method="CP")
#' propCI(x=3, m=7, n=24, ci.method="Blaker")
#' propCI(x=3, m=7, n=24, ci.method="score")
#' propCI(x=3, m=7, n=24, ci.method="soc")
#' 
#' # 3.2) Consider the following situation:
#' #   0 out of 5 groups test positively with groups 
#' #   of size 1 (individual testing), 
#' #   0 out of 5 groups test positively with groups of size 5,
#' #   1 out of 5 groups test positively with groups of size 10, 
#' #   2 out of 5 groups test positively with groups of size 50
#' propCI(x=c(0,0,1,2), m=c(1,5,10,50), n=c(5,5,5,5), 
#'        pt.method="Gart", ci.method="skew-score")
#'        
#' # 4) Estimate a group testing regression model.
#' 
#' # 4.1) Fit a group testing regression model with 
#' #   simple pooling using the "hivsurv" dataset.
#' data(hivsurv)
#' fit1 <- gtReg(type="sp", formula = groupres ~ AGE + EDUC., 
#'               data = hivsurv, groupn = gnum, sens = 0.9, 
#'               spec = 0.9, method = "Xie")
#' summary(fit1)
#' 
#' # 4.2) Simulate data for the halving protocol, and 
#' #   fit a group testing regression model.
#' set.seed(46)
#' gt.data <- gtSim(type="halving", par=c(-6, 0.1), 
#'                  gshape=17, gscale=1.4, size1=1000, 
#'                  size2=5, sens=0.95, spec=0.95)
#' fit2 <- gtReg(type="halving", formula=gres~x, 
#'               data=gt.data, groupn=groupn, subg=subgroup,
#'               retest=retest, sens=0.95, spec=0.95, 
#'               start=c(-6, 0.1), trace=TRUE)
#' summary(fit2)
#' 
#' # This example takes approximately 20 seconds to run.
#' # 4.3) Simulate data in 5x6 array testing form, and 
#' #   fit a group testing regression model.
#' set.seed(9128)
#' array.sim <- gtSim(type="array", par=c(-7, 0.1), 
#'                    size1=c(5,6), size2=c(4,5), sens=0.95, spec=0.95)
#' set1 <- array.sim$dframe
#' \donttest{
#' fit3 <- gtReg(type="array", 
#'               formula=cbind(col.resp, row.resp)~x, 
#'               data=set1, coln=coln, rown=rown, 
#'               arrayn=arrayn, sens=0.95, spec=0.95, 
#'               tol=0.005, n.gibbs=2000, trace=TRUE)
#' summary(fit3)}
#' 
#' @docType package
#' @name binGroup2
"_PACKAGE"
