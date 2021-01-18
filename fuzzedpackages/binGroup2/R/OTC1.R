# Start OTC1() function
###############################################################################
# Updated: Brianna Hitt - 12-06-19 (warning message text)
# Updated: Brianna Hitt - 01-22-20 
#   Removed the printing of the algorithm name
#   Removed the printing of all "messages" other than warnings/errors

#' @title Find the optimal testing configuration for group testing algorithms 
#' that use a single-disease assay
#'
#' @description Find the optimal testing configuration (OTC) using 
#' non-informative and informative hierarchical and array-based group testing 
#' algorithms. Single-disease assays are used at each stage of the algorithms.
#'
#' @param algorithm character string defining the group testing algorithm to be used.
#' Non-informative testing options include two-stage hierarchical ("\kbd{D2}"),
#' three-stage hierarchical ("\kbd{D3}"), square array testing without master
#' pooling ("\kbd{A2}"), and square array testing with master pooling ("\kbd{A2M}").
#' Informative testing options include two-stage hierarchical ("\kbd{ID2}"),
#' three-stage hierarchical ("\kbd{ID3}"), and square array testing without
#' master pooling ("\kbd{IA2}").
#' @param p overall probability of disease that will be used to generate a
#' vector/matrix of individual probabilities. For non-informative algorithms, a
#' homogeneous set of probabilities will be used. For informative algorithms, the
#' \code{\link{expectOrderBeta}} function will be used to generate a heterogeneous
#' set of probabilities. Further details are given under 'Details'. Either 
#' \kbd{p} or \kbd{probabilities} should be specified, but not both.
#' @param probabilities a vector of individual probabilities, which is homogeneous
#' for non-informative testing algorithms and heterogeneous for informative
#' testing algorithms. Either  \kbd{p} or \kbd{probabilities} should be specified,
#' but not both.
#' @param Se a vector of sensitivity values, where one value is given for each 
#' stage of testing (in order). If a single value is provided, sensitivity 
#' values are assumed to be equal to this value for all stages of testing. 
#' Further details are given under 'Details'.
#' @param Sp a vector of specificity values, where one value is given for each 
#' stage of testing (in order). If a single value is provided, specificity 
#' values are assumed to be equal to this value for all stages of testing. 
#' Further details are given under 'Details'.
#' @param group.sz a single group size or range of group sizes for which to
#' calculate operating characteristics and/or find the OTC. The details of group
#' size specification are given under 'Details'.
#' @param obj.fn a list of objective functions which are minimized to find the
#' OTC. The expected number of tests per individual, "\kbd{ET}", will always
#' be calculated. Additional options include "\kbd{MAR}"
#' (the expected number of tests divided by the expected number of correct
#' classifications, described in Malinovsky et al. (2016)), and "\kbd{GR}"
#' (a linear combination of the expected number of tests, the number of
#' misclassified negatives, and the number of misclassified positives,
#' described in Graff & Roeloffs (1972)). See Hitt et al. (2019) for 
#' additional details. The first objective function specified in this list 
#' will be used to determine the results for the top configurations. 
#' Further details are given under 'Details'.
#' @param weights a matrix of up to six sets of weights for the GR function.
#' Each set of weights is specified by a row of the matrix.
#' @param alpha a shape parameter for the beta distribution that specifies the 
#' degree of heterogeneity for the generated probability vector (for 
#' informative testing only).
#' @param trace a logical value indicating whether the progress of 
#' calculations should be printed for each initial group size provided by 
#' the user. The default is \kbd{TRUE}. 
#' @param print.time a logical value indicating whether the length of time 
#' for calculations should be printed. The default is \kbd{TRUE}.
#' @param ... arguments to be passed to the \code{\link{expectOrderBeta}} 
#' function, which generates a vector of probabilities for informative testing 
#' algorithms. Further details are given under 'Details'.
#'
#' @details This function finds the OTC for group testing algorithms 
#' with an assay that tests for one disease and computes the associated 
#' operating characteristics, as described in Hitt et al. (2019).
#' 
#' Available algorithms include two- and three-stage hierarchical testing and
#' array testing with and without master pooling. Both non-informative and informative
#' group testing settings are allowed for each algorithm, except informative
#' array testing with master pooling is unavailable because this method has not
#' appeared in the group testing literature. Operating characteristics calculated are
#' expected number of tests, pooling sensitivity, pooling specificity, pooling
#' positive predictive value, and pooling negative predictive value for each individual.
#' 
#' For informative algorithms where the \kbd{p} argument is specified, the 
#' expected value of order statistics from a beta distribution are found. 
#' These values are used to represent disease risk probabilities for each 
#' individual to be tested. The beta distribution has two parameters: a mean 
#' parameter \kbd{p} (overall disease prevalence) and a shape parameter 
#' \kbd{alpha} (heterogeneity level). Depending on the specified \kbd{p}, 
#' \kbd{alpha}, and overall group size, simulation may be necessary to 
#' generate the vector of individual probabilities. This is done using 
#' \code{\link{expectOrderBeta}} and requires the user to set a seed to 
#' reproduce results. 
#' 
#' Informative two-stage hierarchical (Dorfman) testing is implemented via 
#' the pool-specific optimal Dorfman (PSOD) method described in McMahan et al. 
#' (2012a), where the greedy algorithm proposed for PSOD is replaced by 
#' considering all possible testing configurations. Informative array testing 
#' is implemented via the gradient method (the most efficient array design), 
#' where higher-risk individuals are grouped in the left-most columns of the 
#' array. For additional details on the gradient arrangement method for 
#' informative array testing, see McMahan et al. (2012b).
#' 
#' The sensitivity/specificity values are allowed to vary across stages of 
#' testing. For hierarchical testing, a different sensitivity/specificity value 
#' may be used for each stage of testing. For array testing, a different 
#' sensitivity/specificity value may be used for master pool testing (if included), 
#' row/column testing, and individual testing. The values must be specified 
#' in order of the testing performed. For example, values are specified 
#' as (stage 1, stage 2, stage 3) for three-stage hierarchical testing or 
#' (master pool testing, row/column testing, individual testing) for array 
#' testing with master pooling. A single sensitivity/specificity value may be 
#' specified instead. In this situation, sensitivity/specificity values for all 
#' stages are assumed to be equal.
#'
#' The value(s) specified by \kbd{group.sz} represent the initial (stage 1)
#' group size for hierarchical testing and the row/column size for array 
#' testing. For informative two-stage hierarchical testing,
#' the \kbd{group.sz} specified represents the block size used in the pool-specific
#' optimal Dorfman (PSOD) method, where the initial group (block) is not
#' tested. For more details on informative two-stage hierarchical testing
#' implemented via the PSOD method, see Hitt et al. (2019) and McMahan et al. (2012a).
#'
#' If a single value is provided for \kbd{group.sz} with array testing or
#' non-informative two-stage hierarchical testing, operating
#' characteristics will be calculated and no optimization will be performed.
#' If a single value is provided for \kbd{group.sz} with three-stage hierarchical or
#' informative two-stage hierarchical, the OTC will be
#' found over all possible configurations. If a range of group sizes is specified,
#' the OTC will be found over all group sizes.
#'
#' In addition to the OTC, operating characteristics for some of the other
#' configurations corresponding to each initial group size provided by the user
#' will be displayed. These additional configurations are only determined for whichever
#' objective function ("ET", "MAR", or "GR") is specified first in the
#' function call. If "GR" is the objective function listed first, the
#' first set of corresponding weights will be used.
#' For algorithms where there is only one configuration for each
#' initial group size (non-informative two-stage hierarchical and all array testing
#' algorithms), results for each initial group size are provided. For algorithms where
#' there is more than one possible configuration for each initial group size (informative
#' two-stage hierarchical and all three-stage hierarchical algorithms), two sets of
#' configurations are provided: 1) the best configuration for each initial group size,
#' and 2) the top 10 configurations for each initial group size provided by the user.
#' If a single value is provided for \kbd{group.sz} with array testing or
#' non-informative two-stage hierarchical testing, operating characteristics will
#' not be provided for configurations other than that specified by the user. 
#' Results are sorted by the value of the objective function per individual, \kbd{value}.
#'
#' The displayed overall pooling sensitivity, pooling specificity, pooling positive
#' predictive value, and pooling negative predictive value are weighted
#' averages of the corresponding individual accuracy measures for all
#' individuals within the initial group (or block) for a hierarchical algorithm, or
#' within the entire array for an array-based algorithm.
#' Expressions for these averages are provided in the Supplementary
#' Material for Hitt et al. (2019). These expressions are based on accuracy
#' definitions given by Altman and Bland (1994a, 1994b). Individual 
#' accuracy measures can be calculated using the 
#' \code{\link{operatingCharacteristics1}} (\code{\link{opChar1}}) function.
#'
#' The \kbd{OTC1} function accepts additional arguments, namely \kbd{num.sim}, 
#' to be passed to the \code{\link{expectOrderBeta}} function, which generates 
#' a vector of probabilities for informative group testing algorithms. The 
#' \kbd{num.sim} argument specifies the number of simulations from the beta 
#' distribution when simulation is used. By default, 10,000 simulations are used.
#'
#' @return A list containing:
#' \item{algorithm}{the group testing algorithm used for calculations.}
#' \item{prob}{the probability of disease or the vector of individual 
#' probabilities, as specified by the user.}
#' \item{alpha}{level of heterogeneity for the generated probability vector
#' (for informative testing only).}
#' \item{Se}{the vector of sensitivity values for each stage of testing.}
#' \item{Sp}{the vector of specificity values for each stage of testing.}
#' \item{opt.ET, opt.MAR, opt.GR}{a list of results for each
#' objective function specified by the user, containing:
#' \describe{
#' \item{OTC}{a list specifying elements of the optimal testing configuration,
#' which may include:
#' \describe{
#' \item{Stage1}{group size for the first stage of hierarchical testing, if applicable.}
#' \item{Stage2}{group sizes for the second stage of hierarchical testing, if applicable.}
#' \item{Block.sz}{the block size/initial group size for informative Dorfman testing,
#' which is not tested.}
#' \item{pool.szs}{group sizes for the first stage of testing for informative Dorfman
#' testing.}
#' \item{Array.dim}{the row/column size for array testing.}
#' \item{Array.sz}{the overall array size for array testing (the square of the 
#' row/column size).}}}
#' \item{p.vec}{the sorted vector of individual probabilities, if applicable.}
#' \item{p.mat}{the sorted matrix of individual probabilities in gradient arrangement,
#' if applicable. Further details are given under 'Details'.}
#' \item{ET}{the expected testing expenditure to decode all individuals in the algorithm;
#' this includes all individuals in all groups for hierarchical algorithms or in the
#' entire array for array testing.}
#' \item{value}{the value of the objective function per individual.}
#' \item{Accuracy}{a matrix of overall accuracy measures for the 
#' algorithm. The columns correspond to the pooling sensitivity, 
#' pooling specificity, pooling positive predictive value, and 
#' pooling negative predictive value for the overall algorithm. 
#' Further details are given under 'Details'.}}}
#' \item{Configs}{a data frame containing results for the best configuration 
#' for each initial group size provided by the user. The columns correspond to 
#' the initial group size, configuration (if applicable), overall array size 
#' (if applicable), expected number of tests, value of the objective function 
#' per individual, pooling sensitivity, pooling specificity, pooling positive 
#' predictive value, and pooling negative predictive value. No results are 
#' displayed if a single \kbd{group.sz} is provided. Further details are given 
#' under 'Details'.}
#' \item{Top.Configs}{a data frame containing results for some of the top 
#' configurations for each initial group size provided by the user. The 
#' columns correspond to the initial group size, configuration, 
#' expected number of tests, value of the objective function per individual, 
#' pooling sensitivity, pooling specificity, pooling positive predictive 
#' value, and pooling negative predictive value. No results are displayed for 
#' non-informative two-stage hierarchical testing or for array testing 
#' algorithms. Further details are given under 'Details'.}
#'
#' @section Note: This function returns the pooling positive and negative 
#' predictive values for all individuals even though these measures are 
#' diagnostic specific; e.g., the pooling positive predictive value should 
#' only be considered for those individuals who have tested positive.
#' 
#' Additionally, only stage dependent sensitivity and specificity values are 
#' allowed within the program (no group within stage dependent values are 
#' allowed). See Bilder et al. (2019) for additional information.
#' 
#' @author Brianna D. Hitt
#'
#' @references
#' \insertRef{Altman1994a}{binGroup2}
#'
#' \insertRef{Altman1994b}{binGroup2}
#' 
#' \insertRef{Bilder2019}{binGroup2}
#'
#' \insertRef{Graff1972}{binGroup2}
#'
#' \insertRef{Hitt2019}{binGroup2}
#'
#' \insertRef{Malinovsky2016}{binGroup2}
#'
#' \insertRef{McMahan2012a}{binGroup2}
#'
#' \insertRef{McMahan2012b}{binGroup2}
#'
#' @family OTC functions
#'
#' @examples
#' # Estimated running time for all examples was calculated 
#' #   using a computer with 16 GB of RAM and one core of 
#' #   an Intel i7-6500U processor. Please take this into 
#' #   account when interpreting the run times given.
#' 
#' # Find the OTC for non-informative
#' #   two-stage hierarchical (Dorfman) testing.
#' OTC1(algorithm="D2", p=0.05, Se=0.99, Sp=0.99, 
#'      group.sz=3:100, obj.fn=c("ET", "MAR"), 
#'      trace=TRUE, print.time=TRUE)
#'
#' # Find the OTC for informative two-stage hierarchical 
#' #   (Dorfman) testing.
#' # A vector of individual probabilities is generated using
#' #   the expected value of order statistics from a beta
#' #   distribution with p = 0.01 and a heterogeneity level
#' #   of alpha = 0.5.
#' # This example takes approximately 2.5 minutes to run.
#' \donttest{
#' set.seed(52613)
#' OTC1(algorithm="ID2", p=0.01, Se=0.95, Sp=0.95, group.sz=50,
#'      obj.fn=c("ET", "MAR", "GR"),
#'      weights=matrix(data=c(1, 1, 10, 10, 0.5, 0.5),
#'      nrow=3, ncol=2, byrow=TRUE), alpha=0.5, 
#'      trace=FALSE, print.time=TRUE, num.sim=10000)}
#'
#' # Find the OTC over all possible testing configurations 
#' #   for non-informative three-stage hierarchical testing 
#' #   with a specified group size.
#' OTC1(algorithm="D3", p=0.001, Se=0.95, Sp=0.95, group.sz=18,
#'      obj.fn=c("ET", "MAR", "GR"),
#'      weights=matrix(data=c(1, 1), nrow=1, ncol=2, byrow=TRUE), 
#'      trace=FALSE, print.time=FALSE)
#'
#' # Find the OTC for non-informative three-stage 
#' #   hierarchical testing.
#' # This example takes approximately 20 seconds to run.
#' \donttest{
#' OTC1(algorithm="D3", p=0.06, Se=0.90, Sp=0.90,
#'      group.sz=3:30, obj.fn=c("ET", "MAR", "GR"),
#'      weights=matrix(data=c(1, 1, 10, 10, 100, 100),
#'      nrow=3, ncol=2, byrow=TRUE))}
#'
#' # Find the OTC over all possible configurations
#' #   for informative three-stage hierarchical testing 
#' #   with a specified group size and a heterogeneous 
#' #   vector of probabilities.
#' set.seed(1234)
#' OTC1(algorithm="ID3", 
#'      probabilities=c(0.012, 0.014, 0.011, 0.012, 0.010, 0.015), 
#'      Se=0.99, Sp=0.99, group.sz=6, obj.fn=c("ET","MAR","GR"), 
#'      weights=matrix(data=c(1, 1), nrow=1, ncol=2, byrow=TRUE), 
#'      alpha=0.5, num.sim=5000, trace=FALSE)
#'
#' # Calculate the operating characteristics for 
#' #   non-informative array testing without master pooling 
#' #   with a specified array size.
#' OTC1(algorithm="A2", p=0.005, Se=0.95, Sp=0.95, group.sz=8,
#'      obj.fn=c("ET", "MAR"), trace=FALSE)
#'
#' # Find the OTC for informative array testing without
#' #   master pooling.
#' # A vector of individual probabilities is generated using
#' #   the expected value of order statistics from a beta
#' #   distribution with p = 0.03 and a heterogeneity level
#' #   of alpha = 2. The probabilities are then arranged in
#' #   a matrix using the gradient method.
#' # This example takes approximately 30 seconds to run.
#' \donttest{
#' set.seed(1002)
#' OTC1(algorithm="IA2", p=0.03, Se=0.95, Sp=0.95,
#'      group.sz=3:20, obj.fn=c("ET", "MAR", "GR"),
#'      weights=matrix(data=c(1, 1, 10, 10, 100, 100), 
#'                     nrow=3, ncol=2, byrow=TRUE), alpha=2)}
#'
#' # Find the OTC for non-informative array testing
#' #   with master pooling.
#' # This example takes approximately 20 seconds to run.
#' \donttest{
#' OTC1(algorithm="A2M", p=0.02, Se=0.90, Sp=0.90,
#'      group.sz=3:20, obj.fn=c("ET", "MAR", "GR"),
#'      weights=matrix(data=c(1, 1, 10, 10, 0.5, 0.5, 2, 2, 100, 100, 
#'                            10, 100), nrow=6, ncol=2, byrow=TRUE))}

# Brianna Hitt - 04.02.2020
# Changed cat() to warning()

OTC1 <- function(algorithm, p=NULL, probabilities=NULL, Se=0.99, Sp=0.99, 
                 group.sz, obj.fn=c("ET","MAR"), weights=NULL, alpha=2, 
                 trace=TRUE, print.time=TRUE, ...){
  
  ## make sure that all necessary information is included in the correct format
  if (!(algorithm %in% c("D2", "D3", "A2", "A2M", "ID2", "ID3", "IA2"))) {
    stop("Please specify one of the following algorithms: D2, ID2, D3, ID3, A2, IA2, A2M.")
  }
  
  if(is.null(p) & is.null(probabilities)){
    stop("Please specify an overall probability of disease using the 'p' argument, or specify a vector of individual probabilities using the 'probabilities' argument.")
  } else if(!is.null(p) & !is.null(probabilities)){
    stop("You have specified both an overall probability of disease AND a vector of individual probabilities. Please specify only one option.")
  } else{
    if(!is.null(p) & length(p)>1){
      stop("You have specified a probability vector instead of an overall probability of disease. Please specify an overall probability of disease, and the probability vector will be generated based on the algorithm specified for each group size included in the range.\n")
    }
    if(!is.null(probabilities)){
      if(length(group.sz)==1){
        if((algorithm %in% c("D2", "D3", "ID2", "ID3")) & length(probabilities)!=group.sz){
          stop("The vector of individual probabilities is not the correct length. Please make sure that the length of the probability vector is the same as the specified group size.\n")
        } else if((algorithm %in% c("A2", "A2M", "IA2")) & length(probabilities)!=group.sz^2){
          stop("The vector of individual probabilities is not the correct length. Please make sure that the length of the probability vector is the same as the overall array size (the square of the specified row/column size).\n")
        }
        if((algorithm %in% c("D2", "D3", "A2", "A2M")) & 
           all.equal(probabilities, rep(probabilities[1],length(probabilities)))!=TRUE){
          stop("You have specified a heterogeneous probability vector for a non-informative algorithm. Please specify a homogeneous probability vector using the 'probabilities' argument or specify an overall probability of disease using the 'p' argument.\n")
        }
      } else if(length(group.sz)>1){
        stop("You have specified a probability vector along with a range of group sizes. Please specify a single group size.\n")
      }
    }
  }
  
  Se <- generate.acc(algorithm=algorithm, diseases=1, value=Se, label="sens")
  Sp <- generate.acc(algorithm=algorithm, diseases=1, value=Sp, label="spec")
  
  # check the minimum and maximum group sizes
  if (min(group.sz) < 3) {
    if (algorithm %in% c("D2", "D3", "ID2", "ID3")) {
      stop("Please specify a minimum group size of at least 3.\n")
    } else if (algorithm %in% c("A2", "IA2", "A2M")) {
      stop("Please specify a minimum row/column size of at least 3.\n")
    }
  }
  if(max(group.sz)>=50){
    if(algorithm %in% c("D3", "ID2", "ID3")){
      warning("You have specified a maximum group size of 50 or larger. This function may take a VERY long time to run. Press 'ESC' if you wish to cancel the submitted statements.\n")
    } else if(algorithm %in% c("A2", "A2M", "IA2")){
      warning("You have specified a maximum row/column size of 50 or larger. This function may take a VERY long time to run. Press 'ESC' if you wish to cancel the submitted statements.\n")
    }
  }
  
  if(is.null(obj.fn)){
    stop("Please specify one or more objective functions for which to find the optimal testing configuration.\n")
  }
  
  if (!("ET" %in% obj.fn)) {
    obj.fn <- c(obj.fn, "ET")
  }
  
  if("GR" %in% obj.fn){
    if(is.null(weights)){
      stop("No weights have been specified. The GR function will not be calculated.\n")
    } else if(dim(weights)[2]!=2){
      stop("Please check the dimension of the weights matrix. Each row should specify a set of weights, D1 and D2.\n")
    }
  }
  
  # call function for non-informative two-stage hierarchical (Dorfman) testing
  if(algorithm == "D2"){
    if(!is.null(p)){
      results <- NI.Dorf.OTC1(p=p, Se=Se, Sp=Sp, group.sz=group.sz, 
                              obj.fn=obj.fn, weights=weights, 
                              trace=trace, print.time=print.time, ...)
    } else if(!is.null(probabilities)){
      results <- NI.Dorf.OTC1(p=probabilities, Se=Se, Sp=Sp, group.sz=group.sz, 
                              obj.fn=obj.fn, weights=weights, 
                              trace=trace, print.time=print.time, ...)
    }
  }
  
  # call function for non-informative three-stage hierarchical testing
  if(algorithm == "D3"){
    if(!is.null(p)){
      results <- NI.D3.OTC1(p=p, Se=Se, Sp=Sp, group.sz=group.sz, 
                            obj.fn=obj.fn, weights=weights, 
                            trace=trace, print.time=print.time, ...)
    } else if(!is.null(probabilities)){
      results <- NI.D3.OTC1(p=probabilities, Se=Se, Sp=Sp, group.sz=group.sz, 
                            obj.fn=obj.fn, weights=weights, 
                            trace=trace, print.time=print.time, ...)
    }
  }
  
  # call function for non-informative square array testing without master pooling
  if(algorithm == "A2"){
    if(!is.null(p)){
      results <- NI.Array.OTC1(p=p, Se=Se, Sp=Sp, group.sz=group.sz, 
                               obj.fn=obj.fn, weights=weights, 
                               trace=trace, print.time=print.time, ...)
    } else if(!is.null(probabilities)){
      results <- NI.Array.OTC1(p=probabilities, Se=Se, Sp=Sp, group.sz=group.sz, 
                               obj.fn=obj.fn, weights=weights, 
                               trace=trace, print.time=print.time, ...)
    }
  }
  
  # call function for non-informative square array testing with master pooling
  if(algorithm == "A2M"){
    if(!is.null(p)){
      results <- NI.A2M.OTC1(p=p, Se=Se, Sp=Sp, group.sz=group.sz, 
                             obj.fn=obj.fn, weights=weights, 
                             trace=trace, print.time=print.time, ...)
    } else if(!is.null(probabilities)){
      results <- NI.A2M.OTC1(p=probabilities, Se=Se, Sp=Sp, group.sz=group.sz, 
                             obj.fn=obj.fn, weights=weights, 
                             trace=trace, print.time=print.time, ...)
    }
  }
  
  # call function for informative two-stage hierarchical (Dorfman) testing
  if(algorithm == "ID2"){
    if(!is.null(p)){
      results <- Inf.Dorf.OTC1(p=p, Se=Se, Sp=Sp, group.sz=group.sz, 
                               obj.fn=obj.fn, weights=weights, alpha=alpha, 
                               trace=trace, print.time=print.time, ...)
    } else if(!is.null(probabilities)){
      results <- Inf.Dorf.OTC1(p=probabilities, Se=Se, Sp=Sp, group.sz=group.sz, 
                               obj.fn=obj.fn, weights=weights, alpha=alpha, 
                               trace=trace, print.time=print.time, ...)
    }
  }
  
  # call function for informative three-stage hierarchical testing
  if(algorithm == "ID3"){
    if(!is.null(p)){
      results <- Inf.D3.OTC1(p=p, Se=Se, Sp=Sp, group.sz=group.sz, 
                             obj.fn=obj.fn, weights=weights, alpha=alpha, 
                             trace=trace, print.time=print.time, ...)
    } else if(!is.null(probabilities)){
      results <- Inf.D3.OTC1(p=probabilities, Se=Se, Sp=Sp, group.sz=group.sz, 
                             obj.fn=obj.fn, weights=weights, alpha=alpha, 
                             trace=trace, print.time=print.time, ...)
    }
  }
  
  # call function for informative square array testing without master pooling
  if(algorithm == "IA2"){
    if(!is.null(p)){
      results <- Inf.Array.OTC1(p=p, Se=Se, Sp=Sp, group.sz=group.sz, 
                                obj.fn=obj.fn, weights=weights, alpha=alpha, 
                                trace=trace, print.time=print.time, ...)
    } else if(!is.null(probabilities)){
      results <- Inf.Array.OTC1(p=probabilities, Se=Se, Sp=Sp, group.sz=group.sz, 
                                obj.fn=obj.fn, weights=weights, alpha=alpha, 
                                trace=trace, print.time=print.time, ...)
    }
  }
  
  class(results) <- "OTC"
  results
}



# Summary function for OTC1() and OTC2()
###############################################################################
#' @title Summary method for optimal testing configuration results
#' 
#' @description Produce a summary list for objects of class \kbd{"OTC"} 
#' returned by \code{\link{OTC1}} or \code{\link{OTC2}}.
#' 
#' @param object an object of class \kbd{"OTC"}, providing the optimal testing 
#' configuration and associated operating characteristics for a group testing 
#' algorithm.
#' @param ... currently not used.
#' 
#' @details This function produces a summary list for objects of class 
#' \kbd{"OTC"} returned by \code{\link{OTC1}} or \code{\link{OTC2}}. 
#' It formats the optimal testing configuration, expected number of tests, 
#' expected number of tests per individual, and accuracy measures. 
#' A summary of the results from \code{\link{OTC1}} includes results for all 
#' objective functions specified by the user.
#' 
#' The \kbd{OTC} component of the result gives the optimal testing 
#' configuration, which may include the group sizes for each stage of a 
#' hierarchical testing algorithm or the row/column size and array size for an 
#' array testing algorithm. The \kbd{Tests} component of the result gives the 
#' expected number of tests and the expected number of tests per individual 
#' for the algorithm. 
#' 
#' The \kbd{Accuracy} component gives the overall accuracy measures for the 
#' algorithm. Accuracy measures included are the pooling sensitivity, pooling 
#' specificity, pooling positive predictive value, and pooling negative 
#' predictive value. These values are weighted averages of the corresponding 
#' individual accuracy measures for all individuals in the algorithm. 
#' Expressions for these averages are provided in the Supplementary Material 
#' for Hitt et al. (2019). For more information, see the 'Details' section for 
#' the \code{\link{OTC1}} or \code{\link{OTC2}} function.
#' 
#' @return \kbd{summary.OTC} returns an object of class \kbd{"summary.OTC"}, 
#' a list containing:
#' \item{Algorithm}{character string specifying the name of the group testing 
#' algorithm.}
#' \item{OTC}{matrix detailing the optimal testing configuration from 
#' \kbd{object}. For hierarchical testing, this includes the group sizes for 
#' each stage of testing. For array testing, this includes the array dimension 
#' (row/column size) and the array size (the total number of individuals 
#' in the array).}
#' \item{Tests}{matrix detailing the expected number of tests and expected 
#' number of tests per individual from \kbd{object}}.
#' \item{Accuracy}{matrix detailing the overall accuracy measures for the 
#' algorithm, including the pooling sensitivity, pooling specificity, 
#' pooling positive predictive value, and pooling negative predictive value 
#' for the algorithm from \kbd{object}. Further details are found in the 
#' 'Details' section.}
#' 
#' @author Brianna D. Hitt
#' 
#' @seealso
#' \code{\link{OTC1}} and \code{\link{OTC2}} 
#' for creating an object of class \kbd{"OTC"}.
#' 
#' @examples 
#' # Estimated running time for all examples was calculated 
#' #   using a computer with 16 GB of RAM and one core of 
#' #   an Intel i7-6500U processor. Please take this into 
#' #   account when interpreting the run time given.
#' 
#' # Find the optimal testing configuration for 
#' #   non-informative two-stage hierarchical testing.
#' res1 <- OTC1(algorithm="D2", p=0.01, Se=0.99, Sp=0.99, 
#'              group.sz=3:100, obj.fn=c("ET", "MAR", "GR1"), 
#'              weights=matrix(data=c(1,1), nrow=1, ncol=2))
#' summary(res1)
#' 
#' # Find the optimal testing configuration for 
#' #   informative three-stage hierarchical testing
#' res2 <- OTC1(algorithm="ID3", p=0.025, 
#'              Se=c(0.95, 0.95, 0.99), Sp=c(0.96, 0.96, 0.98), 
#'              group.sz=3:15, obj.fn=c("ET", "MAR"), alpha=2)
#' summary(res2)
#' 
#' # Find the optimal testing configuration for 
#' #   informative array testing without master pooling.
#' # This example takes approximately 30 seconds to run.
#' \donttest{
#' res3 <- OTC1(algorithm="IA2", p=0.05, alpha=2, 
#'              Se=0.90, Sp=0.90, group.sz=3:20, obj.fn="ET")
#' summary(res3)}
#' 
#' # Find the optimal testing configuraiton for 
#' #   informative two-stage hierarchical testing.
#' Se <- matrix(data=c(rep(0.95, 2), rep(0.99, 2)), 
#'              nrow=2, ncol=2, byrow=FALSE)
#' Sp <- matrix(data=c(rep(0.96, 2), rep(0.98, 2)), 
#'              nrow=2, ncol=2, byrow=FALSE)
#' res4 <- OTC2(algorithm="ID2", alpha=c(18.25, 0.75, 0.75, 0.25), 
#'                 Se=Se, Sp=Sp, group.sz=12)
#' summary(res4)
#' 
#' # Find the optimal testing configuration for 
#' #   non-informative three-stage hierarchical testing.
#' # This example takes approximately 3 minutes to run.
#' Se <- matrix(data=c(rep(0.95, 6)), nrow=2, ncol=3)
#' Sp <- matrix(data=c(rep(0.99, 6)), nrow=2, ncol=3)
#' \donttest{
#' res5 <- OTC2(algorithm="D3", p.vec=c(0.95, 0.0275, 0.0175, 0.005), 
#'                 Se=Se, Sp=Sp, group.sz=5:10)
#' summary(res5)}
#' 
#' # Find the optimal testing configuration for
#' #   non-informative array testing with master pooling.
#' # This example takes approximately 10 seconds to run.
#' \donttest{
#' res6 <- OTC2(algorithm="A2M", p.vec=c(0.90, 0.04, 0.04, 0.02), 
#'              Se=rep(0.99, 2), Sp=rep(0.99, 2), group.sz=3:15)
#' summary(res6)}

summary.OTC <- function(object, ...){
  
  # algorithm
  algorithm <- object$algorithm
  cat("\nAlgorithm:", algorithm, "\n\n")
  
  # extract the results for all objective functions
  #   some of these may be NULL if not requested by the user
  opt.ET <- object$opt.ET
  opt.MAR <- object$opt.MAR
  opt.GR1 <- object$opt.GR1
  opt.GR2 <- object$opt.GR2
  opt.GR3 <- object$opt.GR3
  opt.GR4 <- object$opt.GR4
  opt.GR5 <- object$opt.GR5
  opt.GR6 <- object$opt.GR6
  
  all.objfns <- c("ET", "MAR", "GR1", "GR2", "GR3", "GR4", "GR5", "GR6")
  objfn.labels <- all.objfns[paste0("opt.", all.objfns) %in% names(object)]
  
  # create a matrix detailing the optimal configurations for each obj. fn
  stage1 <- c(opt.ET$OTC[[1]], opt.MAR$OTC[[1]], opt.GR1$OTC[[1]], 
              opt.GR2$OTC[[1]], opt.GR3$OTC[[1]], opt.GR4$OTC[[1]], 
              opt.GR5$OTC[[1]], opt.GR6$OTC[[1]])
  
  if (length(opt.ET$OTC) > 1) {
    if (length(opt.ET$OTC[[2]]) > 1) {
      # concatenate pool sizes for ET objective function
      stage2.ET <- opt.ET$OTC[[2]][1]
      for(i in 2:length(opt.ET$OTC[[2]])){
        stage2.ET <- paste(stage2.ET, opt.ET$OTC[[2]][i], sep=",")
      }
      # concatenate pool sizes for MAR objective function
      stage2.MAR <- opt.MAR$OTC[[2]][1]
      for(i in 2:length(opt.MAR$OTC[[2]])){
        stage2.MAR <- paste(stage2.MAR, opt.MAR$OTC[[2]][i], sep=",")
      }
      # concatenate pool sizes for GR1 objective function
      stage2.GR1 <- opt.GR1$OTC[[2]][1]
      for(i in 2:length(opt.GR1$OTC[[2]])){
        stage2.GR1 <- paste(stage2.GR1, opt.GR1$OTC[[2]][i], sep=",")
      }
      # concatenate pool sizes for GR2 objective function
      stage2.GR2 <- opt.GR2$OTC[[2]][1]
      for(i in 2:length(opt.GR2$OTC[[2]])){
        stage2.GR2 <- paste(stage2.GR2, opt.GR2$OTC[[2]][i], sep=",")
      }
      # concatenate pool sizes for GR3 objective function
      stage2.GR3 <- opt.GR3$OTC[[2]][1]
      for(i in 2:length(opt.GR3$OTC[[2]])){
        stage2.GR3 <- paste(stage2.GR3, opt.GR3$OTC[[2]][i], sep=",")
      }
      # concatenate pool sizes for GR4 objective function
      stage2.GR4 <- opt.GR4$OTC[[2]][1]
      for(i in 2:length(opt.GR4$OTC[[2]])){
        stage2.GR4 <- paste(stage2.GR4, opt.GR4$OTC[[2]][i], sep=",")
      }
      # concatenate pool sizes for GR5 objective function
      stage2.GR5 <- opt.GR5$OTC[[2]][1]
      for(i in 2:length(opt.GR5$OTC[[2]])){
        stage2.GR5 <- paste(stage2.GR5, opt.GR5$OTC[[2]][i], sep=",")
      }
      # concatenate pool sizes for GR6 objective function
      stage2.GR6 <- opt.GR6$OTC[[2]][1]
      for(i in 2:length(opt.GR6$OTC[[2]])){
        stage2.GR6 <- paste(stage2.GR6, opt.GR6$OTC[[2]][i], sep=",")
      }
      
      stage2 <- c(stage2.ET, stage2.MAR, stage2.GR1, stage2.GR2, 
                  stage2.GR3, stage2.GR4, stage2.GR5, stage2.GR6)
    } else {
      stage2  <- c(opt.ET$OTC[[2]], opt.MAR$OTC[[2]], opt.GR1$OTC[[2]], 
                   opt.GR2$OTC[[2]], opt.GR3$OTC[[2]], opt.GR4$OTC[[2]], 
                   opt.GR5$OTC[[2]], opt.GR6$OTC[[2]])
    }
  }
  # four-stage hierarchical testing
  # if (length(opt.ET$OTC) > 2) {
  #   if (length(opt.ET$OTC[[3]]) > 1) {
  #     # concatenate pool sizes for ET objective function
  #     stage3.ET <- opt.ET$OTC[[3]][1]
  #     for(i in 2:length(opt.ET$OTC[[3]])){
  #       stage3.ET <- paste(stage3.ET, opt.ET$OTC[[3]][i], sep=",")
  #     }
  #     # concatenate pool sizes for MAR objective function
  #     stage3.MAR <- opt.MAR$OTC[[3]][1]
  #     for(i in 2:length(opt.MAR$OTC[[3]])){
  #       stage3.MAR <- paste(stage3.MAR, opt.MAR$OTC[[3]][i], sep=",")
  #     }
  #     # concatenate pool sizes for GR1 objective function
  #     stage3.GR1 <- opt.GR1$OTC[[3]][1]
  #     for(i in 2:length(opt.GR1$OTC[[3]])){
  #       stage3.GR1 <- paste(stage3.GR1, opt.GR1$OTC[[3]][i], sep=",")
  #     }
  #     # concatenate pool sizes for GR2 objective function
  #     stage3.GR2 <- opt.GR2$OTC[[3]][1]
  #     for(i in 2:length(opt.GR2$OTC[[3]])){
  #       stage3.GR2 <- paste(stage3.GR2, opt.GR2$OTC[[3]][i], sep=",")
  #     }
  #     # concatenate pool sizes for GR3 objective function
  #     stage3.GR3 <- opt.GR3$OTC[[3]][1]
  #     for(i in 2:length(opt.GR3$OTC[[3]])){
  #       stage3.GR3 <- paste(stage3.GR3, opt.GR3$OTC[[3]][i], sep=",")
  #     }
  #     # concatenate pool sizes for GR4 objective function
  #     stage3.GR4 <- opt.GR4$OTC[[3]][1]
  #     for(i in 2:length(opt.GR4$OTC[[3]])){
  #       stage3.GR4 <- paste(stage3.GR4, opt.GR4$OTC[[3]][i], sep=",")
  #     }
  #     # concatenate pool sizes for GR5 objective function
  #     stage3.GR5 <- opt.GR5$OTC[[3]][1]
  #     for(i in 2:length(opt.GR5$OTC[[3]])){
  #       stage3.GR5 <- paste(stage3.GR5, opt.GR5$OTC[[3]][i], sep=",")
  #     }
  #     # concatenate pool sizes for GR6 objective function
  #     stage3.GR6 <- opt.GR6$OTC[[3]][1]
  #     for(i in 2:length(opt.GR6$OTC[[3]])){
  #       stage3.GR6 <- paste(stage3.GR6, opt.GR6$OTC[[3]][i], sep=",")
  #     }
  #     
  #     stage3 <- c(stage3.ET, stage3.MAR, stage3.GR1, stage3.GR2, 
  #                 stage3.GR3, stage3.GR4, stage3.GR5, stage3.GR6)
  #   } else {
  #     stage3  <- c(opt.ET$OTC[[3]], opt.MAR$OTC[[3]], opt.GR1$OTC[[3]], 
  #                  opt.GR2$OTC[[3]], opt.GR3$OTC[[3]], opt.GR4$OTC[[3]], 
  #                  opt.GR5$OTC[[3]], opt.GR6$OTC[[3]])
  #   }
  # }
  
  # columns correspond to objective functions
  # config <- rbind(stage1, stage2)
  # rownames(config) <- stage.labels
  # colnames(config) <- objfn.labels
  
  # rows correspond to objective functions
  if (grepl("Non-informative two-stage", algorithm)) {
    config <- matrix(data = stage1, nrow = length(objfn.labels), ncol = 1)
  } else {
    config <- cbind(stage1, stage2)
  }
  if (grepl("Informative two-stage", algorithm)) {
    stage.labels <- c("Block size", "Group sizes")
  } else if (grepl("hierarchical", algorithm)) {
    # stage.labels <- c("Stage 1", "Stage 2")[1:nrow(config)]
    stage.labels <- c("Stage 1", "Stage 2")[1:ncol(config)]
  } else if (grepl("array", algorithm)) {
    stage.labels <- c("Row/column size", "Array size")
  }
  rownames(config) <- objfn.labels
  colnames(config) <- stage.labels
  
  cat("Optimal testing configuration:\n")
  print(as.data.frame(config))

  # create a matrix detailing the expected number of tests for each obj. fn
  ExpT <- c(opt.ET$ET, opt.MAR$ET, opt.GR1$ET, opt.GR2$ET, 
            opt.GR3$ET, opt.GR4$ET, opt.GR5$ET, opt.GR6$ET)
  value <- c(opt.ET$value, opt.MAR$value, opt.GR1$value, opt.GR2$value, 
             opt.GR3$value, opt.GR4$value, opt.GR5$value, opt.GR6$value)
  
  # columns correspond to objective functions
  # tests <- rbind(ExpT, value)
  # rownames(tests) <- c("Expected number of tests", 
  #                      "Objective function value per individual") 
  # colnames(tests) <- objfn.labels
  
  # rows correspond to objective functions
  tests <- cbind(format(round(ExpT, 2), nsmall = 2), 
                 format(round(value, 4), nsmall = 4))
  rownames(tests) <- objfn.labels
  colnames(tests) <- c("E(T)", "Value")

  cat("\nExpected number of tests:\n")
  print(as.data.frame(tests))
  cat("\nE(T) denotes the expected number of tests.\n")
  cat("Value denotes the objective function value per individual.\n\n")
  
  # create a matrix detailing the accuracy measures for each obj. fn
  if(dim(opt.ET$Accuracy)[1]==1){
    overall.acc <- rbind(opt.ET$Accuracy, opt.MAR$Accuracy, opt.GR1$Accuracy, 
                         opt.GR2$Accuracy, opt.GR3$Accuracy, opt.GR4$Accuracy, 
                         opt.GR5$Accuracy, opt.GR6$Accuracy)
    overall.acc <- format(round(overall.acc, 4), nsmall = 4)
    rownames(overall.acc) <- objfn.labels
    colnames(overall.acc) <- c("PSe", "PSp", "PPPV", "PNPV")
    
    cat("Overall accuracy of the algorithm:\n")
    print(as.data.frame(overall.acc))
  } else if(dim(opt.ET$Accuracy)[1]==2){
    # overall.acc1 <- rbind(opt.ET$Accuracy[1,], opt.MAR$Accuracy[1,], 
    #                       opt.GR1$Accuracy[1,], opt.GR2$Accuracy[1,], 
    #                       opt.GR3$Accuracy[1,], opt.GR4$Accuracy[1,], 
    #                       opt.GR5$Accuracy[1,], opt.GR6$Accuracy[1,])
    # overall.acc1 <- format(round(overall.acc1, 4), nsmall = 4)
    # rownames(overall.acc1) <- objfn.labels
    # colnames(overall.acc1) <- c("PSe", "PSp", "PPPV", "PNPV")
    # 
    # overall.acc2 <- rbind(opt.ET$Accuracy[2,], opt.MAR$Accuracy[2,], 
    #                       opt.GR1$Accuracy[2,], opt.GR2$Accuracy[2,], 
    #                       opt.GR3$Accuracy[2,], opt.GR4$Accuracy[2,], 
    #                       opt.GR5$Accuracy[2,], opt.GR6$Accuracy[2,])
    # overall.acc2 <- format(round(overall.acc2, 4), nsmall = 4)
    # rownames(overall.acc2) <- objfn.labels
    # colnames(overall.acc2) <- c("PSe", "PSp", "PPPV", "PNPV")
    overall.acc <- as.data.frame(format(round(object$opt.ET$Accuracy, 4), nsmall = 4))

    # cat("Overall accuracy of the algorithm for disease 1:\n")
    # print(as.data.frame(overall.acc1))
    # cat("Overall accuracy of the algorithm for disease 2:\n")
    # print(as.data.frame(overall.acc2))
    cat("Overall accuracy of the algorithm:\n")
    print(as.data.frame(overall.acc))
  }
  
  cat("\nPSe denotes the pooling sensitivity.\n")
  cat("PSp denotes the pooling specificity.\n")
  cat("PPPV denotes the pooling positive predictive value.\n")
  cat("PNPV denotes the pooling negative predictive value.\n")
  
  # if (dim(opt.ET$Accuracy)[1] == 1){
  #   res <- list("Algorithm" = algorithm, 
  #               "Configuration" = config, 
  #               "Tests" = tests, 
  #               "Accuracy" = overall.acc)
  # } else {
  #   res <- list("Algorithm" = algorithm, 
  #               "Configuration" = config, 
  #               "Tests" = tests, 
  #               "Accuracy" = list("Disease 1" = overall.acc1,
  #                                 "Disease 2" = overall.acc2))
  # }
  res <- list("Algorithm" = algorithm, 
              "Configuration" = config, 
              "Tests" = tests, 
              "Accuracy" = overall.acc)
  
  class(res) <- "summary.OTC"
  invisible(res)
}




# Supporting functions for OTC1 and the associated calls
###############################################################################

#' @title Determine a vector of probabilities for informative group
#' testing algorithms
#'
#' @description Find the expected value of order statistics from a beta 
#' distribution. This function is used to provide a set of individual 
#' risk probabilities for informative group testing.
#'
#' @param p overall probability of disease that will be used to determine a 
#' vector of individual risk probabilities. This is the expected value of a 
#' random variable with a beta distribution, 
#' \eqn{\frac{\alpha}{\alpha + \beta}}{\alpha/(\alpha + \beta)}.
#' @param alpha a shape parameter for the beta distribution that
#' specifies the degree of heterogeneity for the determined
#' probability vector.
#' @param grp.sz the number of total individuals for which to
#' determine risk probabilities.
#' @param ... arguments to be passed to the \code{beta.dist} function 
#' written by Michael Black for Black et al. (2015).
#'
#' @details This function uses the \code{beta.dist} function from 
#' Black et al. (2015) to determine a vector of individual risk probabilities,
#' ordered from least to greatest. Depending on the specified probability, 
#' \eqn{\alpha} level, and overall group size, simulation may be necessary in 
#' order to determine the probabilities. For this reason, the user should set 
#' a seed in order to reproduce results. The number of simulations can be
#' specified by the user, with 10,000 as the default. The \kbd{expectOrderBeta} 
#' function augments the \code{beta.dist} function by checking whether
#' simulation is needed before attempting to determine the probabilities. 
#' The \kbd{expectOrderBeta} function allows for the number of simulations to 
#' be passed on to the \code{beta.dist} function as an additional argument.
#' See Black et al. (2015) for additional details on the original \kbd{beta.dist} 
#' function.
#'
#' @return A vector of individual risk probabilities.
#'
#' @author Brianna D. Hitt
#'
#' @references
#' \insertRef{Black2015}{binGroup2}
#'
#' @seealso \code{\link{expectOrderBeta}} for generating a vector of individual 
#' risk probabilities and \code{\link{informativeArrayProb}} for arranging a 
#' vector of individual risk probabilities in a matrix for informative array 
#' testing without master pooling.
#'
#' @examples
#' set.seed(8791)
#' expectOrderBeta(p=0.03, alpha=0.5, grp.sz=100)
#'
#' set.seed(52613)
#' expectOrderBeta(p=0.005, alpha=2, grp.sz=40, num.sim=5000)

expectOrderBeta <- function(p, alpha, grp.sz, ...){
  if (is.na(p)) {
    NA
  } else{
    p.try <- suppressWarnings(try(beta.dist2(p = p, alpha = alpha, 
                                             grp.sz = grp.sz, ...), 
                                  silent = TRUE))
    if (class(p.try) == "try-error") {
      beta.dist2(p = p, alpha = alpha, grp.sz = grp.sz, simul = TRUE, ...)
    } else{
      beta.dist2(p = p, alpha = alpha, grp.sz = grp.sz, ...)
    }
  }
}

###############################################################################




# Start MAR.func() function
###############################################################################
#    Brianna Hitt - 4-17-17
#    Purpose: calculates MAR objective function, from Malinovsky, 
#               Albert & Roy (2015)
#      inputs: ET - expected number of tests
#              p.vec - vector of individual probabilities
#              PSe.vec - vector of individual pooling sensitivities
#              PSp.vec - vector of individual pooling specificities
#      Note: The MAR objective function divides ET, the expected number of 
#              tests, by EC, the expected number of correct classifications, 
#              and should be minimized.
#      Note: Malinovsky, Albert, & Roy (2015) maximized the reciprocal, 
#              E(C)/E(T).

MAR.func <- function(ET, p.vec, PSe.vec, PSp.vec){
  EC <- sum(PSe.vec*p.vec + PSp.vec*(1-p.vec))
  ET/EC
}
###############################################################################




# Start GR.func() function
###############################################################################
#    Brianna Hitt - 4-17-17
#    Purpose: calculates GR objective function, from Graff & Roeloffs (1972)
#               M = E(T) + D_1*(# of misclassified negatives) 
#                        + D_2*(# of misclassified positives)
#      inputs: ET - expected number of tests
#              p.vec - vector of individual probabilities
#              PSe.vec - vector of individual pooling sensitivities
#              PSp.vec - vector of individual pooling specificities
#              D1, D2 - weights/costs for misclassification
#      note: this function specifies equal weights of 1 by default

GR.func <- function(ET, p.vec, PSe.vec, PSp.vec, D1=1, D2=1){
  ET + D1*sum((1-PSp.vec)*(1-p.vec)) + D2*sum((1-PSe.vec)*p.vec)
}
###############################################################################




# Start time.it() function
###############################################################################
#    Brianna Hitt - 5-13-17
#    Purpose: calculates the time elapsed
#      inputs: x = object containing the start time

time.it <- function(x) {
  end.time<-proc.time()
  save.time<-end.time-x
  cat("\n Number of minutes running: ", round(save.time[3]/60, 2), "\n \n")
  save.time[3]/60
}
###############################################################################





# Start generate.acc() function
###############################################################################
#    Brianna Hitt - 11-15-19
#    Purpose: generates a vector of sensitivity/specificity values
#             with the appropriate dimensions
#      inputs: algorithm = the group testing algorithm
#              diseases = number of diseases for the diagnostic assay
#              value = the user-specified Se/Sp value for OTC1() and OTC2()
#              label = "sens" for sensitivity, "spec" for specificity

generate.acc <- function(algorithm, diseases, value, label) {
  if (label == "sens") {
    measure <- "sensitivity"
  } else if (label == "spec") {
    measure <- "specificity"
  }
  
  if (diseases == 1) {
    if (length(value) == 1) {
      if (algorithm %in% c("D2", "ID2", "A2", "IA2")) {
        output <- rep(value, 2)
      } else if (algorithm %in% c("D3", "ID3", "A2M")) {
        output <- rep(value, 3)
      } else if (algorithm %in% c("D4", "ID4")) {
        output <- rep(value, 4)
      }
    } else if (algorithm %in% c("D2", "ID2", "A2", "IA2")) {
      if (length(value) == 2) {
        output <- value
      } else{
        stop("Please specify a ", measure, " vector with two values, one for each stage of the testing algorithm.\n")
      }
    } else if (algorithm %in% c("D3", "ID3", "A2M")) {
      if (length(value) == 3) {
        output <- value
      } else{
        stop("Please specify a ", measure, " vector with three values, one for each stage of the testing algorithm.\n")
      }
    } else if (algorithm %in% c("D4", "ID4")) {
      if (length(value) == 4) {
        output <- value
      } else {
        stop("Please specify a ", measure, " vector with four values, one for each stage of the testing algorithm.\n")
      }
    }
  } else if (diseases == 2) {
    if (is.vector(value)) {
      if (length(value) == 2) {
        if (algorithm %in% c("D2", "ID2", "A2", "IA2")) {
          output <- matrix(data = value, nrow = 2, ncol = 2, 
                           dimnames = list(Infection = 1:2, Stage = 1:2))
        } else if (algorithm %in% c("D3", "ID3", "A2M")) {
          output <- matrix(data = value, nrow = 2, ncol = 3, 
                           dimnames = list(Infection = 1:2, Stage = 1:3))
        } else if (algorithm %in% c("D4", "ID4")) {
          output <- matrix(data = value, nrow = 2, ncol = 4, 
                           dimnames = list(Infection = 1:2, Stage = 1:4))
        } else if (algorithm %in% c("D5", "ID5")) {
          output <- matrix(data = value, nrow = 2, ncol = 5, 
                           dimnames = list(Infection = 1:2, Stage = 1:5))
        }
      } else{
        stop("Please specify a matrix of ", measure, " values with the correct dimensions. Each row should correspond to a disease, and each column should correspond to a stage in the algorithm.\n")
      }
    } else if (is.matrix(value)) {
      if (algorithm %in% c("D2", "ID2", "A2", "IA2")) {
        if (dim(value)[1] == 2 & dim(value)[2] == 2) {
          output <- value
        } else{
          stop("Please specify a 2x2 matrix of ", measure, " values with the correct dimensions. Each row should correspond to a disease, and each column should correspond to a stage in the algorithm.\n")
        }
      } else if (algorithm %in% c("D3", "ID3", "A2M")) {
        if (dim(value)[1] == 2 & dim(value)[2] == 3) {
          output <- value
        } else{
          stop("Please specify a 2x3 matrix of ", measure, " values with the correct dimensions. Each row should correspond to a disease, and each column should correspond to a stage in the algorithm.\n")
        }
      } else if (algorithm %in% c("D4", "ID4")) {
        if (dim(value)[1] == 2 & dim(value)[2] == 4) {
          output <- value
        } else {
          stop("Please specify a 2x4 matrix of ", measure, " values with the correct dimensions. Each row should correspond to a disease, and each column should correspond to a stage in the algorithm.\n")
        }
      } else if (algorithm %in% c("D5", "ID5")) {
        if (dim(value)[1] == 2 & dim(value)[2] == 5) {
          output <- value
        } else {
          stop("Please specify a 2x5 matrix of ", measure, " values with the correct dimensions. Each row should correspond to a disease, and each column should correspond to a stage in the algorithm.\n")
        }
      }
    }
  }
}


###############################################################################
