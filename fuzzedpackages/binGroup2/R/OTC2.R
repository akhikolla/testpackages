# Start OTC functions for two diseases
# ##############################################################################
# OTC functions for two diseases - updated 04.23.19
# function to find the OTC for hierarchical testing with two diseases
# Se and Sp are matrices of sensitivity/specificity values
#   corresponding to each disease and each stage of testing



# Start OTC2() function
###############################################################################
#' @title Find the optimal testing configuration for group testing algorithms 
#' that use a multiplex assay for two diseases
#' 
#' @description Find the optimal testing configuration (OTC) using 
#' non-informative and informative hierarchical and array-based group testing 
#' algorithms. Multiplex assays for two diseases are used at each stage of the 
#' algorithms.
#' 
#' @param algorithm character string defining the group testing 
#' algorithm to be used. Non-informative testing options include two-stage 
#' hierarchical ("\kbd{D2}"), three-stage hierarchical ("\kbd{D3}"), 
#' square array testing without master pooling ("\kbd{A2}"), and square array 
#' testing with master pooling ("\kbd{A2M}"). Informative testing options 
#' include two-stage hierarchical ("\kbd{ID2}") and three-stage hierarchical 
#' ("\kbd{ID3}") testing.
#' @param p.vec vector of overall joint probabilities. The joint probabilities
#' are assumed to be equal for all individuals in the algorithm (non-informative 
#' testing only). There are four joint probabilities to 
#' consider: \eqn{p_{00}}{p_00}, the probability that an individual tests 
#' negative for both diseases; \eqn{p_{10}}{p_10}, the probability that an 
#' individual tests positive only for the first disease; \eqn{p_{01}}{p_01}, 
#' the probability that an individual tests positive only for the second 
#' disease; and \eqn{p_{11}}{p_11}, the probability that an individual tests 
#' positive for both diseases. The joint probabilities must sum to 1. 
#' Only one of \kbd{p.vec}, \kbd{probabilities}, or \kbd{alpha} should be specified.
#' @param probabilities matrix of joint probabilities for each individual, 
#' where rows correspond to the four joint probabilities and columns correspond 
#' to each individual in the algorithm. Only one of \kbd{p.vec}, 
#' \kbd{probabilities}, or \kbd{alpha} should be specified.
#' @param alpha vector containing positive shape parameters of the Dirichlet 
#' distribution (for informative testing only). The vector will be used to 
#' generate a heterogeneous matrix of joint probabilities for each individual. 
#' The vector must have length 4. Further details are given under 'Details'. 
#' Only one of \kbd{p.vec}, \kbd{probabilities}, or \kbd{alpha} should be specified. 
#' @param Se matrix of sensitivity values, where one value is given for each 
#' disease (or infection) at each stage of testing. The rows of the matrix correspond 
#' to each disease \eqn{k=1,...,K}, and the columns of the matrix correspond to each 
#' stage of testing \eqn{s=1,...,S}. If a vector of \eqn{K} values is provided, 
#' the sensitivity values associated with disease \eqn{k} are assumed to be equal 
#' to the \eqn{k}th value in the vector for all stages of testing. 
#' Further details are given under 'Details'.
#' @param Sp matrix of specificity values, where one value is given for each 
#' disease (or infection) at each stage of testing. The rows of the matrix correspond 
#' to each disease \eqn{k=1,...,K}, and the columns of the matrix correspond 
#' to each stage of testing \eqn{s=1,...,S}. If a vector of \eqn{K} values is provided, 
#' the specificity values associated with disease \eqn{k} are assumed to be equal 
#' to the \eqn{k}th value in the vector for all stages of testing. 
#' Further details are given under 'Details'.
#' @param ordering matrix detailing the ordering for the binary responses of 
#' the diseases. The columns of the matrix correspond to each disease and the 
#' rows of the matrix correspond to each of the 4 sets of binary responses for 
#' two diseases. This ordering is used with the joint probabilities. The 
#' default ordering is (p_00, p_10, p_01, p_11). 
#' @param group.sz single group size or range of group sizes for which to 
#' calculate operating characteristics and/or find the OTC. The details of 
#' group size specification are given under 'Details'. 
#' @param trace a logical value indicating whether the progress of 
#' calculations should be printed for each initial group size provided by 
#' the user. The default is \kbd{TRUE}. 
#' @param print.time a logical value indicating whether the length of time 
#' for calculations should be printed. The default is \kbd{TRUE}.
#' @param ... additional arguments to be passed to functions for hierarchical 
#' testing with multiplex assays for two diseases.
#' 
#' @details This function finds the OTC for standard group testing algorithms 
#' with a multiplex assay that tests for two diseases and computes the 
#' associated operating characteristics. Calculations for hierarchical group 
#' testing algorithms are performed as described in Bilder et al. (2019) and 
#' calculations for array-based group testing algorithms are performed as 
#' described in Hou et al. (2019).
#' 
#' Available algorithms include two- and three-stage hierarchical testing and 
#' array testing with and without master pooling. Both non-informative and informative
#' group testing settings are allowed for hierarchical algorithms. Only 
#' non-informative group testing settings are allowed for array testing algorithms. 
#' Operating characteristics calculated are expected number of tests, 
#' pooling sensitivity, pooling specificity, pooling positive predictive value, and 
#' pooling negative predictive value for each individual.
#' 
#' For informative algorithms where the \kbd{alpha} argument is specified, a 
#' heterogeneous matrix of joint probabilities for each individual is generated 
#' using the Dirichlet distribution. This is done using 
#' \code{rBeta2009::rdirichlet} and requires the user to set a seed to 
#' reproduce results. See Bilder et al. (2019) for additional details on the 
#' use of the Dirichlet distribution for this purpose.
#' 
#' The sensitivity/specificity values are allowed to vary across stages of 
#' testing. For hierarchical testing, a different sensitivity/specificity value
#' may be used for each stage of testing. For array testing, a different 
#' sensitivity/specificity value may be used for master pool testing (if included), 
#' row/column testing, and individual testing. The values must be specified 
#' in the order of the testing performed. For example, values are specified as 
#' (stage 1, stage 2, stage 3) for three-stage hierarchical testing or 
#' (master pool testing, row/column testing, individual testing) for array 
#' testing with master pooling. A vector of \eqn{K} sensitivity/specificity 
#' values may be specified, and sensitivity/specificity values for all stages 
#' of testing are assumed to be equal. The first value in the vector will be 
#' used at each stage of testing for the first disease, and the second value in 
#' the vector will be used at each stage of testing for the second disease.
#' 
#' The value(s) specified by \kbd{group.sz} represent the initial (stage 1) 
#' group size for hierarchical testing and the row/column size for array 
#' testing. If a single value is provided for \kbd{group.sz} with two-stage 
#' hierarchical or array testing, operating characteristics will be calculated 
#' and no optimization will be performed. If a single value is provided for 
#' \kbd{group.sz} with three-stage hierarchical, the OTC will be found over all 
#' possible configurations with this initial group size. If a range of group 
#' sizes is specified, the OTC will be found over all group sizes.
#' 
#' In addition to the OTC, operating characteristics for some of the other
#' configurations corresponding to each initial group size provided by the user
#' are displayed. For algorithms where there is only one configuration for each
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
#' The displayed overall pooling sensitivity, pooling specificity, pooling 
#' positive predictive value, and pooling negative predictive value are 
#' weighted averages of the corresponding individual accuracy measures for all 
#' individuals within the initial group (or block) for a hierarchical algorithm, or 
#' within the entire array for an array-based algorithm. 
#' Expressions for these averages are provided in the Supplementary Material 
#' for Hitt et al. (2019). These expressions are based on accuracy definitions 
#' given by Altman and Bland (1994a, 1994b). Individual accuracy measures can 
#' be calculated using the \code{\link{operatingCharacteristics2}} 
#' (\code{\link{opChar2}}) function.
#' 
#' @return A list containing:
#' \item{algorithm}{the group testing algorithm used for calculations.}
#' \item{prob.vec}{the vector of joint probabilities provided by the user, 
#' if applicable (for non-informative algorithms only).}
#' \item{joint.p}{the matrix of joint probabilities for each individual 
#' provided by the user, if applicable.}
#' \item{alpha.vec}{the alpha vector provided by the user, if applicable 
#' (for informative algorithms only).}
#' \item{Se}{the matrix of sensitivity values for each disease at each stage of 
#' testing.}
#' \item{Sp}{the matrix of specificity values for each disease at each stage of 
#' testing.}
#' \item{opt.ET}{a list containing:
#' \describe{
#' \item{OTC}{a list specifying elements of the optimal testing configuration, 
#' which may include:
#' \describe{
#' \item{Stage1}{group size for the first stage of hierarchical testing, if 
#' applicable.}
#' \item{Stage2}{group sizes for the second stage of hierarchical testing, if 
#' applicable.}
#' \item{Block.sz}{the block size/initial group size for informative Dorfman 
#' testing, which is not tested.}
#' \item{pool.szs}{group sizes for the first stage of testing for informative 
#' Dorfman testing.}
#' \item{Array.dim}{the row/column size for array testing.}
#' \item{Array.sz}{the overall array size for array testing (the square of the 
#' row/column size).}}}
#' \item{p.mat}{the matrix of joint probabilities for each individual in the 
#' algorithm. Each row corresponds to one of the four joint probabilities. 
#' Each column corresponds to an individual in the testing algorithm.}
#' \item{ET}{the expected testing expenditure for the OTC.}
#' \item{value}{the value of the expected number of tests per individual.}
#' \item{Accuracy}{the matrix of overall accuracy measures for the algorithm. 
#' The rows correspond to each disease. The columns 
#' correspond to the pooling sensitivity, pooling specificity, pooling positive 
#' predictive value, and pooling negative predictive value for the overall 
#' algorithm. Further details are given under 'Details'.}}}
#' \item{Configs}{a data frame containing results for the best configuration 
#' for each initial group size provided by the user. The columns correspond to 
#' the initial group size, configuration (if applicable), overall array size 
#' (if applicable), expected number of tests, value of the objective function 
#' per individual, and accuracy measures for each disease. Accuracy measures 
#' include the pooling sensitivity, pooling specificity, pooling positive 
#' predictive value, and pooling negative predictive value. No results are 
#' displayed if a single \kbd{group.sz} is provided. Further details are given 
#' under 'Details'.}
#' \item{Top.Configs}{a data frame containing results for some of the top 
#' configurations for each initial group size provided by the user. The 
#' columns correspond to the initial group size, configuration, 
#' expected number of tests, value of the objective function per individual, 
#' and accuracy measures for each disease. Accuracy measures include the  
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
#' @author This function was written by Brianna D. Hitt. It calls 
#' \kbd{ET.all.stages.new} and \kbd{PSePSpAllStages}, which were originally 
#' written by Christopher Bilder for Bilder et al. (2019), and \kbd{ARRAY}, 
#' which was originally written by Peijie Hou for Hou et al. (2020). The 
#' functions \kbd{ET.all.stages.new}, \kbd{PSePSpAllStages}, and \kbd{ARRAY} 
#' were obtained from \url{http://chrisbilder.com/grouptesting}. Minor 
#' modifications were made to the functions for inclusion in the binGroup2 
#' package.
#' 
#' @references
#' \insertRef{Altman1994a}{binGroup2}
#' 
#' \insertRef{Altman1994b}{binGroup2}
#' 
#' \insertRef{Bilder2019}{binGroup2}
#' 
#' \insertRef{Hitt2019}{binGroup2}
#' 
#' \insertRef{Hou2019}{binGroup2}
#' 
#' \insertRef{McMahan2012a}{binGroup2}
#' 
#' @family OTC functions
#' @family multiplex testing functions
#' 
#' @examples 
#' # Estimated running time for all examples was calculated 
#' #   using a computer with 16 GB of RAM and one core of 
#' #   an Intel i7-6500U processor. Please take this into 
#' #   account when interpreting the run times given.
#' 
#' # Find the OTC for non-informative two-stage 
#' #   hierarchical (Dorfman) testing
#' Se <- matrix(data = c(0.95, 0.95, 0.99, 0.99), nrow = 2, ncol = 2,
#'              dimnames = list(Infection = 1:2, Stage = 1:2))
#' Sp <- matrix(data = c(0.96, 0.96, 0.98, 0.98), nrow = 2, ncol = 2,
#'              dimnames = list(Infection = 1:2, Stage = 1:2))
#' OTC2(algorithm = "D2", p.vec=c(0.90, 0.04, 0.04, 0.02),
#'      Se = Se, Sp = Sp, group.sz = 3:30)
#'
#' # Find the OTC over all possible testing configurations 
#' #   for informative two-stage hierarchical (Dorfman) 
#' #   testing with a specified group size.
#' # A matrix of joint probabilities for each individual is 
#' #   generated using the Dirichlet distribution.
#' # This examples takes approximately 25 seconds to run.
#' Se <- matrix(data = rep(0.95, 4), nrow = 2, ncol = 2,
#'              dimnames = list(Infection = 1:2, Stage = 1:2))
#' Sp <- matrix(data = rep(0.99, 4), nrow = 2, ncol = 2,
#'              dimnames = list(Infection = 1:2, Stage = 1:2))
#' \donttest{
#' set.seed(1002)
#' OTC2(algorithm = "ID2", alpha=c(18.25, 0.75, 0.75, 0.25),
#'      Se = Se, Sp = Sp, group.sz = 10:20)}
#'      
#' # Find the OTC for non-informative three-stage 
#' #   hierarchical testing.
#' # This example takes approximately 1 minute to run.
#' Se <- matrix(data = rep(0.95, 6), nrow = 2, ncol = 3,
#'              dimnames = list(Infection = 1:2, Stage = 1:3))
#' Sp <- matrix(data = rep(0.99, 6), nrow = 2, ncol = 3,
#'              dimnames = list(Infection = 1:2, Stage = 1:3))
#' \donttest{
#' OTC2(algorithm = "D3", p.vec=c(0.95, 0.02, 0.02, 0.01),
#'      Se = Se, Sp = Sp, group.sz = 3:20)}
#'      
#' # Find the OTC over all possible configurations 
#' #   for informative three-stage hierarchical 
#' #   testing with a specified group size 
#' #   and a heterogeneous matrix of joint 
#' #   probabilities for each individual. 
#' set.seed(8791)
#' Se <- matrix(data = rep(0.95, 6), nrow = 2, ncol = 3,
#'              dimnames = list(Infection = 1:2, Stage = 1:3))
#' Sp <- matrix(data = rep(0.99, 6), nrow = 2, ncol = 3,
#'              dimnames = list(Infection = 1:2, Stage = 1:3))
#' p.unordered <- t(rBeta2009::rdirichlet(n = 12, 
#'                             shape = c(18.25, 0.75, 0.75, 0.25)))
#' p.ordered <- p.unordered[, order(1 - p.unordered[1,])]
#' OTC2(algorithm="ID3", probabilities = p.ordered,
#'          Se=Se, Sp=Sp, group.sz = 12, 
#'          trace=FALSE, print.time=FALSE)
#'                             
#' # Find the OTC for non-informative array testing 
#' #   without master pooling.
#' Se <- matrix(data = rep(0.95, 4), nrow = 2, ncol = 2,
#'              dimnames = list(Infection = 1:2, Stage = 1:2))
#' Sp <- matrix(data = rep(0.99, 4), nrow = 2, ncol = 2,
#'              dimnames = list(Infection = 1:2, Stage = 1:2))
#' OTC2(algorithm = "A2", p.vec=c(0.90, 0.04, 0.04, 0.02),
#'      Se = Se, Sp = Sp, group.sz = 3:12)
#'                   
#' # Find the OTC for non-informative array testing 
#' #   with master pooling.
#' Se <- matrix(data = rep(0.95, 6), nrow = 2, ncol = 3,
#'              dimnames = list(Infection = 1:2, Stage = 1:3))
#' Sp <- matrix(data = rep(0.99, 6), nrow = 2, ncol = 3,
#'              dimnames = list(Infection = 1:2, Stage = 1:3))
#' OTC2(algorithm = "A2M", p.vec=c(0.90, 0.04, 0.04, 0.02),
#'      Se = Se, Sp = Sp, group.sz = 10, 
#'      trace=FALSE, print.time=FALSE)

# Brianna Hitt - 04.02.2020
# Changed cat() to warning()

OTC2 <- function(algorithm, p.vec = NULL, probabilities = NULL, 
                 alpha = NULL, Se, Sp, 
                 ordering = matrix(data = c(0,1,0,1,0,0,1,1), 
                                   nrow = 4, ncol = 2), group.sz, 
                 trace = TRUE, print.time = TRUE, ...){
  
  if (!(algorithm %in% c("D2", "D3", "A2", "A2M", "ID2", "ID3"))) {
    stop("Please specify one of the following algorithms: D2, ID2, D3, ID3, A2, A2M.")
  }
  
  # check that the joint probabilities are specified, either through p.vec, 
  #   probabilities, or alpha arguments
  if (is.null(p.vec) & is.null(probabilities) & is.null(alpha)) {
    if (algorithm %in% c("D2", "D3", "A2", "A2M")) {
      stop("Please specify an overall joint probability vector using the 'p.vec' argument, or specify a matrix of joint probabilities for each individual using the 'probabilities' argument.\n")
    } else if (algorithm %in% c("ID2", "ID3", "IA2")) {
      stop("Please specify a matrix of joint probabilities for each individual using the 'probabilities' argument, or specify a vector of shape parameters for the Dirichlet distribution using the 'alpha' argument.\n")
    }
  } else if (sum(!is.null(p.vec), !is.null(probabilities), !is.null(alpha)) > 1) {
    stop("You have specified more than one of the following arguments: p.vec, probabilities, alpha. Please specify only one option.\n")
  } else {
    if (!is.null(p.vec)) {
      if (length(p.vec) != 4) {
        stop("Please specify an overall joint probability vector of length 4 (one for each of the joint probabilities), and the matrix of individual probabilities will be generated based on the algorithm specified for each group size included in the range.\n")
      }
      if (algorithm %in% c("ID2", "ID3", "IA2")) {
        stop("You have specified an overall joint probability vector for an informative algorithm. Please specify a matrix of joint probabilities for each individual using the 'probabilities' argument, or specify shape parameters for the Dirichlet distribution using the 'alpha' argument.\n")
      }
      if (sum(p.vec) != 1.00) {
        stop("Please specify joint probabilities that sum to 1.")
      }
    } else if (!is.null(probabilities)) {
      if (length(group.sz) == 1) {
        if (dim(probabilities)[1] != 4 | 
            (algorithm %in% c("D2", "D3", "ID2", "ID3") & dim(probabilities)[2] != group.sz) | 
            (algorithm %in% c("A2", "IA2", "A2M") & dim(probabilities)[2] != group.sz^2)) {
          stop("Please specify a matrix of joint probabilities with the correct dimensions. Each row should correspond to one of the four joint probabilities. Each column should correspond to an individual in the algorithm.\n")
        }
        
        if ((algorithm %in% c("D2", "D3", "A2", "A2M")) & 
            (all.equal(probabilities[1,], rep(probabilities[1,1], dim(probabilities)[2])) != TRUE | 
             all.equal(probabilities[2,], rep(probabilities[2,1], dim(probabilities)[2])) != TRUE | 
             all.equal(probabilities[3,], rep(probabilities[3,1], dim(probabilities)[2])) != TRUE | 
             all.equal(probabilities[4,], rep(probabilities[4,1], dim(probabilities)[2])) != TRUE)) {
          stop("You have specified a heterogeneous matrix of joint probabilities for a non-informative algorithm. Please specify a homogeneous matrix of joint probabilities using the 'probabilities' argument or specify an overall joint probability vector using the 'p.vec' argument.\n")
        }
      } else {
        stop("You have specified a matrix of joint probabilities along with a range of group sizes. Please specify a single group size.\n")
      }
    } else if (!is.null(alpha)) {
      if (length(alpha) != 4) {
        stop("Please specify an alpha vector of length 4 (one shape parameter for each of the joint probabilities), and the matrix of individual probabilities will be generated for each group size included in the range.\n")
      }
      if (algorithm %in% c("D2", "D3", "A2", "A2M")) {
        stop("You have specified a vector of shape parameters for the Dirichlet distribution using a non-informative algorithm. Please specify an overall vector of joint probabilities using the 'p.vec' argument, or specify a matrix of joint probabilities for each individual using the 'probabilities' argument.\n")
      }
    }
  }
  
  Se <- generate.acc(algorithm = algorithm, diseases = 2, value = Se, label = "sens")
  Sp <- generate.acc(algorithm = algorithm, diseases = 2, value = Sp, label = "spec")
  
  # check the dimension of the ordering matrix
  if (dim(ordering)[1] != 4 | dim(ordering)[2] != 2) {
    stop("Please specify an ordering matrix with the correct dimensions. The columns should correspond to each disease. The rows should contain each of the four combinations of binary responses for the two diseases.\n")
  }
  
  # check the minimum and maximum group sizes
  if (min(group.sz) < 3) {
    if (algorithm %in% c("D2", "D3", "ID2", "ID3")) {
      stop("Please specify a minimum group size of at least 3.\n")
    } else {
      stop("Please specify a minimum row/column size of at least 3.\n")
    }
  }
  if (max(group.sz) >= 50) {
    if (algorithm %in% c("D3", "ID2", "ID3")) {
      warning("You have specified a maximum group size of 50 or larger. This function may take a VERY long time to run. Press 'ESC' if you wish to cancel the submitted statements.\n")
    } else if (algorithm %in% c("A2", "A2M", "IA2")) {
      warning("You have specified a maximum row/column size of 50 or larger. This function may take a VERY long time to run. Press 'ESC' if you wish to cancel the submitted statements.\n")
    }
  }
  
  # call function for non-informative two-stage hierarchical (Dorfman) testing
  if (algorithm == "D2") {
    if  (!is.null(p.vec)) {
      results <- NI.Dorf.OTC2(p.vec = p.vec, Se = Se, Sp = Sp, 
                              ordering = ordering, group.sz = group.sz, 
                              trace=trace, print.time=print.time, ...)
    } else if (!is.null(probabilities)) {
      results <- NI.Dorf.OTC2(p.vec = probabilities[,1], Se = Se, Sp = Sp, 
                              ordering = ordering, group.sz = group.sz, 
                              trace=trace, print.time=print.time, ...)
    }
  }
  
  # call function for non-informative three-stage hierarchical testing
  if (algorithm == "D3") {
    if (!is.null(p.vec)) {
      results <- NI.D3.OTC2(p.vec = p.vec, Se = Se, Sp = Sp, 
                            ordering = ordering, group.sz = group.sz, 
                            trace=trace, print.time=print.time, ...)
    } else if (!is.null(probabilities)) {
      results <- NI.D3.OTC2(p.vec = probabilities[,1], Se = Se, Sp = Sp, 
                            ordering = ordering, group.sz = group.sz, 
                            trace=trace, print.time=print.time, ...)
    }
  }
  
  # call function for non-informative array testing without master pooling
  if (algorithm == "A2") {
    if (!is.null(p.vec)) {
      results <- NI.Array.OTC2(p.vec = p.vec, Se = Se, Sp = Sp, 
                               group.sz = group.sz, 
                               trace=trace, print.time=print.time, ...)
    } else if (!is.null(probabilities)) {
      results <- NI.Array.OTC2(p.vec = probabilities[,1], Se = Se, Sp = Sp, 
                               group.sz = group.sz, 
                               trace=trace, print.time=print.time, ...)
    }
  }
  
  # call function for non-informative array testing with master pooling
  if (algorithm == "A2M") {
    if (!is.null(p.vec)) {
      results <- NI.A2M.OTC2(p.vec = p.vec, Se = Se, Sp = Sp, 
                             group.sz = group.sz, 
                             trace=trace, print.time=print.time, ...)
    } else if (!is.null(probabilities)) {
      results <- NI.A2M.OTC2(p.vec = probabilities[,1], Se = Se, Sp = Sp, 
                             group.sz = group.sz, 
                             trace=trace, print.time=print.time, ...)
    }
  }
  
  # call function for informative two-stage hierarchical (PSOD) testing
  if (algorithm == "ID2") {
    if  (!is.null(alpha)) {
      results <- Inf.Dorf.OTC2(alpha = alpha, Se = Se, Sp = Sp, 
                               ordering = ordering, group.sz = group.sz, 
                               trace=trace, print.time=print.time, ...)
    } else if (!is.null(probabilities)) {
      results <- Inf.Dorf.OTC2(probabilities = probabilities, Se = Se, Sp = Sp, 
                               ordering = ordering, group.sz = group.sz, 
                               trace=trace, print.time=print.time, ...)
    }
  }
  
  # call function for informative three-stage hierarchical testing
  if (algorithm == "ID3") {
    if  (!is.null(alpha)) {
      results <- Inf.D3.OTC2(alpha = alpha, Se = Se, Sp = Sp, 
                             ordering = ordering, group.sz = group.sz, 
                             trace=trace, print.time=print.time, ...)
    } else if (!is.null(probabilities)) {
      results <- Inf.D3.OTC2(probabilities = probabilities, Se = Se, Sp = Sp, 
                             ordering = ordering, group.sz = group.sz, 
                             trace=trace, print.time=print.time, ...)
    }
  }
  
  # call function for informative array testing not available
  
  class(results) <- "OTC"
  results
}

