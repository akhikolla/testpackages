# Start operatingCharacteristics1() function
###############################################################################
# Brianna Hitt - 12-06-19
# This function is the same as OTC1(), but no longer finds the 
#   optimal testing configuration. It only calculates operating
#   characteristics for a specified testing configuration. It
#   allows the user to specify a configuration using the group.sz 
#   (initial group size) and pool.szs arguments.
# Brianna Hitt - 02-17-20
# Changed the inputs - removed group.sz and pool.szs
# Replaced with hier.config (a group membership matrix) and rowcol.sz
# Added four-stage hierarchical testing for both non-informative and 
#   informative settings

#' @title Calculate operating characteristics for group testing algorithms 
#' that use a single-disease assay
#'
#' @description Calculate operating characteristics, such as the expected 
#' number of tests, for a specified testing configuration using 
#' non-informative and informative hierarchical and array-based group testing 
#' algorithms. Single-disease assays are used at each stage of the algorithms.
#'
#' @param algorithm character string defining the group testing algorithm to be 
#' used. Non-informative testing options include two-stage hierarchical ("\kbd{D2}"),
#' three-stage hierarchical ("\kbd{D3}"), four-stage hierarchical ("\kbd{D4}"), 
#' square array testing without master pooling ("\kbd{A2}"), and square array 
#' testing with master pooling ("\kbd{A2M}"). Informative testing options 
#' include two-stage hierarchical ("\kbd{ID2}"), three-stage hierarchical 
#' ("\kbd{ID3}"), four-stage hierarchical ("\kbd{ID4}"), and square array testing 
#' without master pooling ("\kbd{IA2}").
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
#' @param hier.config a matrix specifying the configuration for a hierarchical 
#' testing algorithm. The rows correspond to the stages of testing, the columns 
#' correspond to each individual to be tested, and the cell values 
#' specify the group number of each individual at each stage. Further details 
#' are given under 'Details'. For array testing algorithms, this argument will be 
#' ignored.
#' @param rowcol.sz the row/column size for array testing algorithms. For hierarchical 
#' testing algorithms, this argument will be ignored.
#' @param alpha a shape parameter for the beta distribution that specifies the degree of
#' heterogeneity for the generated probability vector (for informative testing only).
#' @param a a vector containing indices indicating which individuals to 
#' calculate individual accuracy measures for. If \kbd{NULL}, individual accuracy 
#' measures will be displayed for all individuals in the algorithm. 
#' @param print.time a logical value indicating whether the length of time 
#' for calculations should be printed. The default is \kbd{TRUE}.
#' @param ... arguments to be passed to the \code{\link{expectOrderBeta}} function, which
#' generates a vector of probabilities for informative testing algorithms.
#' Further details are given under 'Details'.
#'
#' @details This function computes the operating characteristics for 
#' group testing algorithms with an assay that tests for one disease, as 
#' described in Hitt et al. (2019).
#'
#' Available algorithms include two-, three-, and four-stage hierarchical testing and
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
#' The matrix specified  by \kbd{hier.config} defines the hierarchical group testing 
#' algorithm for \eqn{I} individuals. The rows of the matrix correspond to the stages 
#' \eqn{s=1,...,S} in the testing algorithm, and the columns correspond to individuals 
#' \eqn{i=1,...I}. The cell values within the matrix represent the group number of 
#' individual \eqn{i} at stage \eqn{s}. For three-stage, four-stage, and non-informative 
#' two-stage hierarchical testing, the first row of the matrix consists of all ones. 
#' This indicates that all individuals in the algorithm are tested together in a single 
#' group in the first stage of testing. For informative two-stage hierarchical testing, 
#' the initial group (block) is not tested. Thus, the first row of the matrix 
#' consists of the group numbers for each individual in the first stage of testing. For 
#' all hierarchical algorithms, the final row of the matrix denotes individual 
#' testing. Individuals who are not tested in a particular stage are represented 
#' by "NA" (e.g., an individual tested in a group of size 1 in the second stage of testing 
#' would not be tested again in a third stage of testing). It is important to note that 
#' this matrix represents the testing that could be performed if each group tests positively 
#' at each stage prior to the last. For more details on this matrix (called a group membership 
#' matrix), see Bilder et al. (2019). 
#' 
#' For array testing without master pooling, the \kbd{rowcol.sz} specified
#' represents the row/column size for initial (stage 1) testing. For array testing
#' with master pooling, the \kbd{rowcol.sz} specified represents the row/column size
#' for stage 2 testing. This is because the master pool size is the overall array 
#' size, given by the square of the row/column size.
#'
#' The displayed overall pooling sensitivity, pooling specificity, pooling positive
#' predictive value, and pooling negative predictive value are weighted
#' averages of the corresponding individual accuracy measures for all
#' individuals within the initial group (or block) for a hierarchical algorithm, or
#' within the entire array for an array-based algorithm.
#' Expressions for these averages are provided in the Supplementary
#' Material for Hitt et al. (2019). These expressions are based on accuracy
#' definitions given by Altman and Bland (1994a, 1994b).
#'
#' The \kbd{operatingCharacteristics1} function accepts additional arguments, 
#' namely \kbd{num.sim}, to be passed to the \code{\link{expectOrderBeta}} function, 
#' which generates a vector of probabilities for informative group testing algorithms. 
#' The \kbd{num.sim} argument specifies the number of simulations from the beta 
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
#' \item{Config}{a list specifying elements of the specified testing configuration, 
#' which may include:
#' \describe{
#' \item{Stage1}{group size for the first stage of hierarchical testing, if applicable.}
#' \item{Stage2}{group sizes for the second stage of hierarchical testing, if applicable.}
#' \item{Stage3}{group sizes for the third stage of hierarchical testing, if applicable.}
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
#' \item{value}{the value of the expected number of tests per individual.}
#' \item{Accuracy}{a list containing:
#' \describe{
#' \item{Individual}{a matrix of accuracy measures for each individual 
#' specified in \kbd{a}. The rows correspond to each unique set of accuracy 
#' measures in the algorithm. Individuals with the same set of accuracy 
#' measures are displayed together in a single row of the matrix. The columns 
#' correspond to the pooling sensitivity, pooling specificity, pooling positive 
#' predictive value, pooling negative predictive value, and the indices for the 
#' individuals in each row of the matrix.}
#' \item{Overall}{a matrix of overall accuracy measures for the algorithm. 
#' The columns correspond to the pooling sensitivity, pooling specificity, 
#' pooling positive predictive value, and pooling negative predictive value 
#' for the overall algorithm. Further details are given under 'Details'.}}}
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
#' \insertRef{Hitt2019}{binGroup2}
#'
#' \insertRef{McMahan2012a}{binGroup2}
#'
#' \insertRef{McMahan2012b}{binGroup2}
#'
#' @family operating characteristic functions
#'
#' @examples
#' # Calculate the operating characteristics for non-informative
#' #   two-stage hierarchical (Dorfman) testing.
#' config.mat <- matrix(data = c(rep(1, 10), 1:10), 
#'                      nrow = 2, ncol = 10, byrow = TRUE)
#' opChar1(algorithm="D2", p=0.05, Se=0.99, Sp=0.99, 
#'         hier.config=config.mat)
#' opChar1(algorithm="D2", p=0.05, Se=0.99, Sp=0.99, 
#'         hier.config=config.mat, a=c(1,4), print.time=FALSE)
#'
#' # Calculate the operating characteristics for informative
#' #   two-stage hierarchical (Dorfman) testing.
#' # A vector of individual probabilities is generated using
#' #   the expected value of order statistics from a beta
#' #   distribution with p = 0.01 and a heterogeneity level
#' #   of alpha = 0.5.
#' config.mat <- matrix(data = c(rep(1:3, each = 10), 1:30), 
#'                      nrow = 2, ncol = 30, byrow = TRUE)
#' set.seed(52613)
#' opChar1(algorithm="ID2", p=0.01, Se=0.95, Sp=0.95, 
#'         hier.config=config.mat, alpha=0.5, num.sim=10000)
#' # Equivalent code using a heterogeneous vector of 
#' #   probabilities
#' set.seed(52613)
#' probs <- expectOrderBeta(p=0.01, alpha=0.5, grp.sz=30)
#' opChar1(algorithm="ID2", probabilities=probs, Se=0.95, Sp=0.95, 
#'         hier.config=config.mat)
#' 
#' # Calculate the operating characteristics for
#' #   non-informative three-stage hierarchical testing.
#' config.mat <- matrix(data = c(rep(1, 18), rep(1:3, each = 5), 
#'                               rep(4, 3), 1:18), 
#'                     nrow = 3, ncol = 18, byrow = TRUE)
#' opChar1(algorithm="D3", p=0.001, Se=0.95, Sp=0.95, 
#'         hier.config=config.mat)
#' opChar1(algorithm="D3", p=0.001, Se=c(0.95, 0.95, 0.99), 
#'         Sp=c(0.96, 0.96, 0.98), hier.config=config.mat)
#' 
#' # Calculate the operating characteristics for 
#' #   informative three-stage hierarchical testing, 
#' #   given a heterogeneous vector of probabilities.
#' config.mat <- matrix(data = c(rep(1, 6), rep(1:2, each = 3), 
#'                               1:6), nrow = 3, ncol = 6, 
#'                      byrow = TRUE)
#' set.seed(52613)
#' opChar1(algorithm="ID3", 
#'          probabilities=c(0.012, 0.014, 0.011, 0.012, 0.010, 0.015), 
#'          Se=0.99, Sp=0.99, hier.config=config.mat, 
#'          alpha=0.5, num.sim=5000)
#'
#' # Calculate the operating characteristics for 
#' #   non-informative four-stage hierarchical testing.
#' config.mat <- matrix(data = c(rep(1, 12), rep(1, 8), 
#'                               rep(2, 2), 3, 4, rep(1, 5), 
#'                               rep(2, 3), 3, 4, rep(NA, 2), 
#'                               1:8, rep(NA, 4)), nrow = 4, 
#'                      ncol = 12, byrow = TRUE)
#' opChar1(algorithm="D4", p=0.041, Se=0.99, Sp=0.90, 
#'         hier.config=config.mat)
#'         
#' # Calculate the operating characteristics for 
#' #   informative four-stage hierarchical testing. 
#' # A vector of individual probabilities is generated using
#' #   the expected value of order statistics from a beta
#' #   distribution with p = 0.041 and a heterogeneity level
#' #   of alpha = 0.5.
#' config.mat <- matrix(data = c(rep(1, 12), rep(1, 8), 
#'                               rep(2, 2), 3, 4, rep(1, 5), 
#'                               rep(2, 3), 3, 4, rep(NA, 2), 
#'                               1:8, rep(NA, 4)), nrow = 4, 
#'                      ncol = 12, byrow = TRUE)
#' set.seed(5678)
#' opChar1(algorithm="ID4", p=0.041, Se=0.99, Sp=0.90,
#'         hier.config=config.mat, alpha=0.5)
#' 
#' # Calculate the operating characteristics for
#' #   non-informative array testing without master pooling.
#' opChar1(algorithm="A2", p=0.005, Se=c(0.95, 0.99), 
#'         Sp=c(0.95, 0.99), rowcol.sz=8, a=1)
#'
#' # Calculate the operating characteristics for 
#' #   informative array testing without master pooling.
#' # A vector of individual probabilities is generated using
#' #   the expected value of order statistics from a beta
#' #   distribution with p = 0.03 and a heterogeneity level
#' #   of alpha = 2.
#' set.seed(1002)
#' opChar1(algorithm="IA2", p=0.03, Se=0.95, Sp=0.95,
#'          rowcol.sz=8, alpha=2, a=1:10)
#'
#' # Calculate the operating characteristics for 
#' #   non-informative array testing with master pooling.
#' opChar1(algorithm="A2M", p=0.02, Se=c(0.95,0.95,0.99), 
#'         Sp=c(0.98,0.98,0.99), rowcol.sz=5)

# Brianna Hitt - 04.02.2020
# Changed cat() to warning()

operatingCharacteristics1 <- function(algorithm, p=NULL, probabilities=NULL, 
                                      Se=0.99, Sp=0.99, hier.config = NULL, 
                                      rowcol.sz = NULL, alpha=2,
                                      a=NULL, print.time = TRUE, ...){
  
  ## make sure that all necessary information is included in the correct format
  if (!(algorithm %in% c("D2", "D3", "D4", "A2", "A2M", "ID2", "ID3", "ID4", "IA2"))) {
    stop("Please specify one of the following algorithms: D2, ID2, D3, ID3, D4, ID4, A2, IA2, A2M.")
  }
  
  if (algorithm %in% c("D2", "D3", "D4", "ID2", "ID3", "ID4")) {
    if (is.null(hier.config) | 
        (algorithm %in% c("D2", "ID2") & nrow(hier.config) != 2) | 
        (algorithm %in% c("D3", "ID3") & nrow(hier.config) != 3) | 
        (algorithm %in% c("D4", "ID4") & nrow(hier.config) != 4)) {
      stop("Please provide a matrix specifying the configuration for the hierarchical algorithm. The rows correspond to the stages of testing, the columns correspond to each individual in the algorithm, and the cell values specify the group number of each individual at each stage.\n")
    } else {
      group.sz <- ncol(hier.config)
      dimnames(hier.config) <- NULL
    }
  } else if (algorithm %in% c("A2", "IA2", "A2M")) {
    if (is.null(rowcol.sz)) {
      stop("Please provide a row/column size for the array testing algorithm.\n")
    } else {
      group.sz <- rowcol.sz
    }
  }
  # if ((algorithm %in% c("D2", "D3", "D4", "ID2", "ID3", "ID4") & !is.null(rowcol.sz)) | 
  #     (algorithm %in% c("A2", "IA2", "A2M") & !is.null(hier.config))) {
  #   stop("Please specify a configuration matrix for a hierarchical testing algorithm, or specify a row/column size for an array testing algorithm.\n")
  # }
  
  # the code below is for when hier.config is a list
  # if (algorithm %in% c("D2", "D3", "D4", "ID2", "ID3", "ID4")) {
  #   if (is.null(hier.config) | 
  #       (algorithm %in% c("D2", "ID2") & length(hier.config) != 1) | 
  #       (algorithm %in% c("D3", "ID3") & length(hier.config) != 2) | 
  #       (algorithm %in% c("D4", "ID4") & length(hier.config) != 3)) {
  #     stop("Please provide a list specifying the configuration for the hierarchical algorithm. Each item in the list corresponds to the group sizes for each stage of testing. The list should not include an item detailing the individual stage of testing.\n")
  #   } else {
  #     group.sz <- hier.config[1]
  #   }
  # } else if (algorithm %in% c("A2", "IA2", "A2M")) {
  #   if (is.null(rowcol.sz)) {
  #     stop("Please provide a row/column size for the array testing algorithm.\n")
  #   } else {
  #     group.sz <- sum(hier.config[1])
  #   }
  # }

  if (algorithm %in% c("D3", "ID3")) {
    stage2 <- get.pools(hier.config[2,])
    stage3 <- NULL
  } else if (algorithm %in% c("D4", "ID4")) {
    stage2 <- get.pools(hier.config[2,])
    stage3 <- get.pools(hier.config[3,])
  } else if (algorithm == "ID2") {
    stage2 <- get.pools(hier.config[1,])
    stage3 <- NULL
  } else if (algorithm %in% c("A2", "IA2", "A2M")) {
    stage2 <- NULL
    stage3 <- NULL
  }
  # if (algorithm %in% c("D3", "ID3")) {
  #   stage2 <- hier.config[2]
  #   stage3 <- NULL
  # } else if (algorithm %in% c("D4", "ID4")) {
  #   stage2 <- hier.config[2]
  #   stage3 <- hier.config[3]
  # } else if (algorithm == "ID2") {
  #   stage2 <- hier.config[1]
  #   stage3 <- NULL
  # }
  
  if (is.null(p) & is.null(probabilities)) {
    stop("Please specify an overall probability of disease using the 'p' argument, or specify a vector of individual probabilities using the 'probabilities' argument.\n")
  } else if (!is.null(p) & !is.null(probabilities)) {
    stop("You have specified both an overall probability of disease AND a vector of individual probabilities. Please specify only one option.\n")
  } else{
    if (!is.null(p) & length(p) > 1) {
      stop("You have specified a probability vector instead of an overall probability of disease. Please specify an overall probability of disease, and the probability vector will be generated based on the algorithm specified for each group size included in the range.\n")
    }
    if (!is.null(probabilities)) {
      if ((algorithm %in% c("D2", "D3", "D4", "ID2", "ID3", "ID4")) & 
          length(probabilities) != group.sz) {
        stop("The vector of individual probabilities is not the correct length. Please make sure that the length of the probability vector is the same as the specified group size.\n")
      } else if ((algorithm %in% c("A2", "A2M", "IA2")) & 
                 length(probabilities) != group.sz^2) {
        stop("The vector of individual probabilities is not the correct length. Please make sure that the length of the probability vector is the same as the overall array size (the square of the specified group size).\n")
      }
      if ((algorithm %in% c("D2", "D3", "D4", "A2", "A2M")) & 
          all.equal(probabilities, rep(probabilities[1], length(probabilities))) != TRUE) {
        stop("You have specified a heterogeneous probability vector for a non-informative algorithm. Please specify a homogeneous probability vector using the 'probabilities' argument or specify an overall probability of disease using the 'p' argument.\n")
      }
    }
  }
  
  Se <- generate.acc(algorithm = algorithm, diseases = 1, value = Se, label = "sens")
  Sp <- generate.acc(algorithm = algorithm, diseases = 1, value = Sp, label = "spec")
  
  # check for an initial group when using hierarchical testing 
  if (algorithm %in% c("D2", "D3", "D4", "ID3", "ID4")) {
    if (all.equal(hier.config[1,], rep(1, ncol(hier.config))) != TRUE) {
      stop("Please specify a configuration where all individuals in the algorithm are tested together in the first stage of testing.\n")
    }
  } else if (algorithm == "ID2") {
    if (all.equal(hier.config[1,], rep(1, ncol(hier.config))) == TRUE) {
      stop("Please specify a configuration for informative two-stage hierarchical (PSOD) testing.\n")
    }
  }
  # the code below is for when hier.config is a list
  # if (algorithm == "D2" & length(hier.config[1]) > 1) {
  #   stop("Please provide a list specifying the configuration of the hierarchical algorithm. For non-informative two-stage hierarchical testing, the list should consist of a single component specifying the initial group size.\n")
  # }
  # if (algorithm %in% c("D3", "ID3") & length(hier.config[1]) > 1) {
  #   stop("Please provide a list specifying the configuration of the hierarchical algorithm. For three-stage hierarchical testing, the list should consist of two components. The first component specifies the initial group size and the second component specifies the group sizes for the second stage of testing.\n")
  # }
  # if (algorithm %in% c("D4", "ID4") & length(hier.config[1]) > 1) {
  #   stop("Please provide a list specifying the configuration of the hierarchical algorithm. For four-stage hierarchical testing, the list should consist of three components. The first component specifies the initial group size. The second and third components specify the group sizes for the second and third stages of testing, respectively.\n")
  # }
  # if (algorithm == "ID2" & length(hier.config[1]) <= 1) {
  #   stop("Please provide a list specifying the configuration of the hierarchical algorithm. For informative two-stage hierarchical testing, the list should consist of a single component specifying the group sizes for the first stage of testing.\n")
  # }
  
  # check the minimum and maximum group sizes
  if (group.sz < 3) {
    if (algorithm %in% c("D2", "D3", "D4", "ID2", "ID3", "ID4")) {
      stop("Please specify a configuration with an initial group size of at least 3.\n")
    } else if (algorithm %in% c("A2", "IA2", "A2M")) {
      stop("Please specify a row/column size of at least 3.\n")
    }
  } 
  # the code below is for when hier.config is a list
  # if ((algorithm %in% c("D2", "D3", "D4", "ID3", "ID4") & hier.config[1] < 3) | 
  #     (algorithm == "ID2" & sum(hier.config[1]) < 3)) {
  #   stop("Please specify a configuration with an initial group size of at least 3.\n")
  # }
  
  if (group.sz >= 50) {
    if (algorithm %in% c("D3", "D4", "ID2", "ID3", "ID4")) {
      warning("You have specified a configuration with a group size of 50 or larger. This function may take a VERY long time to run. Press 'ESC' if you wish to cancel the submitted statements.\n")
    } else if (algorithm %in% c("A2", "IA2", "A2M")) {
      warning("You have specified a row/column size of 50 or larger. This function may take a VERY long time to run. Press 'ESC' if you wish to cancel the submitted statements.\n")
    }
  } 
  # the code below is for when hier.config is a list
  # if ((algorithm %in% c("D3", "D4", "ID3", "ID4") & hier.config[1] >= 50) | 
  #     (algorithm == "ID2" & sum(hier.config[1]) >= 50)) {
  #   warning("You have specified a configuration with a group size of 50 or larger. This function may take a VERY long time to run. Press 'ESC' if you wish to cancel the submitted statements.\n")
  # }

  if (is.null(a)) {
    if (algorithm %in% c("D2", "D3", "D4", "ID2", "ID3", "ID4")) {
      a <- 1:group.sz
    } else if (algorithm %in% c("A2", "A2M", "IA2")) {
      a <- 1:group.sz^2
    }
  }
  
  # call function for non-informative two-stage hierarchical (Dorfman) testing
  if (algorithm == "D2") {
    if (!is.null(p)) {
      results <- NI.Dorf.calc1(p = p, Se = Se, Sp = Sp, 
                               group.sz = group.sz, a = a, 
                               trace = trace, print.time = print.time)
    } else if (!is.null(probabilities)) {
      results <- NI.Dorf.calc1(p = probabilities, Se = Se, Sp = Sp, 
                               group.sz = group.sz, a = a, 
                               trace = trace, print.time = print.time)
    }
  }
  
  # call function for non-informative three-stage hierarchical testing
  if (algorithm == "D3") {
    if (!is.null(p)) {
      results <- NI.D3.calc1(p = p, Se = Se, Sp = Sp, 
                             group.sz = group.sz, pool.szs = stage2, a = a, 
                             trace = trace, print.time = print.time)
    } else if (!is.null(probabilities)) {
      results <- NI.D3.calc1(p = probabilities, Se = Se, Sp = Sp, 
                             group.sz = group.sz, pool.szs = stage2, a = a, 
                             trace = trace, print.time = print.time)
    }
  } 
  
  # call function for non-informative four-stage hierarchical testing
  if (algorithm == "D4") {
    if (!is.null(p)) {
      results <- NI.D4.calc1(p = p, Se = Se, Sp = Sp, 
                             group.sz = group.sz, stage2 = stage2, 
                             stage3 = stage3, a = a, trace = 
                               trace, print.time = print.time)
    } else if (!is.null(probabilities)) {
      results <- NI.D4.calc1(p = probabilities, Se = Se, Sp = Sp, 
                             group.sz = group.sz, stage2 = stage2, 
                             stage3 = stage3, a = a, trace = 
                               trace, print.time = print.time)
    }
  } 
  
  # call function for non-informative square array testing without master pooling
  if (algorithm == "A2") {
    if (!is.null(p)) {
      results <- NI.Array.calc1(p = p, Se = Se, Sp = Sp, 
                                group.sz = group.sz, a = a, 
                                trace = trace, print.time = print.time)
    } else if (!is.null(probabilities)) {
      results <- NI.Array.calc1(p = probabilities, Se = Se, Sp = Sp, 
                                group.sz = group.sz, a = a, 
                                trace = trace, print.time = print.time)
    }
  }
  
  # call function for non-informative square array testing with master pooling
  if (algorithm == "A2M") {
    if (!is.null(p)) {
      results <- NI.A2M.calc1(p = p, Se = Se, Sp = Sp, 
                              group.sz = group.sz, a = a, 
                              trace = trace, print.time = print.time)
    } else if (!is.null(probabilities)) {
      results <- NI.A2M.calc1(p = probabilities, Se = Se, Sp = Sp, 
                              group.sz = group.sz, a = a, 
                              trace = trace, print.time = print.time)
    }
  }
  
  # call function for informative two-stage hierarchical (Dorfman) testing
  if (algorithm == "ID2") {
    if (!is.null(p)) {
      results <- Inf.Dorf.calc1(p = p, Se = Se, Sp = Sp, 
                                group.sz = group.sz, pool.szs = stage2, 
                                alpha = alpha, a = a, trace = 
                                  trace, print.time = print.time, ...)
    } else if (!is.null(probabilities)) {
      results <- Inf.Dorf.calc1(p = probabilities, Se = Se, Sp = Sp, 
                                group.sz = group.sz, pool.szs = stage2, 
                                alpha = alpha, a = a, trace = 
                                  trace, print.time = print.time, ...)
    }
  }
  
  # call function for informative three-stage hierarchical testing
  if (algorithm == "ID3") {
    if (!is.null(p)) {
      results <- Inf.D3.calc1(p = p, Se = Se, Sp = Sp, 
                              group.sz = group.sz, pool.szs = stage2, 
                              alpha = alpha, a = a, 
                              trace = trace, print.time = print.time, ...)
    } else if (!is.null(probabilities)) {
      results <- Inf.D3.calc1(p = probabilities, Se = Se, Sp = Sp, 
                              group.sz = group.sz, pool.szs = stage2, 
                              alpha = alpha, a = a, 
                              trace = trace, print.time = print.time, ...)
    }
  }
  
  # call function for informative four-stage hierarchical testing
  if (algorithm == "ID4") {
    if (!is.null(p)) {
      results <- Inf.D4.calc1(p = p, Se = Se, Sp = Sp, 
                             group.sz = group.sz, stage2 = stage2, 
                             stage3 = stage3, alpha = alpha, 
                             a = a, trace = trace, print.time = print.time, ...)
    } else if (!is.null(probabilities)) {
      results <- Inf.D4.calc1(p = probabilities, Se = Se, Sp = Sp, 
                             group.sz = group.sz, stage2 = stage2, 
                             stage3 = stage3, alpha = alpha, 
                             a = a, trace = trace, print.time = print.time, ...)
    }
  } 
  
  # call function for informative square array testing without master pooling
  if (algorithm == "IA2") {
    if (!is.null(p)) {
      results <- Inf.Array.calc1(p = p, Se = Se, Sp = Sp, 
                                 group.sz = group.sz, alpha = alpha, a = a, 
                                 trace = trace, print.time = print.time, ...)
    } else if (!is.null(probabilities)) {
      results <- Inf.Array.calc1(p = probabilities, Se = Se, Sp = Sp, 
                                 group.sz = group.sz, alpha = alpha, a = a, 
                                 trace = trace, print.time = print.time, ...)
    }
  }
  
  class(results) <- "opChar"
  results
}




###############################################################################
#' @rdname operatingCharacteristics1
opChar1 <- operatingCharacteristics1




# Summary function for operatingCharacteristics1() and 
#   operatingCharacteristics2()
###############################################################################

#' @title Summary method for operating characteristics results
#' 
#' @description Produce a summary list for objects of class 
#' \kbd{"opChar"} returned by \code{\link{operatingCharacteristics1}} 
#' (\kbd{opChar1}) or \code{\link{operatingCharacteristics2}} 
#' (\kbd{opChar2}).
#' 
#' @param object an object of class \kbd{"opChar"}, providing the calculated 
#' operating characteristics for a group testing algorithm.
#' @param ... currently not used.
#' 
#' @details This function produces a summary list for objects of 
#' class \kbd{"opChar"} returned by \code{\link{operatingCharacteristics1}} 
#' (\kbd{opChar1}) or \code{\link{operatingCharacteristics2}} 
#' (\kbd{opChar2}). It formats the testing configuration, expected number 
#' of tests, expected number of tests per individual, and accuracy measures. 
#' 
#' The \kbd{Configuration} component of the result
#' gives the testing configuration, which may include the group sizes for 
#' each stage of a hierarchical testing algorithm or the row/column size and
#' array size for an array testing algorithm. The \kbd{Tests} component 
#' of the result gives the expected number of tests and the expected 
#' number of tests per individual for the algorithm. 
#' 
#' The \kbd{Accuracy} component gives the individual accuracy measures for each 
#' individual in \kbd{object} and the overall accuracy measures for the 
#' algorithm. Accuracy measures included are the pooling sensitivity, pooling 
#' specificity, pooling positive predictive value, and pooling negative 
#' predictive value. The overall accuracy measures displayed are weighted 
#' averages of the corresponding individual accuracy measures for all individuals 
#' in the algorithm. Expressions for these averages are provided in the 
#' Supplementary Material for Hitt et al. (2019). For more information, see the 
#' 'Details' section for the \code{\link{operatingCharacteristics1}} 
#' (\kbd{opChar1}) or \code{\link{operatingCharacteristics2}} (\kbd{opChar2}) 
#' function.
#' 
#' @return \kbd{summary.opChar} returns an object of class \kbd{"summary.opChar"}, 
#' a list containing:
#' \item{Algorithm}{character string specifying the name of the group testing 
#' algorithm.}
#' \item{Configuration}{matrix detailing the configuration from \kbd{object}. 
#' For hierarchical testing, this includes the group sizes for each stage of 
#' testing. For array testing, this includes the array dimension (row/column 
#' size) and the array size (the total number of individuals in the array).}
#' \item{Tests}{matrix detailing the expected number of tests and expected 
#' number of tests per individual from \kbd{object}}.
#' \item{Accuracy}{a list containing:
#' \describe{
#' \item{Individual}{matrix detailing the accuracy measures for each 
#' individual from \kbd{object} (for objects returned by \code{\link{opChar1}}).}
#' \item{Disease 1 Individual}{matrix detailing the accuracy measures 
#' pertaining to disease 1 for each individual from \kbd{object} 
#' (for objects returned by \code{\link{opChar2}}).}
#' \item{Disease 2 Individual}{matrix detailing the accuracy measures 
#' pertaining to disease 2 for each individual from \kbd{object}
#' (for objects returned by \code{\link{opChar2}}).}
#' \item{Overall}{matrix detailing the overall accuracy measures for 
#' the algorithm from \kbd{object}.}}}
#' 
#' @author Brianna D. Hitt
#' 
#' @seealso
#' \code{\link{operatingCharacteristics1}} (\kbd{opChar1}) and 
#' \code{\link{operatingCharacteristics2}} (\kbd{opChar2}) for creating 
#' an object of class \kbd{"opChar"}.
#' 
#' @examples 
#' # Calculate the operating characteristics for 
#' #   non-informative four-stage hierarchical testing.
#' config.mat <- matrix(data=c(rep(1, 24), rep(1, 16), rep(2, 8), 
#'                             rep(1, 8), rep(2, 8), rep(3, 4), 
#'                             rep(4, 2), rep(5, 2), 1:24), 
#'                      nrow=4, ncol=24, byrow=TRUE)
#' calc1 <- opChar1(algorithm="D4", p=0.01, Se=0.99, Sp=0.99, 
#'                  hier.config=config.mat, a=c(1, 9, 17, 21, 23))
#' summary(calc1)
#' 
#' # Calculate the operating characteristics for 
#' #   informative array testing without master pooling.
#' calc2 <- opChar1(algorithm="IA2", p=0.025, alpha=0.5, 
#'                  Se=0.95, Sp=0.99, rowcol.sz=12)
#' summary(calc2)
#' 
#' # Calculate the operating characteristics for 
#' #   informative two-stage hierarchical testing.
#' config.mat <- matrix(data=c(rep(1, 5), rep(2, 4), 1, 1:10), 
#'                      nrow=2, ncol=10, byrow=TRUE)
#' Se <- matrix(data=c(rep(0.95, 2), rep(0.99, 2)), 
#'              nrow=2, ncol=2, byrow=FALSE)
#' Sp <- matrix(data=c(rep(0.96, 2), rep(0.98, 2)), 
#'              nrow=2, ncol=2, byrow=FALSE)
#' calc3 <- opChar2(algorithm="ID2", alpha=c(18.25, 0.75, 0.75, 0.25), 
#'                  Se=Se, Sp=Sp, hier.config=config.mat)
#' summary(calc3)
#' 
#' # Calculate the operating characteristics for 
#' #   non-informative array testing with master pooling.
#' calc4 <- opChar2(algorithm="A2M", p.vec=c(0.92, 0.05, 0.02, 0.01), 
#'                  Se=rep(0.95, 2), Sp=rep(0.99, 2), rowcol.sz=8)
#' summary(calc4)

summary.opChar <- function(object, ...) {
  
  # algorithm
  algorithm <- object$algorithm
  cat("\nAlgorithm:", algorithm, "\n\n")
  
  # configuration
  config <- object$Config
  # create a data frame detailing the configuration
  stage1 <- config[[1]]
  stage2 <- NULL
  stage3 <- NULL
  stage4 <- NULL
  if (length(config) > 1) {
    # three-stage and informative two-stage hierarchical testing
    if (length(config[[2]]) > 1) {
      stage2 <- config[[2]][1]
      for(i in 2:length(config[[2]])) {
        stage2 <- paste(stage2, config[[2]][i], sep=",")
      }
    } else{
      # array testing and non-informative two-stage hierarchical testing
      stage2 <- config[[2]]
    }
  }
  # four-stage hierarchical testing
  if (length(config) > 2) {
    if (length(config[[3]]) > 1) {
      stage3 <- config[[3]][1]
      for(i in 2:length(config[[3]])) {
        stage3 <- paste(stage3, config[[3]][i], sep=",")
      }
    } else{
      stage3 <- config[[3]]
    }
  }
  # five-stage hierarchical testing
  if (length(config) > 3) {
    if (length(config[[4]]) > 1) {
      stage4 <- config[[4]][1]
      for(i in 2:length(config[[4]])) {
        stage4 <- paste(stage4, config[[4]][i], sep=",")
      }
    } else{
      stage4 <- config[[4]]
    }
  }
  config <- rbind(stage1, stage2, stage3, stage4)
  if (grepl("Informative two-stage", algorithm)) {
    dimnames(config) <- list(c("Block size", "Group sizes"), "")
  } else if (grepl("hierarchical", algorithm)) {
    dimnames(config) <- list(c("Stage 1", "Stage 2", "Stage 3", "Stage 4")[1:length(config)], "")
  } else if (grepl("array", algorithm)) {
    dimnames(config) <- list(c("Row/column size", "Array size"), "")
  }
  
  cat("Testing configuration:\n")
  for (i in 1:length(config)) {
    cat(paste0(rownames(config)[i], ": ", config[i], "\n"))
  }
  
  ET <- format(round(object$ET, 2), nsmall = 2)
  value <- format(round(object$value, 4), nsmall = 4)
  # create a matrix detailing the expected number of tests
  tests <- matrix(data = c(ET, value), nrow=2, ncol=1)
  rownames(tests) <- c("Expected number of tests", 
                       "Expected number of tests per individual")
  cat("\nExpected number of tests:", ET)
  cat("\nExpected number of tests per individual:", value)
  
  # create a matrix for accuracy measures
  if (class(object$Accuracy) == "matrix") {
    ind.acc <- NULL
    overall.acc <- as.data.frame(format(round(object$Accuracy, 4), nsmall = 4))
  } else if (length(object$Accuracy) == 2) {
    ind.acc <- object$Accuracy$Individual
    ind.acc[, 1:4] <- format(round(ind.acc[, 1:4], 4), nsmall = 4)
    colnames(ind.acc) <- c("PSe", "PSp", "PPPV", "PNPV", "Individuals")
    
    cat("\n\nAccuracy for individuals:\n")
    print(as.data.frame(ind.acc))
    
    overall.acc <- as.data.frame(format(round(object$Accuracy$Overall, 4), nsmall = 4))
  } else if (length(object$Accuracy) > 2) {
    ind.acc.dis1 <- object$Accuracy$'Disease 1 Individual'
    ind.acc.dis1[, 1:4] <- format(round(ind.acc.dis1[, 1:4], 4), nsmall = 4)
    colnames(ind.acc.dis1) <- c("PSe", "PSp", "PPPV", "PNPV", "Individuals")
    
    ind.acc.dis2 <- object$Accuracy$'Disease 2 Individual'
    ind.acc.dis2[, 1:4] <- format(round(ind.acc.dis2[, 1:4], 4), nsmall = 4)
    colnames(ind.acc.dis2) <- c("PSe", "PSp", "PPPV", "PNPV", "Individuals")
    
    cat("\n\nDisease 1 accuracy for individuals:\n")
    print(as.data.frame(ind.acc.dis1))
    cat("\nDisease 2 accuracy for individuals:\n")
    print(as.data.frame(ind.acc.dis2))
    
    overall.acc <- as.data.frame(format(round(object$Accuracy$Overall, 4), nsmall = 4))
  }
  colnames(overall.acc) <- c("PSe", "PSp", "PPPV", "PNPV")
  
  cat("\nOverall accuracy of the algorithm:\n")
  print(as.data.frame(overall.acc))
  
  cat("\nPSe denotes the pooling sensitivity.\n")
  cat("PSp denotes the pooling specificity.\n")
  cat("PPPV denotes the pooling positive predictive value.\n")
  cat("PNPV denotes the pooling negative predictive value.\n")
  
  if (class(object$Accuracy) == "matrix" | 
      length(object$Accuracy) == 2) {
    res <- list("Algorithm" = algorithm, 
                "Configuration" = config, 
                "Tests" = tests, 
                "Accuracy" = list("Individual" = ind.acc, 
                                  "Overall" = overall.acc))
  } else {
    res <- list("Algorithm" = algorithm, 
                "Configuration" = config, 
                "Tests" = tests, 
                "Accuracy" = list("Disease 1 Individual" = ind.acc.dis1,
                                  "Disease 2 Individual" = ind.acc.dis2,
                                  "Overall" = overall.acc))
  }
  
  class(res) <- "summary.opChar"
  invisible(res)
}




# Start supporting functions for operating characteristics functions
########################################################################
# This function converts a numeric vector to a character string, 
#   where values are separated by commas with no spaces
# This function is to be used when converting a vector of individual
#   indices or group sizes in a group testing algorithm to a character 
#   string to be displayed in a single column in a matrix/data frame

vec2string <- function(vec) {
  res.string <- vec[1]
  if (length(vec) > 1) {
    for(i in 2:length(vec)) {
      res.string <- paste(res.string, vec[i], sep=",")
    }
  }
  res.string
}




# This function takes a matrix of individual accuracy measures, 
#   and finds the indices for each unique row of results
# The output is a matrix of the unique individual accuracy measures
#   with an additional column specifying the indices for all individuals
#   with the same accuracy measure values

get.unique.index <- function(results, col.num, rowlabel = NULL) {
  
  if (is.null(dim(results))) {
    results <- matrix(data = results, nrow=1, ncol=length(results), 
                      dimnames = list(rowlabel, names(results)))
  }
  
  results <- as.data.frame(results)
  rows <- rownames(results)
  
  # account for individuals with same results as the first row
  if (nrow(results) == 1) {
    index <- rows
  } else {
    index <- rows[which(results[,col.num] == results[1,col.num])]
  }
  
  index.string <- vec2string(index)
  
  new <- cbind(results[index[1],], as.character(index.string))
  included <- as.character(index)
  
  # keep going, until all individuals are accounted for
  while(length(included) < nrow(results)) {
    new.start <- as.character(min(as.numeric(suppressWarnings(rows[rows!=included]))))
    index <- rows[which(results[,col.num] == results[new.start,col.num])]
    index.string <- vec2string(index)
    new <- rbind(new, cbind(results[index[1],], as.character(index.string)))
    included <- as.character(sort(as.numeric(c(included, index))))
  }
  
  colnames(new)[length(colnames(new))] <- c("individuals")
  rownames(new) <- NULL
  new
}




# Brianna Hitt - 03.03.2020
#   Check whether all individual have equal accuracy measures
#   If not, list all the indices for each unique row
#   If yes, the individuals column will read "All"

check.all.equal <- function(results, col.num) {

  results <- as.data.frame(results)
  rows <- rownames(results)
  
  # account for individuals with same results as the first row
  index <- rows[which(results[,col.num] == results[1,col.num])]
  
  if (length(index) == nrow(results)) {
    index.string <- "All"
    new <- cbind(results[1,], as.character(index.string))
    colnames(new)[length(colnames(new))] <- c("individuals")
    rownames(new) <- NULL
    new
  }
}




# This function creates a vector of group sizes from a group membership matrix
get.pools <- function(group.nums){
  pools.table <- table(group.nums)
  pool.szs <- pools.table[[1]]
  for (i in 2:dim(pools.table)) {
    pool.szs <- c(pool.szs, pools.table[[i]])
  }
  pool.szs
}

#