# Start Operating Characteristic function for two diseases
# ##############################################################################
# Brianna Hitt - 02.10.2020
# Operating characteristics functions for two diseases
# function to calculate operating characteristics for hierarchical and 
#   array testing with two diseases, for a specified configuration
# Se and Sp are matrices of sensitivity/specificity values
#   corresponding to each disease and each stage of testing
# Brianna Hitt - 02.20.2020
# Added four- and five-stage hierarchical testing as options
# Changed the input from group.sz and pool.szs to hier.config (a group 
#   membership matrix) and rowcol.sz

# Start operatingCharacteristics2() function
###############################################################################
#' @title Calculate operating characteristics for group testing algorithms 
#' that use a multiplex assay for two diseases
#' 
#' @description Calculate operating characteristics, such as the expected 
#' number of tests, for a specified testing configuration using 
#' non-informative and informative hierarchical and array-based group testing 
#' algorithms. Multiplex assays for two diseases are used at each stage of the 
#' algorithms.
#'  
#' @param algorithm character string defining the group testing 
#' algorithm to be used. Non-informative testing options include two-stage 
#' hierarchical ("\kbd{D2}"), three-stage hierarchical ("\kbd{D3}"), four-stage 
#' hierarchical ("\kbd{D4}"), five-stage hierarchical ("\kbd{D5}"), 
#' square array testing without master pooling ("\kbd{A2}"), and square array 
#' testing with master pooling ("\kbd{A2M}"). Informative testing options 
#' include two-stage hierarchical ("\kbd{ID2}"), three-stage hierarchical 
#' ("\kbd{ID3}"), four-stage hierarchical ("\kbd{ID4}"), and five-stage 
#' hierarchical ("\kbd{ID5}") testing.
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
#' @param alpha a vector containing positive shape parameters of the Dirichlet 
#' distribution (for informative testing only). The vector will be used to 
#' generate a heterogeneous matrix of joint probabilities for each individual. 
#' The vector must have length 4. Further details are given under 'Details'.
#' Only one of \kbd{p.vec}, \kbd{probabilities}, or \kbd{alpha} should be specified.
#' @param Se matrix of sensitivity values, where one value is given for each 
#' disease (or infection) at each stage of testing. The rows of the matrix 
#' correspond to each disease \eqn{k=1,...,K}, and the columns of the 
#' matrix correspond to each stage of testing \eqn{s=1,...,S}. If a vector of 
#' \eqn{K} values is provided, the sensitivity values associated with disease 
#' \eqn{k} are assumed to be equal to the \eqn{k}th value in the vector for 
#' all stages of testing. Further details are given under 'Details'.
#' @param Sp a matrix of specificity values, where one value is given for each 
#' disease (or infection) at each stage of testing. The rows of the matrix 
#' correspond to each disease \eqn{k=1,...,K}, and the columns of the 
#' matrix correspond to each stage of testing \eqn{s=1,...,S}. If a vector of 
#' \eqn{K} values is provided, the specificity values associated with disease 
#' \eqn{k} are assumed to be equal to the \eqn{k}th value in the vector for 
#' all stages of testing. Further details are given under 'Details'.
#' @param ordering a matrix detailing the ordering for the binary responses of 
#' the diseases. The columns of the matrix correspond to each disease and the 
#' rows of the matrix correspond to each of the 4 sets of binary responses for 
#' two diseases. This ordering is used with the joint probabilities. The 
#' default ordering is (p_00, p_10, p_01, p_11).
#' @param hier.config a matrix specifying the configuration for a hierarchical 
#' testing algorithm. The rows correspond to the stages of testing, the columns 
#' correspond to each individual to be tested, and the cell values 
#' specify the group number of each individual at each stage. Further details 
#' are given under 'Details'. For array testing algorithms, this argument will be 
#' ignored.
#' @param rowcol.sz the row/column size for array testing algorithms. For 
#' hierarchical testing algorithms, this argument will be ignored.
#' @param a a vector containing indices indicating which individuals to 
#' calculate individual accuracy measures for. If \kbd{NULL}, individual accuracy 
#' measures will be displayed for all individuals in the algorithm.
#' @param print.time a logical value indicating whether the length of time 
#' for calculations should be printed. The default is \kbd{TRUE}.
#' @param ... additional arguments to be passed to functions for hierarchical 
#' testing with multiplex assays for two diseases.
#' 
#' @details This function computes the operating characteristics for standard 
#' group testing algorithms with a multiplex assay that tests for two diseases.  
#' Calculations for hierarchical group testing algorithms are performed as 
#' described in Bilder et al. (2019) and calculations for array-based group 
#' testing algorithms are performed as described in Hou et al. (2019). 
#' 
#' Available algorithms include two-, three-, four-, and five-stage hierarchical 
#' testing and array testing with and without master pooling. Both non-informative 
#' and informative group testing settings are allowed for hierarchical algorithms. 
#' Only non-informative group testing settings are allowed for array testing 
#' algorithms. Operating characteristics calculated are expected number of tests, 
#' pooling sensitivity, pooling specificity, pooling positive predictive value, 
#' and pooling negative predictive value for each individual.
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
#' The matrix specified  by \kbd{hier.config} defines the hierarchical group testing 
#' algorithm for \eqn{I} individuals. The rows of the matrix correspond to the stages 
#' \eqn{s=1,...,S} in the testing algorithm, and the columns correspond to individuals 
#' \eqn{i=1,...I}. The cell values within the matrix represent the group number of 
#' individual \eqn{i} at stage \eqn{s}. For three-stage, four-stage, five-stage, and 
#' non-informative two-stage hierarchical testing, the first row of the matrix 
#' consists of all ones. This indicates that all individuals in the algorithm are 
#' tested together in a single group in the first stage of testing. For informative 
#' two-stage hierarchical testing, the initial group (block) is not tested. Thus, 
#' the first row of the matrix consists of the group numbers for each individual in 
#' the first stage of testing. For all hierarchical algorithms, the final row of the 
#' matrix denotes individual testing. Individuals who are not tested in a particular 
#' stage are represented by "NA" (e.g., an individual tested in a group of size 1 in 
#' the second stage of testing would not be tested again in a third stage of testing). 
#' It is important to note that this matrix represents the testing that could be 
#' performed if each group tests positively at each stage prior to the last. For 
#' more details on this matrix (called a group membership matrix), see Bilder et al. 
#' (2019). 
#' 
#' For array testing without master pooling, the \kbd{rowcol.sz} specified
#' represents the row/column size for initial (stage 1) testing. For array testing
#' with master pooling, the \kbd{rowcol.sz} specified represents the row/column size
#' for stage 2 testing. This is because the master pool size is the overall array 
#' size, given by the square of the row/column size.
#' 
#' The displayed overall pooling sensitivity, pooling specificity, pooling 
#' positive predictive value, and pooling negative predictive value are 
#' weighted averages of the corresponding individual accuracy measures for all 
#' individuals within the initial group (or block) for a hierarchical algorithm, 
#' or within the entire array for an array-based algorithm. 
#' Expressions for these averages are provided in the Supplementary Material 
#' for Hitt et al. (2019). These expressions are based on accuracy definitions 
#' given by Altman and Bland (1994a, 1994b).
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
#' \item{Config}{a list specifying elements of the specified testing configuration, 
#' which may include:
#' \describe{
#' \item{Stage1}{group size for the first stage of hierarchical testing, if 
#' applicable.}
#' \item{Stage2}{group sizes for the second stage of hierarchical testing, if 
#' applicable.}
#' \item{Stage3}{group sizes for the third stage of hierarchical testing, if 
#' applicable.}
#' \item{Stage4}{group sizes for the fourth stage of hierarchical testing, if 
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
#' \item{Accuracy}{a list containing:
#' \describe{
#' \item{Disease 1 Individual}{a matrix of accuracy measures, pertaining to the 
#' first disease, for each individual specified in \kbd{a}. The rows correspond 
#' to each unique set of accuracy measures in the algorithm. Individuals with 
#' the same set of accuracy measures are displayed together in a single row of 
#' the matrix. The columns correspond to the pooling sensitivity, pooling 
#' specificity, pooling positive predictive value, pooling negative predictive 
#' value, and the indices for the individuals in each row of the matrix. 
#' Individual accuracy measures are not displayed for array testing algorithms.}
#' \item{Disease 2 Individual}{a matrix of accuracy measures, pertaining to the 
#' second disease, for each individual specified in \kbd{a}. The rows correspond 
#' to each unique set of accuracy measures in the algorithm. Individuals with 
#' the same set of accuracy measures are displayed together in a single row of 
#' the matrix. The columns correspond to the pooling sensitivity, pooling 
#' specificity, pooling positive predictive value, pooling negative predictive 
#' value, and the indices for the individuals in each row of the matrix. 
#' Individual accuracy measures are not displayed for array testing algorithms.}
#' \item{Overall}{a matrix of overall accuracy measures for the algorithm. 
#' The rows correspond to each disease. The columns correspond to the pooling 
#' sensitivity, pooling specificity, pooling positive predictive value, and 
#' pooling negative predictive value for the overall algorithm. Further 
#' details are given under 'Details'.}}}
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
#' @family operating characteristic functions
#' @family multiplex testing functions
#' 
#' @examples 
#' # Calculate the operating characteristics for 
#' #   non-informative two-stage hierarchical 
#' #   (Dorfman) testing.
#' config.mat <- matrix(data = c(rep(1, 24), 1:24), 
#'                      nrow = 2, ncol = 24, byrow = TRUE)
#' Se <- matrix(data=c(0.95, 0.95, 0.95, 0.95),
#'              nrow=2, ncol=2,
#'              dimnames=list(Infection=1:2, Stage=1:2))
#' Sp <- matrix(data=c(0.99, 0.99, 0.99, 0.99),
#'              nrow=2, ncol=2,
#'              dimnames=list(Infection=1:2, Stage=1:2))
#' opChar2(algorithm="D2", p.vec=c(0.90, 0.04, 0.04, 0.02),
#'          Se=Se, Sp=Sp, hier.config=config.mat)
#' opChar2(algorithm="D2", p.vec=c(0.90, 0.04, 0.04, 0.02),
#'          Se=Se, Sp=Sp, hier.config=config.mat, a=c(1, 13, 24), 
#'          print.time = FALSE)
#'                             
#' # Calculate the operating characteristics for informative 
#' #   two-stage hierarchical (Dorfman) testing.
#' # A matrix of joint probabilities for each individual is 
#' #   generated using the Dirichlet distribution.
#' config.mat <- matrix(data = c(rep(1, 5), rep(2, 4), 3, 1:9, NA), 
#'                      nrow = 2, ncol = 10, byrow = TRUE)
#' Se <- matrix(data=c(0.95, 0.95, 0.99, 0.99),
#'              nrow=2, ncol=2,
#'              dimnames=list(Infection=1:2, Stage=1:2))
#' Sp <- matrix(data=c(0.96, 0.96, 0.98, 0.98),
#'              nrow=2, ncol=2,
#'              dimnames=list(Infection=1:2, Stage=1:2))
#' set.seed(8791)
#' opChar2(algorithm="ID2", alpha=c(18.25, 0.75, 0.75, 0.25),
#'          Se=Se, Sp=Sp, hier.config=config.mat)
#' # Equivalent code using a heterogeneous matrix of joint
#' #   probabilities for each individual
#' set.seed(8791)
#' p.unordered <- t(rBeta2009::rdirichlet(n = 10, 
#'                             shape = c(18.25, 0.75, 0.75, 0.25)))
#' p.ordered <- p.unordered[, order(1 - p.unordered[1,])]
#' opChar2(algorithm="ID2", probabilities=p.ordered,
#'         Se=Se, Sp=Sp, hier.config=config.mat)
#'          
#' # Calculate the operating characteristics for 
#' #   non-informative three-stage hierarchical testing.
#' config.mat <- matrix(data = c(rep(1, 10), rep(1, 5), 
#'                               rep(2, 4), 3, 1:9, NA), 
#'                      nrow = 3, ncol = 10, byrow = TRUE)
#' Se <- matrix(data=rep(0.95, 6), nrow=2, ncol=3,
#'              dimnames=list(Infection=1:2, Stage=1:3))
#' Sp <- matrix(data=rep(0.99, 6), nrow=2, ncol=3,
#'              dimnames=list(Infection=1:2, Stage=1:3))
#' opChar2(algorithm="D3", p.vec=c(0.95, 0.02, 0.02, 0.01),
#'          Se=Se, Sp=Sp, hier.config=config.mat)
#' opChar2(algorithm="D3", p.vec=c(0.95, 0.02, 0.02, 0.01),
#'          Se=Se, Sp=Sp, hier.config=config.mat, a=c(1, 6, 10))
#'     
#' # Calculate the operating characteristics for informative 
#' #   three-stage hierarchical testing. 
#' # A matrix of joint probabilities for each individual is 
#' #   generated using the Dirichlet distribution.
#' config.mat <- matrix(data = c(rep(1, 15), 
#'                               rep(c(1, 2, 3), each = 5), 1:15), 
#'                      nrow = 3, ncol = 15, byrow = TRUE)
#' Se <- matrix(data=rep(0.95, 6), nrow=2, ncol=3,
#'              dimnames=list(Infection=1:2, Stage=1:3))
#' Sp <- matrix(data=rep(0.99, 6), nrow=2, ncol=3,
#'              dimnames=list(Infection=1:2, Stage=1:3))
#' opChar2(algorithm="ID3", alpha=c(18.25, 0.75, 0.75, 0.25),
#'          Se=Se, Sp=Sp, hier.config=config.mat)
#'          
#' # Calculate the operating characteristics for 
#' #   non-informative four-stage hierarchical testing. 
#' config.mat <- matrix(data = c(rep(1, 12), rep(1, 6), rep(2, 6), 
#'                               rep(1, 4), rep(2, 2), rep(3, 3), 
#'                               rep(4, 3), 1:12), 
#'                      nrow = 4, ncol = 12, byrow = TRUE)
#' Se <- matrix(data=rep(0.95, 8), nrow=2, ncol=4,
#'              dimnames=list(Infection=1:2, Stage=1:4))
#' Sp <- matrix(data=rep(0.99, 8), nrow=2, ncol=4,
#'              dimnames=list(Infection=1:2, Stage=1:4))
#' opChar2(algorithm="D4", p.vec=c(0.92, 0.05, 0.02, 0.01),
#'          Se=Se, Sp=Sp, hier.config=config.mat)
#'
#' # Calculate the operating characteristics for informative 
#' #   five-stage hierarchical testing. 
#' # A matrix of joint probabilities for each individual is 
#' #   generated using the Dirichlet distribution.
#' config.mat <- matrix(data = c(rep(1, 20), rep(1, 10), rep(2, 10),
#'                               rep(c(1, 2, 3, 4), each = 5), 
#'                               rep(1, 3), rep(2, 2), rep(3, 3), 
#'                               rep(4, 2), rep(5, 3), rep(6, 2),
#'                               rep(7, 3), rep(8, 2), 1:20), 
#'                      nrow = 5, ncol = 20, byrow = TRUE)
#' Se <- matrix(data=rep(0.95, 10), nrow=2, ncol=5,
#'              dimnames=list(Infection=1:2, Stage=1:5))
#' Sp <- matrix(data=rep(0.99, 10), nrow=2, ncol=5,
#'              dimnames=list(Infection=1:2, Stage=1:5))
#' opChar2(algorithm="ID5", alpha=c(18.25, 0.75, 0.75, 0.25),
#'         Se=Se, Sp=Sp, hier.config=config.mat)
#'
#' # Calculate the operating characteristics for 
#' #   non-informative array testing without master pooling.
#' Se <- matrix(data=rep(0.95, 4), nrow=2, ncol=2,
#'              dimnames=list(Infection=1:2, Stage=1:2))
#' Sp <- matrix(data=rep(0.99, 4), nrow=2, ncol=2,
#'              dimnames=list(Infection=1:2, Stage=1:2))
#' opChar2(algorithm="A2", p.vec=c(0.90, 0.04, 0.04, 0.02),
#'          Se=Se, Sp=Sp, rowcol.sz=12)
#'                   
#' # Calculate the operating characteristics for 
#' #   non-informative array testing with master pooling.
#' Se <- matrix(data=rep(0.95, 6), nrow=2, ncol=3,
#'              dimnames=list(Infection=1:2, Stage=1:3))
#' Sp <- matrix(data=rep(0.99, 6), nrow=2, ncol=3,
#'              dimnames=list(Infection=1:2, Stage=1:3))
#' opChar2(algorithm="A2M", p.vec=c(0.90, 0.04, 0.04, 0.02),
#'          Se=Se, Sp=Sp, rowcol.sz=10)

# Brianna Hitt - 04.02.2020
# Changed cat() to warning()

operatingCharacteristics2 <- 
  function(algorithm, p.vec = NULL, probabilities = NULL, alpha = NULL, 
           Se, Sp, hier.config = NULL, rowcol.sz = NULL, 
           ordering = matrix(data = c(0,1,0,1,0,0,1,1), nrow = 4, ncol = 2), 
           a = NULL, print.time = TRUE, ...){
    
    # make sure that all necessary information is included in the correct format
    if (!(algorithm %in% c("D2", "D3", "D4", "D5", "A2", "A2M", "ID2", "ID3", "ID4", "ID5"))) {
      stop("Please specify one of the following algorithms: D2, ID2, D3, ID3, D4, ID4, D5, ID5, A2, A2M.")
    }
    
    if (algorithm %in% c("D2", "D3", "D4", "D5", "ID2", "ID3", "ID4", "ID5")) {
      if (is.null(hier.config) | 
          (algorithm %in% c("D2", "ID2") & nrow(hier.config) != 2) | 
          (algorithm %in% c("D3", "ID3") & nrow(hier.config) != 3) | 
          (algorithm %in% c("D4", "ID4") & nrow(hier.config) != 4) | 
          (algorithm %in% c("D5", "ID5") & nrow(hier.config) != 5)) {
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
  
    # if (algorithm %in% c("D3", "ID3")) {
    #   stage2 <- get.pools(hier.config[2,])
    #   stage3 <- NULL
    #   stage4 <- NULL
    # } else if (algorithm %in% c("D4", "ID4")) {
    #   stage2 <- get.pools(hier.config[2,])
    #   stage3 <- get.pools(hier.config[3,])
    #   stage4 <- NULL
    # } else if (algorithm %in% c("D5", "ID5")) {
    #   stage2 <- get.pools(hier.config[2,])
    #   stage3 <- get.pools(hier.config[3,])
    #   stage4 <- get.pools(hier.config[4,])
    # } else if (algorithm == "ID2") {
    #   stage2 <- get.pools(hier.config[1,])
    #   stage3 <- NULL
    #   stage4 <- NULL
    # } else if (algorithm %in% c("A2", "IA2", "A2M")) {
    #   stage2 <- NULL
    #   stage3 <- NULL
    #   stage4 <- NULL
    # }

    # check that the joint probabilities are specified, either through p.vec, 
    #   probabilities, or alpha arguments
    if (is.null(p.vec) & is.null(probabilities) & is.null(alpha)) {
      if (algorithm %in% c("D2", "D3", "D4", "D5", "A2", "A2M")) {
        stop("Please specify an overall joint probability vector using the 'p.vec' argument, or specify a matrix of joint probabilities for each individual using the 'probabilities' argument.\n")
      } else if (algorithm %in% c("ID2", "ID3", "ID4", "ID5", "IA2")) {
        stop("Please specify a matrix of joint probabilities for each individual using the 'probabilities' argument, or specify a vector of shape parameters for the Dirichlet distribution using the 'alpha' argument.\n")
      }
    } else if (sum(!is.null(p.vec), !is.null(probabilities), !is.null(alpha)) > 1) {
      stop("You have specified more than one of the following arguments: p.vec, probabilities, alpha. Please specify only one option.\n")
    } else {
      if (!is.null(p.vec)) {
        if (length(p.vec) != 4) {
          stop("Please specify an overall joint probability vector of length 4 (one for each of the joint probabilities), and the matrix of individual probabilities will be generated based on the algorithm specified for each group size included in the range.\n")
        }
        if (algorithm %in% c("ID2", "ID3", "ID4", "ID5", "IA2")) {
          stop("You have specified an overall joint probability vector for an informative algorithm. Please specify a matrix of joint probabilities for each individual using the 'probabilities' argument, or specify shape parameters for the Dirichlet distribution using the 'alpha' argument.\n")
        }
        if (sum(p.vec) != 1.00) {
          stop("Please specify joint probabilities that sum to 1.")
        }
      } else if (!is.null(probabilities)) {
        if (length(group.sz) == 1) {
          if (dim(probabilities)[1] != 4 | 
              (algorithm %in% c("D2", "D3", "D4", "D5", "ID2", "ID3", "ID4", "ID5") & dim(probabilities)[2] != group.sz) | 
              (algorithm %in% c("A2", "IA2", "A2M") & dim(probabilities)[2] != group.sz^2)) {
            stop("Please specify a matrix of joint probabilities with the correct dimensions. Each row should correspond to one of the four joint probabilities. Each column should correspond to an individual in the algorithm.\n")
          }
          
          if ((algorithm %in% c("D2", "D3", "D4", "D5", "A2", "A2M")) & 
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
        if (algorithm %in% c("D2", "D3", "D4", "D5", "A2", "A2M")) {
          stop("You have specified a vector of shape parameters for the Dirichlet distribution using a non-informative algorithm. Please specify an overall vector of joint probabilities using the 'p.vec' argument, or specify a matrix of joint probabilities for each individual using the 'probabilities' argument.\n")
        }
      }
    }
    
    Se <- generate.acc(algorithm = algorithm, diseases = 2, value = Se, label = "sens")
    Sp <- generate.acc(algorithm = algorithm, diseases = 2, value = Sp, label = "spec")
    
    # check for an initial group when using hierarchical testing 
    if (algorithm %in% c("D2", "D3", "D4", "D5", "ID3", "ID4", "ID5")) {
      if (all.equal(hier.config[1,], rep(1, ncol(hier.config))) != TRUE) {
        stop("Please specify a configuration where all individuals in the algorithm are tested together in the first stage of testing.\n")
      }
    } else if (algorithm == "ID2") {
      if (all.equal(hier.config[1,], rep(1, ncol(hier.config))) == TRUE) {
        stop("Please specify a configuration for informative two-stage hierarchical (PSOD) testing.\n")
      }
    }
    
    # check the dimension of the ordering matrix
    if (dim(ordering)[1] != 4 | dim(ordering)[2] != 2) {
      stop("Please specify an ordering matrix with the correct dimensions. The columns should correspond to each disease. The rows should contain each of the four combinations of binary responses for the two diseases.\n")
    }
    
    # check the minimum and maximum group sizes
    if (group.sz < 3) {
      if (algorithm %in% c("D2", "D3", "D4", "D5", "ID2", "ID3", "ID4", "ID5")) {
        stop("Please specify a configuration with an initial group size of at least 3.\n")
      } else if (algorithm %in% c("A2", "IA2", "A2M")) {
        stop("Please specify a row/column size of at least 3.\n")
      }
    }
    if (group.sz >= 50) {
      if (algorithm %in% c("D3", "D4", "D5", "ID2", "ID3", "ID4", "ID5")) {
        warning("You have specified a configuration with a group size of 50 or larger. This function may take a VERY long time to run. Press 'ESC' if you wish to cancel the submitted statements.\n")
      } else if (algorithm %in% c("A2", "A2M", "IA2")) {
        warning("You have specified a row/column size of 50 or larger. This function may take a VERY long time to run. Press 'ESC' if you wish to cancel the submitted statements.\n")
      }
    }
    
    if (is.null(a)) {
      if (algorithm %in% c("D2", "D3", "D4", "D5", "ID2", "ID3", "ID4", "ID5")) {
        a <- 1:group.sz
      } else if (algorithm %in% c("A2", "A2M", "IA2")) {
        a <- 1:group.sz^2
      }
    }
    
    # call function for non-informative two-stage hierarchical (Dorfman) testing
    if (algorithm == "D2") {
      if (!is.null(p.vec)) {
        results <- NI.Dorf.calc2(p.vec = p.vec, Se = Se, Sp = Sp, 
                                 ordering = ordering, group.mem = hier.config, 
                                 a = a, trace = trace, print.time = print.time, ...)
      } else if (!is.null(probabilities)) {
        results <- NI.Dorf.calc2(p.vec = probabilities[,1], Se = Se, Sp = Sp, 
                                 ordering = ordering, group.mem = hier.config,  
                                 a = a, trace = trace, print.time = print.time, ...)
      }
    }
    
    # call function for non-informative three-stage hierarchical testing
    if (algorithm == "D3") {
      if (!is.null(p.vec)) {
        results <- NI.D3.calc2(p.vec = p.vec, Se = Se, Sp = Sp, 
                               ordering = ordering, group.mem = hier.config, 
                               a = a, trace = trace, print.time = print.time, ...)
      } else if (!is.null(probabilities)) {
        results <- NI.D3.calc2(p.vec = probabilities[,1], Se = Se, Sp = Sp, 
                               ordering = ordering, group.mem = hier.config, 
                               a = a, trace = trace, print.time = print.time, ...)
      }
    }
    
    # call function for non-informative three-stage hierarchical testing
    if (algorithm == "D4") {
      if (!is.null(p.vec)) {
        results <- NI.D4.calc2(p.vec = p.vec, Se = Se, Sp = Sp, 
                               ordering = ordering, group.mem = hier.config, 
                               a = a, trace = trace, print.time = print.time, ...)
      } else if (!is.null(probabilities)) {
        results <- NI.D4.calc2(p.vec = probabilities[,1], Se = Se, Sp = Sp, 
                               ordering = ordering, group.mem = hier.config, 
                               a = a, trace = trace, print.time = print.time, ...)
      }
    }
    
    # call function for non-informative three-stage hierarchical testing
    if (algorithm == "D5") {
      if (!is.null(p.vec)) {
        results <- NI.D5.calc2(p.vec = p.vec, Se = Se, Sp = Sp, 
                               ordering = ordering, group.mem = hier.config, 
                               a = a, trace = trace, print.time = print.time, ...)
      } else if (!is.null(probabilities)) {
        results <- NI.D5.calc2(p.vec = probabilities[,1], Se = Se, Sp = Sp, 
                               ordering = ordering, group.mem = hier.config, 
                               a = a, trace = trace, print.time = print.time, ...)
      }
    }
    
    # call function for non-informative array testing without master pooling
    if (algorithm == "A2") {
      if (!is.null(p.vec)) {
        results <- NI.Array.calc2(p.vec = p.vec, Se = Se, Sp = Sp, 
                                  group.sz = rowcol.sz, trace = trace, print.time = print.time)
      } else if (!is.null(probabilities)) {
        results <- NI.Array.calc2(p.vec = probabilities[,1], Se = Se, Sp = Sp, 
                                  group.sz = rowcol.sz, trace = trace, print.time = print.time)
      }
    }
    
    # call function for non-informative array testing with master pooling
    if (algorithm == "A2M") {
      if (!is.null(p.vec)) {
        results <- NI.A2M.calc2(p.vec = p.vec, Se = Se, Sp = Sp, 
                                group.sz = rowcol.sz, trace = trace, print.time = print.time)
      } else if (!is.null(probabilities)) {
        results <- NI.A2M.calc2(p.vec = probabilities[,1], Se = Se, Sp = Sp, 
                                group.sz = rowcol.sz, trace = trace, print.time = print.time)
      }
    }
    
    # call function for informative two-stage hierarchical (PSOD) testing
    if (algorithm == "ID2") {
      if (!is.null(alpha)) {
        results <- Inf.Dorf.calc2(alpha = alpha, Se = Se, Sp = Sp, 
                                  ordering = ordering, group.mem = hier.config, 
                                  a = a, trace = trace, print.time = print.time)
      } else if (!is.null(probabilities)) {
        results <- Inf.Dorf.calc2(probabilities = probabilities, 
                                  Se = Se, Sp = Sp, ordering = ordering, 
                                  group.mem = hier.config, a = a, 
                                  trace = trace, print.time = print.time)
      }
    }
    
    # call function for informative two-stage hierarchical (PSOD) testing
    if (algorithm == "ID3") {
      if (!is.null(alpha)) {
        results <- Inf.D3.calc2(alpha = alpha, Se = Se, Sp = Sp, 
                                ordering = ordering, group.mem = hier.config, 
                                a = a, trace = trace, print.time = print.time)
      } else if (!is.null(probabilities)) {
        results <- Inf.D3.calc2(probabilities = probabilities, 
                                Se = Se, Sp = Sp, ordering = ordering, 
                                group.mem = hier.config, a = a, 
                                trace = trace, print.time = print.time)
      }
    }
    
    # call function for informative two-stage hierarchical (PSOD) testing
    if (algorithm == "ID4") {
      if (!is.null(alpha)) {
        results <- Inf.D4.calc2(alpha = alpha, Se = Se, Sp = Sp, 
                                ordering = ordering, group.mem = hier.config, 
                                a = a, trace = trace, print.time = print.time)
      } else if (!is.null(probabilities)) {
        results <- Inf.D4.calc2(probabilities = probabilities, 
                                Se = Se, Sp = Sp, ordering = ordering, 
                                group.mem = hier.config, a = a, 
                                trace = trace, print.time = print.time)
      }
    }
    
    # call function for informative two-stage hierarchical (PSOD) testing
    if (algorithm == "ID5") {
      if (!is.null(alpha)) {
        results <- Inf.D5.calc2(alpha = alpha, Se = Se, Sp = Sp, 
                                ordering = ordering, group.mem = hier.config, 
                                a = a, trace = trace, print.time = print.time)
      } else if (!is.null(probabilities)) {
        results <- Inf.D5.calc2(probabilities = probabilities, 
                                Se = Se, Sp = Sp, ordering = ordering, 
                                group.mem = hier.config, a = a, 
                                trace = trace, print.time = print.time)
      }
    }
    
    # call function for informative array testing not available
    
    class(results) <- "opChar"
    results
  }



###############################################################################
#' @rdname operatingCharacteristics2
opChar2 <- operatingCharacteristics2


#
