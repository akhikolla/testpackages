#' JOUSBoost: A package for probability estimation
#'
#' JOUSBoost implements under/oversampling with jittering for probability estimation.
#'  Its intent is to be used to improve probability estimates that come from
#' boosting algorithms (such as AdaBoost), but is modular enough to be used with
#' virtually any classification algorithm from machine learning.
#'
#' For more theoretical background, consult Mease (2007).
#'
#' @references Mease, D., Wyner, A. and Buja, A. (2007). Costweighted
#' boosting with jittering and over/under-sampling:
#' JOUS-boost. J. Machine Learning Research 8 409-439.
#'
#' @docType package
#' @name JOUSBoost
NULL

#' Dataset of sonar measurements of rocks and mines
#'
#' A dataset containing sonar measurements used to discriminate rocks from mines.
#'
#' @docType data
#' @keywords datasets
#' @name sonar
#' @usage data(sonar)
#' @format A data frame with 208 observations on 61 variables.  The variables
#'         V1-V60 represent the energy within a certain frequency band, and
#'         are to be used as predictors.  The variable y is a class label, 1
#'         for 'rock' and -1 for 'mine'.
#' @source \url{http://archive.ics.uci.edu/ml/}
#' @references Gorman, R. P., and Sejnowski, T. J. (1988). "Analysis of Hidden
#'  Units in a Layered Network
#'  Trained to Classify Sonar Targets" in Neural Networks, Vol. 1, pp. 75-89.
"sonar"





