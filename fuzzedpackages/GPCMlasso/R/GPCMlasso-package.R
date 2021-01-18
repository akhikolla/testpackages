

#' Find DIF in Generalized Partial Credit Models
#' 
#' Performs GPCMlasso, a method to identify DIF in Generalized Partial Credit Models. 
#' A joint parametric model is set up based on an IRT model chosen by the user. Several variables can be considered simultaneously.
#' For each pair between variable and item, a parametric DIF effect is introduced which
#' indicates DIF if the respective parameter is selected (estimated to be unequal zero). 
#' Parameter selection is done using a lasso-type penalization term.
#' 
#' 
#' @name GPCMlasso-package
#' @docType package
#' @author Gunther Schauberger\cr \email{gunther.schauberger@@tum.de}
#' @references Schauberger, Gunther and Mair, Patrick (2019): A Regularization Approach for the Detection of Differential 
#' Item Functioning in Generalized Partial Credit Models, \emph{Behavior Research Methods}, \url{https://link.springer.com/article/10.3758/s13428-019-01224-2}
#' @seealso \code{\link{GPCMlasso}}
#' @keywords package Partial Credit DIF GPCMlasso GPCM DSF
#' @examples
#' data(tenseness_small)
#' 
#' ## formula for simple model without covariates
#' form.0 <- as.formula(paste("cbind(",paste(colnames(tenseness_small)[1:5],collapse=","),")~0"))
#' 
#' ######
#' ## fit simple RSM where loglikelihood and score function are evaluated parallel on 2 cores
#' rsm.0 <- GPCMlasso(form.0, tenseness_small, model = "RSM", 
#' control= ctrl_GPCMlasso(cores=2))
#' rsm.0
#' 
#' \dontrun{
#' ## formula for model with covariates (and DIF detection)
#' form <- as.formula(paste("cbind(",paste(colnames(tenseness_small)[1:5],collapse=","),")~."))
#' 
#' ######
#' ## fit GPCM model with 10 different tuning parameters
#' gpcm <- GPCMlasso(form, tenseness_small, model = "GPCM", 
#'                   control = ctrl_GPCMlasso(l.lambda = 10))
#' gpcm
#' plot(gpcm)
#' pred.gpcm <- predict(gpcm)
#' trait.gpcm <- trait.posterior(gpcm)
#' 
#' ######
#' ## fit RSM, detect differential step functioning (DSF)
#' rsm.DSF <- GPCMlasso(form, tenseness_small, model = "RSM", DSF = TRUE, 
#'                      control = ctrl_GPCMlasso(l.lambda = 10))
#' rsm.DSF
#' plot(rsm.DSF)
#' 
#' ## create binary data set
#' tenseness_small_binary <- tenseness_small
#' tenseness_small_binary[,1:5][tenseness_small[,1:5]>1] <- 2
#' 
#' ######
#' ## fit and cross-validate Rasch model
#' set.seed(1860)
#' rm.cv <- GPCMlasso(form, tenseness_small_binary, model = "RM", cv = TRUE, 
#'                    control = ctrl_GPCMlasso(l.lambda = 10))
#' rm.cv
#' plot(rm.cv)
#' }
NULL


#' Tenseness data from the Freiburg Complaint Checklist
#' 
#' Data from the Freiburg Complaint Checklist. 
#' The data contain all 8 items corresponding to the scale \emph{Tenseness} for 2042 participants of the 
#' standardization sample of the Freiburg Complaint Checklist. 
#' 
#' @name tenseness
#' @docType data
#' @format A data frame containing data from the Freiburg Complaint Checklist with 1847 observations. 
#' All items refer to the scale \emph{Tenseness} and are measured on a 5-point Likert scale where low numbers 
#' correspond to low frequencies or low intensitites of the respective complaint and vice versa. 
#' \describe{ 
#' \item{Clammy_hands}{Do you have clammy hands?}
#' \item{Sweat_attacks}{Do you have sudden attacks of sweating?}
#' \item{Clumsiness}{Do you notice that you behave clumsy?}
#' \item{Wavering_hands}{Are your hands wavering frequently, e.g. when lightning a cigarette or when holding a cup?}
#' \item{Restless_hands}{Do you notice that your hands are restless?}
#' \item{Restless_feet}{Do you notice that your feet are restless?}
#' \item{Twitching_eyes}{Do you notice unvoluntary twitching of your eyes?}
#' \item{Twitching_mouth}{Do you notice unvoluntary twitching of your mouth?}
#' \item{Gender}{Gender of the person}
#' \item{Household}{Does the person live alone in a household or together with somebody?}
#' \item{Income}{Income, categorized to levels from 1 (low income) to 11(high income). For simplicity,
#' due to the high number of categories income can be treated as a metric variable.}
#' \item{WestEast}{Is the person from East Germany (former GDR)?}
#' \item{Abitur}{Does the person have Abitur (A-levels)?}
#' \item{Age}{Age of the person}
#'  }
#' @source 
#' ZPID (2013). PsychData of the Leibniz Institute for Psychology Information ZPID. Trier: Center for Research Data in Psychology.
#' 
#' Fahrenberg, J. (2010). Freiburg Complaint Checklist [Freiburger Beschwerdenliste (FBL)]. Goettingen, Hogrefe.
#' @keywords datasets
#' @examples
#' 
#' data(tenseness)
#' 
NULL

#' Subset of tenseness data from the Freiburg Complaint Checklist
#' 
#' Data from the Freiburg Complaint Checklist. 
#' The data contain 5 items (out of 8) corresponding to the scale \emph{Tenseness} for a subset of 200 participants of the 
#' standardization sample of the Freiburg Complaint Checklist. 
#' 
#' @name tenseness_small
#' @docType data
#' @format A data frame containing data from the Freiburg Complaint Checklist a subset of 200 observations. 
#' The complete data set with 1847 observations can be found in \code{\link{tenseness}}.
#' All items refer to the scale \emph{Tenseness} and are measured on a 5-point Likert scale where low numbers 
#' correspond to low frequencies or low intensitites of the respective complaint and vice versa. 
#' \describe{ 
#' \item{Clammy_hands}{Do you have clammy hands?}
#' \item{Sweat_attacks}{Do you have sudden attacks of sweating?}
#' \item{Clumsiness}{Do you notice that you behave clumsy?}
#' \item{Wavering_hands}{Are your hands wavering frequently, e.g. when lightning a cigarette or when holding a cup?}
#' \item{Restless_hands}{Do you notice that your hands are restless?}
#' \item{Gender}{Gender of the person}
#' \item{Age}{Age of the person}
#'  }
#' @seealso \code{\link{GPCMlasso}}, \code{\link{ctrl_GPCMlasso}}, \code{\link{trait.posterior}}
#' @source 
#' ZPID (2013). PsychData of the Leibniz Institute for Psychology Information ZPID. Trier: Center for Research Data in Psychology.
#' 
#' Fahrenberg, J. (2010). Freiburg Complaint Checklist [Freiburger Beschwerdenliste (FBL)]. Goettingen, Hogrefe.
#' @keywords GPCMlasso
#' @examples
#' data(tenseness_small)
#' 
#' ## formula for simple model without covariates
#' form.0 <- as.formula(paste("cbind(",paste(colnames(tenseness_small)[1:5],collapse=","),")~0"))
#' 
#' ######
#' ## fit simple RSM where loglikelihood and score function are evaluated parallel on 2 cores
#' rsm.0 <- GPCMlasso(form.0, tenseness_small, model = "RSM", 
#' control= ctrl_GPCMlasso(cores=2))
#' rsm.0
#' 
#' \dontrun{
#' ## formula for model with covariates (and DIF detection)
#' form <- as.formula(paste("cbind(",paste(colnames(tenseness_small)[1:5],collapse=","),")~."))
#' 
#' ######
#' ## fit GPCM model with 10 different tuning parameters
#' gpcm <- GPCMlasso(form, tenseness_small, model = "GPCM", 
#'                   control = ctrl_GPCMlasso(l.lambda = 10))
#' gpcm
#' plot(gpcm)
#' pred.gpcm <- predict(gpcm)
#' trait.gpcm <- trait.posterior(gpcm)
#' 
#' ######
#' ## fit RSM, detect differential step functioning (DSF)
#' rsm.DSF <- GPCMlasso(form, tenseness_small, model = "RSM", DSF = TRUE, 
#'                      control = ctrl_GPCMlasso(l.lambda = 10))
#' rsm.DSF
#' plot(rsm.DSF)
#' 
#' ## create binary data set
#' tenseness_small_binary <- tenseness_small
#' tenseness_small_binary[,1:5][tenseness_small[,1:5]>1] <- 2
#' 
#' ######
#' ## fit and cross-validate Rasch model
#' set.seed(1860)
#' rm.cv <- GPCMlasso(form, tenseness_small_binary, model = "RM", cv = TRUE, 
#'                    control = ctrl_GPCMlasso(l.lambda = 10))
#' rm.cv
#' plot(rm.cv)
#' }
NULL
