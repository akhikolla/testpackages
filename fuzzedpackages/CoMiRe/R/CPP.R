#' Collaborative Perinatal Project data
#'
#' @description
#' A subset of the Collaborative Perinatal Project data set (Klebanoff, 2009)
#' focusing on studying the effect of DDE exposure on pregnancies (Longnecker et al., 2001).
#' The dataset contains the following variables for each pregnant women enrolled in the study:
#'\itemize{
#'   \item hosp, factor denoting the hospital where the woman was hospitalized;
#'   \item dde, Dichlorodiphenyldichloroethylene (DDE) concentration in maternal serum;
#'   \item gestage, gestational age (in weeks);
#'   \item mage, age of the moter (in years);
#'   \item mweight, pre pregnancy weigth of the mother (in lbs);
#'   \item mbmi, pre pregnancy BMI of the mother;
#'   \item sei, socio economic index of the mother;
#'   \item smoke, factor. It takes value 2 if the woman is a smoker, 1 otherwise;
#' }
#'
#' @docType data
#' @keywords dataset internal
#' @name CPP
#'
#' @usage data(CPP)
#' @format A data.frame
#'
#' @examples
#' data(CPP)
#' str(CPP)
#'
#' @references
#'
#' Klebanoff M. A. (2009) The collaborative perinatal project: a 50-year retrospective.
#' Paediatric and perinatal epidemiology, 23, 2.
#'
#' Longnecker, M. P., Klebanof, M. A., Zhou, H., Brock, J. (2001)
#' Association between maternal serum concentration of the DDT metabolite
#' DDE and preterm and small-for-gestational-age babies at birth. The Lancet, 358, 110-114.
#'
"CPP"