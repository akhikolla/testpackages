
#' Blood pressure data from a clinical study
#' 
#' Data from 200 subjects
#' 
#' @format A data frame with 2438 rows and 13 variables:
#' \describe{
#'   \item{ID}{Subject identification number}
#'   \item{BIRTH_WT}{birth weight (lbs)}
#'   \item{WEIGHT}{current weight (lbs)}
#'   \item{HEIGHT}{current height (cm)}
#'   \item{BMI}{current body mass index}
#'   \item{age}{current age (yrs)}
#'   \item{dias}{diastolic blood pressure}
#'   \item{sys}{systolic blood pressure}
#'   \item{SexM}{indicator of sex male}
#'   \item{RaceB}{indicator of race black}
#'   \item{RaceW}{indicator of race white}
#'   \item{PHIGHBP}{indicator that either parent had high blood pressure}
#'   \item{PDIABET}{indicator that either parent had diabetes}
#' }
#' @source Data provided by Wanzhu Tu, Indiana University School of Medicine
#' @references Tu, W., Eckert, G. J., DiMeglio, L. A., Yu, Z., Jung, J., and Pratt, J. H. (2011). \emph{Intensified effect of adiposity on blood pressure in overweight and obese children}. Hypertension, 58(5), 818-824.
"bloodpressure"

#' Coral reef data from survey data on 6 sites
#' 
#' Data from 68 subjects
#' 
#' @format A data frame with 269 rows and 14 variables:
#' \describe{
#'   \item{ZONE}{Management zone}
#'   \item{site}{Name of the habitat site}
#'   \item{complexity}{habitat benthic complexity}
#'   \item{rugosity}{a measurement related to terrain complexity}
#'   \item{LC}{cover of low complexity}
#'   \item{HC}{cover of high complexity}
#'   \item{SCORE1}{PCA score 1 from Wilson, Graham, Polunin}
#'   \item{SCORE2}{PCA score 2 from Wilson, Graham, Polunin}
#'   \item{macro}{indicator of race white}
#'   \item{species}{fish species}
#'   \item{abundance}{fish abundance}
#'   \item{biomass}{fish biomass}
#' }
#' @source Data from supplementary material provided for Fisher, R., Wilson, S. K., Sin, T. M., Lee, A. C., and Langlois, T. J. (2018). \emph{A simple function for full-subsets multiple regression in ecology with R}. Ecology and evolution, 8(12), 6104-6113.
#' @references Wilson, S. K., Graham, N. A. J., and Polunin, N. V. (2007). \emph{Appraisal of visual assessments of habitat complexity and benthic composition on coral reefs}. Marine Biology, 151(3), 1069-1076.
"reef"
