#LocalControl package data descriptions

#' Lindner Center for Research and Education study on Abciximab cost-effectiveness and survival
#' @description
#' The effects of Abciximab use on both survival and cardiac billing.
#' @format A data frame with 996 rows and 10 columns:
#' \describe{
#'  \item{lifepres}{Life years preserved post treatment: 0 (died) vs. 11.6 (survived).}
#'  \item{cardbill}{Cardiac related billing in dollars within 12 months.}
#'  \item{abcix}{Indicates whether the patient received Abciximab treatment: 1=yes 0=no.}
#'  \item{stent}{Was a stent depolyed? 1=yes, 0=no.}
#'  \item{height}{Patient height in centimeters.}
#'  \item{female}{Patient sex: 1=female, 0=male.}
#'  \item{diabetic}{Was the patient diabetic? 1=yes, 0=no.}
#'  \item{acutemi}{Had the patient suffered an acute myocardial infarction witih the last seven days? 1=yes, 0=no.}
#'  \item{ejecfrac}{Left ventricular ejection fraction.}
#'  \item{ves1proc}{Number of vessels involved in the first PCI procedure.}
#'}
#' @name lindner
#' @docType data
#' @references Kereiakes DJ, Obenchain RL, Barber BL, Smith A, McDonald M, Broderick TM, Runyon JP, Shimshak TM, Schneider JF, Hattemer CR, Roth EM, Whang DD, Cocks D, Abbottsmith CW. Abciximab provides cost-effective survival advantage in high-volume interventional practice. Am Heart J. 2000;140(4):603-610.
#' @keywords data
NULL


#' Framingham heart study data extract on smoking and hypertension.
#' @description Data collected over a 24 year study suitable for competing
#'   risks survival analysis of hypertension and death as a function of smoking.
#' @format A data frame with 2316 rows and 11 columns:
#' \describe{
#'  \item{female}{Sex of the patient. 1=female, 0=male.}
#'  \item{totchol}{Total cholesterol of patient at study entry.}
#'  \item{age}{Age of the patient at study entry.}
#'  \item{bmi}{Patient body mass index.}
#'  \item{BPVar}{Average units of systolic and diastolic blood pressure above normal: ((SystolicBP-120)/2) + (DiasystolicBP-80)}
#'  \item{heartrte}{Patient heartrate taken at study entry.}
#'  \item{glucose}{Patient blood glucose level.}
#'  \item{cursmoke}{Whether or not the patient was a smoker at the time of study entry.}
#'  \item{outcome}{Did the patient die, experience hypertension, or leave the study without experiencing either event.}
#'  \item{time_outcome}{The time at which the patient experienced outcome.}
#'  \item{cigpday}{Number of cigarettes smoked per day at time of study entry.}
#' }
#' @name framingham
#' @docType data
#' @references
#' \itemize{
#'     \item Dawber TR, Meadors GF, Moore FE Jr. Epidemiological approaches to heart disease: the Framingham Study. Am J Public Health Nations Health. 1951;41(3):279-281.
#'     \item Teaching Datasets - Public Use Datasets. \url{https://biolincc.nhlbi.nih.gov/teaching/}.
#' }
#' @keywords data
NULL


#' @title Simulated cardiac medication data for survival analysis
#' @description
#' This dataset was created to demonstrate the effects of Local Control on correcting bias within a set of data.
#' @format A data frame with 1000 rows and 6 columns:
#' \describe{
#'  \item{id}{Unique identifier for each row.}
#'  \item{time}{Time in years to the outcome specified by status.}
#'  \item{status}{1 if the patient experienced cardiac arrest. 0 if censored before that.}
#'  \item{drug}{Medication the patient received for cardiac health (drug 1 or drug 0).}
#'  \item{age}{Age of the patient, ranges from 18 to 65 years.}
#'  \item{bmi}{Patient body mass index. Majority of observations fall between 22 and 30.}
#' }
#' @name cardSim
#' @docType data
#' @author  Lauve NR, Lambert CG
#' @keywords data
NULL
