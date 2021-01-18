#' Percentage of Body Fat
#'
#' A data frame that contains the proportion of cord fat for individuals calculated through various body measurements of weight, height and circumferences of 252 men who participated in the study by Dr. A. Garth Fisher, Human Performance Research Center, Brigham Young University.
#'@usage data(bodyfat)
#'@details It is possible to find some errors in the table or strange data:
#'
#'One man (case 42) was measured with over 200 pounds in weight who is less than 3 feet tall, considered that he had a typo when typing 29.5 inches and transformed the data into 69.5 inches;
#'
#'There was a man with a negative percentage of body fat, it was decided to exclude this data from the table.
#'
#'Changes to units of measure:
#'
#'Weight was transformed from lbs to kg / 100 (value 1 corresponds to 100kg);
#'
#'Height has been transformed from inches to meters;
#'
#'All columns that were represented in centimeters were transformed into meters.
#' @format This data frame contains the observations of 252 men:
#' \describe{
#'   \item{case}{Case number.}
#'   \item{brozek}{Percent body fat using Brozek's equation: 457/Density - 414.2}
#'   \item{siri}{Percent body fat using Siri's equation: 495/Density - 450}
#'   \item{density}{Density determined from underwater weighing (gm/cm**3).}
#'   \item{age}{Age (years).}
#'   \item{weight}{Weight (kg/100).}
#'   \item{height}{Height (m).}
#'   \item{neck}{Neck circumference (m).}
#'   \item{chest}{Chest circumference (m).}
#'   \item{abdomen}{Abdomen circumference (m) ''at the umbilicus and level
#'     with the iliac crest''.}
#'   \item{forearm}{Forearm circumference (m).}
#'   \item{hip}{Hip circumference (m).}
#'   \item{thigh}{Thigh circumference (m).}
#'   \item{knee}{Knee circumference (m).}
#'   \item{ankle}{Ankle circumference (m).}
#'   \item{biceps}{Biceps (extended) circumference (m).}
#'   \item{wrist}{Wrist circumference (m) ''distal to the styloid processes''.}
#' }
#'@references
#'  \doi{10.1080/10691898.1996.11910505} Johnson, R. W. (1996). Fitting percentage of body fat to simple body measurements. \emph{Journal of Statistics Education}, \bold{4}(1).
#'@references
#'\doi{10.1249/00005768-198504000-00037} Penrose, K. W., Nelson, A. G., & Fisher, A. G. (1985). Generalized body composition prediction equation for men using simple measurement techniques. \emph{Medicine & Science in Sports & Exercise}, \bold{17}(2), 189.
#'@references
#'\doi{10.1016/j.csda.2006.05.006} Royston, P., & Sauerbrei, W. (2007). Improving the robustness of fractional polynomial models by preliminary covariate transformation: A pragmatic approach. \emph{Computational statistics & data analysis}, \bold{51}(9), 4240-4253.
#' @examples
#' data(bodyfat,package="bayesbr")
#' \dontshow{
#' lines = sample(1:251,15)
#' bodyfat = bodyfat[lines,]
#' }
#'bbr = bayesbr(siri ~ age+wrist*neck+chest+
#'              thigh+wrist| wrist, data = bodyfat,
#'              iter = 100)
#'
#'summary(bbr)
#'
#'\donttest{
#'bbr = bayesbr(siri ~ I(age/100)+heigth+chest+
#'              thigh+wrist| wrist,
#'              data = bodyfat,iter = 1000)
#'
#'}
"bodyfat"

