#' Buhaugetal_2009_JCR
#'
#'Subsetted version of survival database extracted from \href{http://bit.ly/2Q1Igo9}{Buhaug et al. (2009)}.
#'It has precisely dated duration data of internal conflict as well as geographic data.
#'Variables Y, Y0 and C were later added by \href{http://bit.ly/38eDsnG}{Bagozzi et al. (2019)}.
#'It is used to estimate the Bayesian Misclassified Failure (MF) Weibull model
#'presented in \href{http://bit.ly/38eDsnG}{Bagozzi et al. (2019)}.
#'
#'
#'\describe{
#'    \item{lndistx}{log conflict-capital distance.}
#'    \item{confbord}{conflict zone at border.}
#'    \item{borddist}{confbord * lndistx centred.}
#'    \item{figcapdum}{rebel fighting capacity at least moderate.}
#'    \item{lgdp_onset}{gdp capita in onset year.}
#'    \item{sip2l_onset}{Gates et al. (2006) SIP code (1 year lag) for the onset year.}
#'    \item{pcw}{post cold war period, 1989+.}
#'    \item{frst}{percentage of forest in conflict zone.}
#'    \item{mt}{percentage of mountains in conflict zone.}
#'    \item{Y}{conflict duration.}
#'    \item{Y0}{elapsed time since inception to Y (t-1).}
#'    \item{C}{censoring variable.}
#'    \item{coupx}{coup d'etat, except if overlapping with other gov't conflict (PHI 1989).}
#' }
#' @docType data
#' @keywords datasets
#' @name Buhaugetal_2009_JCR
#' @usage data(Buhaugetal_2009_JCR)
#' @format A data frame with 1562 rows and 13 variables
#' @source Buhaug, Halvard, Scott Gates, and Päivi Lujala (2009), Geography, rebel capability, and the duration of civil conflict,  Journal of Conflict Resolution 53(4), 544 - 569.
"Buhaugetal_2009_JCR"



