#' NINDS rt-PA stroke trial data
#'
#' @description The \code{ninds} data frame has 1906 rows and 16 columns.
#'
#' @format This data frame contains the following columns:
#' 
#'   \describe{
#'   
#'   \item{\code{ID}}{patient identifier.}
#'
#'   \item{\code{Y}}{Acute ischemic stroke in an ordinal scale. Coded as \code{1} = no symptoms or no significant disability despite symptoms;
#'   \code{2} = slight disability; \code{3} = moderate disability or moderately severe disability; \code{4} = severe disability or dead.}
#'   
#'   \item{\code{intcpt_RE}}{a column of 1's for model setup as a random intercept covariate.}
#'   
#'   \item{\code{smlves_NP}}{small vessel occlusive disease as a non-proportional odds covariate.}
#'   
#'   \item{\code{lvORcs_NP}}{large vessel atherosclerosis / cardioembolic stroke as a non-proportional odds covariate. }
#'   
#'   \item{\code{group}}{treatment group indicator.}
#'   
#'   \item{\code{time3}}{dummy variable to indicate 3 months of follow-up.}
#'   
#'   \item{\code{time6}}{dummy variable to indicate 6 months of follow-up.}
#'   
#'   \item{\code{time12}}{dummy variable to indicate 12 months of follow-up.}
#'   
#'   \item{\code{mrkprior}}{modified Rankin scale prior stroke onset.}
#'   
#'   \item{\code{smlves}}{small vessel occlusive disease as a proportional odds covariate.}
#'   
#'   \item{\code{lvORcs}}{large vessel atherosclerosis / cardioembolic stroke as a proportional odds covariate.}
#'   
#'   \item{\code{smlves.group}}{interaction term between \code{smlves} and \code{group}.}
#'   
#'   \item{\code{lvORcs.group}}{interaction term between \code{lvORcs} and \code{group}.}
#'   
#'   \item{\code{surv}}{time to drop out or remaining in severe disability.}
#'   
#'   \item{\code{comprisk}}{event indicator. Coded as \code{0} = censored; \code{1} = drop out; 
#'   \code{2} = remaining in severe disability.}
#'   
#'   }
#' @usage data(ninds)
#' 
#' @references 
#' Li, Ning, et al. "Joint modeling of longitudinal ordinal data and competing risks survival 
#' times and analysis of the NINDS rt-PA stroke trial." Statistics in medicine 29.5 (2010): 
#' 546-557.
#' 
"ninds"

#' Scleroderma lung study data
#' 
#' @description The \code{lung} data frame has 715 rows and 13 columns.
#'
#' @format A balanced data set with respect to the times at which observations 
#'   recorded. The data consists of the following variables on each patient:
#'   
#'   \describe{
#'   
#'   \item{\code{ID}}{patient identifier.}
#
#'   \item{\code{FVC}}{forced vital capacity (\%) determined at 3-month intervals from the baseline.}
#'   
#'   \item{\code{time_RE}}{time at visit at 3-month intervals.}
#'   
#'   \item{\code{intercept}}{column of 1's for model setup.}
#'   
#'   \item{\code{time}}{time at visit at 3-month intervals. Same as \code{time_RE}.}
#'   
#'   \item{\code{FVC0}}{forced vital capacity (\%) at baseline.}
#'   
#'   \item{\code{FIB0}}{baseline lung fibrosis.}
#'   
#'   \item{\code{CYC}}{treatment allocation. Coded as \code{1} = oral cyclophosphamide; \code{0} = placebo.}
#'   
#'   \item{\code{FVC0.CYC}}{interaction term between \code{FVC0} and \code{CYC}.}
#'   
#'   \item{\code{FIB0.CYC}}{interaction term between \code{FIB0} and \code{CYC}.}
#'   
#'   \item{\code{time.CYC}}{interaction term between \code{time} and \code{CYC}.}
#'   
#'   \item{\code{time.CYC}}{interaction term between \code{time} and \code{CYC}.}
#'   
#'   \item{\code{surv}}{time to treatment failure or death.}
#'   
#'   \item{\code{failure_type}}{treatment failure/death indicator. Coded as \code{0} = censored; \code{1} = death; 
#'   \code{2} = treatment failure.}
#'   
#'   }
#' @usage data(lung)
#' 
#' @references
#' 
#' Elashoff, Robert M., Gang Li, and Ning Li. "A joint model for longitudinal measurements and 
#' survival data in the presence of multiple failure types." Biometrics 64.3 (2008): 762-771.
"lung"