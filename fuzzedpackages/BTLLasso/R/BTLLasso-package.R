

#' BTLLasso
#' 
#' Performs BTLLasso, a method to model heterogeneity in paired comparison
#' data. Different types of covariates are allowd to have an influence on the
#' attractivity/strength of the objects. Covariates can be subject-specific, 
#' object-specific or subject-object-specific. L1 penalties are used to reduce the 
#' complexity of the model by enforcing clusters of equal effects or by elimination of irrelevant
#' covariates.  Several additional functions are provided, such as
#' cross-validation, bootstrap intervals, and plot functions.
#' 
#' 
#' @name BTLLasso-package
#' @docType package
#' @author Gunther Schauberger\cr \email{gunther.schauberger@@tum.de}
#' @seealso \code{\link{BTLLasso}}, \code{\link{cv.BTLLasso}}
#' @references Schauberger, Gunther and Tutz, Gerhard (2019): BTLLasso - A Common Framework and Software 
#' Package for the Inclusion  and Selection of Covariates in Bradley-Terry Models, \emph{Journal of 
#' Statistical Software}, 88(9), 1-29, \url{https://doi.org/10.18637/jss.v088.i09}
#' 
#' Schauberger, Gunther and Tutz, Gerhard (2017): Subject-specific modelling 
#' of paired comparison data: A lasso-type penalty approach, \emph{Statistical Modelling},
#' 17(3), 223 - 243
#' 
#' Schauberger, Gunther, Groll Andreas and Tutz, Gerhard (2018): 
#' Analysis of the importance of on-field covariates in the German Bundesliga, 
#' \emph{Journal of Applied Statistics}, 45(9), 1561 - 1578
#' @keywords package BTL Bradley-Terry BTLLasso
#' @examples
#' 
#' \dontrun{
#' op <- par(no.readonly = TRUE)
#' 
#' ##############################
#' ##### Example with simulated data set containing X, Z1 and Z2
#' ##############################
#' data(SimData)
#' 
#' ## Specify control argument
#' ## -> allow for object-specific order effects and penalize intercepts
#' ctrl <- ctrl.BTLLasso(penalize.intercepts = TRUE, object.order.effect = TRUE,
#'                       penalize.order.effect.diffs = TRUE)
#' 
#' ## Simple BTLLasso model for tuning parameters lambda
#' m.sim <- BTLLasso(Y = SimData$Y, X = SimData$X, Z1 = SimData$Z1,
#'                   Z2 = SimData$Z2, control = ctrl)
#' m.sim
#' 
#' par(xpd = TRUE)
#' plot(m.sim)
#' 
#' 
#' ## Cross-validate BTLLasso model for tuning parameters lambda
#' set.seed(1860)
#' m.sim.cv <- cv.BTLLasso(Y = SimData$Y, X = SimData$X, Z1 = SimData$Z1,
#'                         Z2 = SimData$Z2, control = ctrl)
#' m.sim.cv
#' coef(m.sim.cv)
#' logLik(m.sim.cv)
#' 
#' head(predict(m.sim.cv, type="response"))
#' head(predict(m.sim.cv, type="trait"))
#' 
#' plot(m.sim.cv, plots_per_page = 4)
#' 
#' 
#' ## Example for bootstrap intervals for illustration only
#' ## Don't calculate bootstrap intervals with B = 20!!!!
#' set.seed(1860)
#' m.sim.boot <- boot.BTLLasso(m.sim.cv, B = 20, cores = 20)
#' m.sim.boot
#' plot(m.sim.boot, plots_per_page = 4)
#' 
#' 
#' ##############################
#' ##### Example with small version from GLES data set
#' ##############################
#' data(GLESsmall)
#' 
#' ## extract data and center covariates for better interpretability
#' Y <- GLESsmall$Y
#' X <- scale(GLESsmall$X, scale = FALSE)
#' Z1 <- scale(GLESsmall$Z1, scale = FALSE)
#' 
#' ## vector of subtitles, containing the coding of the X covariates
#' subs.X <- c('', 'female (1); male (0)')
#' 
#' ## Cross-validate BTLLasso model
#' m.gles.cv <- cv.BTLLasso(Y = Y, X = X, Z1 = Z1)
#' m.gles.cv
#' 
#' coef(m.gles.cv)
#' logLik(m.gles.cv)
#' 
#' head(predict(m.gles.cv, type="response"))
#' head(predict(m.gles.cv, type="trait"))
#' 
#' par(xpd = TRUE, mar = c(5,4,4,6))
#' plot(m.gles.cv, subs.X = subs.X, plots_per_page = 4, which = 2:5)
#' paths(m.gles.cv, y.axis = 'L2')
#' 
#' 
#' ##############################
#' ##### Example with Bundesliga data set
#' ##############################
#' data(Buli1516)
#' 
#' Y <- Buli1516$Y5
#' 
#' Z1 <- scale(Buli1516$Z1, scale = FALSE)
#' 
#' ctrl.buli <- ctrl.BTLLasso(object.order.effect = TRUE, 
#'                            name.order = "Home", 
#'                            penalize.order.effect.diffs = TRUE, 
#'                            penalize.order.effect.absolute = FALSE,
#'                            order.center = TRUE, lambda2 = 1e-2)
#' 
#' set.seed(1860)
#' m.buli <- cv.BTLLasso(Y = Y, Z1 = Z1, control = ctrl.buli)
#' m.buli
#' 
#' par(xpd = TRUE, mar = c(5,4,4,6))
#' plot(m.buli)
#' 
#' 
#' ##############################
#' ##### Example with Topmodel data set
#' ##############################
#' data("Topmodel2007", package = "psychotree")
#' 
#' Y.models <- response.BTLLasso(Topmodel2007$preference)
#' X.models <- scale(model.matrix(preference~., data = Topmodel2007)[,-1])
#' rownames(X.models) <- paste0("Subject",1:nrow(X.models))
#' colnames(X.models) <- c("Gender","Age","KnowShow","WatchShow","WatchFinal")
#' 
#' set.seed(5)
#' m.models <- cv.BTLLasso(Y = Y.models, X = X.models)
#' plot(m.models, plots_per_page = 6)
#' 
#' par(op)
#' }
NULL


#' Bundesliga Data 2015/16 (Buli1516)
#' 
#' Data from the German Bundesliga from the season 2015/16. 
#' The data contain all 306 matches of the season treated as paired comparisons with  5 (Y5) or 3 (Y3) different 
#' response categories. Additionally, different match-specific covariates are given as, for example, 
#' the percentage of ball possession or the total running distance per team and per match.
#' 
#' @name Buli1516
#' @docType data
#' @format A list containing data from the German Bundesliga with 306 observations. 
#' The list contains both information on the response (paired comparisons) and different covariates.
#' \describe{ 
#' \item{Y5}{A response.BTLLasso object with 5 response categories for the Buli1516 data including
#' \itemize{
#' \item{response: Ordinal paired comparison response vector} 
#' \item{first.object: Vector containing the first-named team per paired comparison (home team)}
#' \item{second.object: Vector containing the second-named team per paired comparison (away team)}
#' \item{subject: Vector containing a match-day identifier per paired comparison}
#' \item{with.order} Vector containing information that each match has to be considered including an order effect.
#' }}
#' \item{Y3}{A response.BTLLasso object with 3 response categories for the Buli1516 data including
#' \itemize{
#' \item{response: Ordinal paired comparison response vector} 
#' \item{first.object: Vector containing the first-named team per paired comparison (home team)}
#' \item{second.object: Vector containing the second-named team per paired comparison (away team)}
#' \item{subject: Vector containing a match-day identifier per paired comparison}
#' \item{with.order} Vector containing information that each match has to be considered including an order effect.
#' }}
#' \item{Z1}{Matrix containing all team-match-specific covariates
#' \itemize{
#' \item{Distance: Total amount of km run} 
#' \item{BallPossession: Percentage of ball possession}
#' \item{TacklingRate: Rate of won tacklings}
#' \item{ShotsonGoal: Total number of shots on goal} 
#' \item{CompletionRate: Percentage of passes reaching teammates} 
#' \item{FoulsSuffered: Number of fouls suffered} 
#' \item{Offside: Number of offsides (in attack)}
#' }
#' }
#' \item{Z2}{Matrix containing all the average market values of the teams as a team-specific covariate} 
#' }
#'  @references Schauberger, Gunther and Tutz, Gerhard (2019): BTLLasso - A Common Framework and Software 
#' Package for the Inclusion  and Selection of Covariates in Bradley-Terry Models, \emph{Journal of 
#' Statistical Software}, to appear
#' 
#' Schauberger, Gunther and Tutz, Gerhard (2017): Subject-specific modelling 
#' of paired comparison data: A lasso-type penalty approach, \emph{Statistical Modelling},
#' 17(3), 223 - 243
#' 
#' Schauberger, Gunther, Groll Andreas and Tutz, Gerhard (2018): 
#' Analysis of the importance of on-field covariates in the German Bundesliga, 
#' \emph{Journal of Applied Statistics}, 45(9), 1561 - 1578
#' @source
#' \url{https://www.kicker.de/}
#' @keywords datasets
#' @seealso \code{\link{Buli1415}}, \code{\link{Buli1617}}, \code{\link{Buli1718}}
#' @examples
#' \dontrun{
#' op <- par(no.readonly = TRUE)
#' 
#' data(Buli1516)
#' 
#' Y <- Buli1516$Y5
#' Z1 <- scale(Buli1516$Z1, scale = FALSE)
#' 
#' ctrl.buli <- ctrl.BTLLasso(object.order.effect = TRUE, 
#'                            name.order = "Home", 
#'                            penalize.order.effect.diffs = TRUE, 
#'                            penalize.order.effect.absolute = FALSE,
#'                            order.center = TRUE, lambda2 = 1e-2)
#' 
#' set.seed(1860)
#' m.buli <- cv.BTLLasso(Y = Y, Z1 = Z1, control = ctrl.buli)
#' m.buli
#' 
#' par(xpd = TRUE, mar = c(5,4,4,6))
#' plot(m.buli)
#' 
#' par(op)
#' }
NULL

#' Bundesliga Data 2014/15 (Buli1415)
#' 
#' Data from the German Bundesliga from the season 2014/15. 
#' The data contain all 306 matches of the season treated as paired comparisons with 5 (Y5) or 3 (Y3) different 
#' response categories. Additionally, different match-specific covariates are given as, for example, 
#' the percentage of ball possession or the total running distance per team and per match.
#' 
#' @name Buli1415
#' @docType data
#' @format A list containing data from the German Bundesliga with 306 observations. 
#' The list contains both information on the response (paired comparisons) and different covariates.
#' \describe{ 
#' \item{Y5}{A response.BTLLasso object with 5 response categories for the Buli1516 data including
#' \itemize{
#' \item{response: Ordinal paired comparison response vector} 
#' \item{first.object: Vector containing the first-named team per paired comparison (home team)}
#' \item{second.object: Vector containing the second-named team per paired comparison (away team)}
#' \item{subject: Vector containing a match-day identifier per paired comparison}
#' \item{with.order} Vector containing information that each match has to be considered including an order effect.
#' }}
#' \item{Y3}{A response.BTLLasso object with 3 response categories for the Buli1516 data including
#' \itemize{
#' \item{response: Ordinal paired comparison response vector} 
#' \item{first.object: Vector containing the first-named team per paired comparison (home team)}
#' \item{second.object: Vector containing the second-named team per paired comparison (away team)}
#' \item{subject: Vector containing a match-day identifier per paired comparison}
#' \item{with.order} Vector containing information that each match has to be considered including an order effect.
#' }}
#' \item{Z1}{Matrix containing all team-match-specific covariates
#' \itemize{
#' \item{Distance: Total amount of km run} 
#' \item{BallPossession: Percentage of ball possession}
#' \item{TacklingRate: Rate of won tacklings}
#' \item{ShotsonGoal: Total number of shots on goal} 
#' \item{CompletionRate: Percentage of passes reaching teammates} 
#' \item{FoulsSuffered: Number of fouls suffered} 
#' \item{Offside: Number of offsides (in attack)}
#' }
#' }
#' }
#' @references Schauberger, Gunther and Tutz, Gerhard (2019): BTLLasso - A Common Framework and Software 
#' Package for the Inclusion  and Selection of Covariates in Bradley-Terry Models, \emph{Journal of 
#' Statistical Software}, to appear
#' 
#' Schauberger, Gunther and Tutz, Gerhard (2017): Subject-specific modelling 
#' of paired comparison data: A lasso-type penalty approach, \emph{Statistical Modelling},
#' 17(3), 223 - 243
#' 
#' Schauberger, Gunther, Groll Andreas and Tutz, Gerhard (2018): 
#' Analysis of the importance of on-field covariates in the German Bundesliga, 
#' \emph{Journal of Applied Statistics}, 45(9), 1561 - 1578
#' @source
#' \url{https://www.kicker.de/}
#' @keywords datasets
#' @seealso \code{\link{Buli1516}}, \code{\link{Buli1617}}, \code{\link{Buli1718}}
#' @examples
#' \dontrun{
#' op <- par(no.readonly = TRUE)
#' 
#' data(Buli1415)
#' 
#' Y <- Buli1415$Y5
#' Z1 <- scale(Buli1415$Z1, scale = FALSE)
#' 
#' ctrl.buli <- ctrl.BTLLasso(object.order.effect = TRUE, 
#'                            name.order = "Home", 
#'                            penalize.order.effect.diffs = TRUE, 
#'                            penalize.order.effect.absolute = FALSE,
#'                            order.center = TRUE, lambda2 = 1e-2)
#' 
#' set.seed(1860)
#' m.buli <- cv.BTLLasso(Y = Y, Z1 = Z1, control = ctrl.buli)
#' m.buli
#' 
#' par(xpd = TRUE, mar = c(5,4,4,6))
#' plot(m.buli)
#' 
#' par(op)
#' }
NULL


#' Bundesliga Data 2016/17 (Buli1617)
#' 
#' Data from the German Bundesliga from the season 2016/17. 
#' The data contain all 306 matches of the season treated as paired comparisons with 5 (Y5) or 3 (Y3) different 
#' response categories. Additionally, different match-specific covariates are given as, for example, 
#' the percentage of ball possession or the total running distance per team and per match.
#' 
#' @name Buli1617
#' @docType data
#' @format A list containing data from the German Bundesliga with 306 observations. 
#' The list contains both information on the response (paired comparisons) and different covariates.
#' \describe{ 
#' \item{Y5}{A response.BTLLasso object with 5 response categories for the Buli1516 data including
#' \itemize{
#' \item{response: Ordinal paired comparison response vector} 
#' \item{first.object: Vector containing the first-named team per paired comparison (home team)}
#' \item{second.object: Vector containing the second-named team per paired comparison (away team)}
#' \item{subject: Vector containing a match-day identifier per paired comparison}
#' \item{with.order} Vector containing information that each match has to be considered including an order effect.
#' }}
#' \item{Y3}{A response.BTLLasso object with 3 response categories for the Buli1516 data including
#' \itemize{
#' \item{response: Ordinal paired comparison response vector} 
#' \item{first.object: Vector containing the first-named team per paired comparison (home team)}
#' \item{second.object: Vector containing the second-named team per paired comparison (away team)}
#' \item{subject: Vector containing a match-day identifier per paired comparison}
#' \item{with.order} Vector containing information that each match has to be considered including an order effect.
#' }}
#' \item{Z1}{Matrix containing all team-match-specific covariates
#' \itemize{
#' \item{Distance: Total amount of km run} 
#' \item{BallPossession: Percentage of ball possession}
#' \item{TacklingRate: Rate of won tacklings}
#' \item{ShotsonGoal: Total number of shots on goal} 
#' \item{CompletionRate: Percentage of passes reaching teammates} 
#' \item{FoulsSuffered: Number of fouls suffered} 
#' \item{Offside: Number of offsides (in attack)}
#' \item{Corners: Number of corners (in attack)}
#' }
#' }
#' }
#'  @references Schauberger, Gunther and Tutz, Gerhard (2019): BTLLasso - A Common Framework and Software 
#' Package for the Inclusion  and Selection of Covariates in Bradley-Terry Models, \emph{Journal of 
#' Statistical Software}, to appear
#' 
#' Schauberger, Gunther and Tutz, Gerhard (2017): Subject-specific modelling 
#' of paired comparison data: A lasso-type penalty approach, \emph{Statistical Modelling},
#' 17(3), 223 - 243
#' 
#' Schauberger, Gunther, Groll Andreas and Tutz, Gerhard (2018): 
#' Analysis of the importance of on-field covariates in the German Bundesliga, 
#' \emph{Journal of Applied Statistics}, 45(9), 1561 - 1578
#' @source
#' \url{https://www.kicker.de/}
#' @keywords datasets
#' @seealso \code{\link{Buli1415}}, \code{\link{Buli1516}}, \code{\link{Buli1718}}
#' @examples
#' \dontrun{
#' op <- par(no.readonly = TRUE)
#' 
#' data(Buli1617)
#' 
#' Y <- Buli1617$Y5
#' Z1 <- scale(Buli1617$Z1, scale = FALSE)
#' 
#' ctrl.buli <- ctrl.BTLLasso(object.order.effect = TRUE, 
#'                            name.order = "Home", 
#'                            penalize.order.effect.diffs = TRUE, 
#'                            penalize.order.effect.absolute = FALSE,
#'                            order.center = TRUE, lambda2 = 1e-2)
#' 
#' set.seed(1860)
#' m.buli <- cv.BTLLasso(Y = Y, Z1 = Z1, control = ctrl.buli)
#' m.buli
#' 
#' par(xpd = TRUE, mar = c(5,4,4,6))
#' plot(m.buli)
#' 
#' par(op)
#' }
NULL


#' Bundesliga Data Response Data (BuliResponse)
#' 
#' Data from the German Bundesliga from the season 2015/16. The data contain all 
#' variables from the 306 matches that are necessary to create the respective 
#' \code{response.BTLLasso} object from the data set \code{\link{Buli1516}}.  The purpose
#' of the data set is to provide an example how \code{response.BTLLasso} objects can be created.
#' 
#' @name BuliResponse
#' @docType data
#' @format A data set containing all information that is necessary to create a response object
#' for the Bundesliga data \code{link{Buli1516}}
#' \describe{ 
#' \item{Result}{Ordinal, 5-categorical results from Bundesliga season 2015/16.}
#' \item{TeamHome}{Abbreviation of home team.}
#' \item{TeamAway}{Abbreviation of away team.}
#' \item{Matchday}{Matchdays from 1 to 34.}
#' }
#'  @references Schauberger, Gunther and Tutz, Gerhard (2019): BTLLasso - A Common Framework and Software 
#' Package for the Inclusion  and Selection of Covariates in Bradley-Terry Models, \emph{Journal of 
#' Statistical Software}, to appear
#' 
#' Schauberger, Gunther and Tutz, Gerhard (2017): Subject-specific modelling 
#' of paired comparison data: A lasso-type penalty approach, \emph{Statistical Modelling},
#' 17(3), 223 - 243
#' 
#' Schauberger, Gunther, Groll Andreas and Tutz, Gerhard (2018): 
#' Analysis of the importance of on-field covariates in the German Bundesliga, 
#' \emph{Journal of Applied Statistics}, 45(9), 1561 - 1578
#' @source
#' \url{https://www.kicker.de/}
#' @keywords datasets
#' @examples
#' \dontrun{
#' data(BuliResponse)
#' 
#' Y.Buli <- response.BTLLasso(response = BuliResponse$Result, 
#'                             first.object = BuliResponse$TeamHome,
#'                             second.object = BuliResponse$TeamAway,
#'                             subject = BuliResponse$Matchday)
#' }
NULL


#' German Longitudinal Election Study (GLES)
#' 
#' Data from the German Longitudinal Election Study (GLES), see Rattinger et
#' al. (2014). The GLES is a long-term study of the German electoral process.
#' It collects pre- and post-election data for several federal elections, the
#' data used here originate from the pre-election study for 2013.
#' 
#' @name GLES
#' @docType data
#' @format A list containing data from the German Longitudinal Election Study with 2003 
#' (partly incomplete) observations. 
#' The list contains both information on the response (paired comparisons) and different covariates.
#' \describe{ 
#' \item{Y}{A response.BTLLasso object for the GLES data including
#' \itemize{
#' \item{response: Ordinal paired comparison response vector} 
#' \item{first.object: Vector containing the first-named party per paired comparison}
#' \item{second.object: Vector containing the second-named party per paired comparison}
#' \item{subject: Vector containing a person identifier per paired comparison}
#' \item{with.order} Automatically generated vector containing information on order effect. Irrelevant, because 
#' no order effect needs to be included in the analysis of GLES data. 
#' }}
#' \item{X}{Matrix containing all eight person-specific covariates
#' \itemize{
#' \item{Age: Age in years} 
#' \item{Gender (0: male, 1: female)}
#' \item{EastWest (0: West Germany, 1: East Germany)}
#' \item{PersEcon: Personal economic situation, 1: good or very good,
#' 0: else} 
#' \item{Abitur: School leaving certificate, 1: Abitur/A
#' levels, 0: else} 
#' \item{Unemployment: 1: currently unemployed, 0:
#' else} 
#' \item{Church: Frequency of attendence in a
#' church/synagogue/mosque/..., 1: at least once a month, 0: else}
#' \item{Migration: Are you a migrant / not German since birth? 1: yes,
#' 0: no} 
#' }
#' }
#' \item{Z1}{Matrix containing all four person-party-specific covariates
#' \itemize{
#' \item{Climate: Self-perceived distance of each person to all five parties with respect to 
#' ones attitude towards climate change.}
#' \item{SocioEcon: Self-perceived distance of each person to all five parties with respect to 
#' ones attitude towards socio-economic issues.}
#' \item{Immigration: Self-perceived distance of each person to all five parties with respect to 
#' ones attitude towards immigration.}
#' }
#' }
#' }
#' @references Rattinger, H., S. Rossteutscher, R. Schmitt-Beck, B. Wessels,
#' and C. Wolf (2014): Pre-election cross section (GLES 2013). \emph{GESIS Data
#' Archive, Cologne ZA5700 Data file Version 2.0.0.}
#' 
#' Schauberger, Gunther and Tutz, Gerhard (2019): BTLLasso - A Common Framework and Software 
#' Package for the Inclusion  and Selection of Covariates in Bradley-Terry Models, \emph{Journal of 
#' Statistical Software}, to appear
#' 
#' Schauberger, Gunther and Tutz, Gerhard (2017): Subject-specific modelling 
#' of paired comparison data: A lasso-type penalty approach, \emph{Statistical Modelling},
#' 17(3), 223 - 243
#' @source
#' \url{https://gles-en.eu/}
#' @keywords datasets
#' @examples
#' \dontrun{
#' op <- par(no.readonly = TRUE)
#' 
#' data(GLES)
#' Y <- GLES$Y
#' X <- scale(GLES$X, scale = FALSE)
#' 
#' subs <- c("(in years)","female (1); male (0)","East Germany (1); West Germany (0)",
#'           "(very) good (1); else (0)", "Abitur/A levels (1); else (0)", 
#'           "currently unemployed (1); else (0)","at least once a month (1); else (0)",
#'           "yes (1); no (0)")
#' 
#' set.seed(5)
#' m.gles <- cv.BTLLasso(Y = Y, X = X, control = ctrl.BTLLasso(l.lambda = 50))
#' 
#' par(xpd = TRUE, mar = c(5,4,4,6))
#' plot(m.gles, subs.X = subs)
#' 
#' par(op)
#' }
#' 
NULL





#' Subset of the GLES data set with 200 observations and 4 covariates.
#' 
#' This is a subset of the \code{\link{GLES}} data set from the German
#' Longitudinal Election Study (GLES), see Rattinger et al. (2014). The subset contains 
#' only 200 of the 2003 observations and only  a small part of the covariates. The GLES is
#' a long-term study of the German electoral process. It collects pre- and
#' post-election data for several federal elections, the data used here
#' originate from the pre-election study for 2013.
#' 
#' @name GLESsmall
#' @docType data
#' @format A list containing data from the German Longitudinal Election Study with 200 observations. 
#' The list contains both information on the response (paired comparisons) and different covariates.
#' \describe{ 
#' \item{Y}{A response.BTLLasso object for the GLES data including
#' \itemize{
#' \item{response: Ordinal paired comparison response vector} 
#' \item{first.object: Vector containing the first-named party per paired comparison}
#' \item{second.object: Vector containing the second-named party per paired comparison}
#' \item{subject: Vector containing a person identifier per paired comparison}
#' \item{with.order} Automatically generated vector containing information on order effect. Irrelevant, because 
#' no order effect needs to be included in the analysis of GLES data. 
#' }}
#' \item{X}{Matrix containing all eight person-specific covariates
#' \itemize{
#' \item{Age: Age in years} 
#' \item{Gender (0: male, 1: female)}
#' }
#' }
#' \item{Z1}{Matrix containing all four person-party-specific covariates
#' \itemize{
#' \item{Climate: Self-perceived distance of each person to all five parties with respect to 
#' ones attitude towards climate change.}
#' \item{Immigration: Self-perceived distance of each person to all five parties with respect to 
#' ones attitude towards immigration.}
#' }
#' }
#' }
#' @references Rattinger, H., S. Rossteutscher, R. Schmitt-Beck, B. Wessels,
#' and C. Wolf (2014): Pre-election cross section (GLES 2013). \emph{GESIS Data
#' Archive, Cologne ZA5700 Data file Version 2.0.0.}
#' 
#' Schauberger, Gunther and Tutz, Gerhard (2019): BTLLasso - A Common Framework and Software 
#' Package for the Inclusion  and Selection of Covariates in Bradley-Terry Models, \emph{Journal of 
#' Statistical Software}, to appear
#' 
#' Schauberger, Gunther and Tutz, Gerhard (2017): Subject-specific modelling 
#' of paired comparison data: A lasso-type penalty approach, \emph{Statistical Modelling},
#' 17(3), 223 - 243
#' @source
#' \url{https://gles-en.eu/}
#' @keywords datasets
#' @seealso \code{\link{GLES}}
#' @examples
#' 
#' \dontrun{
#' op <- par(no.readonly = TRUE)
#' 
#' data(GLESsmall)
#' 
#' ## extract data and center covariates for better interpretability
#' Y <- GLESsmall$Y
#' X <- scale(GLESsmall$X, scale = FALSE)
#' Z1 <- scale(GLESsmall$Z1, scale = FALSE)
#' 
#' ## vector of subtitles, containing the coding of the X covariates
#' subs.X <- c('', 'female (1); male (0)')
#' 
#' ## Cross-validate BTLLasso model
#' m.gles.cv <- cv.BTLLasso(Y = Y, X = X, Z1 = Z1)
#' m.gles.cv
#' 
#' coef(m.gles.cv)
#' logLik(m.gles.cv)
#' 
#' head(predict(m.gles.cv, type="response"))
#' head(predict(m.gles.cv, type="trait"))
#' 
#' par(xpd = TRUE, mar = c(5,4,4,6))
#' plot(m.gles.cv, subs.X = subs.X, plots_per_page = 4, which = 2:5)
#' paths(m.gles.cv, y.axis = 'L2')
#' 
#' par(op)
#' }
#' 
NULL

#' Simulated data set for illustration
#' 
#' This data set is a simulated data set including all possible types of covariates (X, Z1 and Z2)
#' and is intended to serve for illustration purpose. The data set contains paired comparisons between
#' four objects with five different response categories from 200 subjects.
#' 
#' @name SimData
#' @docType data
#' @format A list containing  simulated data for 200 observations. 
#' The list contains both information on the response (paired comparisons) and different covariates.
#' \describe{ 
#' \item{Y}{A response.BTLLasso object with simulated responses including
#' \itemize{
#' \item{response: Ordinal paired comparison response vector} 
#' \item{first.object: Vector containing the first-named object per paired comparison}
#' \item{second.object: Vector containing the second-named object per paired comparison}
#' \item{subject: Vector containing a subject identifier per paired comparison}
#' \item{with.order} Automatically generated vector containing information on order effect. Each paired 
#' comparison is associated with an order effect.
#' }}
#' \item{X}{Matrix containing both subject-specific covariates
#' \itemize{
#' \item{X_var1} 
#' \item{X_var2}
#' }
#' }
#' \item{Z1}{Matrix containing both subject-object-specific covariates
#' \itemize{
#' \item{Z1_var1}
#' \item{Z1_var2}
#' }
#' }
#' \item{Z2}{Matrix containing both object-specific covariates
#' \itemize{
#' \item{Z2_var1}
#' \item{Z2_var2}
#' }
#' }
#' }
#' @keywords datasets
#' @examples
#' 
#' \dontrun{
#' op <- par(no.readonly = TRUE)
#' 
#' data(SimData)
#' 
#' ## Specify control argument
#' ## -> allow for object-specific order effects and penalize intercepts
#' ctrl <- ctrl.BTLLasso(penalize.intercepts = TRUE, object.order.effect = TRUE,
#'                       penalize.order.effect.diffs = TRUE)
#' 
#' ## Simple BTLLasso model for tuning parameters lambda
#' m.sim <- BTLLasso(Y = SimData$Y, X = SimData$X, Z1 = SimData$Z1,
#'                   Z2 = SimData$Z2, control = ctrl)
#' m.sim
#' 
#' par(xpd = TRUE)
#' plot(m.sim)
#' 
#' 
#' ## Cross-validate BTLLasso model for tuning parameters lambda
#' set.seed(1860)
#' m.sim.cv <- cv.BTLLasso(Y = SimData$Y, X = SimData$X, Z1 = SimData$Z1,
#'                         Z2 = SimData$Z2, control = ctrl)
#' m.sim.cv
#' coef(m.sim.cv)
#' logLik(m.sim.cv)
#' 
#' head(predict(m.sim.cv, type="response"))
#' head(predict(m.sim.cv, type="trait"))
#' 
#' plot(m.sim.cv, plots_per_page = 4)
#' 
#' 
#' ## Example for bootstrap intervals for illustration only
#' ## Don't calculate bootstrap intervals with B = 20!!!!
#' set.seed(1860)
#' m.sim.boot <- boot.BTLLasso(m.sim.cv, B = 20, cores = 20)
#' m.sim.boot
#' plot(m.sim.boot, plots_per_page = 4)
#' 
#' par(op)
#' }
NULL

#' Bundesliga Data 2017/18 (Buli1718)
#' 
#' Data from the German Bundesliga from the season 2017/18. 
#' The data contain all 306 matches of the season treated as paired comparisons with 5 (Y5) or 3 (Y3) different 
#' response categories. Additionally, different match-specific covariates are given as, for example, 
#' the percentage of ball possession or the total running distance per team and per match.
#' 
#' @name Buli1718
#' @docType data
#' @format A list containing data from the German Bundesliga with 306 observations. 
#' The list contains both information on the response (paired comparisons) and different covariates.
#' \describe{ 
#' \item{Y5}{A response.BTLLasso object with 5 response categories for the Buli1516 data including
#' \itemize{
#' \item{response: Ordinal paired comparison response vector} 
#' \item{first.object: Vector containing the first-named team per paired comparison (home team)}
#' \item{second.object: Vector containing the second-named team per paired comparison (away team)}
#' \item{subject: Vector containing a match-day identifier per paired comparison}
#' \item{with.order} Vector containing information that each match has to be considered including an order effect.
#' }}
#' \item{Y3}{A response.BTLLasso object with 3 response categories for the Buli1516 data including
#' \itemize{
#' \item{response: Ordinal paired comparison response vector} 
#' \item{first.object: Vector containing the first-named team per paired comparison (home team)}
#' \item{second.object: Vector containing the second-named team per paired comparison (away team)}
#' \item{subject: Vector containing a match-day identifier per paired comparison}
#' \item{with.order} Vector containing information that each match has to be considered including an order effect.
#' }}
#' \item{Z1}{Matrix containing all team-match-specific covariates
#' \itemize{
#' \item{Distance: Total amount of km run} 
#' \item{BallPossession: Percentage of ball possession}
#' \item{TacklingRate: Rate of won tacklings}
#' \item{ShotsonGoal: Total number of shots on goal} 
#' \item{CompletionRate: Percentage of passes reaching teammates} 
#' \item{FoulsSuffered: Number of fouls suffered} 
#' \item{Offside: Number of offsides (in attack)}
#' \item{Corners: Number of corners (in attack)}
#' }
#' }
#' }
#'  @references Schauberger, Gunther and Tutz, Gerhard (2019): BTLLasso - A Common Framework and Software 
#' Package for the Inclusion  and Selection of Covariates in Bradley-Terry Models, \emph{Journal of 
#' Statistical Software}, to appear
#' 
#' Schauberger, Gunther and Tutz, Gerhard (2017): Subject-specific modelling 
#' of paired comparison data: A lasso-type penalty approach, \emph{Statistical Modelling},
#' 17(3), 223 - 243
#' 
#' Schauberger, Gunther, Groll Andreas and Tutz, Gerhard (2018): 
#' Analysis of the importance of on-field covariates in the German Bundesliga, 
#' \emph{Journal of Applied Statistics}, 45(9), 1561 - 1578
#' @source
#' \url{https://www.kicker.de/}
#' @keywords datasets
#' @seealso \code{\link{Buli1415}}, \code{\link{Buli1516}}, \code{\link{Buli1617}}
#' @examples
#' \dontrun{
#' op <- par(no.readonly = TRUE)
#' 
#' data(Buli1718)
#' 
#' Y <- Buli1718$Y5
#' Z1 <- scale(Buli1718$Z1, scale = FALSE)
#' 
#' ctrl.buli <- ctrl.BTLLasso(object.order.effect = TRUE, 
#'                            name.order = "Home", 
#'                            penalize.order.effect.diffs = TRUE, 
#'                            penalize.order.effect.absolute = FALSE,
#'                            order.center = TRUE, lambda2 = 1e-2)
#' 
#' set.seed(1860)
#' m.buli <- cv.BTLLasso(Y = Y, Z1 = Z1, control = ctrl.buli)
#' m.buli
#' 
#' par(xpd = TRUE, mar = c(5,4,4,6))
#' plot(m.buli)
#' 
#' par(op)
#' }
NULL


