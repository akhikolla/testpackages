#' A multilevel dataset from Hofmann, Griffin, and Gavin (2000).
#' 
#' A multilevel dataset from Hofmann, Griffin, and Gavin (2000). 
#' @format A data frame with 1,000 observations and 7 variables.
#' \describe{
#'     \item{\code{id}}{a numeric vector of group ids.}
#'     \item{\code{helping}}{a numeric vector of the helping outcome variable construct.}
#'     \item{\code{mood}}{a level 1 mood predictor.}
#'     \item{\code{mood_grp_mn}}{a level 2 variable of the group mean of mood.}
#'     \item{\code{cohesion}}{a level 2 covariate measuring cohesion.}
#'     \item{\code{mood_grp_cent}}{group-mean centered mood predictor.}
#'     \item{\code{mood_grd_cent}}{grand-mean centered mood predictor.}
#' }
#' @source 
#' Hofmann, D.A., Griffin, M.A., & Gavin, M.B. (2000). The application of hierarchical linear modeling to management research. In K.J. Klein, & S.W.J. Kozlowski (Eds.), Multilevel theory, research, and methods in organizations: Foundations, extensions, and new directions (pp. 467-511).  Hoboken, NJ: Jossey-Bass.
#' @references Aguinis, H., & Culpepper, S.A. (2015). An expanded decision
#' making procedure for examining cross-level interaction effects with
#' multilevel modeling. \emph{Organizational Research Methods}.
#' Available at: \url{http://hermanaguinis.com/pubs.html}
#' @seealso \code{\link[lme4]{lmer}}, \code{\link{model.matrix}},
#'          \code{\link[lme4]{VarCorr}}, \code{\link[RLRsim]{LRTSim}},
#'          \code{\link{simICCdata}}
#' @examples
#' \dontrun{
#' 
#' if(requireNamespace("lme4") && requireNamespace("RLRsim")){ 
#' data(Hofmann)
#' library("lme4")
#' 
#' # Random-Intercepts Model
#' lmmHofmann0 = lmer(helping ~ (1|id), data = Hofmann)
#' vy_Hofmann = var(Hofmann[,'helping'])
#' 
#' # Computing icca
#' VarCorr(lmmHofmann0)$id[1,1]/vy_Hofmann
#' 
#' # Estimating Group-Mean Centered Random Slopes Model, no level 2 variables
#' lmmHofmann1  <- lmer(helping ~ mood_grp_cent + (mood_grp_cent |id),
#'                      data = Hofmann, REML = FALSE)
#' X_Hofmann = model.matrix(lmmHofmann1)
#' P = ncol(X_Hofmann)
#' T1_Hofmann  = VarCorr(lmmHofmann1)$id[1:P,1:P]
#' 
#' # Computing iccb
#' icc_beta(X_Hofmann, Hofmann[,'id'], T1_Hofmann, vy_Hofmann)$rho_beta
#' 
#' # Performing LR test
#' # Need to install 'RLRsim' package
#' library("RLRsim")
#' lmmHofmann1a  <- lmer(helping ~ mood_grp_cent + (1 | id),
#'                       data = Hofmann, REML = FALSE)
#' obs.LRT <- 2*(logLik(lmmHofmann1) - logLik(lmmHofmann1a))[1]
#' X <- getME(lmmHofmann1,"X")
#' Z <- t(as.matrix(getME(lmmHofmann1,"Zt")))
#' sim.LRT <- LRTSim(X, Z, 0, diag(ncol(Z)))
#' (pval <- mean(sim.LRT > obs.LRT))
#' } else {
#'  stop("Please install packages `RLRsim` and `lme4` to run the above example.")
#' }
#' }
"Hofmann"

#' Simulated data example from Aguinis and Culpepper (2015).
#'
#' A simulated data example from Aguinis and Culpepper (2015) to demonstrate
#' the \code{\link{icc_beta}} function for computing the proportion of variance
#' in the outcome variable that is attributed to heterogeneity in slopes due to
#' higher-order processes/units. 
#' @format A data frame with 900 observations (i.e., 30 observations nested 
#' within 30 groups) on the following 6 variables.
#' \describe{
#'     \item{\code{l1id}}{A within group ID variable.}
#'     \item{\code{l2id}}{A group ID variable.}
#'     \item{\code{one}}{A column of 1's for the intercept.}
#'     \item{\code{X1}}{A simulated level 1 predictor.}
#'     \item{\code{X2}}{A simulated level 1 predictor.}
#'     \item{\code{Y}}{A simulated outcome variable.}
#' }
#' @details
#' See Aguinis and Culpepper (2015) for the model used to simulate the dataset. 
#' @source Aguinis, H., & Culpepper, S.A. (2015). An expanded decision
#' making procedure for examining cross-level interaction effects with
#' multilevel modeling. \emph{Organizational Research Methods}.
#' Available at: \url{http://www.hermanaguinis.com/pubs.html}
#' @seealso \code{\link[lme4]{lmer}}, \code{\link{model.matrix}},
#'          \code{\link[lme4]{VarCorr}}, \code{\link[RLRsim]{LRTSim}},
#'          \code{\link{Hofmann}}
#' @examples
#' \dontrun{
#' data(simICCdata)
#' if(requireNamespace("lme4")){ 
#' library("lme4")
#' 
#' # computing icca
#' vy <- var(simICCdata$Y)
#' lmm0 <- lmer(Y ~ (1|l2id), data = simICCdata, REML = FALSE)
#' VarCorr(lmm0)$l2id[1,1]/vy
#' 
#' # Create simICCdata2
#' grp_means = aggregate(simICCdata[c('X1','X2')], simICCdata['l2id'],mean)
#' colnames(grp_means)[2:3] = c('m_X1','m_X2')
#' simICCdata2 = merge(simICCdata, grp_means, by='l2id')
#' 
#' # Estimating random slopes model
#' lmm1  <- lmer(Y ~ I(X1-m_X1) + I(X2-m_X2) + (I(X1-m_X1) + I(X2-m_X2) | l2id),
#'               data = simICCdata2, REML = FALSE)
#' X <- model.matrix(lmm1)
#' p <- ncol(X)
#' T1 <- VarCorr(lmm1) $l2id[1:p,1:p]
#' # computing iccb
#' # Notice '+1' because icc_beta assumes l2ids are from 1 to 30.
#' icc_beta(X, simICCdata2$l2id+1, T1, vy)$rho_beta
#' } else {
#'  stop("Please install `lme4` to run the above example.")
#' }
#' }
"simICCdata"