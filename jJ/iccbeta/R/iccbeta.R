#' Intraclass correlation used to assess variability of lower-order
#' relationships across higher-order processes/units.
#'
#' A function and vignettes for computing the intraclass correlation described
#' in Aguinis & Culpepper (2015). iccbeta quantifies the share of variance
#' in an outcome variable that is attributed to heterogeneity in slopes due to
#' higher-order processes/units.
#' 
#' @param x    A `lmer` model object or a design matrix with no missing values.
#' @param ...  Additional parameters...
#' 
#' @return 
#' A `list` with:
#' 
#' - `J`
#' - `means`
#' - `XcpXc`
#' - `Nj`
#' - `rho_beta`
#' 
#' @author 
#' Steven Andrew Culpepper
#' 
#' @seealso 
#' 
#' [lme4::lmer()], [model.matrix()],
#' [lme4::VarCorr()], [RLRsim::LRTSim()],
#' [iccbeta::Hofmann], and [iccbeta::simICCdata]
#'          
#' @references
#' Aguinis, H., & Culpepper, S.A. (2015). An expanded decision making
#' procedure for examining cross-level interaction effects with multilevel
#' modeling. _Organizational Research Methods_. Available at:
#' <http://hermanaguinis.com/pubs.html>
#' 
#' @export
#' @examples
#' \dontrun{
#' 
#' if(requireNamespace("lme4") && requireNamespace("RLRsim")){
#' 
#' ## Example 1: Simulated Data Example from Aguinis & Culpepper (2015) ----
#' data(simICCdata)
#' library("lme4")
#' 
#' # Computing icca
#' vy <- var(simICCdata$Y)
#' lmm0 <- lmer(Y ~ (1 | l2id), data = simICCdata, REML = FALSE)
#' VarCorr(lmm0)$l2id[1, 1]/vy
#' 
#' # Create simICCdata2
#' grp_means = aggregate(simICCdata[c('X1', 'X2')], simICCdata['l2id'], mean)
#' colnames(grp_means)[2:3] = c('m_X1', 'm_X2')
#' simICCdata2 = merge(simICCdata, grp_means, by='l2id')
#' 
#' # Estimating random slopes model
#' lmm1  <- lmer(Y ~ I(X1 - m_X1) + I(X2 - m_X2) + 
#'                  (I(X1 - m_X1) + I(X2 - m_X2) | l2id),
#'               data = simICCdata2, REML = FALSE)
#' 
#' ## iccbeta calculation on `lmer` object
#' icc_beta(lmm1)
#' 
#' ## Manual specification of iccbeta
#' 
#' # Extract components from model.
#' X <- model.matrix(lmm1)
#' p <- ncol(X)
#' T1  <- VarCorr(lmm1)$l2id[1:p,1:p]
#' 
#' # Note: vy was computed under "icca"
#' 
#' # Computing iccb
#' # Notice '+1' because icc_beta assumes l2ids are from 1 to 30.
#' icc_beta(X, simICCdata2$l2id + 1, T1, vy)$rho_beta
#' 
#' ## Example 2: Hofmann et al. (2000)   ----
#' 
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
#' lmmHofmann1 <- lmer(helping ~ mood_grp_cent + (mood_grp_cent |id),
#'                     data = Hofmann, REML = FALSE)
#' 
#' ## Automatic calculation of iccbeta using the lmer model
#' amod = icc_beta(lmmHofmann1)
#' 
#' ## Manual calculation of iccbeta
#' 
#' X_Hofmann <- model.matrix(lmmHofmann1)
#' P <- ncol(X_Hofmann)
#' T1_Hofmann <- VarCorr(lmmHofmann1)$id[1:P,1:P]
#' 
#' # Computing iccb
#' bmod = icc_beta(X_Hofmann, Hofmann[,'id'], T1_Hofmann, vy_Hofmann)$rho_beta
#' 
#' # Performing LR test
#' library("RLRsim")
#' lmmHofmann1a <- lmer(helping ~ mood_grp_cent + (1 |id),
#'                      data = Hofmann, REML = FALSE)
#' obs.LRT <- 2*(logLik(lmmHofmann1) - logLik(lmmHofmann1a))[1]
#' X <- getME(lmmHofmann1,"X")
#' Z <- t(as.matrix(getME(lmmHofmann1,"Zt")))
#' sim.LRT <- LRTSim(X, Z, 0, diag(ncol(Z)))
#' (pval <- mean(sim.LRT > obs.LRT))
#' } else {
#'  stop("Please install packages `RLRsim` and `lme4` to run the above example.") 
#' } 
#' 
#' }
icc_beta = function(x, ...) {
    UseMethod("icc_beta")
}

#' @importFrom lme4 VarCorr
#' @importClassesFrom lme4 lmerMod
#' @importFrom stats model.frame model.matrix model.response var
#' @importFrom methods slot
#' @rdname icc_beta
#' @export
icc_beta.lmerMod = function(x, ...) {
    
    stopifnot(inherits(x, "lmerMod"))
    
    # Retrieve the model frame.
    mf = model.frame(x)
    
    # Extract the design matrix
    X_design <- model.matrix(x)
    
    # Extract the response and calculate variance
    vy <- var(model.response(mf))
    
    # Extract the estimated variance-covariance matrix of a lmer model fit
    p <- ncol(X_design)

    # Figure out the class term:
    term_labels <- attr(mf, "names")
    
    # Retrieve the _LAST_ variable (assume it is a group variable.)
    # To Do: What is the best way to make this approach robust ?
    group_term <- term_labels[length(term_labels)]

    T1  <- VarCorr(x)[[group_term]][1:p,1:p]
    
    # Extract the l2ids
    l2ids <- slot(x, "frame")[[group_term]]
    
    # Enforce IDs greater than 1
    low_id <- min(l2ids)
    
    if(low_id < 0) {
        stop("Please ensure that the L2 ids are not negative.")
    } else if(low_id == 0) {
        message("Added +1 to L2 IDs so that they are greater than 0.")
        
        # Notice '+1' because icc_beta assumes l2ids are from 1 to 30.
        l2ids <- l2ids + 1 
    }
    
    icc_beta.default(X_design, l2ids, T1, vy)
}

#' @param l2id A `vector` that identifies group membership. The vector
#'             must be coded as a sequence of integers from 1 to J,
#'             the number of groups.
#' @param T    A `matrix` of the estimated variance-covariance matrix of
#'             a lmer model fit.
#' @param vy   The variance of the outcome variable.
#' @rdname icc_beta
#' @export
icc_beta.default = function(x, l2id, T, vy, ...) {

    icc_beta_cpp(x, l2id, T, vy)
}


