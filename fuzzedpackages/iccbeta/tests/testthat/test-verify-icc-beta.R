context("test-verify-icc-beta")

test_that("Check S3 works with manual approach", {
    
    data(Hofmann)
    library("lme4")
    
    # Random-Intercepts Model
    lmmHofmann0 = lmer(helping ~ (1|id), data = Hofmann)
    vy_Hofmann = var(Hofmann[,'helping'])
    
    # Computing icca
    VarCorr(lmmHofmann0)$id[1,1]/vy_Hofmann
    
    # Estimating Group-Mean Centered Random Slopes Model, no level 2 variables
    lmmHofmann1 <- lmer(helping ~ mood_grp_cent + (mood_grp_cent |id),
                        data = Hofmann, REML = FALSE)
    
    ## Automatic calculation of iccbeta using the lmer model
    amod = icc_beta(lmmHofmann1)
    
    ## Manual calculation of iccbeta
    
    X_Hofmann <- model.matrix(lmmHofmann1)
    P <- ncol(X_Hofmann)
    T1_Hofmann <- VarCorr(lmmHofmann1)$id[1:P,1:P]
    
    # Computing iccb
    bmod = icc_beta(X_Hofmann, Hofmann[,'id'], T1_Hofmann, vy_Hofmann)
    
    expect_equal(amod, bmod)
})


test_that("Check S3 works with manual approach bad IDs", {
    
    data(simICCdata)
    library("lme4")
    
    # Computing icca
    vy <- var(simICCdata$Y)
    lmm0 <- lmer(Y ~ (1 | l2id), data = simICCdata, REML = FALSE)
    VarCorr(lmm0)$l2id[1, 1]/vy
    
    # Create simICCdata2
    grp_means = aggregate(simICCdata[c('X1', 'X2')], simICCdata['l2id'], mean)
    colnames(grp_means)[2:3] = c('m_X1', 'm_X2')
    simICCdata2 = merge(simICCdata, grp_means, by='l2id')
    
    # Estimating random slopes model
    lmm1  <- lmer(Y ~ I(X1 - m_X1) + I(X2 - m_X2) + 
                      (I(X1 - m_X1) + I(X2 - m_X2) | l2id),
                  data = simICCdata2, REML = FALSE)
    
    ## Automatic iccbeta calculation on `lmer` object
    amod = icc_beta(lmm1)
    
    ## Manual specification of iccbeta
    
    # Extract components from model.
    X <- model.matrix(lmm1)
    p <- ncol(X)
    T1  <- VarCorr(lmm1)$l2id[1:p,1:p]
    
    # Note: vy was computed under "icca"
    
    # Computing iccb
    # Notice '+1' because icc_beta assumes l2ids are from 1 to 30.
    bmod = icc_beta(X, simICCdata2$l2id + 1, T1, vy)
    
    expect_equal(amod, bmod)
})
