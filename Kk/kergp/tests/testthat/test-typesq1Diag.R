context("q1Diag")

## ============================================================================
## check that for this covariance class, 'covMat' gives the
## same result whatever be the class of input used
## "factor", "character" or "integer"
## ============================================================================

SchoolLevels <- c("Bad", "Mean", "Good")

School <- factor(1L:3L, labels = SchoolLevels)
     
## covariance 
myCov <- q1Diag(School, input = "School", cov = "hete",
                intAsChar = FALSE)
coef(myCov) <- c(1.1, 2.2, 3.3)

dfInt <- data.frame(School = sample(3, size = 40, replace = TRUE))
CInt <- covMat(myCov, X = dfInt)

dfChar <- within(dfInt, School <- SchoolLevels[School])  
CChar <- covMat(myCov, X = dfChar)

dfFact <- within(dfChar, School <- factor(School, levels = SchoolLevels))  
CFact <- covMat(myCov, X = dfFact)

err <- max(c(max(abs(CInt - CChar)), max(abs(CInt - CFact))))
test_that(desc = "Identical covariances using int, char or factor: ",
          code = expect_true(err < 1e-10))

err <- max(c("IntChar" = max(abs(attr(CInt, "gradient") -
                                     attr(CChar, "gradient"))),
             "IntFact" = max(abs(attr(CInt, "gradient") -
                                     attr(CFact, "gradient")))))
test_that(desc = "Identical cov. Gradients using int, char or factor: ",
          code = expect_true(err < 1e-10))
