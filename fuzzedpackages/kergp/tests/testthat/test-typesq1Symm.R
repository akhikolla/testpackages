library(kergp)

context("q1Symm")

## ============================================================================
## check that for this covariance class, 'covMat' gives the
## same result whatever be the class of input used
## "factor", "character" or "integer"
## ============================================================================

set.seed(314159)

SchoolLevels <- c("Bad", "Mean", "Good")

School <- factor(1L:3L, labels = SchoolLevels)

## covariance 
myCov <- q1Symm(School, input = "School", cov = "hete",
                intAsChar = FALSE)
coef(myCov) <- c(runif(3, min = 0, max = pi), rexp(3))

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
