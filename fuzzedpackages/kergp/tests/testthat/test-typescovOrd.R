context("covOrd")

## ============================================================================
## check that for this covariance class, 'covMat' gives the
## same result whatever be the class of input used
## "factor", "character" or "integer"
## ============================================================================
LettersLevels <- letters[1:6]
u <- ordered(1:6, labels = letters[1:6])
     
myCov <- covOrd(ordered = u, cov = "homo", intAsChar = FALSE)
coef(myCov) <- c(mean = 0.5, sd = 1, theta = 3, sigma2 = 2)
myCov

dfInt <- data.frame(u = sample(6, size = 40, replace = TRUE))
CInt <- covMat(myCov, X = dfInt)

dfChar <- within(dfInt, u <- LettersLevels[u])  
CChar <- covMat(myCov, X = dfChar)

dfFact <- within(dfChar, u <- factor(u, levels = LettersLevels))  
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
