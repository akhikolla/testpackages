context("q1LowRank")

## ============================================================================
## check that for this covariance class, 'covMat' gives the
## same result whatever be the class of input used
## "factor", "character" or "integer"
## ============================================================================

set.seed(314159)

LettersLevels <- letters[1:8]
myLetters <- factor(LettersLevels)
myCov <- q1LowRank(factor = myLetters, input = "u", rank = 3, cov = "homo",
                   intAsChar = FALSE)
np <- length(coef(myCov))
set.seed(123)
coef(myCov) <- c(runif(np - 1, min = 0, max = pi), rexp(1))

dfInt <- data.frame(u = sample(8, size = 40, replace = TRUE))
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
