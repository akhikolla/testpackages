context("covOrd")

## ============================================================================
## define an ordinal kernel. Choose levels likely to generate problems
## because 'as.factor(as.character(8:12))' puts the levels in the
## order 10, 11, 12, 8, 9
## ============================================================================

kcat <- covOrd(ordered = 8:12,
               k1Fun1 = k1Fun1Cos, 
               warpFun = "norm",
               hasGrad = TRUE, cov = "homo")
## print(levels(kcat))

## as used in the code of the 'covMat' method for the class "covQual"

object <- kcat
C1 <- object@covLevels(object@par, lowerSQRT = FALSE, compGrad = FALSE)

C2 <- object@covLevMat

## as used by an human

C3 <- covMat(kcat)
attr(C3, "gradient") <- NULL

## with a 'X' formal given either with the factor alone or with more
## columns

dfu <- data.frame(u = as.ordered(8:12))
dfux <- data.frame(u = as.ordered(8:12), x = rep(0, 5))

C4 <- covMat(kcat, X = dfu)
C5 <- covMat(kcat, X = dfux)
attr(C4, "gradient") <- attr(C5, "gradient") <- NULL    

test_that(desc = "covMat code", expect_lt(max(abs(C1 - C2)), 1e-8))
test_that(desc = "covMat 'X' missing", expect_lt(max(abs(C1 - C3)), 1e-8))
test_that(desc = "covMat 'X' one column", expect_lt(max(abs(C1 - C4)), 1e-8))
test_that(desc = "covMat 'X' two columns",  expect_lt(max(abs(C1 - C5)), 1e-8))
    
    

