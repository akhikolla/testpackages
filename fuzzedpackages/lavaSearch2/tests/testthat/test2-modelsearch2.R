### test2-modelsearch2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jan 22 2018 (11:45) 
## Version: 
## Last-Updated: jun 27 2019 (14:21) 
##           By: Brice Ozenne
##     Update #: 24
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * header
rm(list = ls())
if(FALSE){ ## already called in test-all.R
    library(testthat)
    library(lavaSearch2)
}

lava.options(symbols = c("~","~~"))

context("lavaSearch2")

## * example with one additional link
n <- 100
m.sim <- lvm(Y~E+0*X1)
m <- lvm(Y~E)
addvar(m) <- ~X1

set.seed(12)
df.sim <- lava::sim(m.sim, n=100, latent = FALSE)
e.base <- estimate(m, data = df.sim)

test_that("Score 1 link",{
    GS.score <- modelsearch(e.base, silent = TRUE)
    index.coef <- which(GS.score$res[,"Index"]=="Y~X1")

    search.holm <- modelsearch2(e.base, method.p.adjust = "holm", trace = 0)
    ## c("E", "hessian", "varS", "outer", "sandwich", "robust", "num"),   "outer"
    
    expect_equal(as.double(GS.score$test[index.coef,"Test Statistic"]),
                 as.double(search.holm$sequenceTest[[1]][1,"statistic"]), tol = 1e-9)
    expect_equal(as.double(GS.score$test[index.coef,"P-value"]),
                 as.double(search.holm$sequenceTest[[1]][1,"p.value"]), tol = 1e-9)


    search.approx <- modelsearch2(e.base, method.p.adjust = "fastmax", trace = 0,
                                  method.maxdist = "approximate")
    search.resampling <- modelsearch2(e.base, method.p.adjust = "fastmax", trace = 0,
                                      method.maxdist = "resampling")
    search.bootstrap <- modelsearch2(e.base, method.p.adjust = "fastmax", trace = 0,
                                      method.maxdist = "bootstrap")
    
    expect_equal(search.approx$sequenceTest[[1]][1,"statistic"],
                 search.holm$sequenceTest[[1]][1,"statistic"])
    expect_equal(search.approx$sequenceTest[[1]][1,"statistic"],
                 search.resampling$sequenceTest[[1]][1,"statistic"], tol = 1e-3)
    expect_equal(search.approx$sequenceTest[[1]][1,"statistic"],
                 search.bootstrap$sequenceTest[[1]][1,"statistic"], tol = 1e-3)

    expect_equal(round(search.resampling$sequenceTest[[1]][1,"p.value"],2),
                 0.24)
    expect_equal(round(search.bootstrap$sequenceTest[[1]][1,"p.value"],2),
                 0.24)

    
})


##----------------------------------------------------------------------
### test2-modelsearch2.R ends here
