context("Check native R Functions")
library(SAGMM)

test_that("Test gainFactors", {
    
    number<-100
    burnin<-10
    result<-gainFactors(number, burnin)
    
    #types correct
    expect_is(result, "numeric")
    #size
    expect_equal(length(result), number)
    #no missing
    expect_identical(result, na.omit(result))
    
})

test_that("Test generateSimData", {
    
    result<-generateSimData(ngroups=10, Dimensions=10, Number=10^2)
   
    #types correct
    expect_is(result, "list")
    expect_is(result[[1]], "matrix")
    expect_is(result[[2]], "integer")
    expect_is(result[[3]], "MixSim")
    
    #no missing
    expect_identical(result, na.omit(result))
    
})

test_that("Test SAGMMFit", {
    sims<-generateSimData(ngroups=10, Dimensions=10, Number=10^3)
    result<-SAGMMFit(sims$X, sims$Y)
    
    # #types correct
    expect_is(result, "list")
    expect_is(result[[1]], "integer")
    expect_is(result[[2]], "logical")
    expect_is(result[[3]], "numeric")
    expect_is(result[[4]], "numeric")
    expect_is(result[[5]], "numeric")
    expect_is(result[[7]], "kmeans")
    expect_is(result[[6]], "numeric")
    expect_is(result[[8]], "numeric")
    expect_is(result[[9]], "matrix")
    expect_is(result[[10]], "list")
    
    #dims
    expect_equal(length(result), 10)
   
    #no missing
    expect_identical(result, na.omit(result))
    
    result<-SAGMMFit(sims$X, ngroups=5)
    # #types correct
    expect_is(result, "list")
    expect_is(result[[1]], "integer")
    expect_is(result[[2]], "logical")
    expect_is(result[[3]], "numeric")
    expect_is(result[[4]], "logical")
    expect_is(result[[5]], "logical")
    expect_is(result[[7]], "kmeans")
    expect_is(result[[6]], "logical")
    expect_is(result[[8]], "numeric")
    expect_is(result[[9]], "matrix")
    expect_is(result[[10]], "list")
    
    #dims
    expect_equal(length(result), 10)
    
})

test_that("Test SAGMMFit Plots", {
    sims<-generateSimData(ngroups=5, Dimensions=2, Number=10^2)
    result<-SAGMMFit(sims$X, ngroups=4, plot=TRUE)
    
    expect_is(result, "list")
    expect_equal(length(result), 10)
    expect_is(result[[2]], "NULL")

})
