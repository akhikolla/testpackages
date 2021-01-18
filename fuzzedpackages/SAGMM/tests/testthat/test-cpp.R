context("Check exported C++ Functions")
library(SAGMM)

test_that("mahalanobis_HD", {
    
    size<-sample(2:10,1)
    x<-runif(size)
    y<-runif(size)
    A <- matrix(runif(size^2)*2-1, ncol=size) 
    Sigma <- t(A) %*% A
    result<-mahalanobis_HD(x,y,Sigma)
    
    #types correct
    expect_is(result, "numeric")
   
    #no missing
    expect_identical(result, na.omit(result))
    
})

test_that("norm_HD", {
    
    size<-sample(2:10,1)
    x<-runif(size)
    y<-runif(size)
    A <- matrix(runif(size^2)*2-1, ncol=size) 
    Sigma <- t(A) %*% A
    result<-norm_HD(x,y,Sigma)
    
    #types correct
    expect_is(result, "numeric")
    
    #no missing
    expect_identical(result, na.omit(result))
    
})

#mainloop tested via SAGMMFit
