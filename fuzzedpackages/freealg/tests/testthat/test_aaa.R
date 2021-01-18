## This file follows the structure of aaa.R in the free group package.

## Define some checker functions, and call them at the end.  They
## should all return TRUE if the package works, and stop with error if
## a test is failed.  Function checker1() has one argument, checker2()
## two, and checker3() has three.  

test_that("Test suite aaa.R",{

checker1 <- function(x){
    expect_true(x==x, info=x)

    expect_true(x == x + constant(0), info=x)
    expect_true(x == x + 0, info=x)
    expect_true(x == (x + 4) -4, info=x)
    expect_true(x == -(-x), info=x)
    expect_true(x == +(+x), info=x)

    expect_true(x+x-x == x, info=x)

    expect_true(is.zero(x-x), info=x)

    expect_true(0*x == constant(0), info=x)
    expect_true(1*x == x, info=x)
    expect_true(2*x == x+x, info=x)
    expect_true(3*x == x+x+x, info=x)
    expect_true(4*x == x+x+x+x, info=x)
    expect_true(5*x == x+x+x+x+x, info=x)
    expect_true(6*x == x+x+x+x+x+x, info=x)

    expect_true(x^0 == constant(1), info=x)
    expect_true(x^1 == x, info=x)
    expect_true(x^2 == x*x, info=x)
    expect_true(x^3 == x*x*x, info=x)
    expect_true(x^4 == x*x*x*x, info=x)
    
    ## check constant() and constant<-():
    ## checks below include 
    y <- x
    expect_true(constant(x) == constant(y), info=x)
    constant(y) <- 4
    expect_true(constant(y) == 4, info=x)
    constant(y) <- 0
    expect_true(constant(y) == 0, info=x)

  

  return(TRUE)
}  # checker1() closes


checker2 <- function(x,y){
  expect_true(x == -y+x+y, info=list(x,y))
  expect_true(x+y == x-(-y), info=list(x,y))

  expect_true(x+y == y+x, info=list(x,y))

  expect_true((-x)*y == -(x*y), info=list(x,y))
  expect_true(x*(-y) == -(x*y), info=list(x,y))

  ##  expect_true(x*y == y*x)  
  return(TRUE)
}

checker3 <- function(x,y,z){
  expect_true(x+(y+z) == (x+y)+z, info=list(x,y,z)) # additive associativity
  expect_true(x*(y*z) == (x*y)*z, info=list(x,y,z)) # multiplicative associativity

  expect_true(x*(y+z) == x*y + x*z, info=list(x,y,z))  # left distributivity
  expect_true((y+z)*x == y*x + z*x, info=list(x,y,z))  # right distributivity
  
  return(TRUE)
} # checker3() closes


for(i in 1:2){
  for(inc in c(TRUE,FALSE)){
    x <- rfalg(5,include.negative=inc)
    y <- rfalg(5,include.negative=inc)
    z <- rfalg(5,include.negative=inc)
    
    checker1(x)
    checker2(x,y)
    checker3(x,y,z)
  }
}

p0 <- as.freealg("0")
p1 <- as.freealg("1")
p2 <- as.freealg("1+2x+3Xy")
p3 <- as.freealg("1+2x+3Xy-4YYYYxxdxd")
p4 <- as.freealg("1+2x+3Xy+4YYYYxxdxd")

checker1(p0)
checker1(p1)
checker1(p2)
checker1(p3)
checker1(p4)

checker2(p1,p2)
checker2(p2,p3)
checker2(p3,p4)
checker2(p4,p1)

checker3(p1,p1,p1)
checker3(p1,p2,p3)
checker3(p2,p3,p4)
checker3(p3,p4,p1)
checker3(p4,p1,p1)

})
