test_that("check dimension of input", {
  # skip_on_cran()
  n = 25; X = matrix(floor(runif(n*3,0,100)), n, 3);
  Y= matrix(floor(runif(n-1,0,10)), n-1, 1)
  expect_error(Data.matrix(X, Y, E=NULL), "Length of Y does not match", ignore.case = TRUE)

  Y= matrix(floor(runif(n,0,10)), n, 1)
  expect_error(Data.matrix(X, Y, E=NULL, clin=1), "clin has a different number of rows.", ignore.case = TRUE)
  expect_error(Data.matrix(X, Y, E=NULL), "E factors must be provided.", ignore.case = TRUE)
  expect_error(Data.matrix(X, Y, E=1), "E has a different number of rows", ignore.case = TRUE)

})


test_that("check design matrix", {
  # skip_on_cran()
  n = 25; env = 2; size = env+1
  X = scale(matrix(floor(runif(n*3,0,100)), n, 3), scale=FALSE)
  E = scale(matrix(floor(runif(n*env,0,20)), n, env), scale=FALSE)
  Y= matrix(floor(runif(n,0,10)), n, 1)
  dat=Data.matrix(X, Y, E)

  expect_equal(dim(dat$xx), c(n, 3*size))
  expect_equal(dat$xx[,1], X[,1])
  expect_equal(dat$xx[,2], X[,1]*E[,1])
  expect_equal(dat$xx[,3], X[,1]*E[,2])
  expect_equal(dat$xx[,size+1], X[,2])
  expect_equal(dat$xx[,size+2], X[,2]*E[,1])
})

