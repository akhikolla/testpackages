
test_that("check_parameters_RBVS-SS_g", {
  # skip_on_cran()
  n = 25; env = 2; s = 3; size = env+1
  X = scale(matrix(rnorm(n*s,0,5), n, s), scale=FALSE)
  E = scale(matrix(rnorm(n*env,0,1), n, env), scale=FALSE)
  Y= 15 + rnorm(n)
  fit=roben(X, Y, E, structure="g")

  expect_equal(fit$iterations, 10000)
  expect_equal(fit$burn.in, 5000)
  expect_equal(dim(fit$posterior$GS.alpha), c(10000, env+1))
  expect_equal(dim(fit$posterior$GS.beta), c(10000, s*size))
  # expect_equal(class(fit), c("roben", "Sparse", "RBVS"))
  expect_equal(length(fit$coefficient$Int), 1)
  expect_equal(length(fit$coefficient$clin), 0)
  expect_equal(length(fit$coefficient$E), env)
  expect_equal(dim(fit$coefficient$GE), c(size, s))
  expect_gt(fit$coefficient$Int, 0)
  expect_equal(sum(fit$coefficient$GE==0), size*s)

})

test_that("check_parameters_RBVS-SS_sg", {
  # skip_on_cran()
  n = 25; env = 2; s = 3; size = env+1
  X = scale(matrix(rnorm(n*s,0,5), n, s), scale=FALSE)
  E = scale(matrix(rnorm(n*env,0,1), n, env), scale=FALSE)
  Y= X%*%c(1, 2, 3)+rnorm(n)
  fit=roben(X, Y, E, iterations=5000)

  expect_equal(fit$iterations, 5000)
  expect_equal(fit$burn.in, 2500)
  expect_equal(dim(fit$posterior$GS.alpha), c(5000, env+1))
  expect_equal(dim(fit$posterior$GS.beta), c(5000, s*size))
  # expect_equal(class(fit), c("roben", "Sparse", "RBVS"))
  expect_equal(length(fit$coefficient$Int), 1)
  expect_equal(length(fit$coefficient$clin), 0)
  expect_equal(length(fit$coefficient$E), env)
  expect_equal(dim(fit$coefficient$GE), c(size, s))
  expect_gt(sum(fit$coefficient$GE==0), s)
  expect_equal(sum(fit$coefficient$GE[1,]!=0), s)

})

test_that("check_parameters_RBVS-SS_i", {
  # skip_on_cran()
  n = 25; env = 2; s = 3; size = env+1
  X = scale(matrix(rnorm(n*s,0,2), n, s), scale=FALSE);
  E = scale(matrix(rnorm(n*env,0,5), n, env), scale=FALSE)
  Y= X[,2]*E[,1]*2+rnorm(n)
  fit=roben(X, Y, E, iterations=5000, structure="i")

  expect_equal(dim(fit$posterior$GS.alpha), c(5000, env+1))
  expect_equal(dim(fit$posterior$GS.beta), c(5000, s*size))
  # expect_equal(class(fit), c("roben", "Sparse", "RBVS"))
  expect_named(fit$coefficient$E)
  expect_equal(dim(fit$coefficient$GE), c(size, s))
  expect_gt(sum(fit$coefficient$GE==0), s)
  expect_true(fit$coefficient$GE[2,2]!=0)

})

test_that("check_parameters_RBVS_sg", {
  # skip_on_cran()
  n = 25; env = 2; s = 4; size = env+1
  X = scale(matrix(runif(n*s,0,100), n, s), scale=FALSE);
  E = scale(matrix(runif(n*env,0,20), n, env), scale=FALSE);
  Y= 15 + rnorm(n)
  fit=roben(X, Y, E, iterations=5000, sparse = FALSE, structure="s")

  expect_equal(dim(fit$posterior$GS.alpha), c(5000, env+1))
  expect_equal(dim(fit$posterior$GS.beta), c(5000, s*size))
  # expect_equal(class(fit), c("roben","RBVS"))
  expect_length(fit$coefficient$Int, 1)
  expect_length(fit$coefficient$clin, 0)
  expect_length(fit$coefficient$E, env)
  expect_equal(dim(fit$coefficient$GE), c(size, s))
  expect_equal(sum(fit$coefficient$GE==0), 0)

})

test_that("check_parameters_RBVS_g", {
  # skip_on_cran()
  n = 25; env = 2; s = 4; size = env+1
  X = scale(matrix(rnorm(n*s,0,1), n, s), scale=FALSE);
  E = scale(matrix(rnorm(n*env,0,1), n, env), scale=FALSE);
  Y= 15 + rnorm(n)
  colnames(X) = c("A", "B", "C", "D")
  colnames(E) = c("act", "gl")
  fit=roben(X, Y, E, iterations=5000, sparse = FALSE, structure="g")

  # expect_equal(class(fit), c("roben","RBVS"))
  expect_equal(names(fit$coefficient$E), c("act", "gl"))
  expect_equal(colnames(fit$coefficient$GE), c("A", "B", "C", "D"))
  expect_equal(dim(fit$coefficient$GE), c(size, s))
  expect_equal(sum(fit$coefficient$GE==0), 0)

})

test_that("check_parameters_RBVS_i", {
  # skip_on_cran()
  n = 25; env = 0; s = 4; size = env+1
  X = scale(matrix(rnorm(n*s,0,2), n, s), scale=FALSE);
  Y= 5+ X%*%c(1, 2, 3, 4)+rnorm(n)
  colnames(X) = c("A", "B", "C", "D")
  fit=roben(X, Y, E=NULL, iterations=3000, sparse = FALSE, structure="i", debugging=TRUE)

  # expect_equal(class(fit), c("roben","RBVS"))
  expect_length(fit$coefficient$E, env)
  expect_equal(dim(fit$coefficient$GE), c(size, s))
  expect_equal(colnames(fit$coefficient$GE), c("A", "B", "C", "D"))
  expect_gt(fit$coefficient$Int, 0)
  expect_equal(sum(fit$coefficient$GE[1,]!=0), s)

})

######### BVS

test_that("check_parameters_BVS-SS_sg", {
  # skip_on_cran()
  n = 25; env = 2; s = 4; nclin=2; size = env+1
  X = scale(matrix(rnorm(n*s,0,3), n, s), scale=FALSE)
  E = scale(matrix(rnorm(n*env,0,2), n, env), scale=FALSE)
  clin = scale(matrix(runif(n*nclin,-5,5), n, nclin), scale=FALSE)
  Y= X%*%c(1, 2, 3, 4)+rnorm(n)
  fit=roben(X, Y, E, clin, iterations=8000, robust=FALSE, structure="s")

  expect_equal(fit$iterations, 8000)
  expect_equal(fit$burn.in, 4000)
  expect_equal(dim(fit$posterior$GS.alpha), c(8000, env+nclin+1))
  expect_equal(dim(fit$posterior$GS.beta), c(8000, s*size))
  # expect_equal(class(fit), c("roben","BVS-SS"))
  expect_equal(length(fit$coefficient$Int), 1)
  expect_equal(length(fit$coefficient$clin), nclin)
  expect_equal(length(fit$coefficient$E), env)
  expect_equal(dim(fit$coefficient$GE), c(size, s))
  expect_gt(sum(fit$coefficient$GE==0), s)
  expect_equal(sum(fit$coefficient$GE[1,]!=0), s)

})


test_that("check_parameters_BVS-SS_g", {
  # skip_on_cran()
  n = 25; env = 2; s = 4; nclin=2; size = env+1
  X = scale(matrix(runif(n*s,0,10), n, s), scale=FALSE)
  E = scale(matrix(runif(n*env,0,5), n, env), scale=FALSE)
  clin = scale(matrix(floor(runif(n*nclin,-5,5)), n, nclin), scale=FALSE)
  Y= 15 + rnorm(n)
  fit=roben(X, Y, E, clin, iterations=8000, robust=FALSE, structure="g")

  expect_named(fit$coefficient$E)
  expect_named(fit$coefficient$clin)
  # expect_equal(class(fit), c("roben","BVS-SS"))
  expect_gt(fit$coefficient$Int, 0)
  expect_equal(sum(fit$coefficient$GE==0), size*s)

})

test_that("check_parameters_BVS-SS_i", {
  # skip_on_cran()
  n = 25; env = 2; s = 4; nclin=2; size = env+1
  X = scale(matrix(runif(n*s,0,10), n, s), scale=FALSE)
  E = scale(matrix(runif(n*env,0,5), n, env), scale=FALSE)
  clin = scale(matrix(floor(runif(n*nclin,-5,5)), n, nclin), scale=FALSE)
  Y= X[,3]*E[,2]*2+rnorm(n)
  colnames(X) = c("A", "B", "C", "D")
  colnames(E) = c("act", "gl")
  colnames(clin) = c("age", "wt")
  fit=roben(X, Y, E, clin, iterations=8000, robust=FALSE, structure="i")

  expect_equal(names(fit$coefficient$E), c("act", "gl"))
  expect_equal(names(fit$coefficient$clin), c("age", "wt"))
  expect_equal(colnames(fit$coefficient$GE), c("A", "B", "C", "D"))
  # expect_equal(class(fit), c("roben","BVS-SS"))
  expect_gt(sum(fit$coefficient$GE==0), s)
  expect_true(fit$coefficient$GE[3,3]!=0)

})

test_that("check_parameters_BVS_sg", {
  # skip_on_cran()
  n = 25; env = 3; s = 2; nclin=2; size = env+1
  X = scale(matrix(runif(n*s,0,5), n, s), scale=FALSE)
  E = scale(matrix(runif(n*env,0,5), n, env), scale=FALSE)
  clin = scale(matrix(floor(runif(n*nclin,-5,5)), n, nclin), scale=FALSE)
  Y= 15 + rnorm(n)
  fit=roben(X, Y, E, clin, iterations=5000, robust=FALSE, sparse = FALSE)

  expect_equal(dim(fit$posterior$GS.alpha), c(5000, env+nclin+1))
  expect_equal(dim(fit$posterior$GS.beta), c(5000, s*size))
  # expect_equal(class(fit), c("roben","BVS"))
  expect_equal(dim(fit$coefficient$GE), c(size, s))
  expect_equal(sum(fit$coefficient$GE==0), 0)

})

test_that("check_parameters_BVS_g", {
  # skip_on_cran()
  n = 25; env = 0; s = 4; nclin=2; size = env+1
  X = scale(matrix(runif(n*s,0,5), n, s), scale=FALSE)
  clin = scale(matrix(floor(runif(n*nclin,-5,5)), n, nclin), scale=FALSE)
  Y= 5+ X%*%c(1, 2, 3, 4)+rnorm(n)
  fit=roben(X, Y, E=NULL, clin, iterations=3000, robust=FALSE, sparse = FALSE, structure="g", debugging = TRUE)

  expect_equal(dim(fit$posterior$GS.alpha), c(3000, env+nclin+1))
  expect_equal(dim(fit$posterior$GS.beta), c(3000, s*size))
  # expect_equal(class(fit), c("roben","BVS"))
  expect_length(fit$coefficient$clin, nclin)
  expect_equal(dim(fit$coefficient$GE), c(size, s))
  expect_gt(fit$coefficient$Int, 0)
})
