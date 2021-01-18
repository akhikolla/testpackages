if (requireNamespace("RUnit", quietly=TRUE) && requireNamespace("RcppXPtrUtils", quietly=TRUE) && requireNamespace("r2sundials", quietly=TRUE) && requireNamespace("RcppArmadillo", quietly=TRUE)) {
   library(RUnit)
   library(RcppXPtrUtils)
   library(r2sundials)
   library(RcppArmadillo)

   testSuite <- defineTestSuite(
      name = "r2sundials unit tests",
      dirs = system.file("unitTests", package = "r2sundials"),
      testFuncRegexp = "^[Tt]est.+"
   )
   Sys.setenv("R_TESTS"="")
   tests <- runTestSuite(testSuite)

   printTextProtocol(tests)

   if (getErrors(tests)$nFail > 0) stop("RUnit test failure")
   if (getErrors(tests)$nErr > 0) stop("Errors in RUnit tests")
}
