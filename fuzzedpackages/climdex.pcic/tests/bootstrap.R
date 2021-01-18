load("exemplar_data.rda")
library(PCICt)

if(require("RUnit", quietly=TRUE)) {
  ## Run all the tests
  testsuite <- defineTestSuite("climdex.pcic", dirs=".", testFileRegexp = "^test_.+.R$", testFuncRegexp = "^climdex.pcic.test.+")
  climdex.pcic.test.result <- runTestSuite(testsuite, useOwnErrorHandler=F)
  printTextProtocol(climdex.pcic.test.result)
  stopifnot(climdex.pcic.test.result$climdex.pcic$nFail == 0 && climdex.pcic.test.result$climdex.pcic$nErr == 0)    
}
