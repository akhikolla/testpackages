# -----------------------------------------------------------------------------
# Test if phasePortrait produces the correct numerical output.
#
# Uses reference files created with the R-script in CreateTestCases/.
# These reference files are in the project's subdirectory tests/testthat/,
# which automatically becomes the working directory when test_that() is run.
# -----------------------------------------------------------------------------


# function cleanUp
# Deletes all files from the R session's temporary directory which match the
# naming conventions of phasePortrait's temporary files. Is called before and
# after each test run. Thus, the R session's temporary directory, is always
# absolutely clean and unambiguous as far as phasePortrait output is concerned.
cleanUp <- function() {
  existZWMat <- dir(tempdir(), pattern = "^\\d+(z|w)mat\\d{10}\\.RData$")
  unlink(paste0(tempdir(), "/", existZWMat))
}


# function loadActual
# After a call to cleanup and a subsequent call of phasePortrait, only output
# files from the last phasePortrait run can possibly reside in the R session's
# temporary directory (if we explicitly told phasePortrait not to delete them
# after use). The test cases are so slim, that there is only one zmat and one
# wmat file. We load the "wmat" file for further inspection.
loadActual <- function() {
  existWMat <- dir(tempdir(), pattern = "^\\d+wmat\\d{10}\\.RData$")
  get(load(paste0(tempdir(), "/", existWMat)))
}


# Test cases defined as functions
# -------------------------------

# Test cases for phasePortrait


# Test case 1:
# A rational function given as single string
testCase1 <- function() {
  cleanUp()
  phasePortrait("(2-z)^2*(-1i+z)^3*(4-3i-z)/((2+2i+z)^4)",
                xlim = c(-4, 4), ylim = c(-4, 4),
                invertFlip = FALSE,
                blockSizePx = 2250000,
                res = 150,
                tempDir = NULL,
                deleteTempFiles = FALSE,
                noScreenDevice = TRUE,
                nCores = 2,
                verbose = FALSE)
  referenceWmat <- get(load("1wmatCase001.RData"))
  actualWmat    <- loadActual()
  cleanUp()
  rslt          <- all.equal(referenceWmat, actualWmat)
  rm(referenceWmat, actualWmat)
  return(rslt)
}


# Test case 2:
# A rational function given as single string, but with invertFlip = TRUE
testCase2 <- function() {
  cleanUp()
  phasePortrait("(2-z)^2*(-1i+z)^3*(4-3i-z)/((2+2i+z)^4)",
                xlim = c(-4, 4), ylim = c(-4, 4),
                invertFlip = TRUE,
                blockSizePx = 2250000,
                res = 150,
                tempDir = NULL,
                deleteTempFiles = FALSE,
                noScreenDevice = TRUE,
                nCores = 2,
                verbose = FALSE)
  referenceWmat <- get(load("1wmatCase002.RData"))
  actualWmat    <- loadActual()
  cleanUp()
  rslt          <- all.equal(referenceWmat, actualWmat)
  rm(referenceWmat, actualWmat)
  return(rslt)
}


# Test case 3
# User function with additional default arguments which are _not_ specified
# in the call to phasePortrait
testCase3 <- function() {

  jacobiTheta_1 <- function(z, tau = 1i, nIter = 30) {
    k <- c(1:nIter)
    q <- exp(pi*1i*tau)
    g <- exp(2*pi*1i*z)
    return(1 + sum(q^(k^2)*g^k + q^(k^2)*(1/g)^k))
  }

  cleanUp()
  phasePortrait(jacobiTheta_1,
                xlim = c(-2, 2), ylim = c(-2, 2),
                invertFlip = FALSE,
                blockSizePx = 2250000,
                res = 150,
                tempDir = NULL,
                deleteTempFiles = FALSE,
                noScreenDevice = TRUE,
                nCores = 2,
                verbose = FALSE)
  referenceWmat <- get(load("1wmatCase003.RData"))
  actualWmat    <- loadActual()
  cleanUp()
  rslt          <- all.equal(referenceWmat, actualWmat)
  rm(referenceWmat, actualWmat)
  return(rslt)
}


# Test case 4
# User function with additional default arguments which are specified
# in the call to phasePortrait
testCase4 <- function() {

  jacobiTheta_1 <- function(z, tau = 1i, nIter = 30) {
    k <- c(1:nIter)
    q <- exp(pi*1i*tau)
    g <- exp(2*pi*1i*z)
    return(1 + sum(q^(k^2)*g^k + q^(k^2)*(1/g)^k))
  }

  cleanUp()
  phasePortrait(jacobiTheta_1,
                moreArgs = list(tau = 1i/2 - 1/4, nIter = 30),
                xlim = c(-2, 2), ylim = c(-2, 2),
                invertFlip = FALSE,
                blockSizePx = 2250000,
                res = 150,
                tempDir = NULL,
                deleteTempFiles = FALSE,
                noScreenDevice = TRUE,
                nCores = 2,
                verbose = FALSE,
                autoDereg = TRUE) # Register sequential backend after
                                  # the last phase portrait
  referenceWmat <- get(load("1wmatCase004.RData"))
  actualWmat    <- loadActual()
  cleanUp()
  rslt          <- all.equal(referenceWmat, actualWmat)
  rm(referenceWmat, actualWmat)
  return(rslt)
}


# Test cases for phasePortrait

# Test case 5
# Same as case 1, but with phasePortraitBw
testCase5 <- function() {
  cleanUp()
  phasePortraitBw("(2-z)^2*(-1i+z)^3*(4-3i-z)/((2+2i+z)^4)",
                  xlim = c(-4, 4), ylim = c(-4, 4),
                  invertFlip = FALSE,
                  blockSizePx = 2250000,
                  res = 150,
                  tempDir = NULL,
                  deleteTempFiles = FALSE,
                  noScreenDevice = TRUE,
                  nCores = 2,
                  verbose = FALSE)
  referenceWmat <- get(load("1wmatCase005.RData"))
  actualWmat    <- loadActual()
  cleanUp()
  rslt          <- all.equal(referenceWmat, actualWmat)
  rm(referenceWmat, actualWmat)
  return(rslt)
}


# Test case 6:
# A rational function given as single string, but with invertFlip = TRUE
# Same as case 2, but with phasePortraitBw
testCase6 <- function() {
  cleanUp()
  phasePortraitBw("(2-z)^2*(-1i+z)^3*(4-3i-z)/((2+2i+z)^4)",
                  xlim = c(-4, 4), ylim = c(-4, 4),
                  invertFlip = TRUE,
                  blockSizePx = 2250000,
                  res = 150,
                  tempDir = NULL,
                  deleteTempFiles = FALSE,
                  noScreenDevice = TRUE,
                  nCores = 2,
                  verbose = FALSE)
  referenceWmat <- get(load("1wmatCase006.RData"))
  actualWmat    <- loadActual()
  cleanUp()
  rslt          <- all.equal(referenceWmat, actualWmat)
  rm(referenceWmat, actualWmat)
  return(rslt)
}


# Test case 7
# User function with additional default arguments which are _not_ specified
# in the call to phasePortraitBw
# Same as case 3, but with phasePortraitBw
testCase7 <- function() {

  jacobiTheta_1 <- function(z, tau = 1i, nIter = 30) {
    k <- c(1:nIter)
    q <- exp(pi*1i*tau)
    g <- exp(2*pi*1i*z)
    return(1 + sum(q^(k^2)*g^k + q^(k^2)*(1/g)^k))
  }

  cleanUp()
  phasePortraitBw(jacobiTheta_1,
                  xlim = c(-2, 2), ylim = c(-2, 2),
                  invertFlip = FALSE,
                  blockSizePx = 2250000,
                  res = 150,
                  tempDir = NULL,
                  deleteTempFiles = FALSE,
                  noScreenDevice = TRUE,
                  nCores = 2,
                  verbose = FALSE)
  referenceWmat <- get(load("1wmatCase007.RData"))
  actualWmat    <- loadActual()
  cleanUp()
  rslt          <- all.equal(referenceWmat, actualWmat)
  rm(referenceWmat, actualWmat)
  return(rslt)
}


# Test case 8
# User function with additional default arguments which are specified
# in the call to phasePortraitBw
# Same as case 4, but with phasePortraitBw
testCase8 <- function() {

  jacobiTheta_1 <- function(z, tau = 1i, nIter = 30) {
    k <- c(1:nIter)
    q <- exp(pi*1i*tau)
    g <- exp(2*pi*1i*z)
    return(1 + sum(q^(k^2)*g^k + q^(k^2)*(1/g)^k))
  }

  cleanUp()
  phasePortrait(jacobiTheta_1,
                moreArgs = list(tau = 1i/2 - 1/4, nIter = 30),
                xlim = c(-2, 2), ylim = c(-2, 2),
                invertFlip = FALSE,
                blockSizePx = 2250000,
                res = 150,
                tempDir = NULL,
                deleteTempFiles = FALSE,
                noScreenDevice = TRUE,
                nCores = 2,
                verbose = FALSE,
                autoDereg = TRUE) # Register sequential backend after
  # the last phase portrait
  referenceWmat <- get(load("1wmatCase008.RData"))
  actualWmat    <- loadActual()
  cleanUp()
  rslt          <- all.equal(referenceWmat, actualWmat)
  rm(referenceWmat, actualWmat)
  return(rslt)
}




# The actual tests
test_that("phasePortrait produces correct numerical output", {
  expect_true(testCase1())
  expect_true(testCase2())
  expect_true(testCase3())
  expect_true(testCase4())
})

test_that("phasePortraitBw produces correct numerical output", {
  expect_true(testCase5())
  expect_true(testCase6())
  expect_true(testCase7())
  expect_true(testCase8())
})


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


