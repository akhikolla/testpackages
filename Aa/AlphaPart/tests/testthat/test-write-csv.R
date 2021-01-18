context("test-write-csv")

test_that("Check writing input for write.csv.AlphaPart", {

  ## Making sure we accept only the right class and one file name!
  expect_error(write.csv.AlphaPart(data.frame(x=1:10)))
  expect_error(write.csv.AlphaPart(data.frame(x=1:10), file=c("a", "b")))
})

test_that("Check writing process for write.csv.AlphaPart", {

  ## Partition additive genetic values
  res <- AlphaPart(x=AlphaPart.ped, colPath="country", colBV="bv1")

  ## Write summary on the disk and collect saved file names
  dirT <- tempdir()
  fileName <- file.path(dirT, "AlphaPart")
  retF <- write.csv(x=res, file=fileName)


  ## Check content of files
  tmp <- read.csv2(file=retF[1])
  expect_equal(tmp$agv1,   res$agv1$agv1)
  expect_equal(tmp$agv1_1, res$agv1$agv1_1)

  ## Clean up
  files <- dir(path=dirT, pattern="AlphaPart*")
  unlink(x=files)
})

###############################################################
###############################################################
###############################################################
###############################################################
#write.csv.summaryAlphaPart

test_that("Check writing input for write.csv.summaryAlphaPart", {

  ## Making sure we accept only the right class and one file name!
  expect_error(write.csv.summaryAlphaPart(data.frame(x=1:10)))
  expect_error(write.csv.summaryAlphaPart(data.frame(x=1:10), file=c("a", "b")))
})

test_that("Check writing process for write.csv.summaryAlphaPart", {

  ## Partition additive genetic values
  res <- AlphaPart(x=AlphaPart.ped, colPath="country", colBV=c("bv1"))

  ## Summarize population by generation (=trend)
  ret <- summary(res, by="gen")

  ## Write summary on the disk and collect saved file names
  dirT <- tempdir()
  fileName <- file.path(dirT, "AlphaPart")
  retF <- write.csv(x=ret, file=fileName)


  ## Check content of files
  col <- c("gen", "N", "Sum", "1", "2")
  tmp <- read.csv2(file=retF[1])
  expect_equal(tmp, ret$bv1)

  ## Clean up
  files <- dir(path=dirT, pattern="AlphaPart*")
  unlink(x=files)
})

