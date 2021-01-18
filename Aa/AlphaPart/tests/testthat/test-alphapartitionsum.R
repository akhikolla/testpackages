context("test-alphapartsum")

test_that("Test input for AlphaPartSum", {
  ## Small pedigree with additive genetic (=breeding) values
  ped <- data.frame(  id=c(  1,   2,   3,   4,   5,   6),
                      fid=c(  0,   0,   2,   0,   4,   0),
                      mid=c(  0,   0,   1,   0,   3,   3),
                      loc=c("A", "B", "A", "B", "A", "A"),
                      gen=c(  1,   1,   2,   2,   3,   3),
                      trt1=c(100, 120, 115, 130, 125, 125),
                      trt2=c(100, 110, 105, 100,  85, 110))

  ## Partition additive genetic values
  tmp <- AlphaPart(x=ped, colBV=c("trt1", "trt2"))

  ## Test that we only accept objects of class AlphaPart and summaryAlphaPart
  expect_error(AlphaPartSum(x=ped))
  ## ... and that map is a list
  expect_error(AlphaPartSum(x=tmp, map="A"))
  ## ... defined paths exist when zeroPath=FALSE
  expect_error(AlphaPartSum(x=tmp, map=list(c("X", "A", "B", "Z")), zeroPath=FALSE))
  expect_error(AlphaPartSum(x=tmp, map=list(c("Z", "A", "B", "Z")), zeroPath=FALSE))
  expect_error(AlphaPartSum(x=tmp, map=list(c("X", "A", "B"), c("Y", "A", "Y")), zeroPath=FALSE))
})

test_that("Test output of AlphaPartSum - summing by paths", {
  ped <- data.frame(  id=c(  1,   2,   3,   4,   5,   6),
                      fid=c(  0,   0,   2,   0,   4,   0),
                      mid=c(  0,   0,   1,   0,   3,   3),
                      loc=c("A", "B", "A", "B", "A", "A"),
                      gen=c(  1,   1,   2,   2,   3,   3),
                      trt1=c(100, 120, 115, 130, 125, 125),
                      trt2=c(100, 110, 105, 100,  85, 110))

  ## Partition additive genetic values
  tmp <- AlphaPart(x=ped, colBV=c("trt1", "trt2"))
  ## Sum some partitions (working on object of class AlphaPart)
  tmp2 <- AlphaPartSum(x=tmp, map=list(c("X", "A", "B"), c("A", "B"), c("C", "X", "A")))

  ## Test that the creation of new path is done properly (is sum correct)
  expect_equal(tmp2$trt1$trt1_X, c(100, 120, 115, 130, 125, 125))
  expect_equal(tmp2$trt2$trt2_X, c(100, 110, 105, 100,  85, 110))
  expect_equal(tmp2$trt2$trt2_C, tmp2$trt2$trt2_X + tmp$trt2$trt2_B) ## A is first overwritten with B
  ## Test that overwritting existing path works (A with B in this case)
  expect_equal(tmp2$trt1$trt1_A, tmp$trt1$trt1_B)
  expect_equal(tmp2$trt2$trt2_A, tmp$trt2$trt2_B)
  ## Test that non target/new paths are removed (B in this case)
  expect_true(!("trt2_B" %in% colnames(tmp2$trt2)))
  ## Test that meta info slot is updated properly
  expect_equal(tmp2$info$nP, 3)
  expect_equal(tmp2$info$lP, c("A", "X", "C"))

  ## Test that unexisitng path are set to zero
  tmpE1 <- AlphaPartSum(x=tmp, map=list(c("X", "A", "B", "Z")))
  tmpE2 <- AlphaPartSum(x=tmp, map=list(c("X", "A", "B")))
  expect_equal(tmpE1, tmpE2)
  tmpE1 <- AlphaPartSum(x=tmp, map=list(c("Z", "A", "B", "Z")))
  tmpE2 <- AlphaPartSum(x=tmp, map=list(c("Z", "A", "B")))
  expect_equal(tmpE1, tmpE2)
})

test_that("Test the output of AlphaPartSum - summing by generation", {
  ped <- data.frame(  id=c(  1,   2,   3,   4,   5,   6),
                      fid=c(  0,   0,   2,   0,   4,   0),
                      mid=c(  0,   0,   1,   0,   3,   3),
                      loc=c("A", "B", "A", "B", "A", "A"),
                      gen=c(  1,   1,   2,   2,   3,   3),
                      trt1=c(100, 120, 115, 130, 125, 125),
                      trt2=c(100, 110, 105, 100,  85, 110))

  ## Partition additive genetic values
  tmp <- AlphaPart(x=ped, colBV=c("trt1", "trt2"))
  ## Summarize by generation
  tmpS <- summary(tmp, by="gen")

  ## Sum some partitions (working on object of class summaryAlphaPart)
  tmpS2 <- AlphaPartSum(x=tmpS, map=list(c("X", "A", "B"), c("A", "B")))

  ## Test that the creation of new path is is done properly (is sum correct)
  expect_equal(tmpS2$trt1$X, c(110, 122.5, 125))
  expect_equal(tmpS2$trt2$X, c(105, 102.5, 97.5))
  ## Test that overwritting existing path works (A with B in this case)
  expect_equal(tmpS2$trt1$A, tmpS$trt1$B)
  expect_equal(tmpS2$trt2$A, tmpS$trt2$B)
  ## Test that non target/new paths are removed (B in this case)
  expect_true(!("B" %in% colnames(tmpS2$trt2)))
  ## Test that meta info slot is updated properly
  expect_equal(tmpS2$info$nP, 2)
  expect_equal(tmpS2$info$lP, c("A", "X"))

  ## ... must be equal to
  tmpS3 <- summary(AlphaPartSum(x=tmp, map=list(c("X", "A", "B"), c("A", "B"))), by="gen")

  ## Test that we get the same output if we do AlphaPartSum(summary(AlphaPart())) or summary(AlphaPartSum(AlphaPart()))
  expect_equal(tmpS2, tmpS3)

  }
)
