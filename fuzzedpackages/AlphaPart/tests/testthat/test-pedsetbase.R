context("test-pedsetbase")

test_that("Test whether base population is set properly", {

  ped <- data.frame(      id=1:10,
                         fid=c(0, 0, 0, 1, 1, 1, 3, 3, 3, 5),
                         mid=c(0, 0, 0, 2, 0, 2, 2, 2, 5, 0),
                    birth_dt=c(0, 0, 1, 2, 3, 3, 3, 4, 4, 5) + 2000)

  ped1  <- pedSetBase(x=ped, keep=ped$birth_dt > 2002, unknown=0, report=FALSE)


  expect_equal(nrow(ped1), 6)
  expect_equal(ped1$id,  5:10)
  expect_equal(ped1$fid, c(rep(0, times=5), 5))
  expect_equal(ped1$mid, c(rep(0, times=4), 5, 0))

})


test_that("Test whether base population is set properly with NAs for missing", {
    ped <- data.frame(      id=1:10,
                     fid=c(0, 0, 0, 1, 1, 1, 3, 3, 3, 5),
                     mid=c(0, 0, 0, 2, 0, 2, 2, 2, 5, 0),
                birth_dt=c(0, 0, 1, 2, 3, 3, 3, 4, 4, 5) + 2000)
    ped   <- unknownToNA(x=ped, unknown=0)
    ped2  <- pedSetBase(x=ped, keep=ped$birth_dt > 2002,            report=FALSE)

    expect_equal(ped2$fid, c(rep(NA, times=5), 5))
    expect_equal(ped2$mid, c(rep(NA, times=4), 5, NA))
})
