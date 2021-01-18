context("test-pedfixbirthyear")

test_that("Test pedFixBirthYear", {

  ped0 <- data.frame(     id=c( 1, 2, 3,  4, 5, 6, 7,  8, 9, 10, 11, 12, 13, 14),
                          fid=c( 0, 0, 0,  1, 1, 1, 3,  3, 3,  5,  4,  0,  0, 12),
                          mid=c( 0, 0, 0,  2, 0, 2, 2,  2, 5,  0,  0,  0,  0, 13),
                          birth_dt=c(NA, 0, 1, NA, 3, 3, 3, 3, 4, 4, 5, NA, 6, 6) + 2000)



  ped00 <- data.frame(     id=c( 1, 2, 3,  4, 5, 6, 7,  8, 9, 10, 11, 12, 13, 14),
                          fid=c( 0, 0, 0,  1, 1, 1, 3,  3, 3,  5,  0,  0,  0, 12),
                          mid=c( 0, 0, 0,  2, 0, 2, 2,  2, 5,  0,  0,  0,  0, 13),
                          birth_dt=c(NA, 0, 1, NA, 3, 3, 3, NA, 4, NA, NA, NA, NA, NA) + 2000)


  # The test fails on animals with no data and no offspring!!!
  ped1 <- pedFixBirthYear(x=ped0, interval=1,            report=FALSE)
  ped2 <- pedFixBirthYear(x=ped1, interval=1, down=TRUE, report=FALSE)
  ped3 <- pedFixBirthYear(x=ped2, interval=1,            report=FALSE)

  expect_equal(ped1[1, 4], 2002)
  expect_equal(ped2[c(4, 8, 10), 4], c(2004, 2003, 2004))
  expect_false(is.na(ped3[12, 4]))

  expect_warning(pedFixBirthYear(x=ped0[12:14, ], interval=1, report=FALSE, down=TRUE))
  ## two warnings:
  ##  - no information found for ...
  ##  - missing value for ...
  #

})
