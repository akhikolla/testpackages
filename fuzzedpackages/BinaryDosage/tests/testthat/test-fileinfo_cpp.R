test_that("badfilename", {
  expect_error(GetLineLocations("junk.txt"), "Unable to open file")
})
