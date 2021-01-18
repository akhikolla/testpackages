test_that("decoding works", {
  expect_identical(
    decode_ids(c(
      "32013030-2d30-3033-3338-3733fa30c4fa",
      "NA",
      "00-0033873",
      "NA",
      "32013030-2d30-3032-3739-3434d4d3846d"
    )),
    c("00-0033873", NA_character_, "00-0033873", NA_character_, "00-0027944")
  )
})
