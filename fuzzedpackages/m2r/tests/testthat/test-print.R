context("print methods")










test_that("rings print properly", {

  r <- structure(NA,
    class = c("m2_polynomialring", "m2"),
    m2_name = "m2rintring00000000",
    m2_meta = list(
      vars = list("x", "y"),
      coefring = "CC",
      order = "grevlex"
    )
  )

  # check that print returns the input object
  expect_identical(r, print(r))

  # check that the printed format is correct
  expect_equal(
    capture.output(print(r)),
    "M2 Ring: CC[x,y], grevlex order"
  )

})
