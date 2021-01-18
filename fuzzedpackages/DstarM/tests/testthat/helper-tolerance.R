expect_equal2 <- function (object, expected, tol = 1e-5, ..., info = NULL,
                           label = NULL, expected.label = NULL) {

  expect_equal(
    object = object, expected = expected,
    tol = tol, ..., info = info, label = label,
    expected.label = expected.label
  )
}
