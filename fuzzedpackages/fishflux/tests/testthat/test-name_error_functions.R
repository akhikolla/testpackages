test_that("simple corner case", {
  needs_api()
  check_api()
  #check_name_fishbase
  expect_error(fishflux::check_name_fishbase("name"))
  expect_null(fishflux::check_name_fishbase("Zebrasoma scopas"))
  #name_errors
  expect_message(fishflux::name_errors("name"))
  expect_message(fishflux::name_errors("name"), "Inaccurate species names found:")
  expect_message(expect_equal(fishflux::name_errors("name"), "name"))
  expect_message(expect_null(fishflux::name_errors("Zebrasoma scopas")))
  expect_message(fishflux::name_errors("Zebrasoma scopas"))
  expect_message(fishflux::name_errors("Zebrasoma scopas"), "All species names are correct")
})
