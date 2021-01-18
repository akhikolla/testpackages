context("test_auxilar_code.R")


test_that("online site works", {
  expect_that(chunkR_devel(), prints_text("Opening link"))
})

test_that("get_matrix2dataframe return warning", {
  expect_that(get_matrix2dataframe(), gives_warning("deprecated in chunkR"))
})

test_that("get_matrix2dataframe returns warning", {
  expect_that(get_matrix2dataframe(), gives_warning("deprecated in chunkR"))
})


test_that("get_dataframe returns warning", {
  expect_that(get_dataframe(), gives_warning("deprecated in chunkR"))
})

test_that("get_matrix returns warning", {
  expect_that(get_matrix(), gives_warning("deprecated in chunkR"))
})


test_that("reader returns warning", {
  expect_that(reader(), gives_warning("deprecated in chunkR"))
})

