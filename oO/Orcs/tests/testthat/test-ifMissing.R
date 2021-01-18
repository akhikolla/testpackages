context("ifMissing")

## raster files
logo = system.file("external/rlogo.grd", package = "raster")

ofl = file.path(tmp <- tempdir(), "rlogo.tif")
if (file.exists(ofl)) jnk = file.remove(ofl)

test_that("nonexisting raster file is created", {
  
  s = ifMissing(logo)
  expect_is(s, "RasterBrick")

  s2 = ifMissing(ofl, arg1 = "filename", x = s, datatype = "INT1U")  
  expect_is(s, "RasterBrick")
  
  s3 = ifMissing(ofl)
  expect_is(s, "RasterBrick")
  
  jnk = file.remove(ofl)
})

## text files
data(iris)

ofl = file.path(tmp, "iris.csv")
if (file.exists(ofl)) jnk = file.remove(ofl)

test_that("nonexisting text file is created", {
  
  jnk = ifMissing(ofl, fun1 = write.csv, x = iris, file = ofl, row.names = FALSE)
  expect_null(jnk)
  
  dat = ifMissing(ofl, fun0 = function(x) read.csv(x, stringsAsFactors = TRUE))
  expect_equal(dat, iris)
  
  jnk = file.remove(ofl)
  expect_true(jnk)
  
  fun = function(x, file = "", ...) {
    write.csv(x, file, ...)
    read.csv(file, stringsAsFactors = TRUE)
  }
  dat = ifMissing(ofl, fun1 = fun, arg1 = "file", x = iris, quote = FALSE, row.names = FALSE)
  expect_equal(dat, iris)
  
  jnk = file.remove(ofl)
})
