test_that("decoding works", {
  infile="testdata/mesh.draco"
  expect_is(mesh3d <- draco_decode(infile, mesh3d = T), 'mesh3d')
  expect_known_value(mesh3d, 'testdata/mesh3d.rds')

  expect_known_value(raw <- draco_decode(infile, mesh3d = F),
                     'testdata/meshraw.rds')
  # check that this matches our previous raw format
  raw$faces=raw$faces+1
  expect_known_value(raw, 'testdata/mesh.rds')

  # error when we don't get draco data
  expect_error(draco_decode('testdata/mesh.ply'), 'Bad input')

  # trunc=readBin('tests/testthat/testdata/mesh.draco', what=raw(), n = 1000)
  # writeBin(trunc, 'mesh.draco.truncated')
  expect_error(draco_decode('testdata/mesh.draco.truncated'),
               "Unable to decode triangular mesh")
})

test_that("decoding URL works", {
  skip_if_offline()
  carurl='https://github.com/google/draco/blob/master/testdata/car.drc?raw=true'
  expect_is(car <- draco_decode(carurl, quiet=TRUE),'mesh3d')
})

# (base) Gregs-MBP-2:build jefferis$ ./draco_encoder -i /Users/jefferis/dev/R/dracor/tests/testthat/testdata/mesh.ply -o /Users/jefferis/dev/R/dracor/tests/testthat/testdata/mesh.ply.draco
# Encoder options:
#   Compression level = 7
# Positions: Quantization = 11 bits
#
# Encoded mesh saved to /Users/jefferis/dev/R/dracor/tests/testthat/testdata/mesh.ply.draco (7 ms to encode).
#
# Encoded size = 6074 bytes
#
# For better compression, increase the compression level up to '-cl 10' .
#
# (base) Gregs-MBP-2:build jefferis$ ./draco_encoder -i /Users/jefferis/dev/R/dracor/tests/testthat/testdata/mesh.obj -o /Users/jefferis/dev/R/dracor/tests/testthat/testdata/mesh.obj.draco
# Encoder options:
#   Compression level = 7
# Positions: Quantization = 11 bits
#
# Encoded mesh saved to /Users/jefferis/dev/R/dracor/tests/testthat/testdata/mesh.obj.draco (2 ms to encode).
#
# Encoded size = 6074 bytes
#
# For better compression, increase the compression level up to '-cl 10' .
test_that("draco round trip works", {
  origfile="testdata/mesh.draco"
  rtply="testdata/mesh.ply.draco"
})
