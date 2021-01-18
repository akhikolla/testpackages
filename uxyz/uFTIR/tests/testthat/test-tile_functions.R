test_that("The example is loaded correctly", {
  x <- tile_read(system.file("extdata", "tile.bsp", package="uFTIR"))
  y <- tile_base_corr(x)
  y <- wavealign(y, primpke)
  y <- tile_sam(y)
  y <- get_profile_tile(y, x, 
                        dst_cluster = 5, 
                        clusternames = primpke@clusternames, 
                        plotpol = FALSE, plotpt = FALSE)
  expect_true(nrow(y) == 56)
})

test_that("Wavealign return a SpectralPack object", {
  x <- tile_read(system.file("extdata", "tile.bsp", package="uFTIR"))
  y <- wavealign(x, primpke)
  expect_s4_class(y, "SpectralPack")
})

test_that("Wavealign gives back data.x resampled as data.y", {
  x <- tile_read(system.file("extdata", "tile.bsp", package="uFTIR"))
  y <- wavealign(x, primpke)
  
  wn <- min(primpke@wavenumbers) < x@wavenumbers & x@wavenumbers < max(primpke@wavenumbers)
  wn[c(1, length(wn))] <- FALSE #because we interpolate we loose the borders
  
  expect_true(sum(y@Reference@wavenumbers == x@wavenumbers[wn]) == length(wn)-2)
})

test_that("Wavealign SpectralPack has slots with object having the same wavenumbers", {
  x <- tile_read(system.file("extdata", "tile.bsp", package="uFTIR"))
  x <- wavealign(x, primpke)
  expect_equal(x@Readings@wavenumbers, x@Reference@wavenumbers)
})

test_that("All derivative methods work", {
  x <- tile_read(system.file("extdata", "tile.bsp", package="uFTIR"))
  x <- tile_base_corr(x)
  x <- wavealign(x, primpke)
  expect_s4_class(tile_sam(x, derivative = NULL), "SAM")
  expect_s4_class(tile_sam(x, derivative = 1), "SAM")
  expect_s4_class(tile_sam(x, derivative = 2), "SAM")
})

test_that("Expect clipper to return S3 clipper", {
  x <- tile_read(system.file("extdata", "tile.bsp", package="uFTIR"))
  x <- tile_base_corr(x)
  x <- wavealign(x, primpke)
  x <- tile_sam(x)

  clip <- toClip(rad = 8, segments = 20, centre = c(10,10))
  expect_s3_class(clipper(x, clip@centre, clip@rad, slice = 1), "clipper")
})

test_that("Summary methods are equivalent for Smooth and SAM objects", {
  x <- tile_read(system.file("extdata", "tile.bsp", package="uFTIR"))
  x <- tile_base_corr(x)
  x <- wavealign(x, primpke)
  x <- tile_sam(x)
  
  out1 <- summary_sam(x, mask = NULL, 
                      smooth = TRUE, window = 3, 
                      slice = 1, clusternames = primpke@clusternames)
  x <- smooth_sam(x, as.integer(length(primpke@clusternames)), window = 3, nslices = 1)
  out2 <- summary_sam(x, clusternames = primpke@clusternames)
  expect_equal(out1, out2)
})

test_that("Summary methods are equivalent for Smooth, SAM and clipper objects", {
  x <- tile_read(system.file("extdata", "tile.bsp", package="uFTIR"))
  x <- tile_base_corr(x)
  x <- wavealign(x, primpke)
  x <- tile_sam(x)
  y <- smooth_sam(x, as.integer(length(primpke@clusternames)), window = 3, nslices = 1)
  clip <- toClip(rad = 8, segments = 20, centre = c(10,10))
  z <- clipper(y, centre = clip@centre, rad = clip@rad, slice = 1)
  
  out1 <- summary_sam(x, mask = clip, smooth = TRUE, 
                      window = 3, slice = 1, clusternames = primpke@clusternames, temporal = TRUE)
  out2 <- summary_sam(y, mask = clip, clusternames = primpke@clusternames, temporal = TRUE)
  out3 <- summary_sam(z, clusternames = primpke@clusternames, temporal = TRUE)
  
  expect_true(all.equal(out1, out2) & all.equal(out2, out3))
  
  out4 <- summary_sam(z, clusternames = NULL, temporal = TRUE)
  expect_equal(ncol(out4), 2)
})

test_that("get_profile can plot as side effect", {
  x <- tile_read(system.file("extdata", "tile.bsp", package="uFTIR"))
  y <- tile_base_corr(x)
  y <- wavealign(y, primpke)
  y <- tile_sam(y)
  y <- get_profile_tile(y, x, 
                        dst_cluster = 5, 
                        clusternames = primpke@clusternames, 
                        plotpol = TRUE, plotpt = TRUE)
  expect_true(nrow(y) == 56)
})

test_that("Preprocess is not changing the output whimsically", {
  x <- tile_read(system.file("extdata", "tile.bsp", package="uFTIR"))
  y <- preprocess(x, function(x){x})
  expect_true(sum(y@Spectra == x@Spectra) == prod(dim(x@Spectra)))
})

if(requireNamespace("signal", quietly = TRUE)){
  test_that("Preprocess for SpectralPack object is reshaping the wavenumbers", {
    x <- tile_read(system.file("extdata", "tile.bsp", package="uFTIR"))
    x <- tile_base_corr(x)
    x <- wavealign(x, primpke)
    preprocess(x, function(x){signal::sgolayfilt(x)})
    
    expect_true(dim(x@Readings@Spectra)[3] == length(x@Readings@wavenumbers))
    expect_true(ncol(x@Reference@Spectra) == length(x@Reference@wavenumbers))
    expect_true(length(x@Readings@wavenumbers) == length(x@Reference@wavenumbers))
  })
}
