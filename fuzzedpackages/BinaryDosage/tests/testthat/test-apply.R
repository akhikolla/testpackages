test_that("bdapply", {
  vcf1abdfile <- system.file("extdata", "set1a.vcf", package = "BinaryDosage")
  vcfinfo <- getvcfinfo(vcffiles = vcf1abdfile,
                        index = FALSE)
  vcf1abdfile <- system.file("extdata", "vcf1a.bdose", package = "BinaryDosage")
  expect_error(bdinfo <- getbdinfo(bdfiles = vcf1abdfile), NA)
  func <- getaaf
  aaf1afile <- system.file("extdata", "aaf1a.rds", package = "BinaryDosage")
  aaf1a <- readRDS(aaf1afile)

  expect_error(bdapply(),
               "No binary dosage file information specified")
  expect_error(bdapply(bdinfo = 1,
                       func = func),
               "bdinfo does not contain information about a binary dosage file")
  expect_error(bdapply(bdinfo = vcfinfo,
                       func = func),
               "bdinfo does not contain information about a binary dosage file")

  expect_error(bdapply(bdinfo = bdinfo),
               "No function specified")
  expect_error(bdapply(bdinfo = bdinfo, func = 1),
               "func is not a function")

  expect_equal(unlist(bdapply(bdinfo = bdinfo,
                              func = func)),
               aaf1a,
               tolerance = 4e-5)
})

test_that("vcfapply", {
  vcf1abdfile <- system.file("extdata", "set1a.vcf", package = "BinaryDosage")
  vcfinfo <- getvcfinfo(vcffiles = vcf1abdfile)
  vcf1abdfile <- system.file("extdata", "vcf1a.bdose", package = "BinaryDosage")
  bdinfo <- getbdinfo(bdfiles = vcf1abdfile)
  func <- getaaf
  aaf1afile <- system.file("extdata", "aaf1a.rds", package = "BinaryDosage")
  aaf1a <- readRDS(aaf1afile)

  expect_error(vcfapply(),
               "No vcf file information specified")
  expect_error(vcfapply(vcfinfo = 1,
                        func = func),
               "vcfinfo does not appear to contain information about a vcf file")
  expect_error(vcfapply(vcfinfo = bdinfo,
                        func = func),
               "vcfinfo does not appear to contain information about a vcf file")

  expect_error(vcfapply(vcfinfo = vcfinfo),
               "No function specified")
  expect_error(vcfapply(vcfinfo = vcfinfo, func = 1),
               "func is not a function")

  expect_equal(unlist(vcfapply(vcfinfo = vcfinfo,
                               func = func)),
               aaf1a,
               tolerance = 4e-5)

  vcf1andfile <- system.file("extdata", "set1a_nd.vcf.gz", package = "BinaryDosage")
  expect_error(vcfinfo <- getvcfinfo(vcf1andfile,
                                     gz = TRUE,
                                     index = FALSE),
               NA)
  expect_equal(unlist(vcfapply(vcfinfo = vcfinfo,
                               func = getaaf)),
               aaf1a,
               tolerance = 4e-5)
})

test_that("genapply", {
  gen1afile <- system.file("extdata", "set1a.imp", package = "BinaryDosage")
  geninfo <- getgeninfo(genfiles = gen1afile,
                        snpcolumns = c(1L, 3L, 2L, 4L, 5L),
                        header = TRUE)
  vcf1abdfile <- system.file("extdata", "vcf1a.bdose", package = "BinaryDosage")
  bdinfo <- getbdinfo(bdfiles = vcf1abdfile)
  func <- getaaf
  aaf1afile <- system.file("extdata", "aaf1a.rds", package = "BinaryDosage")
  aaf1a <- readRDS(aaf1afile)

  expect_error(genapply(),
               "No gen file information specified")
  expect_error(genapply(geninfo = 1,
                        func = func),
               "geninfo does not appear to contain information about a gen file")
  expect_error(genapply(geninfo = bdinfo,
                        func = func),
               "geninfo does not appear to contain information about a gen file")

  expect_error(genapply(geninfo = geninfo),
               "No function specified")
  expect_error(genapply(geninfo = geninfo, func = 1),
               "func is not a function")

  expect_equal(unlist(genapply(geninfo = geninfo,
                               func = func)),
               aaf1a,
               tolerance = 4e-5)

  gen4afile <- system.file("extdata", "set4a.imp.gz", package = "BinaryDosage")
  gen4asample <- system.file("extdata", "set4a.sample", package = "BinaryDosage")
  expect_error(geninfo <- getgeninfo(genfiles = c(gen4afile, gen4asample),
                                     snpcolumns = c(1L, 2L, 4L:6L),
                                     startcolumn = 7L,
                                     impformat = 2L,
                                     gz = TRUE,
                                     index = FALSE),
               NA)
  expect_equal(unlist(genapply(geninfo = geninfo,
                               func = func)),
               aaf1a,
               tolerance = 4e-5)

})
