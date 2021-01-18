test_that("subformat2", {
  aaffile <- system.file("extdata", "aaf1a.rds", package = "BinaryDosage")
  aaf <- readRDS(aaffile)
  # Testing subformat 2
  vcf1afile <- system.file("extdata", "set1a.vcf", package = "BinaryDosage")

  bdfile12 <- tempfile()
  famfile12 <- tempfile()
  mapfile12 <- tempfile()
  expect_error(vcftobd(vcffiles = vcf1afile,
                       bdfiles = c(bdfile12, famfile12, mapfile12),
                       format = 1L,
                       subformat = 2L),
               NA)
  expect_error(getbdinfo(bdfile12),
               "Binary dosage file format 1, 2, and 3 require family and map files")
  expect_error(bdinfo12 <- getbdinfo(c(bdfile12, famfile12, mapfile12)),
               NA)
  expect_equal(unlist(bdapply(bdinfo12, getaaf)),
               aaf,
               tolerance = 4e-5)

  bdfile22 <- tempfile()
  famfile22 <- tempfile()
  mapfile22 <- tempfile()
  expect_error(vcftobd(vcffiles = vcf1afile,
                       bdfiles = c(bdfile22, famfile22, mapfile22),
                       format = 2L,
                       subformat = 2L),
               NA)
  expect_error(bdinfo22 <- getbdinfo(c(bdfile22, famfile22, mapfile22)),
               NA)
  expect_equal(unlist(bdapply(bdinfo22, getaaf)),
               aaf,
               tolerance = 4e-5)

  bdfile32 <- tempfile()
  famfile32 <- tempfile()
  mapfile32 <- tempfile()
  expect_error(vcftobd(vcffiles = vcf1afile,
                       bdfiles = c(bdfile32, famfile32, mapfile32),
                       format = 3L,
                       subformat = 2L),
               NA)
  expect_error(bdinfo32 <- getbdinfo(c(bdfile32, famfile32, mapfile32)),
               NA)
  expect_equal(unlist(bdapply(bdinfo32, getaaf)),
               aaf,
               tolerance = 4e-5)

  bdfile34 <- tempfile()
  famfile34 <- tempfile()
  mapfile34 <- tempfile()
  expect_error(vcftobd(vcffiles = vcf1afile,
                       bdfiles = c(bdfile34, famfile34, mapfile34),
                       format = 3L,
                       subformat = 4L),
               NA)
  expect_error(bdinfo34 <- getbdinfo(c(bdfile34, famfile34, mapfile34)),
               NA)
  expect_equal(unlist(bdapply(bdinfo34, getaaf)),
               aaf,
               tolerance = 4e-5)

  bdfile42 <- tempfile()
  expect_error(vcftobd(vcffiles = vcf1afile,
                       bdfiles = bdfile42,
                       format = 4L,
                       subformat = 2L),
               NA)
  expect_error(bdinfo42 <- getbdinfo(bdfile42),
               NA)
  expect_equal(unlist(bdapply(bdinfo42, getaaf)),
               aaf,
               tolerance = 4e-5)

  bdfile44 <- tempfile()
  expect_error(vcftobd(vcffiles = vcf1afile,
                       bdfiles = bdfile44,
                       format = 4L,
                       subformat = 4L,
                       bdoptions = c("aaf", "maf", "rsq")),
               NA)
  expect_error(bdinfo44 <- getbdinfo(bdfile44),
               NA)
  expect_equal(unlist(bdapply(bdinfo44, getaaf)),
               aaf,
               tolerance = 4e-5)

})

test_that("subformat1", {
  aaffile <- system.file("extdata", "aaf1a.rds", package = "BinaryDosage")
  aaf <- readRDS(aaffile)
  # Testing subformat 1
  vcf2afile <- system.file("extdata", "set2a.vcf", package = "BinaryDosage")

  bdfile11 <- tempfile()
  famfile11 <- tempfile()
  mapfile11 <- tempfile()
  expect_error(vcftobd(vcffiles = vcf2afile,
                       bdfiles = c(bdfile11, famfile11, mapfile11),
                       format = 1L,
                       subformat = 1L),
               NA)
  expect_error(bdinfo11 <- getbdinfo(c(bdfile11, famfile11, mapfile11)),
               NA)
  expect_equal(unlist(bdapply(bdinfo11, getaaf)),
               aaf,
               tolerance = 4e-5)

  bdfile21 <- tempfile()
  famfile21 <- tempfile()
  mapfile21 <- tempfile()
  expect_error(vcftobd(vcffiles = vcf2afile,
                       bdfiles = c(bdfile21, famfile21, mapfile21),
                       format = 2L,
                       subformat = 1L),
               NA)
  expect_error(bdinfo21 <- getbdinfo(c(bdfile21, famfile21, mapfile21)),
               NA)
  expect_equal(unlist(bdapply(bdinfo21, getaaf)),
               aaf,
               tolerance = 4e-5)

  bdfile31 <- tempfile()
  famfile31 <- tempfile()
  mapfile31 <- tempfile()
  expect_error(vcftobd(vcffiles = vcf2afile,
                       bdfiles = c(bdfile31, famfile31, mapfile31),
                       format = 3L,
                       subformat = 1L),
               NA)
  expect_error(bdinfo31 <- getbdinfo(c(bdfile31, famfile31, mapfile31)),
               NA)
  expect_equal(unlist(bdapply(bdinfo31, getaaf)),
               aaf,
               tolerance = 4e-5)

  bdfile33 <- tempfile()
  famfile33 <- tempfile()
  mapfile33 <- tempfile()
  expect_error(vcftobd(vcffiles = vcf2afile,
                       bdfiles = c(bdfile33, famfile33, mapfile33),
                       format = 3L,
                       subformat = 3L),
               NA)
  expect_error(bdinfo33 <- getbdinfo(c(bdfile33, famfile33, mapfile33)),
               NA)
  expect_equal(unlist(bdapply(bdinfo33, getaaf)),
               aaf,
               tolerance = 4e-5)

  bdfile41 <- tempfile()
  expect_error(vcftobd(vcffiles = vcf2afile,
                       bdfiles = bdfile41),
               NA)
  expect_error(bdinfo41 <- getbdinfo(bdfiles = bdfile41),
               NA)
  expect_equal(unlist(bdapply(bdinfo41, getaaf)),
               aaf,
               tolerance = 4e-5)

  bdfile43 <- tempfile()
  expect_error(vcftobd(vcffiles = vcf2afile,
                       bdfiles = bdfile43,
                       format = 4L,
                       subformat = 3L),
               NA)
  expect_error(bdinfo43 <- getbdinfo(bdfiles = bdfile43),
               NA)
  expect_equal(unlist(bdapply(bdinfo43, getaaf)),
               aaf,
               tolerance = 4e-5)
})
