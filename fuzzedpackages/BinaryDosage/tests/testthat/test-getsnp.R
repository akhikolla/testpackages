test_that("getsnp", {
  vcf1afile <- system.file("extdata", "set1a.vcf", package = "BinaryDosage")
  vcfinfo <- getvcfinfo(vcffiles = vcf1afile)
  vcf1abdfile <- system.file("extdata", "vcf1a.bdose", package = "BinaryDosage")
  bdinfo <- getbdinfo(bdfiles = vcf1abdfile)

  expect_error(getsnp(snp = 1),
               "bdinfo missing")
  expect_error(getsnp(bdinfo = 1,
                      snp = 1),
               "bdinfo is not of class genetic-info")
  expect_error(getsnp(bdinfo = vcfinfo,
                      snp = 1),
               "bdinfo does not contain information about a binary dosage file")

  expect_error(getsnp(bdinfo = bdinfo),
               "No SNP specified")
  expect_error(getsnp(bdinfo = bdinfo,
                      snp = 1:3),
               "snp must be of length one")
  expect_error(getsnp(bdinfo = bdinfo,
                      snp = 1.3),
               "snp must be a character or integer value")
  expect_error(getsnp(bdinfo = bdinfo,
                      snp = TRUE),
               "snp must be a character or integer value")
  expect_error(getsnp(bdinfo = bdinfo,
                      snp = "xyz"),
               "Cannot find SNP in bdinfo")
  expect_error(getsnp(bdinfo = bdinfo,
                      snp = 100),
               "snp value out or range")

  expect_error(getsnp(bdinfo = bdinfo,
                      snp = 1,
                      FALSE),
               NA)
  expect_error(getsnp(bdinfo = bdinfo,
                      snp = "1:11000:T:C",
                      TRUE),
               NA)
})
