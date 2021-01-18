test_that("readerrors", {
  vcf1a <- system.file("extdata", "set1a.vcf", package = "BinaryDosage")
  bdbadversion <- system.file("extdata", "badversion.bdose", package = "BinaryDosage")

  expect_error(ReadBinaryDosageHeader("junk"),
               "Unable to open binary dosage file")
  expect_error(ReadBinaryDosageHeader(vcf1a),
               "File does not appear to be a binary dosage file")
  expect_error(ReadBinaryDosageHeader(bdbadversion),
               "Unknown binary dosage file fromat")
})
