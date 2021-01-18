test_that("filemismatch12", {
  gen3afile <- system.file("extdata", "set3a.imp", package = "BinaryDosage")
  gen3asample <- system.file("extdata", "set3a.sample", package = "BinaryDosage")
  bdfile3a <- tempfile()
  bdfam3a <- tempfile()
  bdmap3a <- tempfile()
  gentobd(genfiles = c(gen3afile, gen3asample),
          snpcolumns = c(0L, 2L:5L),
          bdfiles = c(bdfile3a, bdfam3a, bdmap3a),
          format = 3)

  vcf1basnpfile <- system.file("extdata", "set1b_asnp.vcf", package = "BinaryDosage")
  bdfile1basnp <- tempfile()
  bdfam1basnp <- tempfile()
  bdmap1basnp <- tempfile()
  vcftobd(vcffiles = vcf1basnpfile,
          bdfiles = c(bdfile1basnp, bdfam1basnp, bdmap1basnp),
          format = 3)
  expect_error(getbdinfo(bdfiles = c(bdfile1basnp, bdfam3a, bdmap3a)),
               "Subject file does not line up with binary dosage file")

  vcf2bfile <- system.file("extdata", "set2b.vcf", package = "BinaryDosage")
  bdfile2b <- tempfile()
  bdfam2b <- tempfile()
  bdmap2b <- tempfile()
  vcftobd(vcffiles = vcf2bfile,
          bdfiles = c(bdfile2b, bdfam2b, bdmap2b),
          format = 3)
  expect_error(getbdinfo(bdfiles = c(bdfile2b, bdfam3a, bdmap3a)),
               "Subject file does not line up with binary dosage file")
})

test_that("filemismatch34", {
  gen3afile <- system.file("extdata", "set3a.imp", package = "BinaryDosage")
  gen3asample <- system.file("extdata", "set3a.sample", package = "BinaryDosage")
  bdfile3a <- tempfile()
  bdfam3a <- tempfile()
  bdmap3a <- tempfile()
  gentobd(genfiles = c(gen3afile, gen3asample),
          snpcolumns = c(0L, 2L:5L),
          bdfiles = c(bdfile3a, bdfam3a, bdmap3a),
          format = 3)

  vcf1basnpfile <- system.file("extdata", "set1b_asnp.vcf", package = "BinaryDosage")
  bdfile1basnp <- tempfile()
  bdfam1basnp <- tempfile()
  bdmap1basnp <- tempfile()
  vcftobd(vcffiles = vcf1basnpfile,
          bdfiles = c(bdfile1basnp, bdfam1basnp, bdmap1basnp),
          format = 3,
          subformat = 4)
  expect_error(getbdinfo(bdfiles = c(bdfile1basnp, bdfam3a, bdmap3a)),
               "Subject file does not line up with binary dosage file")
  expect_error(getbdinfo(bdfiles = c(bdfile1basnp, bdfam1basnp, bdmap3a)),
               "Map file does not line up with binary dosage file")

  vcf2bfile <- system.file("extdata", "set2b.vcf", package = "BinaryDosage")
  bdfile2b <- tempfile()
  bdfam2b <- tempfile()
  bdmap2b <- tempfile()
  vcftobd(vcffiles = vcf2bfile,
          bdfiles = c(bdfile2b, bdfam2b, bdmap2b),
          format = 3,
          subformat = 3)
  expect_error(getbdinfo(bdfiles = c(bdfile2b, bdfam3a, bdmap3a)),
               "Subject file does not line up with binary dosage file")
  expect_error(getbdinfo(bdfiles = c(bdfile2b, bdfam2b, bdmap3a)),
               "Map file does not line up with binary dosage file")
})
