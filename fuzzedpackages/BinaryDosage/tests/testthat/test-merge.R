test_that("merge", {
  expect_error(bdmerge(bdfiles = c("file1", "file2")),
               "No output files specified")
  expect_error(bdmerge(mergefiles = 1L,
                       bdfiles = c("file1", "file2")),
               "Output file names must be a character values")
  expect_error(bdmerge(mergefiles = c("file1", "file2", "file3"),
                       bdfiles = c("file1", "file2")),
               "Only one file name is needed when using format 4")
  expect_error(bdmerge(mergefiles = "file1",
                       format = 3,
                       bdfiles = c("file1", "file2")),
               "Three file names are required when using formats 1, 2, and 3")
  expect_error(bdmerge(mergefiles = "",
                       bdfiles = c("file1", "file2")),
               "Output file names cannot be blank")

  expect_error(bdmerge(mergefiles = "file1"),
               "No files specified")
  expect_error(bdmerge(mergefiles = "file1",
                       bdfiles = 1),
               "bdfiles must be a character vector")
  expect_error(bdmerge(mergefiles = "file1",
                       bdfiles = "file1"),
               "At least two binary dosage files must be specified")

  expect_error(bdmerge(mergefiles = "file1",
                       bdfiles = c("file1", "file2"),
                       famfiles = 1,
                       mapfiles = c("file1", "file2")),
               "famfiles must be a character vector")
  expect_error(bdmerge(mergefiles = "file1",
                       bdfiles = c("file1", "file2"),
                       famfiles = c("file1", "file2"),
                       mapfiles = 1),
               "mapfiles must be a character vector")
  expect_error(bdmerge(mergefiles = "file1",
                       bdfiles = c("file1", "file2"),
                       famfiles = c("file1", "file2")),
               "If famfiles is specified, mapfiles must be specified")
  expect_error(bdmerge(mergefiles = "file1",
                       bdfiles = c("file1", "file2"),
                       mapfiles = c("file1", "file2")),
               "If mapfiles is specified, famfiles must be specified")
  expect_error(bdmerge(mergefiles = "file1",
                       bdfiles = c("file1", "file2"),
                       famfiles = c("file1", "file2"),
                       mapfiles = "file1"),
               "If famfiles and mapfiles are specified they must have the same length as bdfiles")

  expect_error(bdmerge(mergefiles = "file1",
                       format = "a",
                       bdfiles = c("file1", "file2")),
               "format must be an integer value")
  expect_error(bdmerge(mergefiles = "file1",
                       format = 1:2,
                       bdfiles = c("file1", "file2")),
               "format must be an integer vector of length 1")
  expect_error(bdmerge(mergefiles = "file1",
                       format = 1.2,
                       bdfiles = c("file1", "file2")),
               "format must be an integer")
  expect_error(bdmerge(mergefiles = "file1",
                       format = 5,
                       bdfiles = c("file1", "file2")),
               "format must be an integer value from 1 to 4")

  expect_error(bdmerge(mergefiles = "file1",
                       subformat = "a",
                       bdfiles = c("file1", "file2")),
               "subformat must be an integer value")
  expect_error(bdmerge(mergefiles = "file1",
                       subformat = 1:2,
                       bdfiles = c("file1", "file2")),
               "subformat must be an integer vector of length 1")
  expect_error(bdmerge(mergefiles = "file1",
                       subformat = 1.2,
                       bdfiles = c("file1", "file2")),
               "subformat must be an integer")
  expect_error(bdmerge(mergefiles = "file1",
                       subformat = 5,
                       bdfiles = c("file1", "file2")),
               "subformat must be an integer value from 0 to 4")
  expect_error(bdmerge(mergefiles = "file1",
                       format = 2,
                       subformat = 3,
                       bdfiles = c("file1", "file2")),
               "subformat must be an integer value from 0 to 2 for formats 1 and 2")

  expect_error(bdmerge(mergefiles = "file1",
                       bdfiles = c("file1", "file2"),
                       onegroup = 1),
               "onegroup must be logical value")
  expect_error(bdmerge(mergefiles = "file1",
                       bdfiles = c("file1", "file2"),
                       onegroup = c(TRUE, TRUE)),
               "onegroup must be a logical vector of length 1")

  expect_error(bdmerge(mergefiles = "file1",
                       bdfiles = c("file1", "file2"),
                       bdoptions = 1),
               "bdoptions must be a character vector")
  expect_error(bdmerge(mergefiles = "file1",
                       bdfiles = c("file1", "file2"),
                       onegroup = FALSE,
                       bdoptions = "aaf"),
               "bdoptions can only be used if onegroup is TRUE")
  expect_error(bdmerge(mergefiles = c("file1", "file2", "file3"),
                       bdfiles = c("file1", "file2"),
                       format = 3,
                       bdoptions = "aaf"),
               "bdoptions can only be used if format = 4")
  expect_error(bdmerge(mergefiles = "file1",
                       bdfiles = c("file1", "file2"),
                       bdoptions = "abc"),
               "Only valid bdoptions are aaf, maf, and rsq")

  expect_error(bdmerge(mergefiles = "file1",
                       bdfiles = c("file1", "file2"),
                       snpjoin = "abc"),
               "snpjoin must have a value of either \"inner\" or \"outer\"")

  bdvcf1afile <- system.file("extdata", "vcf1a.bdose", package = "BinaryDosage")
  bdvcf1bfile <- system.file("extdata", "vcf1b.bdose", package = "BinaryDosage")
  mergefiles <- tempfile()

  expect_error(bdmerge(mergefiles = mergefiles,
                       bdfiles = c(bdvcf1afile, bdvcf1bfile),
                       bdoptions = c("aaf", "maf", "rsq")),
               NA)
  expect_error(bdinfo <- getbdinfo(mergefiles), NA)

  vcf1afile <- system.file("extdata", "set1a.vcf", package = "BinaryDosage")
  vcf1ainfo <- system.file("extdata", "set1a.info", package = "BinaryDosage")
  vcf1bfile <- system.file("extdata", "set1b.vcf", package = "BinaryDosage")
  vcf1binfo <- system.file("extdata", "set1b.info", package = "BinaryDosage")
  bdfile1a <- tempfile()
  bdfile1b <- tempfile()
  vcftobd(vcffiles = c(vcf1afile, vcf1ainfo),
          bdfiles = bdfile1a)
  vcftobd(vcffiles = c(vcf1bfile, vcf1binfo),
          bdfiles = bdfile1b)

  bdinfo1a <- getbdinfo(bdfile1a)
  bdinfo1b <- getbdinfo(bdfile1b)

  mergefiles1g <- tempfile()
  expect_error(bdmerge(mergefiles = mergefiles1g,
                       bdfiles = c(bdfile1a, bdfile1b),
                       onegroup = FALSE),
               NA)
  expect_error(getbdinfo(mergefiles1g), NA)

  vcf2afile <- system.file("extdata", "set2a.vcf", package = "BinaryDosage")
  vcf2bfile <- system.file("extdata", "set2b.vcf", package = "BinaryDosage")

  bdfile2a <- tempfile()
  bdfam2a <- tempfile()
  bdmap2a <- tempfile()
  bdfile2b <- tempfile()
  bdfam2b <- tempfile()
  bdmap2b <- tempfile()
  vcftobd(vcffiles = vcf2afile,
          bdfiles = c(bdfile2a, bdfam2a, bdmap2a),
          format = 3L)
  vcftobd(vcffiles = vcf2bfile,
          bdfiles = c(bdfile2b, bdfam2b, bdmap2b),
          format = 3L)

  bdinfo2a <- getbdinfo(c(bdfile2a, bdfam2a, bdmap2a))
  bdinfo2b <- getbdinfo(c(bdfile2b, bdfam2b, bdmap2b))

  mergefiles2g <- tempfile()
  expect_error(bdmerge(mergefiles = mergefiles2g,
                       bdfiles = c(bdfile2a, bdfile2b),
                       famfiles = c(bdfam2a, bdfam2b),
                       mapfiles = c(bdmap2a, bdmap2b)),
               NA)
  expect_error(getbdinfo(mergefiles2g), NA)

  vcf1brfile <- system.file("extdata", "set1b_rsub.vcf", package = "BinaryDosage")
  bdfile1br <- tempfile()
  vcftobd(vcffiles = vcf1brfile,
          bdfiles = bdfile1br)
  bdinfo1br <- getbdinfo(bdfile1br)

  mergefile3 <- tempfile()
  expect_error(bdmerge(mergefiles = mergefile3,
                       bdfiles = c(bdfile1a, bdfile1br)),
               "There are duplicate samples in the files to merge")

  gen3bfile <- system.file("extdata", "set3b.imp", package = "BinaryDosage")
  gen3bsample <- system.file("extdata", "set3b.sample", package = "BinaryDosage")
  bdfile3b <- tempfile()
  gentobd(genfiles = c(gen3bfile, gen3bsample),
          snpcolumns = c(0L, 2L:5L),
          bdfiles = bdfile3b)
  bdinfo3b <- getbdinfo(bdfiles = bdfile3b)
  mergefile3b <- tempfile()
  expect_error(bdmerge(mergefiles = mergefile3b,
                       bdfiles = c(bdfile1a, bdfile3b)),
               "Some files use FID and others do not")

  bdvcf1afile <- system.file("extdata", "vcf1a.bdose", package = "BinaryDosage")
  vcf1brs <- system.file("extdata", "set1b_rssnp.vcf", package = "BinaryDosage")
  bdinfo1a <- getbdinfo(bdvcf1afile)
  bdfile1brs <- tempfile()
  vcftobd(vcffiles = vcf1brs,
          bdfiles = bdfile1brs)
  bdinfo1brs <- getbdinfo(bdfiles = bdfile1brs)
  merge1brs <- tempfile()
  expect_error(bdmerge(mergefiles = merge1brs,
                       bdfiles = c(bdvcf1afile, bdfile1brs)),
               NA)
  expect_error(getbdinfo(merge1brs), NA)
})
