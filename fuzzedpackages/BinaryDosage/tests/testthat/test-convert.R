test_that("vcftobd", {
  expect_error(vcftobd(bdfiles = "file1"),
               "No VCF file specified")
  expect_error(vcftobd(vcffiles = "file1"),
               "No output files specified")
  expect_error(vcftobd(vcffiles = "file1",
                       bdfiles = 1),
               "bdfiles must be a vector of characters")

  expect_error(vcftobd(vcffiles = "file1",
                       bdfiles = "file2",
                       format = ""),
               "format must be an integer value")
  expect_error(vcftobd(vcffiles = "file1",
                       bdfiles = "file2",
                       format = 1:2),
               "format must be an integer vector of length 1")
  expect_error(vcftobd(vcffiles = "file1",
                       bdfiles = "file2",
                       format = 1.2),
               "format must be an integer")
  expect_error(vcftobd(vcffiles = "file1",
                       bdfiles = "file2",
                       format = 5),
               "format must be an integer value from 1 to 4")

  expect_error(vcftobd(vcffiles = "file1",
                       bdfiles = "file2",
                       subformat = ""),
               "subformat must be an integer value")
  expect_error(vcftobd(vcffiles = "file1",
                       bdfiles = "file2",
                       subformat = 1:2),
               "subformat must be an integer vector of length 1")
  expect_error(vcftobd(vcffiles = "file1",
                       bdfiles = "file2",
                       subformat = 1.2),
               "subformat must be an integer")
  expect_error(vcftobd(vcffiles = "file1",
                       bdfiles = "file2",
                       subformat = 5),
               "subformat must be an integer value from 0 to 4")
  expect_error(vcftobd(vcffiles = "file1",
                       bdfiles = "file2",
                       format = 2,
                       subformat = 3),
               "subformat must be an integer value from 0 to 2 for formats 1 and 2")

  expect_error(vcftobd(vcffiles = "file1",
                       bdfiles = c("file2", "file3"),
                       format = 4),
               "Only one output file name is needed when using format 4")
  expect_error(vcftobd(vcffiles = "file1",
                       bdfiles = c("file2", "file3"),
                       format = 3),
               "Three output file names are required when using formats 1, 2, and 3")
  expect_error(vcftobd(vcffiles = "file1",
                       bdfiles = c("file2", "file3", ""),
                       format = 3),
               "Output file names cannot be blank")

  expect_error(vcftobd(vcffiles = "file1",
                       bdfiles = "file2",
                       snpidformat = "A"),
               "snpidformat must be an integer value")
  expect_error(vcftobd(vcffiles = "file1",
                       bdfiles = "file2",
                       snpidformat = 1:2),
               "snpidformat must be an integer vector of length 1")
  expect_error(vcftobd(vcffiles = "file1",
                       bdfiles = "file2",
                       snpidformat = 1.2),
               "snpidformat must be an integer")
  expect_error(vcftobd(vcffiles = "file1",
                       bdfiles = "file2",
                       snpidformat = 4),
               "snpidformat must be and integer from -1 to 3")

  expect_error(vcftobd(vcffiles = "file1",
                       bdfiles = "file2",
                       bdoptions = 4),
               "bdoptions must be a character array")
  expect_error(vcftobd(vcffiles = "file1",
                       bdfiles = c("file2", "file3", "file4"),
                       format = 3,
                       bdoptions = "aaf"),
               "bdoptions can only be used with format 4")
  expect_error(vcftobd(vcffiles = "file1",
                       bdfiles = "file2",
                       bdoptions = "abc"),
               "Only valid bdoptions are aaf, maf, and rsq")

  vcf1afile = system.file("extdata", "set1a.vcf", package = "BinaryDosage")
  vcf1ainfo <- system.file("extdata", "set1a.info", package = "BinaryDosage")
  bdfiles <- tempfile()
  expect_error(vcftobd(vcffiles = c(vcf1afile, vcf1ainfo),
                       bdfiles = bdfiles,
                       bdoptions = c("aaf", "maf", "rsq")),
               NA)
  expect_error(bdinfo <- getbdinfo(bdfiles = bdfiles), NA)

  vcf2afile = system.file("extdata", "set2a.vcf", package = "BinaryDosage")
  bdfiles2 <- tempfile()
  expect_error(vcftobd(vcffiles = c(vcf2afile),
                       bdfiles = bdfiles2,
                       snpidformat = -1),
               NA)
  expect_error(bdinfo <- getbdinfo(bdfiles = bdfiles2), NA)

})

test_that("gentobd", {
  expect_error(gentobd(bdfiles = "file1"),
               "No gen file specified")
  expect_error(gentobd(genfiles = "file1"),
               "No output files specified")

  gen3afile <- system.file("extdata", "set3a.imp", package = "BinaryDosage")
  gen3asample <- system.file("extdata", "set3a.sample", package = "BinaryDosage")
  bdfile3a <- tempfile()
  bdfam3a <- tempfile()
  bdmap3a <- tempfile()
  expect_error(getgeninfo(genfiles = c(gen3afile, gen3asample),
                          snpcolumns = c(0L, 2L:5L),
                          chromosome = "",
                          index = TRUE),
               NA)
  expect_error(gentobd(genfiles = c(gen3afile, gen3asample),
                       snpcolumns = c(0L, 2L:5L),
                       bdfiles = c(bdfile3a, bdfam3a, bdmap3a),
                       chromosome = "",
                       format = 3),
               NA)
  expect_error(getbdinfo(bdfiles = c(bdfile3a, bdfam3a, bdmap3a)), NA)

  gen2afile <- system.file("extdata", "set2a.imp", package = "BinaryDosage")
  gen2asample <- system.file("extdata", "set2a.sample", package = "BinaryDosage")
  bdfile2a <- tempfile()
  expect_error(gentobd(genfiles = c(gen2afile, gen2asample),
                       snpcolumns = c(1L, 3L, 2L, 4L, 5L),
                       impformat = 1L,
                       bdfiles = bdfile2a,
                       snpidformat = -1L,
                       bdoptions = c("aaf", "maf", "rsq")),
               NA)
  expect_error(getbdinfo(bdfiles = bdfile2a), NA)

  gen4afile <- system.file("extdata", "set4a.imp.gz", package = "BinaryDosage")
  gen4asample <- system.file("extdata", "set4a.sample", package = "BinaryDosage")
  bdfile4a <- tempfile()
  expect_error(gentobd(genfiles = c(gen4afile, gen4asample),
                       snpcolumns = c(1L:2L, 4L:6L),
                       startcolumn = 7L,
                       impformat = 2L,
                       gz = TRUE,
                       bdfile = bdfile4a),
               NA)
  expect_error(getbdinfo(bdfiles = bdfile4a), NA)
})

test_that("writecpp", {
  expect_error(WriteBinaryDosageBaseHeader(filename = "",
                                           format = 1L,
                                           subformat = 2L),
               "Unable to create output file")
  expect_error(WriteBinaryDosageHeader4A(filename = "",
                                         headerEntries = 1L,
                                         numSubjects = 1L,
                                         numSNPs = 1L,
                                         groups = 1L,
                                         fid = "",
                                         sid = "",
                                         snpid = "",
                                         chromosome = "",
                                         location = 1L,
                                         reference = "",
                                         alternate = "",
                                         aaf = 1,
                                         maf = 1,
                                         avgCall = 1,
                                         rsq = 1,
                                         offsets = 1L,
                                         numIndices = 1L),
               "Unable to open file for read/write")
  expect_error(WriteBinaryDosageIndicesC(filename = "",
                                         headersize = 0L,
                                         datasize = 1L),
               "Unable to open file for read/write")
  expect_error(updatesnpinfo(filename = "",
                             offset = 0L,
                             value = 1),
               "Unable to open file for read/write")
  expect_error(WriteBinaryDosageHeader3A(filename = "",
                                         numSubjects = 1L),
               "Unable to open file for appending")
  expect_error(WriteBinaryDosageHeader3B(filename = "",
                                         md5samples = "",
                                         md5SNPs = "",
                                         numIndices = 1L),
               "Unable to open file for appending")
  expect_error(WriteBinaryDosageDataC(filename = "",
                                      dosage = 1,
                                      us = 1L,
                                      base = 1L),
               "Unable to open file for appending")
  expect_error(WriteBinaryP1P2Data(filename = "",
                                   p1 = 1,
                                   p2 = 1,
                                   us = 1L,
                                   base = 1L),
               "Unable to open file for appending")
  expect_error(WriteBinaryCompressed(filename = "",
                                     dosage = 1,
                                     p0 = 0,
                                     p1 = 1,
                                     p2 = 0,
                                     snpnumber = -1L,
                                     datasize = 1L,
                                     us = 1L),
               "Unable to open file for appending")
  bdfile <- tempfile()
  dosage <- as.numeric(c(NA, 1))
  p0 <- as.numeric(c(NA, NA))
  p1 <- as.numeric(c(NA, NA))
  p2 <- as.numeric(c(NA, NA))
  snpnumber = 0L
  datasize = 0L
  us <- integer(4)
  expect_error(WriteBinaryCompressed(filename = bdfile,
                                     dosage = dosage,
                                     p0 = p0,
                                     p1 = p1,
                                     p2 = p2,
                                     snpnumber = snpnumber,
                                     datasize = datasize,
                                     us = us),
               NA)
  rdosage <- numeric(2)
  rp0 <- numeric(2)
  rp1 <- numeric(2)
  rp2 <- numeric(2)
  rus <- integer(4)
  expect_error(ReadBinaryDosageDataCompressed(filename = bdfile,
                                              index = 0L,
                                              datasize = datasize,
                                              numsub = 2L,
                                              dosage = rdosage,
                                              p0 = rp0,
                                              p1 = rp1,
                                              p2 = rp2,
                                              us = rus),
               NA)
  expect_true(is.na(rdosage[1]))
  expect_equal(rdosage[2], 1)
  expect_true(all(is.na(rp0)))
  expect_true(all(is.na(rp1)))
  expect_true(all(is.na(rp2)))

  bdfile <- tempfile()
  p1 <- as.numeric(c(NA, 0.0025))
  p2 <- as.numeric(c(NA, 0.999))
  us <- integer(4)
  base <- 3L
  expect_error(WriteBinaryP1P2Data(filename = bdfile,
                                   p1 = p1,
                                   p2 = p2,
                                   us = us,
                                   base = base),
               NA)

  rdosage <- numeric(2)
  rp0 <- numeric(2)
  rp1 <- numeric(2)
  rp2 <- numeric(2)
  rus <- integer(4)
  expect_error(ReadBinaryDosageDataP1P2(filename = bdfile,
                                        headersize = 0L,
                                        numsub = 2L,
                                        snp = 1L,
                                        dosage = rdosage,
                                        p0 = rp0,
                                        p1 = rp1,
                                        p2 = rp2,
                                        us = rus,
                                        base = base),
               NA)
  expect_true(is.na(rdosage[1]))
  expect_equal(rdosage[2], 2)
  expect_true(is.na(rp0[1]))
  expect_true(is.na(rp1[1]))
  expect_true(is.na(rp2[1]))
  expect_equal(rp0[2], 0)
  expect_equal(rp1[2], 0.0025)
  expect_equal(rp2[2], 0.999)
})
