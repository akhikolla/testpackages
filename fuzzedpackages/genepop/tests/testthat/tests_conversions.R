cat("test_conversions.R:\n")

locinfile <- genepopExample('sample.txt')
outfile <- conversion(locinfile, format="Fstat", "sample.txt.DAT")

outfile <- diploidize(inputFile = locinfile,outputFile="Dsample.txt")
testthat::expect_equal(readLines(outfile)[30]," last pop, 0101 002001 0101 0401 0807 0202 ")
