# Satisfy R CMD check
options("mc.cores" = 2)

# StoxAcoustic
# NMD Echosounder
example <- system.file("testresources","libas_ListUserFile20__L40.0-2259.9_small.xml", package="RstoxData")

context("test-StoxAcoustic: DOM echosounder")
defaultParseEchosounder <- readXmlFile(example, stream = F)
sa1 <- StoxAcoustic(list(defaultParseEchosounder))
expect_equal(nrow(sa1$Beam), 2)

context("test-StoxAcoustic: stream parse echosounder")
streamParseEchosounder <- readXmlFile(example, stream = T)
sa2 <- StoxAcoustic(list(streamParseEchosounder))
expect_equal(nrow(sa2$Beam), 2)

context("test-StoxAcoustic: all.equal")
expect_true(all.equal(sa1, sa2))


## ICES Acoustic
icesFiles <- c("ICES_Acoustic_1.xml", "ICES_Acoustic_2.xml")
exampleDir <- system.file("testresources","", package="RstoxData")

for(item in icesFiles) {
	context(paste("test-StoxAcoustic: ICES acoustic data", item, "to StoxAcoustic: DOM"))
    icesDataA <- StoxAcoustic(list(readXmlFile(paste0(exampleDir, "/", item), stream = F)))
	expect_equal(nrow(icesDataA$Log), 2)

	context(paste("test-StoxAcoustic: ICES acoustic data", item, "to StoxAcoustic: Stream"))
    icesDataB <- StoxAcoustic(list(readXmlFile(paste0(exampleDir, "/", item), stream = T)))
	expect_equal(nrow(icesDataB$Log), 2)

	context(paste("test-StoxAcoustic: ICES acoustic data", item, "to StoxAcoustic DOM == stream"))
    expect_true(all.equal(icesDataA, icesDataB))
}
