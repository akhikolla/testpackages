# Satisfy R CMD check
options("mc.cores" = 2)

# StoxBiotic
## NMD Biotic v3.1
example <- system.file("testresources","biotic3.1_example.xml", package="RstoxData")

context("test-StoxBiotic: using DOM biotic v3.1")
defaultParseBiotic <- readXmlFile(example, stream = F)
sb1 <- StoxBiotic(list(defaultParseBiotic))
expect_equal(nrow(sb1$Haul), 2)

context("test-StoxBiotic: using stream parse biotic v3.1")
streamParseBiotic <- readXmlFile(example, stream = T)
sb2 <- StoxBiotic(list(streamParseBiotic))
expect_equal(nrow(sb2$Haul), 2)

context("test-StoxBiotic: all.equal")
expect_true(all.equal(sb1, sb2))

## NMD Biotic v3
example <- system.file("testresources","biotic_v3_example.xml", package="RstoxData")

context("test-StoxBiotic: using DOM biotic v3")
defaultParseBiotic <- readXmlFile(example, stream = F)
sb1 <- StoxBiotic(list(defaultParseBiotic))
expect_equal(nrow(sb1$Haul), 2)

context("test-StoxBiotic: using stream parse biotic v3.1")
streamParseBiotic <- readXmlFile(example, stream = T)
sb2 <- StoxBiotic(list(streamParseBiotic))
expect_equal(nrow(sb2$Haul), 2)

context("test-StoxBiotic: all.equal")
expect_true(all.equal(sb1, sb2))

## ICES Biotic
icesFiles <- c("ICES_Biotic_1.xml", "ICES_Biotic_2.xml")
exampleDir <- system.file("testresources","", package="RstoxData")

for(item in icesFiles) {
	context(paste("test-StoxBiotic: ICES biotic data", item, "to StoxBiotic: DOM"))
        icesDataA <- StoxBiotic(list(readXmlFile(paste0(exampleDir, "/", item), stream = F)))
	expect_equal(nrow(icesDataA$Individual), 4)

	context(paste("test-StoxBiotic: ICES biotic data", item, "to StoxBiotic: Stream"))
        icesDataB <- StoxBiotic(list(readXmlFile(paste0(exampleDir, "/", item), stream = T)))
	expect_equal(nrow(icesDataB$Individual), 4)

	context(paste("test-StoxBiotic: ICES biotic data", item, "to StoxBiotic DOM == stream"))
        expect_true(all.equal(icesDataA, icesDataB))
}

## ICES Biotic (missing individual)
item <- "ICES_Biotic_1_missingind.xml"
exampleDir <- system.file("testresources","", package="RstoxData")

context(paste("test-StoxBiotic: ICES biotic data (with missing individual)", item, "to StoxBiotic: DOM"))
icesDataA <- StoxBiotic(list(readXmlFile(paste0(exampleDir, "/", item), stream = F)))
expect_equal(nrow(icesDataA$Individual), 6)
expect_equal(icesDataA$Sample[1, ]$SampleCount, 6)


context(paste("test-StoxBiotic: ICES biotic data (with missing individual)", item, "to StoxBiotic: Stream"))
icesDataB <- StoxBiotic(list(readXmlFile(paste0(exampleDir, "/", item), stream = T)))
expect_equal(nrow(icesDataB$Individual), 6)
expect_equal(icesDataB$Sample[1, ]$SampleCount, 6)

context(paste("test-StoxBiotic: ICES biotic data (with missing individual)", item, "to StoxBiotic DOM == stream"))
expect_true(all.equal(icesDataA, icesDataB))
