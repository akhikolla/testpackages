# Satisfy R CMD check
options("mc.cores" = 2)


context("test-StoxExport: DATRAS export")
example <- system.file("testresources", "biotic_v3_example.xml", package="RstoxData")	
data <- ReadBiotic(example)

data[[1]]$fishstation[, stationstartdate := stationstopdate]
datras1 <- RstoxData::WriteICESDatras(data, Combine = FALSE)
datras2 <- RstoxData::WriteICESDatras(data, Combine = TRUE)
expect_equal(nrow(datras1[[1]]$HH), 2)

context("test-StoxExport: ICES biotic export")
example <- system.file("testresources", "biotic_v3_example.xml", package="RstoxData")
data <- ReadBiotic(example)

data[[1]]$fishstation[, stationstartdate := stationstopdate]
icesbiotic1 <- RstoxData::WriteICESBiotic(data, Combine = FALSE)
icesbiotic2 <- RstoxData::WriteICESBiotic(data, Combine = TRUE)
expect_equal(nrow(icesbiotic1[[1]]$Haul), 2)


context("test-StoxExport: ICES acoustic export #1")
example <- system.file("testresources", "ICES_Acoustic_1.xml", package="RstoxData")
data <- ReadAcoustic(example)
icesacoustic1 <- RstoxData::WriteICESAcoustic(data, Combine = FALSE)
icesacoustic2 <- RstoxData::WriteICESAcoustic(data, Combine = TRUE)

expect_equal(nrow(icesacoustic1[[1]]$Data), 11)

context("test-StoxExport: ICES acoustic export #2")
example <- system.file("testresources", "ICES_Acoustic_2.xml", package="RstoxData")
data <- ReadAcoustic(example)
icesacoustic1 <- RstoxData::WriteICESAcoustic(data, Combine = FALSE)
icesacoustic2 <- RstoxData::WriteICESAcoustic(data, Combine = TRUE)

expect_equal(nrow(icesacoustic1[[1]]$Data), 12)
