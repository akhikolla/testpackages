
context("test-stoxLanding")
landingXML <- readXmlFile(system.file("testresources", "landing.xml", package="RstoxData"), stream = T)
flatSL <- StoxLanding(landingXML)
expected_colums <- c("Species",
                     "Year",
                     "CatchDate",
                     "Gear",
                     "Area",
                     "SubArea",
                     "Coastal",
                     "N62Code",
                     "VesselLengthGroup",
                     "CountryVessel",
                     "LandingSite",
                     "CountryLanding",
                     "Usage",
                     "RoundWeightKilogram"
                     )
expect_equivalent(expected_colums, names(flatSL))
expect_true(is.numeric(flatSL$RoundWeightKilogram))
expect_true(is.numeric(flatSL$Year))
expect_true(is.character(flatSL$CountryVessel))
expect_true(length(flatSL$CatchDate) > 1 & "POSIXct" %in% class(flatSL$CatchDate))

context("test-stoxLanding missing values in aggColumns")
weightPre <- sum(flatSL$RoundWeightKilogram)
landingXML$Mottaker$Mottaksstasjon[2] <- NA
flatSL <- StoxLanding(landingXML)
expect_equal(sum(is.na(flatSL$LandingSite)), 1)
weightPost <- sum(flatSL$RoundWeightKilogram)
expect_equal(weightPre, weightPost)

context("test-stoxLanding is.StoxLandingData")
landingXML <- readXmlFile(system.file("testresources", "landing.xml", package="RstoxData"), stream = T)
flatSL <- StoxLanding(landingXML)
expect_true(is.StoxLandingData(flatSL))
expect_false(is.StoxLandingData(landingXML))

expect_false(is.LandingData(flatSL))
expect_true(is.LandingData(landingXML))

context("test-stoxLanding append columns")
#regular append
landingXML <- readXmlFile(system.file("testresources", "landing.xml", package="RstoxData"), stream = T)
flatSL <- StoxLanding(landingXML, appendColumns = "Lokasjon_kode")
expected_colums <- c("Species",
                     "Year",
                     "CatchDate",
                     "Gear",
                     "Area",
                     "SubArea",
                     "Coastal",
                     "N62Code",
                     "VesselLengthGroup",
                     "CountryVessel",
                     "LandingSite",
                     "CountryLanding",
                     "Usage",
                     "Lokasjon_kode",
                     "RoundWeightKilogram"
)
expect_equivalent(expected_colums, names(flatSL))

#append and rename
landingXML <- readXmlFile(system.file("testresources", "landing.xml", package="RstoxData"), stream = T)
flatSL <- StoxLanding(landingXML, appendColumns = "Lokasjon_kode", appendColumnsNames = c("Location"))
expected_colums <- c("Species",
                     "Year",
                     "CatchDate",
                     "Gear",
                     "Area",
                     "SubArea",
                     "Coastal",
                     "N62Code",
                     "VesselLengthGroup",
                     "CountryVessel",
                     "LandingSite",
                     "CountryLanding",
                     "Usage",
                     "Location",
                     "RoundWeightKilogram"
)
expect_equivalent(expected_colums, names(flatSL))

#append and non-existing
landingXML <- readXmlFile(system.file("testresources", "landing.xml", package="RstoxData"), stream = T)
flatSL <- StoxLanding(landingXML, appendColumns = "Non-exisiting", appendColumnsNames = c("Location"))
expected_colums <- c("Species",
                     "Year",
                     "CatchDate",
                     "Gear",
                     "Area",
                     "SubArea",
                     "Coastal",
                     "N62Code",
                     "VesselLengthGroup",
                     "CountryVessel",
                     "LandingSite",
                     "CountryLanding",
                     "Usage",
                     "Location",
                     "RoundWeightKilogram"
)
expect_equivalent(expected_colums, names(flatSL))