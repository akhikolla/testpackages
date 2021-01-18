#' Tools to Read and Manipulate Fisheries Data
#'
#' Set of tools to read and manipulate various data formats for fisheries. Mainly catered towards scientific trawl survey sampling ('biotic') data, acoustic trawl data, and commercial fishing catch ('landings') data. Among the supported data formats are the data products from the Norwegian Institute Marine Research ('IMR') and the International Council for the Exploration of the Sea (ICES).
#'
#' The RstoxData package contains functions for reading, filtering and writing biotic, acoustic and landing data as XML files. Filtering can be done by R syntax such as longitude > 10, or by pre defined functions such as inside(). On computers that return errors when trying to run the Rtools through RStudio (most institutional Windows machines), install the binary directly from https://github.com/StoXProject/RstoxData/releases. Download the newest RstoxData zip file, click the "Packages" tab -> "Install" -> "Install from:" "Package Archive File" -> "Install". If the installer does not complain, the package is installed correctly.
#' @docType package
#' @name RstoxData
#'
"_PACKAGE"

# Global variables
utils::globalVariables(c(
	 "RstoxDataEnv", "xsdObjects", ".", "..allDuplicated", "..colAgg", "..colList", "..columns",
	 "..digits", "..keep", "..key", "..signifDigits", "..sourceColumns",
	 "..targetAndSourceVariables", "..varToExtract", "..x", "AcousticCategory", "Addition", "age",
	 "agingstructure", "ap", "aphia", "BeamKey", "bottomdepthstart", "bottomdepthstop",
	 "catCatchWgt", "catchcount", "catchpartnumber", "catchproducttype",
	 "CatchSpeciesCategoryNumber", "CatchSpeciesCategoryWeight", "CatchSpeciesCode",
	 "CatchSubsampledNumber", "CatchSubsampleWeight", "catchweight", "cc", "Constant", "Country",
	 "cruise", "Cruise", "CruiseKey", "cw", "DateTime", "direction", "DoorType", "EchoType", "EDSU",
	 "FishID", "fishingdepthmax", "fishingdepthmin", "freq", "g", "gear", "Gear", "gearcondition",
	 "GearExp", "gearflow", "HaulNo", "HaulNumber", "HaulVal", "HaulValidity", "hv",
	 "inapplicableFormats", "individualweight", "isCrustacean", "isHerringOrSprat",
	 "isHerringOrSpratOrMackerel", "latitudeend", "latitudestart", "LengthClass", "LengthCode",
	 "lengthmeasurement", "lengthsamplecount", "lengthsampleweight", "lenInterval", "level",
	 "lngtClass", "lngtCode", "LocalID", "LogDuration", "LogKey", "LogOrigin", "longitudeend",
	 "longitudestart", "lsc", "lsCountTot", "maturationstage", "maturity", "meanW",
	 "MiddleDateTime", "missionstartdate", "missionstopdate", "ms", "N", "nation", "nInd", "noMeas",
	 "NumberAtLength", "nWithWeight", "parasite", "platformname", "preferredagereading", "Quarter",
	 "readability", "ReplaceBy", "reportInMM", "res", "rowIndex", "s", "SaCategory", "sampleFac",
	 "samplequality", "sampletype", "Scaling", "serialnumber", "sex", "Ship", "sp", "specialstage",
	 "specimenid", "SpecVal", "start_time", "StartDateTime", "startyear", "station",
	 "stationstartdate", "stationstarttime", "stationstopdate", "stationstoptime", "stationtype",
	 "StatRec", "stomach", "StopDateTime", "stoxBioticObject", "subFactor", "SubsampledNumber",
	 "subWeight", "suffixes", "Survey", "SweepLngt", "target", "Time", "tissuesample", "totalNo",
	 "totWeight", "transceiver", "trawldoorarea", "trawldoorspread", "trawldoortype",
	 "trawldoorweight", "VariableName", "verticaltrawlopening", "WeightMeasurement",
	 "winddirection", "windspeed", "wingspread", "wiredensity", "wirediameter", "wirelength"))

.onLoad <- function(libname, pkgname) {
	# Initiate the RstoxData environment:
	initiateRstoxData()
} 

# Try to unload dynamic library
.onUnload <- function (libpath) {
	library.dynam.unload("RstoxData", libpath)
} 

## usethis namespace: start
#' @useDynLib RstoxData, .registration = TRUE
## usethis namespace: end
NULL

## usethis namespace: start
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL

