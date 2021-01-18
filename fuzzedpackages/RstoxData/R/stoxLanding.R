#' Check if argument is LandingData
#' @description 
#'  Checks if argument conforms to specification for \code{\link[RstoxData]{LandingData}}
#' @param LandingData argument to be checked for data conformity
#' @return logical, TRUE if argument conformed to specification for \code{\link[RstoxData]{LandingData}}
#' @name is.LandingData
#' @export
is.LandingData <- function(LandingData){

  if (!is.list(LandingData)){
    return(FALSE)
  }
  if (!("Seddellinje" %in% names(LandingData))){
    return(FALSE)
  }
  if (!data.table::is.data.table(LandingData$Seddellinje)){
    return(FALSE)
  }
  if (!all(c("Dokumentnummer", 
             "Linjenummer", 
             "Art_kode", 
             "Registreringsmerke_seddel", 
             "SisteFangstdato", 
             "Redskap_kode")
              %in% names(LandingData$Seddellinje))){
    return(FALSE)
  }
  
  return(TRUE)
}

#' Extracts aggregated Landings from NMD - landings (namespace: http://www.imr.no/formats/landinger/v2)
#' @noRd
extractNMDlandingsV2 <- function(LandingData, appendColumns=character(), appendColumnsNames=appendColumns){
  flatLandings <- LandingData$Seddellinje
  for (part in names(LandingData)[!(names(LandingData) %in% c("Landingsdata", "Seddellinje", "metadata"))]){
    keys <- names(LandingData[[part]])[names(LandingData[[part]]) %in% names(flatLandings)]
    flatLandings <- merge(LandingData[[part]], flatLandings, by=keys)
  }
  
  #
  # Note: if non-character columns are added to aggColumns. Handle accoridngly in NA-aggregation below
  #
  sourceColumns <- c("Art_kode", 
                     "Fangst\u00E5r", 
                     "SisteFangstdato", 
                     "Redskap_kode", 
                     "Hovedomr\u00E5de_kode",
                     "Lokasjon_kode",
                     "KystHav_kode", 
                     "NordS\u00F8rFor62GraderNord", 
                     "Lengdegruppe_kode", 
                     "Fart\u00F8ynasjonalitet_kode",
                     "Mottaksstasjon",
                     "Mottakernasjonalitet_kode",
                     "HovedgruppeAnvendelse_kode")
  sourceColumns <- c(sourceColumns, appendColumns)
  
  outputColumns <- c( "Species",
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
                      "Usage")
  outputColumns <- c(outputColumns, appendColumnsNames)
  
  # add NAs for missing columns
  # this is done because the underlaying format is supposed to be generalized to hetergenous formats like StoxBiotic.
  for (col in appendColumns[!(appendColumns %in% names(flatLandings))]){
    flatLandings[[col]] <- NA
  }
  
  flatLandings <- flatLandings[,c(sourceColumns, "Rundvekt"), with=F]
  
  aggList <- list()
  for (i in 1:length(sourceColumns)){
    sourcename <- sourceColumns[i]
    outputname <- outputColumns[i]
    flatLandings[[outputname]] <- flatLandings[[sourcename]]
    flatLandings[[outputname]][is.na(flatLandings[[sourcename]])] <- "<NA>" #set NAs to text-string for aggregation
    aggList[[outputname]] <- flatLandings[[outputname]]
  }
  names(aggList) <- outputColumns
  aggLandings <- stats::aggregate(list(Rundvekt=flatLandings$Rundvekt), by=aggList, FUN=function(x){sum(x, na.rm=T)})
  aggLandings <- aggLandings[,c(outputColumns, "Rundvekt")]
  
  #reset NAs
  for (aggC in outputColumns){
    aggLandings[[aggC]][aggLandings[[aggC]] == "<NA>"] <- NA
  }
  aggLandings$RoundWeightKilogram <- aggLandings$Rundvekt
  outputColumns <- c(outputColumns, "RoundWeightKilogram")
  
  # format conversions
  cd <- as.POSIXct(aggLandings$CatchDate, format="%d.%m.%Y")
  attributes(cd)$tzone <- "UTC"
  aggLandings$CatchDate <- as.POSIXct(substr(as.character(cd),1,10), format="%Y-%m-%d", tzone="UTC")
  
  aggLandings$Year <- as.integer(aggLandings$Year)
  
  return(data.table::as.data.table(aggLandings[,outputColumns]))
}

#' Convert landing data
#' @description
#'  StoX function
#'  Convert landing data to the aggregated format \code{\link[RstoxData]{StoxLandingData}}
#' @details 
#'  All columns that are not the ones aggregated (weight), will be used as aggregation variables.
#'  This includes any columns added with 'appendColumns' and may not make much sense for continuous variables.
#'  
#'  If 'LandingData' does not contain columns identified by 'appendColumns'. NA columns will be added.
#' 
#'  Correspondences indicate which field a value is derived from, not necessarily verbatim copied.
#' 
#'  Correspondence to LandingData (http://www.imr.no/formats/landinger/v2):
#'  \describe{
#'   \item{Species}{Art_kode}
#'   \item{Year}{Fangstår}
#'   \item{CatchDate}{SisteFangstdato}
#'   \item{Gear}{Redskap_kode}
#'   \item{Area}{Hovedområde_kode}
#'   \item{SubArea}{Lokasjon_kode}
#'   \item{Coastal}{KystHav_kode}
#'   \item{N62Code}{NordSørFor62GraderNord}
#'   \item{VesselLengthGroup}{Lengdegruppe_kode}
#'   \item{CountryVessel}{Fartøynasjonalitet_kode}
#'   \item{LandingSite}{Mottaksstasjon}
#'   \item{CountryLanding}{Landingsnasjon_kode}
#'   \item{Usage}{HovedgruppeAnvendelse_kode}
#'   \item{RoundWeightKilogram}{Rundvekt}
#'  }
#'  
#' @param LandingData Sales-notes data. See \code{\link[RstoxData]{LandingData}}
#' @param appendColumns character() vector that identifies additional columns in \code{\link[RstoxData]{LandingData}} to append to \code{\link[RstoxData]{StoxLandingData}}.
#' @param appendColumnsNames character() vector that defines the names of the columns in 'appendColumns' in the output.
#' @return \code{\link[RstoxData]{StoxLandingData}}, aggregated landings data.
#' @name StoxLanding
#' @export
StoxLanding <- function(LandingData, appendColumns=character(), appendColumnsNames=appendColumns){
  
  if (length(appendColumns) != length(appendColumnsNames)){
    stop("elements in appendColumnNames must correspond to elements in appendColumns")
  }
  
  return(extractNMDlandingsV2(LandingData, appendColumns, appendColumnsNames))
  
}

#' Check if argument is StoxLandingData
#' @description 
#'  Checks if argument conforms to specification for \code{\link[RstoxData]{StoxLandingData}}
#' @param StoxLandingData argument to be checked for data conformity
#' @return logical, TRUE if argument conformed to specification for \code{\link[RstoxData]{StoxLandingData}}
#' @name is.StoxLandingData
#' @export
is.StoxLandingData <- function(StoxLandingData){
  
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
  
  if (!data.table::is.data.table(StoxLandingData)){
    return(FALSE)
  }
  
  if (!all(expected_colums %in% names(StoxLandingData))){
    return(FALSE)
  }
  
  return(TRUE)
}
