#' Parses landings (sales notes)
#' @description
#'  Parses sales notes data from the Norwegian Directorate of Fisheries (FDIR) on the LSS format
#' @details
#'  The LSS format is a pipe-separated format encoding landings (sales-notes).
#'  It is provided to IMR on a regular basis from FDIR.
#'  Column headers are in Norwegian.
#'
#'  Historically, columns in the landings provided from FDIR has been adapted for each data delivery
#'  Lately data deliveries has become standardized, but in order to support variants
#'  adherence to the standardization is not enforced by this function, unless option 'strict' is selected.
#'  If column names does not match specification, but data is otherwise parse-able, a warning will be issued.
#'  
#'  If 'strict' is not selection, data types may be inferred from data.
#'  The parameter 'guessMax' limits how many lines are inspected for data type inference (passed to \code{\link[readr]{read_delim}})
#' @param file path to file with LSS landings
#' @param encoding encoding for 'file'
#' @param guessMax passed to \code{\link[readr]{read_delim}}, unless 'strict' is true
#' @param strict enforce strict adherence to data format.
#' @return data.table with LSS landings
#' @importFrom data.table as.data.table
#' @importFrom readr default_locale read_delim col_character col_date col_datetime col_double col_integer col_number cols
#' @export
readLssFile <- function(file, encoding="latin1", guessMax = 100000, strict=T){
  
  spec_land <- cols(
    Dokumentnummer = col_character(),
    `Dokumenttype (kode)` = col_character(),
    Dokumenttype = col_character(),
    `Dokument versjonsnummer` = col_character(),
    `Dokument salgsdato` = col_character(),
    `Dokument versjonstidspunkt` = col_character(),
    `Salgslag ID` = col_character(),
    `Salgslag (kode)` = col_character(),
    Salgslag = col_character(),
    `Mottakernasjonalitet (kode)` = col_character(),
    Mottakernasjonalitet = col_character(),
    Mottaksstasjon = col_character(),
    `Landingskommune (kode)` = col_character(),
    Landingskommune = col_character(),
    `Landingsfylke (kode)` = col_character(),
    Landingsfylke = col_character(),
    `Landingsnasjon (kode)` = col_character(),
    Landingsnasjon = col_character(),
    Produksjonsanlegg = col_character(),
    `Produksjonskommune (kode)` = col_character(),
    Produksjonskommune = col_character(),
    `Fiskerkommune (kode)` = col_character(),
    Fiskerkommune = col_character(),
    `Fiskernasjonalitet (kode)` = col_character(),
    Fiskernasjonalitet = col_character(),
    Fartoynavn = col_character(),
    `Fartoy ID` = col_character(),
    `Registreringsmerke (seddel)` = col_character(),
    `Radiokallesignal (seddel)` = col_character(),
    `Storste lengde` = col_double(),
    `Lengdegruppe (kode)` = col_character(),
    Lengdegruppe = col_character(),
    `Bruttotonnasje 1969` = col_double(),
    `Bruttotonnasje annen` = col_double(),
    Byggear = col_integer(),
    Ombyggingsar = col_integer(),
    Motorkraft = col_double(),
    Motorbyggear = col_integer(),
    `Fartoy gjelder fra dato` = col_character(),
    `Fartoy gjelder til dato` = col_character(),
    `Fartoytype (kode)` = col_character(),
    Fartoytype = col_character(),
    `Kvotefartoy reg.merke` = col_character(),
    `Fartoykommune (kode)` = col_character(),
    Fartoykommune = col_character(),
    `Fartoyfylke (kode)` = col_character(),
    Fartoyfylke = col_character(),
    `Fartoynasjonalitet (kode)` = col_character(),
    Fartoynasjonalitet = col_character(),
    `Mottakende fartoy reg.merke` = col_character(),
    `Mottakende fartoy rkal` = col_character(),
    `Mottakende fartoytype (kode)` = col_character(),
    `Mottakende fart.type` = col_character(),
    `Mottakende fartoynasj. (kode)` = col_character(),
    `Mottakende fart.nasj` = col_character(),
    Fangstar = col_integer(),
    `Siste fangstdato` = col_date(format="%d.%m.%Y"),
    `Kvotetype (kode)` = col_character(),
    Kvotetype = col_character(),
    `Redskap (kode)` = col_character(),
    Redskap = col_character(),
    `Redskap - hovedgruppe (kode)` = col_character(),
    `Redskap - hovedgruppe` = col_character(),
    `Fangstfelt (kode)` = col_character(),
    `Kyst/hav (kode)` = col_character(),
    `Hovedomrade (kode)` = col_character(),
    Hovedomrade = col_character(),
    `Lokasjon (kode)` = col_character(),
    `Sone (kode)` = col_character(),
    Sone = col_character(),
    Omradegruppering = col_character(),
    `Hovedomrade FAO (kode)` = col_character(),
    `Hovedomrade FAO` = col_character(),
    `Nord/sor for 62 grader nord` = col_character(),
    `Fangstdagbok (nummer)` = col_character(),
    `Fangstdagbok (turnummer)` = col_character(),
    Landingsdato = col_character(),
    Landingsklokkeslett = col_character(),
    `Dellanding (signal)` = col_character(),
    `Neste mottaksstasjon` = col_character(),
    `Forrige mottakstasjon` = col_character(),
    Linjenummer = col_integer(),
    `Art - FDIR (kode)` = col_character(),
    `Art - FDIR` = col_character(),
    `Art - gruppe (kode)` = col_character(),
    `Art - gruppe` = col_character(),
    `Art - hovedgruppe (kode)` = col_character(),
    `Art - hovedgruppe` = col_character(),
    `Art FAO (kode)` = col_character(),
    `Art FAO` = col_character(),
    `Produkttilstand (kode)` = col_character(),
    Produkttilstand = col_character(),
    `Konserveringsmate (kode)` = col_character(),
    Konserveringsmate = col_character(),
    `Landingsmate (kode)` = col_character(),
    Landingsmate = col_character(),
    `Kvalitet (kode)` = col_character(),
    Kvalitet = col_character(),
    `Storrelsesgruppering (kode)` = col_character(),
    `Anvendelse (kode)` = col_character(),
    Anvendelse = col_character(),
    `Anvendelse hovedgruppe (kode)` = col_character(),
    `Anvendelse hovedgruppe` = col_character(),
    `Antall stykk` = col_integer(),
    Bruttovekt = col_double(),
    Produktvekt = col_double(),
    Rundvekt = col_double()
  )
  names(spec_land$cols)[26] <- "Fart\u00F8ynavn"
  names(spec_land$cols)[27] <- "Fart\u00F8y ID"
  names(spec_land$cols)[30] <- "St\u00F8rste lengde"
  names(spec_land$cols)[35] <- "Bygge\u00E5r"
  names(spec_land$cols)[36] <- "Ombyggings\u00E5r"
  names(spec_land$cols)[38] <- "Motorbygge\u00E5r"
  names(spec_land$cols)[39] <- "Fart\u00F8y gjelder fra dato"
  names(spec_land$cols)[40] <- "Fart\u00F8y gjelder til dato"
  names(spec_land$cols)[41] <- "Fart\u00F8ytype (kode)"
  names(spec_land$cols)[42] <- "Fart\u00F8ytype"
  names(spec_land$cols)[43] <- "Kvotefart\u00F8y reg.merke"
  names(spec_land$cols)[44] <- "Fart\u00F8ykommune (kode)"
  names(spec_land$cols)[45] <- "Fart\u00F8ykommune"
  names(spec_land$cols)[46] <- "Fart\u00F8yfylke (kode)"
  names(spec_land$cols)[47] <- "Fart\u00F8yfylke"
  names(spec_land$cols)[48] <- "Fart\u00F8ynasjonalitet (kode)"
  names(spec_land$cols)[49] <- "Fart\u00F8ynasjonalitet"
  names(spec_land$cols)[50] <- "Mottakende fart\u00F8y reg.merke"
  names(spec_land$cols)[51] <- "Mottakende fart\u00F8y rkal"
  names(spec_land$cols)[52] <- "Mottakende fart\u00F8ytype (kode)"
  names(spec_land$cols)[54] <- "Mottakende fart\u00F8ynasj. (kode)"
  names(spec_land$cols)[56] <- "Fangst\u00E5r"
  names(spec_land$cols)[66] <- "Hovedomr\u00E5de (kode)"
  names(spec_land$cols)[67] <- "Hovedomr\u00E5de"
  names(spec_land$cols)[71] <- "Omr\u00E5degruppering"
  names(spec_land$cols)[72] <- "Hovedomr\u00E5de FAO (kode)"
  names(spec_land$cols)[73] <- "Hovedomr\u00E5de FAO"
  names(spec_land$cols)[74] <- "Nord/s\u00F8r for 62 grader nord"
  names(spec_land$cols)[93] <- "Konserveringsm\u00E5te (kode)"
  names(spec_land$cols)[94] <- "Konserveringsm\u00E5te"
  names(spec_land$cols)[95] <- "Landingsm\u00E5te (kode)"
  names(spec_land$cols)[96] <- "Landingsm\u00E5te"
  names(spec_land$cols)[99] <- "St\u00F8rrelsesgruppering (kode)"
  
  loc <- readr::default_locale()
  loc$decimal_mark <- ","
  loc$encoding <- encoding
  if (strict){
    headers <- names(readr::read_delim(file, delim="|", col_names=T, col_types=paste(rep("c",107), collapse=""), trim_ws=TRUE, na=c("", "na", "NA"), locale=loc, n_max = 1))
    
    if (length(headers) != length(spec_land$cols)){
      stop("Number of columns in file does not match specification.")
    }
    if (!all(headers == names(spec_land$cols))){
      differences <- sum(headers != names(spec_land$cols))
      warning(paste("StoX: Header names does not match specification,", differences, "column names differ."))
    }
      
    db <- readr::read_delim(file, delim="|", col_names=names(spec_land$cols), trim_ws=TRUE, na=c("", "na", "NA"), locale=loc, col_types = spec_land, skip = 1) 
    db <- data.table::as.data.table(db) 
    db$`Siste fangstdato` <- as.POSIXct(db$`Siste fangstdato`)
  }
  else{
    db <- readr::read_delim(file, delim="|", col_names=T, trim_ws=TRUE, na=c("", "na", "NA"), locale=loc, guess_max = guessMax)    
    db <- data.table::as.data.table(db) 
  }
  return(db)
}


#' Read pipe separated file with specified columns
#' @noRd
read_psv <- function(file, encoding, col_types){
  loc <- readr::default_locale()
  loc$decimal_mark <- ","
  loc$encoding <- encoding
  db <- readr::read_delim(file, delim="|", col_names=T, trim_ws=TRUE, na=c("", "na", "NA"),locale=loc, col_types = col_types)
  return(db)
}

#' Parses logbooks (ERS) 
#' @description 
#'  Parses electronic logbooks (ERS) from tabular format delivered by Directorate of Fisheries (FDIR)
#' @details 
#'  The format is a pipe-separated format encoding aggregated ERS records (logbooks).
#'  It is provided to IMR on a regular basis from FDIR.
#'  Column headers are in Norwegian.
#' @param file path to file
#' @param encoding encoding for 'file'
#' @return data.table() with logbooks
#' @export
readErsFile <- function(file, encoding="latin1"){
  
  spec_log <- cols(
    RC = col_character(),
    REGM = col_character(),
    STORSTE_LENGDE = col_double(),
    BRUTTOTONNASJE = col_integer(),
    MOTORKRAFT = col_integer(),
    TM1 = col_character(),
    AKTIVITET_KODE = col_character(),
    AKTIVITET = col_character(),
    PUMPET_FRA = col_character(),
    FANGSTAR = col_integer(),
    STARTTIDSPUNKT = col_datetime(format = "%Y-%m-%d %H:%M:%S"),
    START_LT = col_double(),
    START_LG = col_double(),
    SONE = col_character(),
    KVOTETYPE_KODE = col_character(),
    KVOTETYPE = col_character(),
    REDSKAP_FAO = col_character(),
    REDSKAP_NS = col_character(),
    REDSKAP = col_character(),
    REDSKAPSSPESIFIKASJON_KODE = col_character(),
    REDSKAPSSPESIFIKASJON = col_character(),
    MASKEVIDDE = col_integer(),
    REDSKAP_PROBLEMER_KODE = col_character(),
    REDSKAP_PROBLEMER = col_character(),
    STOPPTIDSPUNKT = col_datetime(format = "%Y-%m-%d %H:%M:%S"),
    STOPP_LT = col_double(),
    STOPP_LG = col_double(),
    VARIGHET = col_integer(),
    INNSATS = col_number(),
    SILD_BESTAND_KODE = col_character(),
    SILD_BESTAND_NS = col_character(),
    SILD_BESTAND = col_character(),
    HOVEDART_FAO = col_character(),
    HOVEDART_NS = col_character(),
    HOVEDART = col_character(),
    INT_OMR_GML_START = col_character(),
    INT_OMR_NY_START = col_character(),
    INT_OMR_GML_STOPP = col_character(),
    INT_OMR_NY_STOPP = col_character(),
    HAV_DYBDE_START = col_number(),
    HAV_DYBDE_STOPP = col_number(),
    LOKASJON_START = col_character(),
    LOKASJON_STOPP = col_character(),
    TREKK_AVSTAND_METER = col_integer(),
    FANGSTART_FAO = col_character(),
    FANGSTART_NS = col_character(),
    FANGSTART = col_character(),
    RUNDVEKT = col_double()
  )
  names(spec_log$cols) <- c(names(spec_log$cols)[1:2], "ST\u00D8RSTE_LENGDE", names(spec_log$cols)[4:9], "FANGST\u00C5R", names(spec_log$cols)[11:length(spec_log$cols)])
  
  logb <- read_psv(file, encoding, col_types=spec_log)
  
  return(data.table::as.data.table(logb))
}
