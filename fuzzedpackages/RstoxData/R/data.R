#' @title xsdObjects
#' @description Pre-processed XSD file objects
#' @format A list with 4 elements
#' \describe{
#'   \item{\code{landingerv2.xsd}}{List Landing Format v2}
#'   \item{\code{nmdbioticv1.xsd}}{List NMD Biotic Format v1}
#'   \item{\code{nmdbioticv1.1.xsd}}{List NMD Biotic Format v1.1}
#'   \item{\code{nmdbioticv1.2.xsd}}{List NMD Biotic Format v1.2}
#'   \item{\code{nmdbioticv1.3.xsd}}{List NMD Biotic Format v1.3}
#'   \item{\code{nmdbioticv1.4.xsd}}{List NMD Biotic Format v1.4}
#'   \item{\code{nmdbioticv3.xsd}}{List NMD Biotic Format v3}
#'   \item{\code{nmdbioticv3.1.xsd}}{List NMD Biotic Format v3.1}
#'   \item{\code{nmdechosounderv1.xsd}}{List NMD Echosounder Format v1}
#' }
#' @source \url{https://www.imr.no/formats}
"xsdObjects"

#' @title stoxBioticObject
#' @description Pre-processed objects for raw XML data to StoXBiotic format
"stoxBioticObject"


##################################################
##################################################
#' General parameters of RstoxData.
#' 
#' All functions referring to a project, a model, a process or an output table use the same parameters, listed here.
#' 
#' @param processData The current data produced by a previous instance of the function.
#' @param UseProcessData Logical: If TRUE use the existing function output in the process. 
#' @param NumberOfCores The number of cores to use (defaulted to 1), truncated to the number of avaliable cores.
#' 
#' @name general_arguments
#' 
NULL


##################################################
##################################################
#' StoX data types of the RstoxData package
#' 
#' StoX data types are the data types used to transfer data and information between processes in a StoX estimation model. The data types are divided into two types, the \code{\link{ModelData}} and \code{\link{ProcessData}}.
#' 
#' @name DataTypes
#' 
NULL

##################################################
##################################################
#' StoX data types of the RstoxData package
#' 
#' StoX data types are the data types used to transfer data and information between processes in a StoX estimation model.
#' 
#' @details
#' This RstoxData package produces the following StoX data types:
#' \itemize{
#' \item{\code{\link{BioticData}}}
#' \item{\code{\link{StoxBioticData}}}
#' \item{\code{\link{MergeStoxBioticData}}}
#' \item{\code{\link{AcousticData}}}
#' \item{\code{\link{StoxAcousticData}}}
#' \item{\code{\link{MergeStoxAcousticData}}}
#' \item{\code{\link{LandingData}}}
#' \item{\code{\link{StoxLandingData}}}
#' }
#' 
#' @param BioticData \code{\link{BioticData}}.
#' @param StoxBioticData \code{\link{StoxBioticData}}.
#' @param AcousticData \code{\link{AcousticData}}.
#' @param StoxAcousticData \code{\link{StoxAcousticData}}.
#'
#' @seealso \href{https://github.com/StoXProject/RstoxBase}{RstoxBase} and \href{https://github.com/StoXProject/RstoxFDA}{RstoxFDA} for a list of all StoX data types produced by the other official StoX function packages.
#' 
#' @name ModelData
#' 
NULL


##################################################
##################################################
#' Process data used in estimation models in StoX
#' 
#' The process data of the RstoxData package. 
#' 
#' @details
#' \itemize{
#' \item{\code{\link{StoxBioticTranslation}}}
#' }
#' 
#' @name ProcessData
#' 
#' @seealso \code{\link{ModelData}} for model data types and \code{\link{DataTypes}} for all data types produced by \code{\link{RstoxData}}.
#' 
NULL


##################################################
##################################################
#' StoX data type BioticData
#' 
#' Biotic data read from biotic xml files.
#' 
#' @details
#' This StoX data type is produced by \code{\link{ReadBiotic}}, and contains one list per input biotic file holding the tables read from each file, added a table named "metadata" holding the input file path and format. Currently supported are NMDBiotic1.4 (\url{https://www.imr.no/formats/nmdbiotic/v1.4/}), NMDBiotic3.0 (\url{https://www.imr.no/formats/nmdbiotic/v3/}), and ICESBiotic (\url{https://ices.dk/data/data-portals/Pages/acoustic.aspx}, click on "Acoustic data format" to download the format description).
#' 
#' @seealso \code{\link{DataTypes}} for a list of all StoX data types produced by \code{\link{RstoxData}}
#' 
#' @name BioticData
NULL


##################################################
##################################################
#' StoX data type StoxBioticData
#' 
#' Biotic data stored in the StoxBiotic format, which contains the variables needed for most estimation models used by StoX.
#' 
#' @details
#' This StoX data type is produced by \code{\link{StoxBiotic}}, and contains the tables Cruise, Station, Haul, SpeciesCategory, Sample and Individual in that hierarchical order.
#' 
#' @seealso \code{\link{DataTypes}} for a list of all StoX data types produced by \code{\link{RstoxData}}
#' 
#' @name StoxBioticData
#' 
NULL


##################################################
##################################################
#' StoX data type MergeStoxBioticData
#' 
#' Merged \code{\link{StoxBioticData}}.
#' 
#' @details
#' This StoX data type is produced by \code{\link{MergeStoxBiotic}}, and contains one merged table of \code{\link{StoxBioticData}}.
#' 
#' @seealso \code{\link{DataTypes}} for a list of all StoX data types produced by \code{\link{RstoxData}}
#' 
#' @name MergeStoxBioticData
#' 
NULL


##################################################
##################################################
#' StoX data type AcousticData
#' 
#' Biotic data read from biotic xml files.
#' 
#' @details
#' This StoX data type is produced by \code{\link{ReadAcoustic}}, and contains one list per input acoustic file holding the tables read from each file, added a table named "metadata" holding the input file path and format. Currently supported are NMDEchosounder1 (\url{https://www.imr.no/formats/nmdechosounder/v1/}), and ICESAcoustic (\url{https://ices.dk/data/data-portals/Pages/acoustic.aspx}, click on "Acoustic data format" to download the format description). 
#' 
#' @seealso \code{\link{DataTypes}} for a list of all StoX data types produced by \code{\link{RstoxData}}
#' 
#' @name AcousticData
#' 
NULL


##################################################
##################################################
#' StoX data type StoxAcousticData
#' 
#' Acoustic data stored in the StoxAcoustic format, which contains the variables needed for most estimation models used by StoX.
#' 
#' @details
#' This StoX data type is produced by \code{\link{StoxAcoustic}}, and contains the tables Cruise, Log, Beam, AcousticCategory, ChannelReference and NASC in that hierarchical order.
#' 
#' @seealso \code{\link{DataTypes}} for a list of all StoX data types produced by \code{\link{RstoxData}}
#' 
#' @name StoxAcousticData
#' 
NULL


##################################################
##################################################
#' StoX data type MergeStoxAcousticData
#' 
#' Merged \code{\link{StoxAcousticData}}.
#' 
#' @details
#' This StoX data type is produced by \code{\link{MergeStoxAcoustic}}, and contains one merged table of \code{\link{StoxAcousticData}}.
#' 
#' @seealso \code{\link{DataTypes}} for a list of all StoX data types produced by \code{\link{RstoxData}}
#' 
#' @name MergeStoxAcousticData
#' 
NULL


#' LandingData
#' 
#' @section Data:
#' One entry 'Seddellinje' is one line of a sales-note or landing-note. 
#' These are issued as fish is landed, and a complete set of these for a period
#' can be considered a census of all first hand sale of fish sold from Norwegian vessels.
#' 
#' @section Format:
#' list() of \code{\link[data.table]{data.table}} 
#' representing the different complexTypes in namespace http://www.imr.no/formats/landinger/v2
#' For ease of merging: all top level attributes are repeated for all tables. And all line-identifying variables are included as top-level attributes.
#' 
#' @seealso \code{\link{DataTypes}} for a list of all StoX data types produced by \code{\link{RstoxData}}
#' 
#' @name LandingData
#' 
NULL


#' StoxLandingData
#'
#' Table (\code{\link[data.table]{data.table}}) with aggregated weight of landings from landing records.
#'
#' @section Column definitions:
#'  \describe{
#'   \item{Species}{character() code for species category (species identified by market or regulation standards. Several codes may code the same species or stock, and some catch may be recorded only at higher taxonomic classifications)}
#'   \item{Year}{integer() Year of catch}
#'   \item{CatchDate}{POSIXct() Date of catch (last catch on trip) in UTC}
#'   \item{Gear}{character() Code for gear used for catch (dominant gear for trip)}
#'   \item{Area}{character() Area code for the position where the catch was caught (dominant area for trip)}
#'   \item{SubArea}{character() Subdivision of area code for the position where the catch was caught (dominant area for trip)}
#'   \item{Coastal}{character() Code indicating whether catch was taken within coastal delimitation line (dominant side for trip)}
#'   \item{N62Code}{character() Code indicating whether catch was taken north or south of 62 deg. Lat. (dominant side for trip)}
#'   \item{VesselLengthGroup}{character() Length group for vessel}
#'   \item{CountryVessel}{character() Country of the vessel that caught the catch}
#'   \item{LandingSite}{character() Code identifying landing site (buyer of catch)}
#'   \item{CountryLanding}{character() Country where catch was landed}
#'   \item{Usage}{character() Code for market usage of catch.}
#'   \item{RoundWeightKilogram}{numeric() Weight of round catch in kg.}
#'  }
#'
#' @seealso \code{\link{DataTypes}} for a list of all StoX data types produced by \code{\link{RstoxData}}
#' 
#' @name StoxLandingData
#'
NULL


##################################################
##################################################
#' Translation definition (from file) for \code{\link{StoxBioticData}}.
#' 
#' @details
#' This StoX data type is produced by \code{\link{DefineStoxBioticTranslation}}, and contains the columns VariableName, Value and NewValue.
#' 
#' @seealso \code{\link{DataTypes}} for a list of all StoX data types produced by \code{\link{RstoxData}}
#' 
#' @name StoxBioticTranslation
#' 
NULL


