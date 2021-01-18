# Function to compare with ices vocabulary of allowed values
compareICES <- function(url, field) {
  pg <- tryCatch(
        {
         read_xml(url)
        },
        error = function(e){
                warning(paste("StoX: Url", url, "is not exist or no internet connection available."))
                return(NA)
        }
  )
  recs <- xml_find_all(pg, "//Key")
  vals <- trimws(xml_text(recs))
  for(x in field){
    if(!x %in% vals){
      warning(paste0("StoX: ", x, " not defined in ", url))
    }
  }
}     

# Get quarter representation from a date
getQuarter <- function(stationstartdate) {
    x <- format(as.Date(stationstartdate, format="%Y-%m-%dZ"), "%m")
    return(floor((as.numeric(x) - 1) / 3 + 1))
}

# Get maturity indicator for a species
getDATRASMaturity <- function(quarter, aphia, specialstage, maturationstage) {

          temp <-  as.data.table(cbind(q=quarter, ap=aphia, sp=specialstage, ms=maturationstage, res=NA))

          temp[, `:=`(sp = as.numeric(sp), ms = as.numeric(ms), res = as.numeric(res) )]

          temp[, isHerringOrSpratOrMackerel := ifelse(ap %in% c("126417", "126425", "127023"), TRUE, FALSE)]

          temp[!is.na(sp) & isHerringOrSpratOrMackerel == TRUE,  res := ifelse(sp <= 2, 61, ifelse(sp <= 5, 62, 60 + sp - 3))]
          temp[!is.na(sp) & isHerringOrSpratOrMackerel == FALSE, res := 60 + sp]

          temp[is.na(sp) & !is.na(ms), res := ifelse(ms == 5 & q == "3", NA, 60 + ms)]

          return(temp$res)
}

# Convert gear number to sweep length
getGOVSweepByEquipment <- function(gear) {
    cnvTbl <- c("3120" = NA,
                "3190" = 60,
                "3191" = 60,
                "3192" = 60,
                "3193" = 110,
                "3194" = 110,
                "3195" = 110,
                "3196" = 60,
                "3197" = 110,
                "3198" = 60,
                "3199" = 60)

    x <- cnvTbl[as.character(gear)]
    x[is.null(x)] <- NA
    return(x)
}

# Get Haul validity
getHaulVal <- function(gearcondition, samplequality) {
    temp <-  as.data.table(cbind(g=gearcondition, s=samplequality))
    temp[, res:="I"]
    temp[(is.na(g) | g %in% c("1","2")) &
        (is.na(s) | s %in% c("0", "1")), res:="V"]

    return(temp$res)
}

# Generate ICES rectangle from a coordinate
# Stolen from: https://github.com/cran/mapplots/blob/master/R/ices.rect.R
getICESrect <- function(lat, lng){
    x <- floor(lng+60)+1000
    y <- floor(lat*2)-71+100
    num1<- substr(y,2,3)
    lett <- LETTERS[as.numeric(substr(x,2,3))]
    num2 <- substr(x,4,4)
    paste(num1,lett,num2,sep='')
}

# Get distance in meters between two coordinates
getDistanceMeter <- function(lat1, lon1, lat2, lon2) {
    x <-  acos( sin(lat1*pi/180)*sin(lat2*pi/180) + cos(lat1*pi/180)*cos(lat2*pi/180)*cos(lon2*pi/180-lon1*pi/180) ) * 6371000
    return(x)
}

# Calculate time diff
getTimeDiff <- function(stationstartdate, stationstarttime, stationstopdate, stationstoptime) {

    t0 <- ifelse(is.na(stationstartdate) | is.na(stationstarttime), NA, gsub("Z", " ", paste0(stationstartdate, stationstarttime)))
    t1 <- ifelse(is.na(stationstopdate) | is.na(stationstoptime), NA, gsub("Z", " ", paste0(stationstopdate, stationstoptime)))
    
    start <- as.POSIXct(t0)
    end <- as.POSIXct(t1)

    return(round(difftime(end, start, units = "mins")))
}

# Get ICES ship data
#' @importFrom xml2 xml_ns_strip xml_find_all xml_text
getICESShipCode <- function(platformname) {

    construct <- function(shipName) {
        # We have to remove "."," " and use uppercase
        shipName <- toupper(gsub("[[:space:][:punct:]]", "", shipName))

        # Replace the nordic character with AA
        shipName <- gsub("\u00C5", "AA", shipName)

        data <- tryCatch(
            {
                read_xml("https://vocab.ices.dk/services/pox/GetCodeList/SHIPC")
            },
                error = function(e){return(NA)}
        )

        # Can't download from ICES
        if (is.na(data))
            return(NA)

        xml_ns_strip(data)
        nodes <- xml_find_all(data, paste0("//Code[contains(translate(Description[normalize-space()],'abcdefghijklmnopqrstuvwxyz. ','ABCDEFGHIJKLMNOPQRSTUVWXYZ'), \"", shipName, "\")]/Key"))

        # Ship not found
        if (length(nodes) < 1) {
          return(NA)
        }

        # Get the latest matching ships
        shipCode <- xml_text(tail(nodes,1))

        return(shipCode)
    }

    nm <- unique(platformname)
    y <- unlist(lapply(nm, construct))
    names(y) <- nm

    x <- y[as.character(platformname)]
    x[is.null(x)] <- NA

    return(x)
}


#' @importFrom utils write.table
exportCSV <- function(filename, data, combine = FALSE, overwrite = FALSE, na = "") {

  ovw <- function(fn, ow) {
    if (file.exists(fn)) {
      if (ow) {
        file.remove(fn)
      } else {
        file.rename(fn, paste0(fn, ".old"))
      }
    }
  }

  message("Resulting CSV file is saved as:")

  if (combine) {
    ovw(filename, overwrite)
    suppressWarnings(lapply(data, write.table, file = filename, append = TRUE, na = na, row.names = FALSE, quote = FALSE, sep = ","))
    message(filename)
  } else {
    subname <- names(data)

    if (length(subname) < length(data)) {
      subname <- seq_len(length(data))
    }
    for (i in seq_len(length(data))) {
      subfilename <- paste0(tools::file_path_sans_ext(filename), "_", subname[i], ".", tools::file_ext(filename))
      ovw(subfilename, overwrite)
      suppressWarnings(write.table(data[[i]], file = subfilename, append = FALSE, na = na, row.names = FALSE, quote = FALSE, sep = ","))
      message(subfilename)
    }
  }
  return(TRUE)
}





#' Write ICES acoustic CSV format file
#'
#' Given an \code{AcousticData} object, this function will write an ICES acoustic CSV file. Note that this function only supports
#' \code{AcousticData} object that is created from reading an ICES acoustic XML file.
#'
#' @param AcousticData A \code{AcousticData} object from an ICES acoustic XML format file.
#' @param Combine Logical: If TRUE stack the output tables for each acosutic file.
#'
#' @return List of data.table objects in the ICES acoustic CSV format.
#'
#' @export
WriteICESAcoustic <- function(
	AcousticData, 
	Combine = FALSE
){
	
	out <- lapply(
		AcousticData, 
		prepareICESAcoustic, 
		Combine = Combine
	)
	
	return(out)
}


prepareICESAcoustic <- function(
	AcousticDataOne, 
	Combine = FALSE
){
	
		if(AcousticDataOne$metadata$useXsd=='icesAcoustic'){
			
			#Remove echo type
			AcousticDataOne$Data$EchoType<-NULL      
			
			#Fix notation of metadata
			translate <- function(xx) {
				res <- AcousticDataOne$vocabulary$value[match(xx, AcousticDataOne$vocabulary$id)]
				# Fix NAs
				res[which(is.na(res))] <- xx[which(is.na(res))]
				return(res)
			}
			
			AcousticDataOne$Instrument$TransducerLocation <- translate(AcousticDataOne$Instrument$TransducerLocation)
			AcousticDataOne$Instrument$TransducerBeamType <- translate(AcousticDataOne$Instrument$TransducerBeamType)
			AcousticDataOne$Calibration$AcquisitionMethod <- translate(AcousticDataOne$Calibration$AcquisitionMethod)
			AcousticDataOne$Calibration$ProcessingMethod <- translate(AcousticDataOne$Calibration$ProcessingMethod)
			AcousticDataOne$DataAcquisition$SoftwareName <- translate(AcousticDataOne$DataAcquisition$SoftwareName)
			AcousticDataOne$DataAcquisition$StoredDataFormat <- translate(AcousticDataOne$DataAcquisition$StoredDataFormat)
			AcousticDataOne$DataProcessing$SoftwareName <- translate(AcousticDataOne$DataProcessing$SoftwareName)
			AcousticDataOne$DataProcessing$TriwaveCorrection <- translate(AcousticDataOne$DataProcessing$TriwaveCorrection)
			AcousticDataOne$DataProcessing$OnAxisGainUnit <- translate(AcousticDataOne$DataProcessing$OnAxisGainUnit)
			AcousticDataOne$Cruise$Country <- translate(AcousticDataOne$Cruise$Country)
			AcousticDataOne$Cruise$Platform <- translate(AcousticDataOne$Cruise$Platform)
			AcousticDataOne$Cruise$Organisation <- translate(AcousticDataOne$Cruise$Organisation)
			# Take the last survey code
			AcousticDataOne$Survey$Code <- translate(tail(AcousticDataOne$Survey$Code, 1))
			AcousticDataOne$Log$Origin <- translate(AcousticDataOne$Log$Origin)
			AcousticDataOne$Log$Validity <- translate(AcousticDataOne$Log$Validity)
			AcousticDataOne$Sample$PingAxisIntervalType <- translate(AcousticDataOne$Sample$PingAxisIntervalType)
			AcousticDataOne$Sample$PingAxisIntervalUnit <- translate(AcousticDataOne$Sample$PingAxisIntervalUnit)
			AcousticDataOne$Sample$PingAxisIntervalOrigin <- translate(AcousticDataOne$Sample$PingAxisIntervalOrigin)
			AcousticDataOne$Data$SaCategory <- translate(AcousticDataOne$Data$SaCategory)
			AcousticDataOne$Data$SaCategory <- translate(AcousticDataOne$Data$SaCategory)
			AcousticDataOne$Data$Type <- translate(AcousticDataOne$Data$Type)
			AcousticDataOne$Data$Unit <- translate(AcousticDataOne$Data$Unit)
			
			
			#Check metadata towards ices definitions
			compareICES('https://acoustic.ices.dk/Services/Schema/XML/AC_TransducerLocation.xml',unique(AcousticDataOne$Instrument$TransducerLocation))
			compareICES('https://acoustic.ices.dk/Services/Schema/XML/AC_TransducerBeamType.xml',unique(AcousticDataOne$Instrument$TransducerBeamType))
			compareICES('https://acoustic.ices.dk/Services/Schema/XML/AC_AcquisitionMethod.xml',unique(AcousticDataOne$Calibration$AcquisitionMethod))
			compareICES('https://acoustic.ices.dk/Services/Schema/XML/AC_ProcessingMethod.xml',unique(AcousticDataOne$Calibration$ProcessingMethod))
			compareICES('https://acoustic.ices.dk/Services/Schema/XML/AC_DataAcquisitionSoftwareName.xml',unique(AcousticDataOne$DataAcquisition$SoftwareName))
			compareICES('https://acoustic.ices.dk/Services/Schema/XML/AC_StoredDataFormat.xml',unique(AcousticDataOne$DataAcquisition$StoredDataFormat))
			compareICES('https://acoustic.ices.dk/Services/Schema/XML/AC_DataAcquisitionSoftwareName.xml',unique(AcousticDataOne$DataAcquisition$SoftwareName))
			compareICES('https://acoustic.ices.dk/Services/Schema/XML/AC_DataProcessingSoftwareName.xml',unique(AcousticDataOne$DataProcessing$SoftwareName))
			compareICES('https://acoustic.ices.dk/Services/Schema/XML/AC_TriwaveCorrection.xml',unique(AcousticDataOne$DataProcessing$TriwaveCorrection))
			compareICES('https://acoustic.ices.dk/Services/Schema/XML/AC_OnAxisGainUnit.xml',unique(AcousticDataOne$DataProcessing$OnAxisGainUnit))
			compareICES('https://acoustic.ices.dk/Services/Schema/XML/ISO_3166.xml',unique(AcousticDataOne$Cruise$Country))
			compareICES('https://acoustic.ices.dk/Services/Schema/XML/SHIPC.xml',unique(AcousticDataOne$Cruise$Platform))
			compareICES('https://acoustic.ices.dk/Services/Schema/XML/EDMO.xml',unique(AcousticDataOne$Cruise$Organisation))
			compareICES('https://acoustic.ices.dk/Services/Schema/XML/AC_Survey.xml',unique(AcousticDataOne$Survey$Code))
			compareICES('https://acoustic.ices.dk/Services/Schema/XML/AC_LogOrigin.xml',unique(AcousticDataOne$Log$Origin))
			compareICES('https://acoustic.ices.dk/Services/Schema/XML/AC_LogValidity.xml',unique(AcousticDataOne$Log$Validity))
			compareICES('https://acoustic.ices.dk/Services/Schema/XML/AC_PingAxisIntervalType.xml',unique(AcousticDataOne$Sample$PingAxisIntervalType))
			compareICES('https://acoustic.ices.dk/Services/Schema/XML/AC_PingAxisIntervalUnit.xml',unique(AcousticDataOne$Sample$PingAxisIntervalUnit))
			compareICES('https://acoustic.ices.dk/Services/Schema/XML/AC_PingAxisIntervalOrigin.xml',unique(AcousticDataOne$Sample$PingAxisIntervalOrigin))
			compareICES('https://acoustic.ices.dk/Services/Schema/XML/AC_SaCategory.xml',unique(AcousticDataOne$Data$SaCategory))
			compareICES('https://acoustic.ices.dk/Services/Schema/XML/AC_AcousticDataType.xml',unique(AcousticDataOne$Data$Type))
			compareICES('https://acoustic.ices.dk/Services/Schema/XML/AC_DataUnit.xml',unique(AcousticDataOne$Data$Unit))
			
			
			
			# #Test to see if the log distance is always increasing
			# if(is.unsorted((AcousticDataOne$Log$Distance))){warning('The log distance in Log field is not chronological with time')}
			# if(is.unsorted((AcousticDataOne$Sample$Distance))){warning('The log distance in Sample field is not chronological with time')}
			# if(is.unsorted((AcousticDataOne$Data$Distance))){warning('The log distance in Data field is not chronological with time')}
			
			
			#Make the Instrument field
			hl<-c()
			hl$Instrument <- 'Instrument'
			hl$Header <- 'Record'
			#tmp <- data.frame(AcousticDataOne$Instrument)
			tmp <- data.table::data.table(AcousticDataOne$Instrument)
			names(tmp)<-paste0('Instrument',names(tmp))
			hl <- cbind(hl,tmp)
			#HInst <- format(hl, trim=TRUE, width=0)
			HInst <- hl
			
			
			#Make the Calibration field
			hl<-c()
			hl$Calibration <- 'Calibration'
			hl$Header <- 'Record'
			#tmp <- data.frame(AcousticDataOne$Calibration)
			tmp <- data.table::data.table(AcousticDataOne$Calibration)
			names(tmp)<-paste0('Calibration',names(tmp))
			hl <- cbind(hl,tmp)
			#HCal <- format(hl, trim=TRUE, width=0)
			HCal <- hl
			
			
			#Make the DataAcquisition field
			hl<-c()
			hl$DataAcquisition <- 'DataAcquisition'
			hl$Header <- 'Record'
			#tmp <- data.frame(AcousticDataOne$DataAcquisition)
			tmp <- data.table::data.table(AcousticDataOne$DataAcquisition)
			names(tmp)<-paste0('DataAcquisition',names(tmp))
			hl <- cbind(hl,tmp)
			#HDatA <- format(hl, trim=TRUE, width=0)
			HDatA <- hl
			
			
			#Make the DataProcessing Field
			hl<-c()
			hl$DataProcessing <- 'DataProcessing'
			hl$Header <- 'Record'
			#tmp <- data.frame(AcousticDataOne$DataProcessing)
			tmp <- data.table::data.table(AcousticDataOne$DataProcessing)
			names(tmp)<-paste0('DataProcessing',names(tmp))
			hl <- cbind(hl,tmp)
			#HDatP <- format(hl, trim=TRUE, width=0)
			HDatP <- hl
			
			
			#Make the Cruise Field
			hl<-c()
			hl$Cruise <- 'Cruise'
			hl$Header <- 'Record'
			ttt <- unique(AcousticDataOne$Survey$Code)
			hl$CruiseSurvey <- ttt[[1]]
			#tmp <- data.frame(AcousticDataOne$Cruise)
			tmp <- data.table::data.table(AcousticDataOne$Cruise)
			names(tmp)<-paste0('Cruise',names(tmp))
			hl <- cbind(hl,tmp)
			#HCru <- format(hl, trim=TRUE, width=0)
			HCru <- hl
			
			
			#Make the Data Field
			hl<-c()
			hl$Data <- 'Data'
			hl$Header <- 'Record'
			
			
			tmp_log <- unique(AcousticDataOne$Log)
			tmp_log <- lapply(tmp_log, unlist)
			# tmp_log$Origin<-unlist(tmp_log$Origin)
			# tmp_log$Validity<-unlist(tmp_log$Validity)
			names(tmp_log)[names(tmp_log)=='LocalID']<-'CruiseLocalID'
			names(tmp_log)[names(tmp_log)!='CruiseLocalID']<-paste0('Log',names(tmp_log)[names(tmp_log)!='CruiseLocalID'])
			
			tmp_sample <- lapply(AcousticDataOne$Sample, unlist)
			# tmp_sample$SamplePingAxisIntervalType<-unlist(tmp_sample$SamplePingAxisIntervalType)
			# tmp_sample$SamplePingAxisIntervalUnit<-unlist(tmp_sample$SamplePingAxisIntervalUnit)
			# tmp_sample$SamplePingAxisIntervalUnit<-unlist(tmp_sample$SamplePingAxisIntervalUnit)
			names(tmp_sample)[names(tmp_sample)=='LocalID']<-'CruiseLocalID'
			names(tmp_sample)[names(tmp_sample)=='Distance']<-'LogDistance'
			names(tmp_sample)[names(tmp_sample)=='Instrument']<-'InstrumentID'
			names(tmp_sample)[names(tmp_sample)=='Calibration']<-'CalibrationID'
			names(tmp_sample)[names(tmp_sample)=='DataAcquisition']<-'DataAcquisitionID'
			names(tmp_sample)[names(tmp_sample)=='DataProcessing']<-'DataProcessingID'
			names(tmp_sample)[!names(tmp_sample)%in%c('CruiseLocalID','LogDistance','InstrumentID','CalibrationID','DataAcquisitionID','DataProcessingID')] <-paste0('Sample',names(tmp_sample)[!names(tmp_sample)%in%c('CruiseLocalID','LogDistance','InstrumentID','CalibrationID','DataAcquisitionID','DataProcessingID')])
			
			
			tmp_data <- lapply(AcousticDataOne$Data, unlist)
			
			# tmp_data$SaCategory<-unlist(tmp_data$SaCategory)
			# tmp_data$Type<-unlist(tmp_data$Type)
			# tmp_data$Unit<-unlist(tmp_data$Unit)
			# tmp_data$Value<-formatC(tmp_data$Value,format='fg')
			
			names(tmp_data)[names(tmp_data)=='LocalID']<-'CruiseLocalID'
			names(tmp_data)[names(tmp_data)=='Distance']<-'LogDistance'
			names(tmp_data)[names(tmp_data)=='ChannelDepthUpper']<-'SampleChannelDepthUpper'
			names(tmp_data)[!names(tmp_data)%in%c('CruiseLocalID','LogDistance','SampleChannelDepthUpper')]<-paste0('Data',names(tmp_data)[!names(tmp_data)%in%c('CruiseLocalID','LogDistance','SampleChannelDepthUpper')])
			
			
			
			tmp_sub <- merge(as.data.table(tmp_log), as.data.table(tmp_sample), by=intersect(names(tmp_log), names(tmp_sample)))
			tmp_sub <- merge(as.data.table(tmp_data), as.data.table(tmp_sub), by=intersect(names(tmp_data), names(tmp_sub)))
			
			HDat <- cbind(as.data.table(hl), tmp_sub)
			
			tmp <- list(Instrument=HInst,
						Calibration=HCal,
						DataAcquisition=HDatA,
						DataProcessing=HDatP,
						Cruise=HCru,
						Data=HDat
			)
			
			if(Combine) {
				tmp <- data.table::rbindlist(tmp, fill = TRUE)
			}
			
		}
		else{
			stop('StoX: only ices acoustic format is allowed')
		}
	
	return(tmp)
}



#' Write ICES biotic CSV format file
#'
#' Given an \code{BioticData} object, this function will write an ICES biotic CSV file. Note that this function only supports
#' \code{BioticData} object that is created from reading an NMD biotic version 3 XML file.
#'
#' @param BioticData a \code{BioticData} object from an XML file with NMD biotic version 3 format.
#' @param CruiseSurvey A string naming the survey. Must be one of the names listed on \url{https://vocab.ices.dk/?ref=1453} or NONE (the default).
#' @param CruiseOrganisation An integer code representing the organization running the cruise. See \url{https://vocab.ices.dk/?ref=1398} for a list of possible codes (e.g., Institute of Marine Research: 612).
#' @param AllowRemoveSpecies ICES submission will not allow the resulting CSV file to be uploaded if the file contains species not listed in
#'        https://acoustic.ices.dk/Services/Schema/XML/SpecWoRMS.xml . Setting this parameter to TRUE will remove the unlisted species records.
#' @param Combine Logical: If TRUE stack the output tables for each acosutic file.
#'        
#' @return List of data.table objects in the ICES acoustic CSV format.
#'
#' @export
WriteICESBiotic <- function(
	BioticData, 
	CruiseSurvey = "NONE", 
	CruiseOrganisation = integer(), 
	AllowRemoveSpecies = TRUE, 
	Combine = FALSE
) {

	out <- lapply(
		BioticData, 
		prepareICESBiotic, 
		CruiseSurvey = CruiseSurvey, 
		CruiseOrganisation = CruiseOrganisation, 
		AllowRemoveSpecies = AllowRemoveSpecies, 
		Combine = Combine
	)
	
	return(out)
} 



prepareICESBiotic <- function(
	BioticDataOne, 
	CruiseSurvey = "NONE", 
	CruiseOrganisation = integer(), 
	AllowRemoveSpecies = TRUE, 
	Combine = FALSE
) {
	
	if(!(BioticDataOne$metadata$useXsd %in% c("nmdbioticv3", "nmdbioticv3.1"))) {
		warning("StoX: writeICESBiotic: Only NMD Biotic version 3 and 3.1 data accepted as input!")
		return(NA)
	}
	
	cruiseRaw <- BioticDataOne$mission
	
	Cruise <- cruiseRaw[, .(
		Cruise = "Cruise",
		Header = "Record",
		CruiseSurvey = CruiseSurvey,
		CruiseCountry = "NO",
		CruiseOrganisation = CruiseOrganisation,
		CruisePlatform = getICESShipCode(platformname),
		CruiseStartDate = gsub("Z", "", missionstartdate),
		CruiseEndDate = gsub("Z", "", missionstopdate),
		CruiseLocalID = cruise
	)]
	
	CruiseLocalID <- Cruise$CruiseLocalID
	
	haulRaw <- merge(cruiseRaw, BioticDataOne$fishstation)
	
	Haul <- haulRaw[, .(
		Haul = "Haul",
		Header = "Record",
		CruiseLocalID = cruise,
		HaulGear = gear,
		HaulNumber = serialnumber,
		HaulStationName = station,
		HaulStartTime = ifelse(is.na(stationstartdate) | is.na(stationstarttime), NA, gsub("Z", " ", paste0(stationstartdate, substr(stationstarttime, 1, 5)))),
		HaulDuration = getTimeDiff(stationstartdate, stationstarttime, stationstopdate, stationstoptime),
		HaulValidity = getHaulVal(gearcondition, samplequality),
		HaulStartLatitude = latitudestart,
		HaulStartLongitude = longitudestart,
		HaulStopLatitude = latitudeend,
		HaulStopLongitude = longitudeend,
		HaulStatisticalRectangle = getICESrect(latitudestart, longitudestart),
		HaulMinTrawlDepth = ifelse(is.na(fishingdepthmin), fishingdepthmax, fishingdepthmin),
		HaulMaxTrawlDepth = fishingdepthmax,
		HaulBottomDepth = ifelse(bottomdepthstop > fishingdepthmax, bottomdepthstop, NA),
		HaulDistance = round(getDistanceMeter(latitudestart, longitudestart, latitudeend, longitudeend)),
		HaulNetopening = verticaltrawlopening,
		HaulCodendMesh = NA,
		HaulSweepLength = getGOVSweepByEquipment(gear),
		HaulGearExceptions = NA,
		HaulDoorType = trawldoortype,
		HaulWarpLength = wirelength,
		HaulWarpDiameter = wirediameter,
		HaulWarpDensity = wiredensity,
		HaulDoorSurface = trawldoorarea,
		HaulDoorWeight = trawldoorweight,
		HaulDoorSpread = trawldoorspread,
		HaulWingSpread = wingspread,
		HaulBuoyancy = NA,
		HaulKiteArea = NA,
		HaulGroundRopeWeight = NA,
		HaulRigging = NA,
		HaulTickler = NA,
		HaulHydrographicStationID = NA,
		HaulTowDirection = direction,
		HaulSpeedGround = NA,
		HaulSpeedWater = gearflow,
		HaulWindDirection = winddirection,
		HaulWindSpeed = windspeed,
		HaulSwellDirection = NA,
		HaulSwellHeight = NA,
		HaulLogDistance = NA,
		HaulStratum = NA
	)]
	
	catchRaw <- merge(BioticDataOne$catchsample, haulRaw, by = intersect(names(BioticDataOne$catchsample), names(haulRaw)))
	
	# We must filter records with aphia == NA
	catchRaw <- catchRaw[!is.na(aphia)]
	
	Catch <- catchRaw[, .(
		Catch = "Catch",
		Header = "Record",
		CruiseLocalID = cruise,
		HaulGear = gear,
		HaulNumber = serialnumber,
		CatchDataType = "R",
		CatchSpeciesCode = aphia,
		CatchSpeciesValidity = ifelse(is.na(catchproducttype), 0, catchproducttype),
		CatchSpeciesCategory = catchpartnumber,
		CatchSpeciesCategoryNumber = catchcount,
		CatchWeightUnit = "kg",
		CatchSpeciesCategoryWeight = catchweight,
		CatchSpeciesSex = NA,
		CatchSubsampledNumber = lengthsamplecount,
		CatchSubsamplingFactor = catchcount / lengthsamplecount,
		CatchSubsampleWeight = lengthsampleweight,
		CatchLengthCode = NA,
		CatchLengthClass = NA,
		CatchLengthType = "1",
		CatchNumberAtLength = lengthsamplecount,
		CatchWeightAtLength = NA
	)]
	
	# Logic for missing important records
	Catch[is.na(CatchSpeciesCategoryNumber) & is.na(CatchSpeciesCategoryWeight) & !is.na(CatchSubsampledNumber), CatchSpeciesCategoryNumber := CatchSubsampledNumber]
	Catch[is.na(CatchSpeciesCategoryNumber) & is.na(CatchSpeciesCategoryWeight) & !is.na(CatchSubsampleWeight), CatchSpeciesCategoryWeight := CatchSubsampleWeight]
	
	# NA means that nothing is subsampled
	Catch[!is.na(CatchSpeciesCategoryWeight) & is.na(CatchSubsampleWeight), CatchSubsampleWeight := 0]
	
	# Set Haul without any catch as invalid hauls
	'%ni%' <- Negate('%in%')
	Haul[HaulNumber %ni% unique(Catch$HaulNumber), HaulValidity := "I"]
	
	# Combine required tables for the Biology level
	indRaw <- BioticDataOne$individual
	indRaw[is.na(preferredagereading), preferredagereading := 1]
	
	baseAge <- intersect(names(indRaw), names(BioticDataOne$agedetermination))
	indRaw <- merge(indRaw, BioticDataOne$agedetermination, by.x=c(baseAge, "preferredagereading"), by.y= c(baseAge, "agedeterminationid"), all.x = TRUE)
	indRaw <- merge(catchRaw, indRaw, by = intersect(names(catchRaw), names(indRaw)))
	
	Biology <- indRaw[, .(
		Biology = "Biology",
		Header = "Record",
		CruiseLocalID = cruise,
		HaulGear = gear,
		HaulNumber = serialnumber,
		CatchSpeciesCode = aphia,
		CatchSpeciesCategory = catchpartnumber,
		BiologyStockCode = NA,
		BiologyFishID = specimenid,
		BiologyLengthCode = "mm",
		BiologyLengthClass = length * 1000,
		BiologyWeightUnit = 'gr',
		BiologyIndividualWeight = individualweight * 1000,
		BiologyIndividualSex = ifelse(is.na(sex), NA, ifelse(sex == "1", "F", "M")),
		BiologyIndividualMaturity = getDATRASMaturity(getQuarter(stationstartdate), aphia, specialstage, maturationstage),
		BiologyMaturityScale = "M6",
		BiologyIndividualAge = age,
		BiologyAgePlusGroup = NA,
		BiologyAgeSource = "Otolith",
		BiologyGeneticSamplingFlag = NA,
		BiologyStomachSamplingFlag = NA,
		BiologyParasiteSamplingFlag = NA,
		BiologyIndividualVertebraeCount = NA
	)]
	
	if(AllowRemoveSpecies) {
		message("All species that is not listed in https://acoustic.ices.dk/Services/Schema/XML/SpecWoRMS.xml are automatically removed. Set AllowRemoveSpecies = FALSE to prevent this.")
		# Check for valid aphias, mark other as invalid
		xmlRaw <- read_xml("https://acoustic.ices.dk/Services/Schema/XML/SpecWoRMS.xml")
		validCodes <- xml_text(xml_find_all(xmlRaw, "//Code//Key"))
		Catch <- Catch[CatchSpeciesCode %in% validCodes, ]
		Biology <- Biology[CatchSpeciesCode %in% validCodes, ]
	} else {
		message("AllowRemoveSpecies is set to FALSE. Will only give warning for records with species that is not accepted by the ICES system.")
		compareICES("https://acoustic.ices.dk/Services/Schema/XML/SpecWoRMS.xml", unique(Catch$CatchSpeciesCode))
	}
	
	bioticOutput <- list(Cruise = Cruise, Haul = Haul, Catch = Catch, Biology = Biology)
	
	if(Combine) {
		bioticOutput <- data.table::rbindlist(bioticOutput, fill = TRUE)
	}
	
	return(bioticOutput)
}



#' Write ICES DATRAS (NS-IBTS) format file
#'
#' Given an \code{BioticData} object, this function will write an ICES DATRAS (NS-IBTS) file. Note that this function only supports
#' \code{BioticData} object that is created from reading an NMD biotic version 3 XML file.
#'
#' @param BioticData a \code{BioticData} object from an XML file with NMD biotic version 3 format.
#' @param AddStationType additional StationType to be included. By default only fish stations with StationType == NA are included.
#' @param Combine Logical: If TRUE stack the output tables for each acosutic file.
#'
#' @return List of data.table objects in the ICES DATRAS CSV format.
#'
#' @importFrom stats aggregate
#' @importFrom data.table copy
#' @export
WriteICESDatras <- function(
	BioticData, 
	AddStationType = NA_character_, 
	Combine = FALSE
) {

  out <- lapply(
  	BioticData, 
  	prepareICESDatras, 
  	AddStationType = AddStationType, 
  	Combine = Combine
  )

  return(out)
}


prepareICESDatras <- function(
	BioticDataOne, 
	AddStationType, 
	Combine = FALSE
) {
	# Check input is a NMD Biotic v3 data
	if(!(BioticDataOne$metadata$useXsd %in% c("nmdbioticv3", "nmdbioticv3.1"))) {
		warning("StoX: writeICESBiotic: Only NMD Biotic version 3 and 3.1 data accepted as input!")
		return(NA)
	}
	
	# Construct user-defined additional stations (default to NA stationtypes)
	if(length(AddStationType) > 0 && !all(is.na(AddStationType))) {
		targetStationType <- AddStationType
	} else {
		targetStationType <- c()
	}
	
	## 1. HH ##
	getTSCountryByIOC <- function(nation) {
		cnvTbl <- c("58" = "NO")
		
		x <- cnvTbl[as.character(nation)]
		x[is.null(x)] <- NA
		return(x)
	}
	
	'%ni%' <- Negate('%in%')
	
	getGearExp <- function(sweep, year, serialnumber, depth) {
		
		temp <-  as.data.table(cbind(sweep, year, serialnumber, depth))
		temp[, res:= "S"]
		
		temp[year == 2011 & serialnumber > 24362 & depth >= 70 | year == 2012
			 | year == 2011 & serialnumber >= 24135 & depth >= 70, res:="ST"]
		
		return (temp$res)
	}
	
	getYear <- function(stationstartdate) {
		format(as.Date(stationstartdate, format="%Y-%m-%dZ"), "%Y")
	}
	
	getMonth <- function(stationstartdate) {
		format(as.Date(stationstartdate, format="%Y-%m-%dZ"), "%m")
	}
	
	getDay <- function(stationstartdate) {
		format(as.Date(stationstartdate, format="%Y-%m-%dZ"), "%d")
	}
	
	getTimeShot <- function(stationstarttime) {
		
		timeshot <- function(y) {
			if(length(y) == 3) {
				return(paste0(y[1], y[2]))
			} else {
				return(NA)
			}
		}
		
		x <- strsplit(stationstarttime, ":")
		
		return(unlist(lapply(x, timeshot)))
	}
	
	getQuarter <- function(stationstartdate) {
		x <- format(as.Date(stationstartdate, format="%Y-%m-%dZ"), "%m")
		return(floor((as.numeric(x) - 1) / 3 + 1))
	}
	
	# Adopted from: https://www.mathworks.com/matlabcentral/fileexchange/62180-sunriseset-lat-lng-utcoff-date-plot
	getDayNight <- function(stationstartdate, stationstarttime, latitudestart, longitudestart, UTCoff = 0) {
		
		deg2rad <- function(val) {
			return(val * (pi / 180))
		}
		
		rad2deg <- function(val) {
			return(val * (180 / pi))
		}
		
		datetime0 <- as.POSIXct("1990-12-30", tz = "GMT")
		
		uniqueDates <- unique(stationstartdate)
		
		nDaysA = as.numeric(difftime(uniqueDates, datetime0, units = "days")) # Number of days since 01/01
		
		nTimes = 24*3600;                       # Number of seconds in the day
		tArray = seq(0, 1, length = nTimes);
		
		ssTab <- list()
		
		for(idx in seq_len(length(nDaysA))) {
			
			nDays <- nDaysA[idx]
			lat <- latitudestart[idx]
			lng <- longitudestart[idx]
			localdate <- as.POSIXct(uniqueDates[idx], tz = "GMT")
			
			# Compute
			# Letters correspond to colums in the NOAA Excel
			E = tArray;
			F = nDays+2415018.5+E-UTCoff/24;
			G = (F-2451545)/36525;
			I = (280.46646+G * (36000.76983+G*0.0003032)) %% 360;
			J = 357.52911+G * (35999.05029-0.0001537*G);
			K = 0.016708634-G * (0.000042037+0.0000001267*G);
			L = sin(deg2rad(J)) * (1.914602-G * (0.004817+0.000014*G))+sin(deg2rad(2*J)) * (0.019993-0.000101*G)+sin(deg2rad(3*J))*0.000289;
			M = I+L;
			P = M-0.00569-0.00478*sin(deg2rad(125.04-1934.136*G));
			Q = 23+(26+((21.448-G * (46.815+G * (0.00059-G*0.001813))))/60)/60;
			R = Q+0.00256*cos(deg2rad(125.04-1934.136*G));
			T = rad2deg(asin(sin(deg2rad(R)) * sin(deg2rad(P))));
			U = tan(deg2rad(R/2)) * tan(deg2rad(R/2));
			V = 4*rad2deg(U * sin(2*deg2rad(I))-2*K * sin(deg2rad(J))+4*K * U * sin(deg2rad(J)) *  cos(2*deg2rad(I))-0.5 * U * U * sin(4*deg2rad(I))-1.25 * K * K * sin(2 * deg2rad(J)));
			AB = (E*1440+V+4*lng-60*UTCoff) %% 1440
			
			AC = ifelse (AB/4 < 0, AB/4+180, AB/4-180)
			
			AD = rad2deg(acos(sin(deg2rad(lat))*sin(deg2rad(T))+cos(deg2rad(lat))*cos(deg2rad(T)) * cos(deg2rad(AC))));
			W = rad2deg(acos(cos(deg2rad(90.833)) / (cos(deg2rad(lat))*cos(deg2rad(T)))-tan(deg2rad(lat))*tan(deg2rad(T))));
			X = (720-4*lng-V+UTCoff*60)*60;
			
			sunrise = which.min(abs(X-round(W*4*60) - nTimes*tArray));
			sunrisetime = localdate + sunrise
			sunset = which.min(abs(X+round(W*4*60) - nTimes*tArray));
			sunsettime = localdate + sunset
			
			ssTab[[uniqueDates[idx]]] <- list(sunrise = sunrisetime, sunset = sunsettime)
		}
		
		getDN <- function(x, ssTab) {
			
			y <- ssTab[[format(x, "%Y-%m-%dZ")]]
			
			if(x < y$sunrise || x >= y$sunset) {
				return("N")
			} else {
				return("D")
			}
		}
		
		datetime <- as.POSIXct(gsub("Z", " ", paste0(stationstartdate, stationstarttime)), tz = "GMT")
		
		return(unlist(lapply(datetime, getDN, ssTab)))
	}
	
	finalHH <- merge(BioticDataOne$mission, BioticDataOne$fishstation)
	
	# Make HH records and filter only valid stationtype
	finalHH[, `:=`(
		"RecordType" = "HH",
		"Quarter" = getQuarter(stationstartdate),
		"Country" = getTSCountryByIOC(nation),
		"Ship" = getICESShipCode(platformname),
		"Gear" = "GOV",
		"SweepLngt" = getGOVSweepByEquipment(gear),
		"GearExp" = getGearExp(getGOVSweepByEquipment(gear), startyear, serialnumber, bottomdepthstart),
		"DoorType" = "P",
		"StNo" = serialnumber,
		"HaulNo" = station,
		"Year" = getYear(stationstartdate),
		"Month" = getMonth(stationstartdate),
		"Day" = getDay(stationstartdate),
		"TimeShot" = getTimeShot(stationstarttime),
		"Stratum" = NA,
		"HaulDur" = as.numeric(getTimeDiff(stationstartdate, stationstarttime, stationstopdate, stationstoptime)),
		"DayNight" = getDayNight(stationstartdate, stationstarttime, latitudestart, longitudestart),
		"ShootLat" = round(latitudestart, 4),
		"ShootLong" = round(longitudestart, 4),
		"HaulLat" = round(latitudeend, 4),
		"HaulLong" = round(longitudeend, 4),
		"StatRec" = getICESrect(latitudestart, longitudestart),
		"Depth" = round(bottomdepthstart),
		"HaulVal" = getHaulVal(gearcondition, samplequality),
		"HydroStNo" = NA,
		"StdSpecRecCode" = 1,
		"BycSpecRecCode" = 1,
		"DataType" = "R",
		"Netopening"= round(verticaltrawlopening, 1),
		"Rigging" = NA,
		"Tickler" = NA,
		"Distance" = round(getDistanceMeter(latitudestart, longitudestart, latitudeend, longitudeend)),
		"Warpingt" = round(wirelength),
		"Warpdia" = NA,
		"WarpDen" = NA,
		"DoorSurface" = 4.5,
		"DoorWgt" = 1075,
		"DoorSpread" = ifelse(!is.na(trawldoorspread), round(trawldoorspread, 1), NA),
		"WingSpread" = NA,
		"Buoyancy" = NA,
		"KiteDim" = 0.8,
		"WgtGroundRope" = NA,
		"TowDir" = ifelse(!is.na(direction), round(direction), NA),
		"GroundSpeed" = round(gearflow, 1),
		"SpeedWater" = NA,
		"SurCurDir" = NA,
		"SurCurSpeed" = NA,
		"BotCurDir" = NA,
		"BotCurSpeed" = NA,
		"WindDir" = NA,
		"WindSpeed" = NA,
		"SwellDir" = NA,
		"SwellHeight" = NA,
		"SurTemp" = NA,
		"BotTemp" = NA,
		"SurSal" = NA,
		"BotSal" = NA,
		"ThermoCline" = NA,
		"ThClineDepth" = NA,
		"CodendMesh" = NA ,
		"SecchiDepth" = NA,
		"Turbidity" = NA,
		"TidePhase" = NA,
		"TideSpeed" = NA,
		"PelSampType" = NA,
		"MinTrawlDepth" = NA,
		"MaxTrawlDepth" = NA
	)]
	
	HHraw <- data.table::copy(finalHH[is.na(stationtype) | stationtype %in% targetStationType, c("RecordType", "Quarter", "Country", "Ship", "Gear",
																					 "SweepLngt", "GearExp", "DoorType", "StNo", "HaulNo", "Year", "Month", "Day",
																					 "TimeShot", "Stratum", "HaulDur", "DayNight", "ShootLat", "ShootLong", "HaulLat", "HaulLong",
																					 "StatRec", "Depth", "HaulVal", "HydroStNo", "StdSpecRecCode", "BycSpecRecCode", "DataType", "Netopening",
																					 "Rigging", "Tickler", "Distance", "Warpingt", "Warpdia", "WarpDen", "DoorSurface", "DoorWgt",
																					 "DoorSpread", "WingSpread", "Buoyancy", "KiteDim", "WgtGroundRope", "TowDir", "GroundSpeed",
																					 "SpeedWater", "SurCurDir", "SurCurSpeed", "BotCurDir", "BotCurSpeed", "WindDir", "WindSpeed",
																					 "SwellDir", "SwellHeight", "SurTemp", "BotTemp", "SurSal", "BotSal", "ThermoCline", "ThClineDepth",
																					 "CodendMesh", "SecchiDepth", "Turbidity", "TidePhase", "TideSpeed", "PelSampType", "MinTrawlDepth", "MaxTrawlDepth")]
	)
	
	## 2. HL ##
	
	mergedHL <- merge(BioticDataOne$catchsample, finalHH, by=intersect(names(BioticDataOne$catchsample), names(finalHH)))
	
	groupCA <- c("missiontype", "startyear", "platform", "missionnumber", "serialnumber", "aphia", "sex")
	groupHL <- c(groupCA, "catchpartnumber")
	
	# Remove rows with empty aphia
	mergedHL <- mergedHL[!is.na(aphia)]
	
	getSpecVal <- function(HaulVal, catchcount, lengthsamplecount, catchweight){
		temp <-  as.data.table(cbind(hv=HaulVal, cc=catchcount, lsc=lengthsamplecount, cw=catchweight))
		
		# Default is invalid
		temp[, res := "0"]
		
		temp[!is.na(cc) & !is.na(lsc) & !is.na(cw), res:="1"]
		temp[!is.na(cc) &  is.na(lsc) &  is.na(cw), res:="4"]
		temp[ is.na(cc) &  is.na(lsc) & !is.na(cw), res:="6"]
		temp[!is.na(cc) &  is.na(lsc) & !is.na(cw), res:="7"]
		temp[ is.na(cc) &  is.na(lsc) &  is.na(cw), res:="5"]
		temp[!is.na(cc) & !is.na(lsc) &  is.na(cw), res:="0"]
		
		temp[hv == "I", res:="0"]
		
		return(temp$res)
	}
	
	# Convert Length Measurement Type
	# http://tomcat7.imr.no:8080/apis/nmdapi/reference/v2/dataset/lengthmeasurement?version=2.0
	# http://vocab.ices.dk/?ref=1392
	convLenMeasType <- function(LenMeasType) {
		# Convert table
		ct <- c("B" = 5,
				"C" = 6,
				"E" = 1,
				"F" = 8,
				"G" = 4,
				"H" = 3,
				"J" = 2,
				"L" = 7,
				"S" = 9)
		return(ct[LenMeasType])
	}
	
	mergedHL[, SpecVal := getSpecVal(HaulVal, catchcount, lengthsamplecount, catchweight)]
	
	# Get herring or sprat
	mergedHL[,`:=`( isHerringOrSprat = ifelse(aphia %in% c("126417", "126425"), TRUE, FALSE),
					isCrustacean = ifelse(aphia %in% c("107275", "107276", "107369", "107253", "107703", "107704", "107350", "107254", "107205", "140712", "140687", "140658"), TRUE, FALSE))]
	
	# Calculate lngtCode
	mergedHL[,lngtCode := "1"]
	mergedHL[is.na(sampletype), lngtCode := NA]
	mergedHL[isCrustacean == TRUE, lngtCode := "."]
	mergedHL[isHerringOrSprat == TRUE, lngtCode := "0"]
	
	# lenInterval, and reportInMM
	mergedHL[,`:=`(lenInterval = ifelse(lngtCode=="0", 5, 1), reportInMM = ifelse(lngtCode %ni% c("1", NA), TRUE, FALSE))]
	mergedHL[is.na(lenInterval), lenInterval := 1]
	
	# catCatchWgt & subWeight
	mergedHL[!is.na(catchweight), catCatchWgt := ceiling(catchweight * 1000)]
	mergedHL[!is.na(lengthsampleweight), subWeight := ceiling(lengthsampleweight * 1000)]
	
	# get sampleFac
	mergedHL[, sampleFac := catchweight / lengthsampleweight]
	
	# Merge with individual
	mergedHL <- merge(mergedHL, BioticDataOne$individual, by = intersect(names(mergedHL), names(BioticDataOne$individual)), all.x = TRUE)
	
	# Get count
	mergedHL[, N := sum(!is.na(specimenid)), by = groupHL]
	
	# For the record with empty individual data
	mergedHL[N == 0, `:=`(lngtClass = as.integer(NA), sex = as.character(NA))]
	
	# Get Individual length
	mergedHL[, length := length * 100]
	
	# Some species have very small length in cm, use mm instead
	mergedHL[length < 1, `:=`(lngtCode = ".", lenInterval = 1, reportInMM = TRUE)]
	
	# Process MM length
	mergedHL[reportInMM == TRUE, length := length * 10]
	
	# Get sex
	mergedHL[, sex := ifelse(is.na(sex), as.character(NA), ifelse(sex == "1", "F", "M"))]
	
	# Get lngtClass
	for(interval in unique(mergedHL$lenInterval)) {
		intVec <- seq(0, max(mergedHL$length, na.rm = T), by = interval)
		mergedHL[lenInterval == interval, lngtClass := intVec[findInterval(length, intVec)]]
	}
	
	# Count measured individual
	mergedHL[!is.na(length), lsCountTot := 1]
	
	# Aggregate hlNoAtLngth and lsCountTot
	finalHL <- mergedHL[, .(N, lsCountTot = sum(lsCountTot)), by = c(groupHL,  "lngtClass",
																	 "Quarter",
																	 "Country",
																	 "Ship",
																	 "Gear",
																	 "SweepLngt", "GearExp", "DoorType", "HaulNo", "SpecVal", "catCatchWgt", "sampleFac", "subWeight", "lngtCode", "stationtype", "lengthmeasurement")]
	
	finalHL <- finalHL[!duplicated(finalHL)]
	finalHL[,`:=`(noMeas = sum(lsCountTot)), by = groupHL]
	finalHL[,`:=`(totalNo = noMeas * sampleFac, subFactor = sampleFac)]
	
	HLraw <- data.table::copy(finalHL[is.na(stationtype) | stationtype %in% targetStationType, .("RecordType" = "HL",
																					 "Quarter" = Quarter,
																					 "Country" = Country,
																					 "Ship" = Ship,
																					 "Gear" = Gear,
																					 "SweepLngt" = SweepLngt,
																					 "GearExp" = GearExp,
																					 "DoorType" = DoorType,
																					 "StNo" = serialnumber,
																					 "HaulNo" = HaulNo,
																					 "Year" = startyear,
																					 "SpecCodeType" = "W",
																					 "SpecCode" = aphia,
																					 "SpecVal" = SpecVal,
																					 "Sex" = sex,
																					 "TotalNo" = round(totalNo, 2),
																					 "CatIdentifier" = catchpartnumber,
																					 "NoMeas" = noMeas,
																					 "SubFactor" = round(subFactor, 4),
																					 "SubWgt" = round(subWeight),
																					 "CatCatchWgt" = round(catCatchWgt),
																					 "LngtCode" = lngtCode,
																					 "LngtClass" = lngtClass,
																					 "HLNoAtLngt" = round(lsCountTot, 2),
																					 "DevStage" = NA,
																					 "LenMeasType" = convLenMeasType(lengthmeasurement)
	)]
	)
	
	
	## 3. CA ##
	
	mergedHL[is.na(preferredagereading), preferredagereading := 1]
	baseAge <- intersect(names(mergedHL), names(BioticDataOne$agedetermination))
	mergedCA <- merge(mergedHL, BioticDataOne$agedetermination, by.x=c(baseAge, "preferredagereading"), by.y= c(baseAge, "agedeterminationid"), all.x = TRUE)
	
	# Remove empty individual
	mergedCA <- mergedCA[!is.na(specimenid)]
	
	# Get maturity
	mergedCA[, maturity:=getDATRASMaturity(Quarter, aphia, specialstage, maturationstage)]
	
	# Aggregate count
	mergedCA[!is.na(individualweight), `:=`(nWithWeight =.N, totWeight = sum(individualweight)), by = c(groupCA,  "lngtClass", "maturity", "age")]
	
	finalCA <- mergedCA[, .(nInd =.N), by = c(groupCA,  "lngtClass", "maturity", "age",
											  "Quarter",
											  "Country",
											  "Ship",
											  "Gear",
											  "SweepLngt", "GearExp", "DoorType", "HaulNo", "SpecVal", "StatRec", "lngtCode", "stationtype", "nWithWeight", "totWeight",
											  "specimenid", "tissuesample", "stomach", "agingstructure", "readability", "parasite")]
	finalCA[!is.na(nWithWeight),  meanW := totWeight / nWithWeight]
	
	# Convert aging structure source
	# http://tomcat7.imr.no:8080/apis/nmdapi/reference/v2/dataset/agingstructure?version=2.0
	# http://vocab.ices.dk/?ref=1507
	convAgeSource <- function(AgeSource) {
		# Convert table
		ct <- c("1" = "scale",
				"2" = "otolith",
				"4" = "df-spine",
				"6" = "spine",
				"7" = "vertebra",
				"8" = "caudal-thorn")
		return(ct[AgeSource])
	}
	
	CAraw <- data.table::copy(finalCA[is.na(stationtype) | stationtype %in% targetStationType, .("RecordType" = "CA",
																					 "Quarter" = Quarter,
																					 "Country" = Country,
																					 "Ship" = Ship,
																					 "Gear" = Gear,
																					 "SweepLngt" = SweepLngt,
																					 "GearExp" = GearExp,
																					 "DoorType" = DoorType,
																					 "StNo" = serialnumber,
																					 "HaulNo" = HaulNo,
																					 "Year" = startyear,
																					 "SpecCodeType" = "W",
																					 "SpecCode" = aphia,
																					 "AreaType" = "0",
																					 "AreaCode" = StatRec,
																					 "LngtCode" = lngtCode,
																					 "LngtClass" = lngtClass,
																					 "Sex" = sex,
																					 "Maturity" = maturity,
																					 "PlusGr" = as.character(NA),
																					 "AgeRings" = ifelse(!is.na(age), age, NA),
																					 "CANoAtLngt" = nInd,
																					 "IndWgt" = ifelse(!is.na(meanW), round(meanW * 1000, 1), NA),
																					 "MaturityScale" = "M6",
																					 "FishID" = specimenid,
																					 "GenSamp" = ifelse(!is.na(tissuesample), "Y", "N"),
																					 "StomSamp" = ifelse(!is.na(stomach), "Y", "N"),
																					 "AgeSource" = convAgeSource(agingstructure),
																					 "AgePrepMet" = NA,
																					 "OtGrading" = ifelse(readability %in% as.character(c(1:4)), readability, NA),  # From http://tomcat7.imr.no:8080/apis/nmdapi/reference/v2/dataset/otolithreadability?version=2.0 and http://vocab.ices.dk/?ref=1395
																					 "ParSamp" = ifelse(!is.na(parasite), "Y", "N")
	)]
	)
	
	
	## Prepare for cleaning (mostly from the old Rstox. Will change as we improve) ##
	hh <- HHraw
	hl <- HLraw
	ca <- CAraw
	
	## WARN #0:
	# It's possible to have two same aphia (but different species, e.g. SILD05) catch sampes in a haul.
	# We need to combine them if we have two different TotalNo and catcatchwgt.
	
	# Find duplicate species in a haul
	dupl <- aggregate(catchcategory ~ aphia + serialnumber, BioticDataOne$catchsample, FUN = function(x) length(unique(x)))
	dupl <- dupl[dupl$catchcategory > 1, ]
	
	# Find the above in our DATRAS HL
	if(nrow(dupl)) {
		found <- aggregate(CatCatchWgt ~ StNo + SpecCode + Sex + CatIdentifier, hl[(hl$SpecCode %in% dupl$aphia & hl$StNo %in% dupl$serialnumber),], FUN = function(x) length(unique(x)))
		found <- found[found$CatCatchWgt > 1, ]
		for(iz in seq_len(nrow(found))) {
			tmpHL <- hl[hl$StNo==found[iz, "StNo"] & hl$SpecCode==found[iz, "SpecCode"] & hl$Sex==found[iz, "Sex"] & hl$CatIdentifier==found[iz, "CatIdentifier"], ]
			combinedCatCatchWgt <- tmpHL
			# Fix CatCatchWgt
			hl[hl$StNo==found[iz, "StNo"] & hl$SpecCode==found[iz, "SpecCode"] & hl$Sex==found[iz, "Sex"] & hl$CatIdentifier==found[iz, "CatIdentifier"], "CatCatchWgt"] <- round(mean(tmpHL$CatCatchWgt))
			# Fix CatCatchWgt
			hl[hl$StNo==found[iz, "StNo"] & hl$SpecCode==found[iz, "SpecCode"] & hl$Sex==found[iz, "Sex"] & hl$CatIdentifier==found[iz, "CatIdentifier"], "SubWgt"] <- round(mean(tmpHL$SubWgt))
			# Fix totalNo
			hl[hl$StNo==found[iz, "StNo"] & hl$SpecCode==found[iz, "SpecCode"] & hl$Sex==found[iz, "Sex"] & hl$CatIdentifier==found[iz, "CatIdentifier"], "TotalNo"] <- sum(unique(tmpHL$TotalNo))
			# Fix noMeas
			hl[hl$StNo==found[iz, "StNo"] & hl$SpecCode==found[iz, "SpecCode"] & hl$Sex==found[iz, "Sex"] & hl$CatIdentifier==found[iz, "CatIdentifier"], "NoMeas"] <- sum(tmpHL$HLNoAtLngt)
			# Finally, fix SubFactor
			hl[hl$StNo==found[iz, "StNo"] & hl$SpecCode==found[iz, "SpecCode"] & hl$Sex==found[iz, "Sex"] & hl$CatIdentifier==found[iz, "CatIdentifier"], "SubFactor"] <- sum(unique(tmpHL$TotalNo))/sum(tmpHL$HLNoAtLngt)
		}
	}
	
	## WARN #1:
	# Find species with different SpecVal, if any of them have SpecVal == 1, delete any other records with different SpecVal
	# otherwise, use the lowest SpecVal value for all
	
	tmp <- aggregate(SpecVal ~ SpecCode + StNo, hl, FUN = function(x) length(unique(x)))
	tmp <- tmp[tmp$SpecVal>1, ]
	
	for( rownum in seq_len(nrow(tmp)) ) {
		tmpSpecs <- hl[(hl$StNo==tmp$StNo[rownum] & hl$SpecCode==tmp$SpecCode[rownum]),]$SpecVal
		if(any(tmpSpecs == 1))
			hl <- hl[!(hl$StNo==tmp$StNo[rownum] & hl$SpecCode==tmp$SpecCode[rownum] & hl$SpecVal!=1),]
		else
			hl[(hl$StNo==tmp$StNo[rownum] & hl$SpecCode==tmp$SpecCode[rownum]), c("SpecVal")] <- min(tmpSpecs)
	}
	
	## SpecVal Conditionals
	hl[hl$SpecVal==0, c("Sex", "TotalNo", "CatIdentifier", "NoMeas", "SubFactor", "SubWgt", "CatCatchWgt", "LngtCode", "LngtClass", "HLNoAtLngt")] <- NA
	
	hl[hl$SpecVal==4, c("NoMeas", "SubWgt", "CatCatchWgt", "LngtCode", "LngtClass", "HLNoAtLngt")] <- NA
	hl[hl$SpecVal==4, c("SubFactor")] <- 1
	
	hl[hl$SpecVal==5, c("TotalNo", "NoMeas", "SubWgt", "CatCatchWgt", "LngtCode", "LngtClass", "HLNoAtLngt")] <- NA
	hl[hl$SpecVal==5, c("SubFactor")] <- 1
	
	hl[hl$SpecVal==6, c("TotalNo", "NoMeas", "LngtCode", "LngtClass", "HLNoAtLngt")] <- NA
	
	hl[hl$SpecVal==7, c("NoMeas", "LngtCode", "LngtClass", "HLNoAtLngt")] <- NA
	
	hl[hl$SpecVal==10, c("CatCatchWgt")] <- NA
	
	## WARN #2:
	## will now get errors in DATRAS upload for duplicate records
	hl <- hl[!duplicated(hl),]
	
	## hl and ca contain 0-tow info - must throw these out
	hl <- hl[hl$StNo %in% hh$StNo,]
	ca <- ca[ca$StNo %in% hh$StNo,]
	# throw out ca records for Invalid hauls
	ca <- ca[!ca$StNo %in% hh$StNo[hh$HaulVal=='I'],]
	
	##########################################
	## Removing some benthos - this won't be needed in the future
	## keep 11725 138139 138482 138483 140600 140621 140624 140625 141443 141444 141449 153083 153131-- these are cephaolopods
	## required benthos: 107205
	hl <- hl[!hl$SpecCode %in% c(230,558,830,883,1302,1839,100635,100706,100930,103929,106048,106087,106204,106733,106791,
								 106854,106928,107044,107218,107230,107240,107273,107292,107318,107330,107346,107397,107398,107551,
								 107616,107643,111374,111597,111604,116986,117302,117809,117815,117890,123117,123867,123920,123970,
								 123987,124319,124418,124913,124929,124934,125128,125131,125134,129196,129229,130464,130867,132072,
								 132480,135144,135302,137704,137732,138223,138239,138760,138899,139004,139488,140299,140627,141753,
								 144129,150642,178639,181228,23986719494,21263,100817,100982,106738,107160,107232,107277,107322,
								 107323,107327,107387,107531,107552,107564,107649,107651,111367,123080,123083,123084,123776,123813,
								 124043,124154,124160,124287,124535,125166,125333,128517,129840,138802,138878,138920,140467,140717,
								 143755,145541,145546,145548,532031,589677,1762,123082,149),]
	
	ca <- ca[!ca$SpecCode %in% c(230,558,830,883,1302,1839,100635,100706,100930,103929,106048,106087,106204,106733,106791,
								 106854,106928,107044,107218,107230,107240,107273,107292,107318,107330,107346,107397,107398,107551,
								 107616,107643,111374,111597,111604,116986,117302,117809,117815,117890,123117,123867,123920,123970,
								 123987,124319,124418,124913,124929,124934,125128,125131,125134,129196,129229,130464,130867,132072,
								 132480,135144,135302,137704,137732,138223,138239,138760,138899,139004,139488,140299,140627,141753,
								 144129,150642,178639,181228,23986719494,21263,100817,100982,106738,107160,107232,107277,107322,
								 107323,107327,107387,107531,107552,107564,107649,107651,111367,123080,123083,123084,123776,123813,
								 124043,124154,124160,124287,124535,125166,125333,128517,129840,138802,138878,138920,140467,140717,
								 143755,145541,145546,145548,532031,589677,1762,123082,149),]
	
	#more benthods 10216 = skate egg case
	hl <- hl[!hl$SpecCode %in% c(443,938,1131,1292,1337,1360,19494,22988,100751,100757,100790,101054,103484,104062,
								 106122,106669,107011,107052,107148,107239,107388,107563,110690,110911,110956,111411,117136,
								 117258,123260,123276,123321,123335,123574,123593 ,123851,123922,123985,124085,125158,125269,
								 128506,130467,130987,131779,134591,137683,141872,146142 ,149864,445590,510534,105,175,927,1107,
								 1135,1267,100793),]
	hl <- hl[!hl$SpecCode %in% c(105,175,927,1107,1135,1267,100793,103443,103692,106057,106835,106903,107558,110908,111361,
								 117940,122348,123160,123426,124257,125027,125284,131495,135294,135301,135306,138992,140528,140687,
								 167882,178527,239867,291396,106763,137656,117225,100653,125125,100698,131774,134366,123386,117228,
								 117994,138923,123127,137701,123320,131629 ,152391,1363,214,103543,106994,103450,129400,140143,
								 146420,141905,22496,988,103717,107163,982,985,123622,102145,1082,10216,103483),]
	
	ca <- ca[!ca$SpecCode %in% c(443,938,1131,1292,1337,1360,19494,22988,100751,100757,100790,101054,103484,104062,
								 106122,106669,107011,107052,107148,107239,107388,107563,110690,110911,110956,111411,117136,
								 117258,123260,123276,123321,123335,123574,123593 ,123851,123922,123985,124085,125158,125269,
								 128506,130467,130987,131779,134591,137683,141872,146142 ,149864,445590,510534,105,175,927,1107,
								 1135,1267,100793),]
	ca <- ca[!ca$SpecCode %in% c(105,175,927,1107,1135,1267,100793,103443,103692,106057,106835,106903,107558,110908,111361,
								 117940,122348,123160,123426,124257,125027,125284,131495,135294,135301,135306,138992,140528,140687,
								 167882,178527,239867,291396,106763,137656,117225,100653,125125,100698,131774,134366,123386,117228,
								 117994,138923,123127,137701,123320,131629 ,152391,1363,214,103543,106994,103450,129400,140143,
								 146420,141905,22496,988,103717,107163,982,985,123622,102145,1082,10216,103483),]
	
	hl <- hl[!hl$SpecCode %in% c(NA, 101,106769,106782,107010,107726,122478,123506,12437,124951,128539,129402,196221,205077,124373, 123187, 124710),]
	ca <- ca[!ca$SpecCode %in% c(NA, 101,106769,106782,107010,107726,122478,123506,12437,124951,128539,129402,196221,205077,124373, 123187, 124710),]
	
	## IU: Filter out additional benthos:
	benthosSpecCodes <- c(104,956,966,1128,1296,1367,1608,11707,100782,100839,100854,103439,103732,104040,105865,106041,106673,106702,106789,106834,107152,
						  107205,107264,110749,110916,110993,111152,111355,111365,117093,117195,118445,122626,123204,123255,123613,124147,124151,124324,124670,
						  128490,128503,129563,130057,134691,136025,137710,138018,138068,138477,138631,138749,138938,140166,140173,140480,140625,141904,141929,
						  149854,152997,532035,816800)
	
	hl <- hl[!hl$SpecCode %in% benthosSpecCodes,]
	ca <- ca[!ca$SpecCode %in% benthosSpecCodes,]
	
	
	## WARN #3:
	## ca records with no HL records
	## these records are because there is no catch weight
	## DATRAS does not accept length info without catch weight
	## so create a line in the HL for each, but give SpecValue=4 and delete ca record
	
	#IU: Improved cleaning#
	# Use join to find missing value in HL
	if (nrow(ca) > 0) {
		#testca <- unique(data.frame(StNo=ca$StNo, SpecCode=ca$SpecCode, ca=TRUE))
		#testhl <- unique(data.frame(StNo=hl$StNo, SpecCode=hl$SpecCode, hl=TRUE))
		testca <- unique(data.table::data.table(StNo=ca$StNo, SpecCode=ca$SpecCode, ca=TRUE))
		testhl <- unique(data.table::data.table(StNo=hl$StNo, SpecCode=hl$SpecCode, hl=TRUE))
		tt <- merge(testca, testhl, by = c("StNo","SpecCode"), all=TRUE)
		missingHL <- tt[is.na(tt$hl),]
		
		# Populate missing value in HL
		for(idxHL in seq_len(nrow(missingHL))) {
			r <- missingHL[idxHL,]
			tmp <- hl[hl$StNo==r$StNo,][1,]
			tmp$SpecCode <- r$SpecCode
			tmp$SpecVal <- 4
			tmp$TotalNo <- c(hh$HaulDur[hh$StNo==r$StNo])
			tmp$CatCatchWgt <- NA
			hl <- rbind(hl,tmp)
		}
	}
	## WARN #4:
	# Use plus group for herring and mackerel individuals with age ring above 15
	ca[ which((ca$SpecCode==127023 | ca$SpecCode==126417) & ca$AgeRings >= 15), c("PlusGr", "AgeRings")] <- list("+", 15)
	
	# Order HL
	hl <- hl[order(hl$StNo),]
	
	#
	datrasOutput <- list(HH=hh, HL=hl, CA=ca)
	
	if(Combine) {
		datrasOutput <- data.table::rbindlist(datrasOutput, fill = TRUE)
	}
	
	return(datrasOutput)
}

