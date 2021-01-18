processXSD <- function(doc, path = NULL) {

	getNameType <- function(x) {
		y <- xml_attrs(x)
		return(c(y[["name"]], y[["type"]]))
	}

	getNameTypeExt <- function(x, base) {
		y <- xml_attrs(x)
		return(c(base, y[["base"]]))
	}

	getRecNameType <- function(x, rootName = "", recEnv) {
		# Extract name
		sName <- x[2]

		# Get rootname
		if(rootName == "")
			rootName <- x[1]

		# Clean sName
		sName <- gsub(paste0(rootName, ":"), "", sName)

		# Extract elements
		y <- xml_find_all(doc, paste0("//", defNS, "complexType[@name=\"", sName, "\"]//", defNS, "element"))

		# This is needed for echosounder v1 SA records
		extension <- xml_find_all(doc, paste0("//", defNS, "complexType[@name=\"", sName, "\"]//", defNS, "extension"))

		# If no children and "KeyType" (This might be specific to NMDBioticv3) or "*IDREFType*" (specific to ICES XSDs)
		if(length(y) > 0 && (!grepl("KeyType", sName) || !grepl("IDREFType", sName))) {

			if(grepl("IDREFType", sName)) message("IDREFType\n")

			z <- lapply(lapply(y, getNameType), getRecNameType, rootName, recEnv)

			# ICES data have extension and children
			if(length(extension) > 0) {
				ext <- lapply(lapply(extension, getNameTypeExt, x[1]), getRecNameType, rootName, recEnv)
				z <- c(z, ext[[1]]$members)
			}

			# Prepare flat
			recEnv$flat[[x[1]]] <- z
			recEnv$flatAttr[[x[1]]] <- sapply(xml_find_all(doc, paste0("//", defNS, "complexType[@name=\"", sName, "\"]//", defNS, "attribute/@name")), function(xx) xml_text(xx))

			# Remove nested elements
			recEnv$flat[[x[1]]] <- lapply(recEnv$flat[[x[1]]], function(xx){ if(is.list(xx)) return(xx[[1]][1]) else return(xx) })
			return(list(x, members=z))
		# Below is specific for Echosounder v1's SA records (with XSD extension) and NMDBioticv1.x ("StringDescriptionType")
		} else if (length(extension) > 0 && !grepl("StringDescriptionType", sName))  {
			z <- lapply(extension, getNameTypeExt, x[1])

			# Prepare flat
			recEnv$flat[[x[1]]] <- z
			recEnv$flatAttr[[x[1]]] <- sapply(xml_find_all(doc, paste0("//", defNS, "complexType[@name=\"", sName, "\"]//", defNS, "attribute/@name")), function(xx) xml_text(xx))

			# Remove nested elements
			recEnv$flat[[x[1]]] <- lapply(recEnv$flat[[x[1]]], function(xx){ if(is.list(xx)) return(xx[[1]][1]) else return(xx) })
			return(list(x, members=z))
		} else {
			return(x)
		}
	}

	# Get the default namespace
	defNS <- names(xml_ns(doc))[[grep("XMLSchema", as.list(xml_ns(doc)))]]

	if(length(defNS) > 0)
		defNS <- paste0(defNS[[1]], ":")
	else
		defNS <- ""

	# See if we need to include more file(s) (include schemaLocation)
	extraXSD <- xml_find_all(doc, paste0("//", defNS, "include"))
	if(length(extraXSD) > 0) {
		message("We have extra XSDs to be included!\n")
		exFiles <- xml_attr(extraXSD, "schemaLocation")
		exObj <- lapply(paste0(path, "/include/", exFiles), read_xml)
		lapply(exObj, function(x) lapply(xml_children(x), function(y) xml_add_child(doc, y)))
	}

	rootInfo <- getNameType(xml_find_all(doc, paste0("/", defNS, "schema/", defNS, "element"))[[1]])

	r_e <- new.env()
	r_e$flat <- list()
	r_e$flatAttr <- list()

	# start the recursive search
	getRecNameType(rootInfo, recEnv = r_e)

	return(list(flat = r_e$flat, flatAttr = r_e$flatAttr, rootInfo = rootInfo, doc = doc))

}

processMetadata <- function(flat, flatAttr, rootInfo, xsdFile, xsdDoc) {

	# Recursively try to find the column types
	traceType <- function(xx) {
		if(grepl("xs[:]|xsd[:]", xx))
			return(xx)
		else {
			# Try to search the root "type"
			findstring <- paste0('//*[@name="', xx, '"]')
			el <- xml_find_all(xsdDoc, findstring)
			elu <- unique(xml_attr(el, "type"))
			res <- unlist(lapply(elu, traceType))

			# Now try to search the root "base"
			if(is.null(res)) {
				findstring <- paste0('//*[@name="', xx, '"]//*[@base]')
				el <- xml_find_all(xsdDoc, findstring)
				elu <- unique(xml_attr(el, "base"))
				res <- unlist(lapply(elu, traceType))
			}
			return(res)
		}
	}

	# Process headers and dimensions
	getAttrTable <- function(root, headAttr = c(), xpath="", recEnv) {

		rootStart <- root[1]

		state <- flat[[rootStart]]
		attrs <- flatAttr[[rootStart]]

		# Combine with attributs
		tableElements <- c(headAttr, attrs)
		tableTypes <- rep("attr", length(tableElements))

		# List prefix for next children	
		prefix <- c(headAttr, attrs)

		# Get length of header
		headLen <- length(headAttr)

		# Get number of elements
		tempXpath <- paste(xpath, rootStart, sep='/')
		elemNumXpath <- paste0("count(", tempXpath, ")")
		xpath <- tempXpath

		for(s in 1:length(state)) {
			if(length(state[[s]]) > 1) {
				tableElements <- c(tableElements, state[[s]][1])
				tableTypes <- c(tableTypes, state[[s]][2])
			} else {
				getAttrTable(state[[s]], prefix, xpath, recEnv = recEnv)
			}
		}

		recEnv$tableHeaders[[rootStart]] <- unlist(tableElements)
		recEnv$tablePrefix[[rootStart]] <- unlist(prefix)
		recEnv$tableTypes[[rootStart]] <- unlist(tableTypes)
		recEnv$levelDims[[rootStart]] <- elemNumXpath
	}

	# Meta data before go to C++

	r_e <- new.env()
	r_e$tableHeaders <- list()
	r_e$tablePrefix <- list()
	r_e$levelDims <- list()

	# For element types
	r_e$tableTypes <- list()

	# Defining root
	root <- rootInfo[1]

	# Start recursive
	getAttrTable(rootInfo, recEnv = r_e)

	# Fill in missing values
	missingSets <- setdiff(names(flatAttr), names(r_e$tableHeaders))
	lapply(missingSets, function (x) {
			r_e$tablePrefix[[x]] <- character(0)
			r_e$tableHeaders[[x]] <- character(0)
	})
	r_e$tablePrefix[[root]] <- character(0)

	# Function to get information about children nodes
	getChildren <- function(root, flat) {

		x <- flat[[root]]

		children <- c()
		for(it in 1:length(x))
			if(length(x[[it]]) == 1)
				children <- c(children, x[[it]])
		return (children)
	} 

	# Get tree information
	treeStruct <- lapply(names(flat), getChildren, flat)
	names(treeStruct) <- names(flat)

	# Get prefix length information
	prefixLens <- lapply(r_e$tablePrefix, function(x) length(x))
	names(prefixLens) <- names(r_e$tablePrefix)

	# Unlist couple of metadata
	levelDims <- unlist(r_e$levelDims)
	prefixLens <- unlist(prefixLens)

	# Process types that are still unresolved
	for(nm in names(r_e$tableHeaders)) {
		hd <- r_e$tableHeaders[[nm]]
		tp <- r_e$tableTypes[[nm]]
		if(length(hd)) {
			for(it in 1:length(hd)){
				bla <- NULL
				if(!grepl("xs[:]|xsd[:]", tp[it])) {
					# Exclude obvious types
					if (tp[it] == "IDREF" || tp[it] == "IDREFType" || tp[it] == "ID") {
						bla <- "xs:string"
					} else {
						# Try the recursive search
						bla <- traceType(hd[it])
					}
					if(is.null(bla[1]) || is.na(bla[1])) stop("Error in determining types!")
					r_e$tableTypes[[nm]][it] <- bla[1]
				}
			}
		}
	}

	xsdObject <- list() 
	xsdObject[["root"]] <- root
	xsdObject[["treeStruct"]] <- treeStruct;
	xsdObject[["tableTypes"]] <- r_e$tableTypes;
	xsdObject[["tableHeaders"]] <- r_e$tableHeaders;
	xsdObject[["prefixLens"]] <- prefixLens;
	xsdObject[["levelDims"]] <- levelDims;

	return(xsdObject)
}

#' @importFrom xml2 xml_attrs read_xml xml_find_all xml_find_num xml_ns xml_text xml_add_child xml_attr xml_children
#' @importFrom utils tail
createXsdObject <- function(xsdFile) {

	# Check if XSD exists
	message(paste("Using:", xsdFile))

	# If not exists
	if(!file.exists(xsdFile)) {
		# Get path from local environment
		#fpath <- get("fpath", envir = localEnv)
		fpath <- getRstoxDataDefinitions("fpath")
		xsdFilePath <- paste0(fpath, "/", basename(xsdFile))
		if(!file.exists(xsdFilePath)) {
			message(paste("It seems that", xsdFile, "does not exist or the format is not supported."))
			return(NULL)
		}
	} else {
		xsdFilePath <- xsdFile
	}

	# Parse XSD
	xsdObj <- read_xml(xsdFilePath)

	# Get metadata based on XSD
	metaData <- processXSD(xsdObj, dirname(xsdFile))

	# Process XML file
	ret <- processMetadata(metaData$flat, metaData$flatAttr, metaData$rootInfo, xsdFile, metaData$doc)

	return(ret)
}


#' @importFrom xml2 xml_child read_html xml_find_all
autodetectXml <- function(xmlFile, xsdObjects, verbose) {

	# Read first 500 characters
	tmpText <- readChar(xmlFile, 500)
	bits <- read_html(tmpText)

	# Getting encoding
	tmpG1 <- regexpr('encoding="\\K[^"]*', tmpText, perl=T)
	if(tmpG1 > -1) {
		tmpG2 <- tmpG1 + attr(tmpG1, "match.length") - 1
		xmlEnc <- substr(tmpText, tmpG1, tmpG2)
	} else {
		xmlEnc <- NULL
	}

	# Getting XSD information
	# Peek xmlns
	if(verbose)
		message("Try to use XML namespace")

	tmpG1 <- regexpr('xmlns="\\K[^"]*', tmpText, perl=T)
	if(tmpG1 > -1) {
		tmpG2 <- tmpG1 + attr(tmpG1, "match.length") - 1
		xmlXsd <- substr(tmpText, tmpG1, tmpG2)
		xmlXsd <- paste0(tail(unlist(strsplit(xmlXsd, "/")), 2), collapse = "")
        } else {
		xmlXsd <- NULL
	}

	if(paste0(xmlXsd, ".xsd") %in% names(xsdObjects))
		return(list(xsd = xmlXsd, encoding = xmlEnc))

	if(verbose)
		message("Do manual detection")

	# Do manual detection
	# Later We need to distinguish Biotic v3&v3.1, Biotic v1.4&earlier
	if( length(xml_find_all(bits, "//mission[@startyear]")) )
		xmlXsd <- "nmdbioticv3"
	else if( length(xml_find_all(bits, "//mission[@year]")) )
		xmlXsd <- "nmdbioticv1.4"
	else if( length(xml_find_all(bits, "//biotic")) )
                xmlXsd <- "icesBiotic"
	else if( length(xml_find_all(bits, "//echosounder_dataset")) )
		xmlXsd <- "nmdechosounderv1"
	else if( length(xml_find_all(bits, "//acoustic")) )
		xmlXsd <- "icesAcoustic"
	else if( length(xml_find_all(bits, "//seddellinje")) )
		xmlXsd <- "landingerv2"
	else
		xmlXsd <- NULL

	return(list(xsd = xmlXsd, encoding = xmlEnc))
}

#' @importFrom data.table rbindlist setnames
#' @importFrom xml2 as_list read_xml xml_find_all
getIcesVocabulary <- function(xmlFile) {

	# Read only first 10e4 character
	tmpText <- readChar(xmlFile, 10e4)
	xmlObj <- read_html(tmpText)

	# Apply transformation to get the vocabulary translation table
	ret <- rbindlist(lapply(as_list(xml_find_all(xmlObj, "//vocabulary/*/code")), function(x) if(length(x)>0) return(list(attr(x, "id"), ifelse(length(x) > 0, unlist(x), NA), attr(x, "codetype")))))
	setnames(ret, c("id", "value", "codetype"))
	return(ret)
}
