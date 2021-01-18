
if (!isGeneric("crop<-")) { setGeneric("crop<-", function(x, value) standardGeneric("crop<-")) }	
if (!isGeneric("soil<-")) { setGeneric("soil<-", function(x, value) standardGeneric("soil<-")) }	
if (!isGeneric("biom<-")) { setGeneric("biom<-", function(x, value) standardGeneric("biom<-")) }	
if (!isGeneric("fert<-")) { setGeneric("fert<-", function(x, value) standardGeneric("fert<-")) }	

if (!isGeneric("run")) { setGeneric("run", function(x, ...) standardGeneric("run")) }	

#if (!isGeneric("quefts")) { setGeneric("quefts", function(soil, crop, fert, biom) standardGeneric("quefts")) }	



setMethod("soil<-", signature("Rcpp_QueftsModel", "list"), 
	function(x, value) {
		parameters <- c("K_base_supply", "K_recovery", "N_base_supply", "N_recovery", "P_base_supply", "P_recovery", "UptakeAdjust")
		nms <- names(value)
		if (!all(parameters %in% nms)) stop(paste("parameters missing:", paste(parameters[!(parameters %in% nms)], collapse=", ")))

		value <- value[parameters]
		nms <- names(value)
		
		lapply(1:length(value), function(i) eval(parse(text = paste0("x$soil$", nms[i], " <- ", value[i]))))
		return(x)
	}
)

setMethod("crop<-", signature("Rcpp_QueftsModel", "list"), 
	function(x, value) {
		parameters <- c("KmaxStore", "KmaxVeg", "KminStore", "KminVeg", "Nfix", "NmaxStore", "NmaxVeg", "NminStore", "NminVeg", "PmaxStore", "PmaxVeg", "PminStore", "PminVeg", "Yzero")
		nms <- names(value)
		
		if (!all(parameters %in% nms)) stop(paste("parameters missing:", paste(parameters[!(parameters %in% nms)], collapse=", ")))
		value <- value[parameters]
		nms <- names(value)
		lapply(1:length(value), function(i) eval(parse(text = paste0("x$crop$", nms[i], " <- ", value[i]))))
		return(x)
	}
)


setMethod("fert<-", signature("Rcpp_QueftsModel", "list"), 
	function(x, value) {
		parameters <- c("N", "P", "K")
		value <- value[parameters]
		nms <- names(value)
		lapply(1:length(value), function(i) eval(parse(text = paste0("x$", nms[i], " <- ", value[i]))))
		return(x)
	}
)

setMethod("biom<-", signature("Rcpp_QueftsModel", "list"), 
	function(x, value) {
		parameters <- c("leaf_att", "stem_att", "store_att", "SeasonLength")
		nms <- names(value)
		if (!all(parameters %in% nms)) stop(paste("parameters missing:", paste(parameters[!(parameters %in% nms)], collapse=", ")))
		value <- value[parameters]
		nms <- names(value)
		lapply(1:length(value), function(i) eval(parse(text = paste0("x$", nms[i], " <- ", value[i]))))
		return(x)
	}
)

