

quefts_fert <- function() {
	list(N=60, P=10, K=60)
}


quefts_biom <- function() {
	list(leaf_att=2500, stem_att=2500, store_att=5000, SeasonLength=120);
}


quefts_soil <- function() {
	list(N_base_supply=60, P_base_supply=10, K_base_supply=60,
	     N_recovery=0.5, P_recovery=0.1, K_recovery=0.5,
		 #N_fertilizer=0, P_fertilizer=0, K_fertilizer=0,
		 UptakeAdjust = matrix(c(0, 0, 40, 0.4, 80, 0.7, 120, 1, 240, 1.6, 360, 2, 1000, 2), ncol=2, byrow=TRUE)
		)
}


quefts_crop <- function(name="") {

	if (name == "") {
		d <- list(	
			NminStore=0.0095, NminVeg=0.004,  NmaxStore=0.022,  NmaxVeg=0.0125, 
			PminStore=0.0017, PminVeg=0.0004, PmaxStore=0.0075, PmaxVeg=0.003, 
			KminStore=0.002,  KminVeg=0.005,  KmaxStore=0.006,  KmaxVeg=0.02, 
			Yzero=400, Nfix=0)	
	} else {
		f <- system.file("extdata/quefts_crop_pars.csv", package="Rquefts")
		x <- utils::read.csv(f, stringsAsFactors=FALSE)
		if (name %in% x$crop) {
			d <- as.list(x[x$crop == name, ])
		} else {
			m <- matrix(x$crop, nrow=7)
			m <- apply(m, 2, function(i) paste(i, collapse=", "))
			m <- paste(m, collapse="\n")
			stop(name, " not found. Choose from:\n", m)
		}
	}	
	return(d)
}
				

setMethod ("show" , "Rcpp_QueftsModel", 
	function(object) {
		utils::str(object)
	}
)	


quefts <- function(soil, crop, fert, biom) {
	if (missing(soil)) { soil <- quefts_soil() }
	if (missing(crop)) { crop <- quefts_crop() }
	if (missing(fert)) { fert <- quefts_fert() }
	if (missing(biom)) { biom <- quefts_biom() }	
	x <- QueftsModel$new()
	crop(x) <- crop
	soil(x) <- soil
	fert(x) <- fert
	biom(x) <- biom
	x
}


setMethod("run", signature("Rcpp_QueftsModel"), 
	function(x, soil=NULL, yield=NULL, ...) {
		if (is.null(soil) | is.null(yield)) {
			x$run()
		} else {
			if (is.matrix(soil)) {
				x$runbatch(soil[, "Ns"], soil[, "Ps"], soil[, "Ks"], yield[])		
			} else {
				x$runbatch(soil[["Ns"]][], soil[["Ps"]][], soil[["Ks"]][], yield[])		
			}
		}
	}
)



setMethod("[", c("Rcpp_QueftsModel", "character", "missing"),
function(x, i, j, ... ,drop=TRUE) {
	sapply(i, function(y) eval(parse(text = paste0("x$", y))))
})

setMethod("[", c("Rcpp_QueftsModel", "character", "character"),
function(x, i, j, ... ,drop=TRUE) {
	if (i=="model") {
		sapply(j, function(y) eval(parse(text = paste0("x$", y))))
	} else {
		j <- paste0(i, "$", j)
		sapply(j, function(y) eval(parse(text = paste0("x$", y))))
	}
})

setMethod("[", c("Rcpp_QueftsSoil", "character", "missing"),
function(x, i, j, ... ,drop=TRUE) {
	sapply(i, function(y) eval(parse(text = paste0("x$", y))))
})


setMethod("[", c("Rcpp_QueftsCrop", "character", "missing"),
function(x, i, j, ... ,drop=TRUE) {
	sapply(i, function(y) eval(parse(text = paste0("x$", y))))
})


setReplaceMethod("[", c("Rcpp_QueftsModel", "character", "missing"),
	function(x, i, j, value) {
		sapply(i, function(y) eval(parse(text =  paste0("x$", y, " <- ", value))))
		return(x)
	}
)

setReplaceMethod("[", c("Rcpp_QueftsCrop", "character", "missing"),
	function(x, i, j, value) {
		eval(parse(text = paste0("x$crop", i, " <- ", value)))
		return(x)
	}
)

setReplaceMethod("[", c("Rcpp_QueftsSoil", "character", "missing"),
	function(x, i, j, value) {
		eval(parse(text = paste0("x$soil", i, " <- ", value)))
		return(x)
	}
)

setReplaceMethod("[", c("Rcpp_QueftsModel", "character", "character"),
	function(x, i, j, value) {
		if (i=="model") {	
			eval(parse(text = paste0("x$", j, " <- ", value)))
		} else {
			eval(parse(text = paste0("x$", i, "$", j, " <- ", value)))
		}
		return(x)
	}
)

