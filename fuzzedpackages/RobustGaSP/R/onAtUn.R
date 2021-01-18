##########################################################################
## start-up and clean-up functions
##
## This software is distributed under the terms of the GNU GENERAL
## PUBLIC LICENSE Version 2, June 2016.
##
## Copyright (C) 2016-present Mengyang Gu, Jesus Palomo, James O. Berger
##    
##########################################################################

.onAttach <- function(...) {
	
	date <- date()
	x <- regexpr("[0-9]{4}", date)
	this.year <- substr(date, x[1], x[1] + attr(x, "match.length") - 1)
	
	# echo output to screen
	packageStartupMessage("#########")
	packageStartupMessage("##\n## Robust Gaussian Stochastic Process, RobustGaSP Package")
	packageStartupMessage("## Copyright (C) 2016-", this.year,
			" Mengyang Gu, Jesus Palomo and James O. Berger", sep="")
	packageStartupMessage("#########")
}

.onUnload <- function(libpath) {
	library.dynam.unload("RobustGaSP", libpath)
}