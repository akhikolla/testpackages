# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# The bigMap Package for R.

# Copyright (c) 2018, Joan Garriga <jgarriga@ceab.csic.es>, Frederic Bartumeus <fbartu@ceab.csic.es> (Blanes Centre for Advanced Studies, CEAB-CSIC).

# bigMap is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

# bigMap is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

# You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses.
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Define package environment
bigMap.Env <- new.env()

# Set package environment vars
# default output-data path
mybdm <- '~/bigMap/mybdm/'
assign(mybdm, mybdm, envir = bigMap.Env)
# default local machine name or IP address
local <- 'xxx.xxx.xxx.xxx'
assign(local, local, envir = bigMap.Env)

# Locked package variables
# example-data path
# exbdm <- paste(system.file('extdata', package = 'bigMap'), '/', sep = '')


#' Set/get default path for \var{mybdm}
#'
#' @param path Path to \var{mybdm}.
#'
#' @return The current path value to \var{mybdm}
#'
#' @examples
#'
#' # --- set default path for \var{mybdm}
#' bdm.mybdm('~/mybdm')

bdm.mybdm <- function(path = NULL)
{
	if (!is.null(path)){
		if (Sys.info()['sysname'] != 'Windows') {
			if (substr(path, nchar(path), nchar(path)) != '/'){
				path <- paste(path, '/', sep = '')
			}
		}
		assign(mybdm, path, envir = bigMap.Env)
	}
	return(get(mybdm, envir = bigMap.Env))
}


#' Set/get default local machine name or IP address
#'
#' @param dest Name or IP address of the local machine.
#'
#' @return The current value of \var{local}
#'
#' @examples
#'
#' # --- set default value of \var{local}
#' bdm.local('xxx.255.0.0')
#' bdm.local('mymachine.mydomain.cat')

bdm.local <- function(dest = NULL)
{
	if (!is.null(dest)){
		assign(local, dest, envir = bigMap.Env)
	}
	return(get(local, envir = bigMap.Env))
}


# ------------------------------------------------------------------------------
# +++ Auxiliar: get temp filename (for Linux based systems !!!)
# ------------------------------------------------------------------------------

tName.get <- function(pattern){
	t <- strsplit(tempfile(pattern = pattern), '/', fixed = T)[[1]]
	fBin <- paste(t[length(t)], '.bin', sep = '')
	fDesc <- paste(t[length(t)], '.desc', sep = '')
	# fPath <- paste(t[1:(length(t)-1)], collapse = '/')
	fPath <- '~/'
	list(path = fPath, bin = fBin, desc = fDesc)
}
