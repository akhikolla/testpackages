##
# GUTS R Definitions.
# soeren.vogel@uzh.ch, carlo.albert@eawag.ch, oliver.jakoby@rifcon.de, alexander.singer@rifcon.de, dirk.nickisch@rifcon.de
# License GPL-2
# 2019-05-24


##
# Function guts_setup(...).
#
guts_setup <- function(C, Ct, y, yt, dist = 'lognormal', model = 'Proper', N = 1000, M = 10000) {

	#
	# Check missing arguments and arguments types (numeric, character).
	#
	args_num_names  <- c('C', 'Ct', 'y', 'yt', 'N', 'M')
	args_char_names <- c('dist', 'model')
	args_num_type   <- c(is.numeric(C), is.numeric(Ct), is.numeric(y), is.numeric(yt), is.numeric(N), is.numeric(M))
	args_char_type  <- c(is.character(dist), is.character(model))
	if ( any( !args_num_type ) ) {
		i <- which(!args_num_type)[1]
		stop( paste( "Argument ", args_num_names[i], " must be numeric.", sep='' ) )
	} else if ( any( !args_char_type ) ) {
		i <- which(!args_num_type)[1]
		stop( paste( "Argument ", args_char_names[i], " must be character.", sep='' ) )
	}


	#
	# Check length of single value arguments.
	#
	args_sin_names  <- c('dist', 'model', 'N', 'M')
	args_sin_len    <- c(length(dist), length(model), length(N), length(M))
	for ( i in seq_along(args_sin_len) ) {
		if ( args_sin_len[i] > 1 ) {
			warning( paste( "Argument ", args_sin_names[i], " must be of length 1, only first element used.", sep='' ) )
			assign( args_sin_names[i], get(args_sin_names[i])[1] )
		}
	}


	#
	# Check concentrations and survivors.
	#
	if ( length(C) < 2 ) {
		stop( 'Vector C must be longer than 1.' )
	} else if ( length(Ct) < 2 ) {
		stop( 'Vector Ct must be longer than 1.' )
	} else if ( length(C) != length(Ct) ) {
		stop( 'Vectors C and Ct must have the same length.' )
	} else if ( Ct[1] != 0 ) {
		stop( 'Vector Ct must start at 0.0.' )
	} else if ( any(diff(Ct) <= 0) ) {
		stop( 'Vector Ct must contain unique values in ascending order.' )
	} else if ( length(y) < 2 ) {
		stop( 'Vector y must be longer than 1.' )
	} else if ( length(yt) < 2 ) {
		stop( 'Vector yt must be longer than 1.' )
	} else if ( length(y) != length(yt) ) {
		stop( 'Vectors y and yt must have the same length.' )
	} else if ( yt[1] != 0 ) {
		stop( 'Vector yt must start at 0.0.' )
	} else if ( any(diff(yt) <= 0) ) {
		stop( 'Vector yt must contain unique values in ascending order.' )
	} else if ( any(diff(y) > 0) ) {
		stop( 'Values in vector y must not ascend.' )
	} else if ( min(c(C, Ct, y, yt)) < 0 ) {
		stop( 'Vectors C, Ct, y, yt must contain non-negative values.' )
	}


	#
	# Check Ct and yt length and, if needed, truncate y and yt.
	#
	Ct.last <- Ct[length(Ct)]
	if ( Ct.last < yt[length(yt)] ) {
		i <- which(yt <= Ct.last)
		y <- y[i]
		yt <- yt[i]
		warning( 'Survivor information at time points later than the latest concentration time point are disregarded.' )
	}

  #
	# Check dist and model.
	# Set experiment code.
	# Set par_pos.
	# Set par.
	# Set wpar.
	#
	experiment <- 11 # Default
	mdist  <- tolower(dist)
	mmodel <- tolower(model)
	par_pos <- 1:5
	
	if (mmodel == "sd") {
	  mdist <- "delta"
	  mmodel <- "proper"
	}
	  
	
	if ( mdist == 'lognormal' ) {
		if ( mmodel == 'proper' ) {
			experiment <- 11
			par_pos <- 1:5
		} else if ( mmodel == 'it' ) {
			experiment <- 12
			par_pos <- c(1:2, 4:5)
		} else {
			stop( 'Model must be either "Proper" or "IT".' )
		}
	} else if ( mdist == "delta" ) {
	  if ( mmodel == 'proper' ) {
	    experiment <- 21
	    par_pos <- 1:4
	  } else if ( mmodel == 'it' ) {
	    experiment <- 22
	    par_pos <- c(1:2, 4)
	  } else {
	    stop( 'Model must be either "Proper" or "IT".' )
	  }
	} else if ( mdist == "loglogistic" ) { 
	  if ( mmodel == 'proper' ) {
	    experiment <- 31
	    par_pos <- 1:5
	  } else if ( mmodel == 'it' ) {
	    experiment <- 32
	    par_pos <- c(1:2, 4,5)
	  } else {
	    stop( 'Model must be either "Proper" or "IT".' )
	  }
	} else if ( mdist == "external" ) {
	  if ( mmodel == 'proper' ) {
	    experiment <- 41
	    par_pos <- 1:3
	  } else if ( mmodel == 'it' ) {
	    experiment <- 42
	    par_pos <- c(1:2)
	  } else {
	    stop( 'Model must be either "Proper" or "IT".' )
	  }
	} else {
	  stop( 'Distribution must be either "lognormal", "loglogistic", "external" or "delta".' )
	}
	par <- rep(NA, length(par_pos))
	wpar <- c(0, 0, .Machine$double.xmax, 0, 0)


	#
	# Check sample length and time grid points.
	#
	if ( N < 3 ) {
		stop( 'N must be greater than 2.' )
	} else if ( M < 2 ) {
		stop( 'M must be greater than 1.' )
	}


	#
	# Build GUTS object for return.
	#
	ret <- structure(
		list(
			'C'     = C,
			'Ct'    = Ct,
			'y'     = y,
			'yt'    = yt,
			'dist'  = dist,
			'model' = model,
			'N'     = N,
			'M'     = M,
			'par'   = par,
			'S'     = rep(NA, length(yt)),
			'D'     = rep(NA, M),
			'LL'    = NA,
			'SPPE'  = NA,
			'squares' = NA
		),
		class      = "GUTS",
		experiment = experiment,
		wpar       = wpar,
		par_pos    = par_pos
	)
	invisible( return( ret ) )

} # End of guts_setup()




##
# Function guts_calc_loglikelihood(...).
guts_calc_loglikelihood <- function(gobj, par, external_dist = NULL) {
	invisible(.Call('_GUTS_guts_engine', PACKAGE = 'GUTS', gobj, par, z_dist = external_dist))
	return(gobj[['LL']])
}






##
# Function guts_calc_survivalprobs(...).
guts_calc_survivalprobs <- function(gobj, par, external_dist = NULL) {
	invisible(.Call('_GUTS_guts_engine', PACKAGE = 'GUTS', gobj, par, z_dist = external_dist))
	return(gobj[['S']])
}


##
# Function guts_report_damage(...).
guts_report_damage <- function(gobj) {
  return(
    data.frame(
      time = seq(min(gobj[['yt']]), max(gobj[['yt']]), length.out = gobj[['M']]),
      damage = gobj[['D']]
      )
  )
}

##
# Function guts_report_sppe(...).
guts_report_sppe <- function(gobj) {
  return(gobj[['SPPE']])
}

##
# Function guts_report_squares(...).
guts_report_squares <- function(gobj) {
  return(gobj[['squares']])
}




##
# Printing.
#

# A small helper for printing and wrapping.
.g_print_help <- function( x, width, digits, exdent = 6, prefix = NULL ) {
	y <- paste( round(x, digits = digits), sep = "", collapse = ", " )
	z <- strwrap( y, width = (width-6), indent = 0, exdent = exdent, initial = prefix )
	return( z )
}

# The actual print function.
.g_print <- function( object, width = getOption('width'), digits = getOption('digits') ) {

	# Header
	cat(
		"\n",
		"GUTS object:\n",
		"============\n",
		sep=""
	)

	# Distribution, Model
	cat( "Distribution: ", object$dist, ", model: ", object$model, ".\n", sep="" )

	# Concentrations, Survivors
	cat( "Concentrations (n=", length(object$C), "), survivors (n=", length(object$y), ")", sep="" )
	if ( length(object$C) > 0 ) {
		cat( ":", sep="\n" )
		cat( .g_print_help(object$Ct, width, digits, prefix="  Ct: "), sep="\n" )
		cat( .g_print_help(object$C,  width, digits, prefix="   C: "), sep="\n" )
	}
	if ( length(object$y) > 0 ) {
		cat( .g_print_help(object$yt, width, digits, prefix="  yt: "), sep="\n" )
		cat( .g_print_help(object$y,  width, digits, prefix="   y: "), sep="\n" )
	} else {
		cat( "\n", sep="" )
	}

	# Sample length, Time grid points
	cat( "Sample length: ", object$N, ", Time grid points: ", object$M, ".\n", sep="" )

	# Parameters
	prf <- paste("Parameters (n=", length(object$par), ")", sep="")
	if ( length(object$par) > 0 ) {
		prf <- paste(prf, ": ", sep="")
		cat( .g_print_help(object$par, width, digits, prefix=prf), sep="\n" )
	} else {
		cat( "\n", sep="" )
	}

	# Survival probabilities
	prf <- paste("Survival probabilities (n=", length(object$S), ")", sep="")
	if ( length(object$S) > 0 ) {
		prf <- paste(prf, ": ", sep="")
		cat( .g_print_help(object$S, width, digits, exdent = 2, prefix = prf), sep="\n" )
	} else {
		cat( "\n", sep="" )
	}

	# Damage
	prf <- paste("Damage (n=", length(object$D), "; subset: 30 regularly spaced values)", sep="")
	if ( length(object$D) > 0 ) {
		prf <- paste(prf, ": ", sep="")
		cat( .g_print_help(object$D[seq(1, length(object$D), length.out = 30)], width, digits, exdent = 2, prefix = prf), sep="\n" )
	} else {
		cat( "\n", sep="" )
	}

	#Sum of squares
	cat( "Sum of squares: ", object$squares, "\n", sep="" )
	
	# Loglikelihood
	cat( "Loglikelihood: ", object$LL, "\n", sep="" )

	#SPPE
	cat( "SPPE: ", object$SPPE, "\n", sep="" )

	# Footer
	cat( "\n", sep="" )
}

# Do print.
print.GUTS <- function(x, ...) {
	out <- .g_print(x)
	return(invisible(out))
}





##
# Setters of Fields.
#
#"[[<-" <- function(x, field, value) {
#	UseMethod("[[<-", x)
#}
"[[<-.GUTS" <- function(x, field, value) {
	stop( "Use function `guts_setup()` for changing fields." )
}
#"$<-" <- function(x, field, value) {
#	UseMethod("$<-", x)
#}
"$<-.GUTS" <- function(x, field, value) {
	stop( "Use function `guts_setup()` for changing fields." )
}





##
# Attributes.
#
#'attr<-' <- function(x, which, value) {
#	UseMethod('attr<-', x)
#}
"attr<-.GUTS" <- function(x, which, value) {
	stop( "Use function `guts_setup()` for changing GUTS object fields." )
}
#'attributes<-' <- function(x, which, value) {
#	UseMethod('attributes<-', x)
#}
"attributes<-.GUTS" <- function(x, which, value) {
	stop( "Use function `guts_setup()` for changing attributes of a GUTS object." )
}
#'mostattributes<-' <- function(x, which, value) {
#	UseMethod('mostattributes<-', x)
#}
"mostattributes<-.GUTS" <- function(x, which, value) {
	stop( "Use function `guts_setup()` for changing attributes of a GUTS object." )
}
