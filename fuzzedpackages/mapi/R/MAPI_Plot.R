#' @title Function MAPI_Plot
#' @export
#' 
#' @description Plot a MAPI analysis result
#' 
#' @param resu A spatial object of class 'sf' resulting from a MAPI analysis done using 
#'   \code{\link{MAPI_RunAuto}} or \code{\link{MAPI_RunOnGrid}}.
#' @param tails An optional spatial object of class 'sf' resulting from the post-process with 
#'   \code{\link{MAPI_Tails}} of a MAPI analysis result. Default = NULL (no tails shown).
#' @param samples A data.frame with names and geographical coordinates of samples. 
#'   Column names must be: 'ind', 'x', 'y'.  
#'   Optional column 'errRad' with an error radius for sample locations (eg. GPS uncertainty). 
#'   Coordinates must be projected (not latitude/longitude).
#' @param pal A color ramp, eg. from \pkg{RColorBrewer} (default: orange > light gray > blue)
#' @param shades Number of breaks for the color ramp (default 20)
#' @param main Plot title (none by default)
#' @param upper If TRUE and tails is not NULL, upper-tail significant areas are plotted. TRUE by default.
#' @param lower If TRUE and tails is not NULL, lower-tail significant areas are plotted. TRUE by default.
#' @param upper.border Border color of the upper-tail significant area. "black" by default.
#' @param lower.border Border color of the lower-tail significant area. "gray" by default.
#' 
#' @return Returns the "trellis" object.
#' 
#' @examples
#' \dontrun{
#' data("metric")
#' data("samples")
#' resu <- MAPI_RunAuto(samples, metric, crs=3857, nbPermuts = 1000)
#' tails <- MAPI_Tails(resu)
#' pl <- MAPI_Plot(resu, tails=tails, samples=samples)
#' # Open png driver
#' png("mapiPlotOutput.png", width=1000, type="cairo-png")
#' print(pl) # Do plot in file
#' dev.off() # Close driver
#' }


MAPI_Plot <- function(resu, tails=NULL, samples=NULL, pal=c("#994000", "#CC5800", "#FF8F33", "#FFAD66", "#FFCA99", "#FFE6CC", "#FBFBFB", "#CCFDFF", "#99F8FF", "#66F0FF", "#33E4FF", "#00AACC", "#007A99"), shades=20, main=NA, upper=TRUE, lower=TRUE, upper.border="black", lower.border="gray") {
	
	# DEPRECATED: this function is deprecated, use MAPI_Plot2 instead
	.Deprecated("MAPI_Plot2")
	
	# Check dependencies for plotting
	if (!requireNamespace("grDevices", quietly = TRUE)) { stop("Package \"grDevices\" needed for this function to work. Please install it.", call. = FALSE) }
	if (!requireNamespace("sp", quietly = TRUE)) { stop("Package \"sp\" needed for this function to work. Please install it.", call. = FALSE) }
	if (!requireNamespace("latticeExtra", quietly = TRUE)) { stop("Package \"latticeExtra\" needed for this function to work. Please install it.", call. = FALSE) }
	
	# Drop pertutations column before conversion
	resu$permuts <- NULL
	# Convert to Spatial object
	my.rsp <- sf::as_Spatial(resu)
	
	# Build color ramp for a given number of classes
	my.ramp <- grDevices::colorRampPalette(pal)
	# Add a color number column
	my.rsp$col <- cut(my.rsp$avg_value, breaks=shades, include.lowest=TRUE, include.highest=TRUE)
	
	# Plot MAPI values
	pl <- sp::spplot(my.rsp, "col", col="transparent", col.regions=my.ramp(shades), main=main)
	
	# Add upper tail, if any
	if(any(class(tails)=="sf") && upper==TRUE){
		uTail0 <- tails[tails$tail=="upper", ]
		if (exists("uTail0") && nrow(uTail0) > 0) {
			uTail <- sf::as_Spatial(uTail0)
			pl <- pl + latticeExtra::layer(sp::sp.polygons(uTail,  fill=c("transparent"), lwd=2, col=upper.border), data=list(upper.border=upper.border, lower.border=lower.border))
		}
	}
	
	# Add lower tail, if any
	if(any(class(tails)=="sf") && lower==TRUE){
		lTail0 <- tails[tails$tail=="lower", ]
		if (exists("lTail0") && nrow(lTail0) > 0) {
			lTail <- sf::as_Spatial(lTail0)
			pl <- pl + latticeExtra::layer(sp::sp.polygons(lTail,  fill=c("transparent"), lwd=2, col=lower.border), data=list(upper.border=upper.border, lower.border=lower.border))
		}
	}
	
	# Add samples, if any
	if(!is.null(samples)) {
		pl <- pl + latticeExtra::layer(sp::sp.points(sf::as_Spatial(sf::st_as_sf(samples, coords=c("x", "y"), crs=sf::st_crs(resu))), col='black', pch=15, cex=0.5))
	}
	
	return(pl)
}


