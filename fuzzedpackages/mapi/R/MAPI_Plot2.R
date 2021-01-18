#' @title Function MAPI_Plot2
#' @export
#' 
#' @description Plot a MAPI analysis result with ggplot
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
#' @return Returns the ggplot object.
#' 
#' @examples
#' \dontrun{
#' library(ggplot2)
#' data("metric")
#' data("samples")
#' resu <- MAPI_RunAuto(samples, metric, crs=3857, nbPermuts = 1000)
#' tails <- MAPI_Tails(resu)
#' pl <- MAPI_Plot2(resu, tails=tails, samples=samples)
#' # Save to image
#' ggsave("mapiPlotOutput.png", plot=pl)
#' }


MAPI_Plot2 <- function(resu, tails=NULL, samples=NULL, pal=c("#994000", "#CC5800", "#FF8F33", "#FFAD66", "#FFCA99", "#FFE6CC", "#FBFBFB", "#CCFDFF", "#99F8FF", "#66F0FF", "#33E4FF", "#00AACC", "#007A99"), shades=20, main="", upper=TRUE, lower=TRUE, upper.border="black", lower.border="gray") {

	# Check dependencies for plotting
	if (!requireNamespace("sf", quietly=TRUE)) { stop("Package \"sp\" needed for this function to work. Please install it.", call. = FALSE) }
	if (!requireNamespace("ggplot2", quietly=TRUE)) { stop("Package \"ggplot2\" needed for this function to work. Please install it.", call. = FALSE) }
	
	crs <- sf::st_crs(resu)

	myMin <- min(resu$avg_value, na.rm=TRUE)
	myMax <- max(resu$avg_value, na.rm=TRUE)
	myBreaks <- signif(seq(myMin, myMax, (myMax-myMin)/shades), digits=3)
	
	# Plot MAPI values
	pl <- ggplot2::ggplot() + 
		ggplot2::geom_sf(data=resu, ggplot2::aes_string(fill="avg_value"), colour=NA) + 
		ggplot2::scale_fill_gradientn(colors=pal, name="MAPI values", breaks=myBreaks, guide=ggplot2::guide_colourbar(order=3, barheight=shades)) 
	pl <- pl + ggplot2::labs(title=main, x="Longitude", y="Latitude")

	# Add tails, if any...
	if (any(class(tails)=="sf" & (upper | lower))) {
		tails$tail <- factor(tails$tail, levels=c("lower", "upper"))
		if (!lower) {tails <- tails[tails$tail!="lower", ]}
		if (!upper) {tails <- tails[tails$tail!="upper", ]}
		if (nrow(tails)>0) {
			pl <- pl + ggplot2::geom_sf(data=tails, ggplot2::aes_string(color="tail"), size=1.0, fill=NA, shape=NA, show.legend="polygon", inherit.aes=FALSE) +
				ggplot2::scale_color_manual(name="Significant areas", values=c(lower.border, upper.border), breaks=c("lower", "upper"), labels=c("lower tail", "upper tail"), drop=TRUE, guide=ggplot2::guide_legend(order=2, override.aes=list(shape=NA)))
		}
	}
	
	# Add samples, if any...
	if(any(class(samples)=="data.frame")) {
		samp0 <- sf::st_as_sf(samples, coords=c("x","y"), crs=crs, remove=FALSE)
		if (exists("samp0") && (nrow(samp0) > 0)) {
			samp0$col = as.factor("x")
			pl <- pl + ggplot2::geom_sf(data=samp0, mapping=ggplot2::aes_string(shape="col"), show.legend="point", inherit.aes=FALSE) +
				ggplot2::scale_shape_discrete(name=NULL, breaks=c("x"), labels=c("Sampling points"), guide=ggplot2::guide_legend(order=1, override.aes=list(linetype="blank", color="black", shape=16)))
		}
	}
	
	return(pl)
}
