# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# The bigMap Package for R.

# Copyright (c) 2018, Joan Garriga <jgarriga@ceab.csic.es>, Frederic Bartumeus <fbartu@ceab.csic.es> (Blanes Centre for Advanced Studies, CEAB-CSIC).

# bigMap is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

# bigMap is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

# You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses.
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# ------------------------------------------------------------------------------
# +++ Control graphic parameters
# ------------------------------------------------------------------------------

# Set default graphic parameters
parbdm.pdef <- par(no.readonly = TRUE)

parbdm.set <- function(oma = c(1.5,2,0.5,2), mar = c(0.25,0.5,0.25,0.5), mgp = c(1.5,0.4,0), cex = 1.0, cex.main = 1.0, cex.lab = 0.8, cex.axis = 0.8, bg = '#FFFFFF')
{
	par(oma = oma)
	par(mar = mar)
	par(mgp = mgp)
	par(cex = cex)
	par(cex.main = cex.main)
	par(cex.lab = cex.lab)
	par(cex.axis = cex.axis)
	par(bg = bg)
}

parbdm.def <- function() par(parbdm.pdef)


# ------------------------------------------------------------------------------
# +++ Color palettes
# ------------------------------------------------------------------------------

# <environment: namespace:colorspace>
cs.palette <- function (n, h = c(0, 90), c. = c(80, 30), l = c(30, 90), power = c(0.2,2), fixup = TRUE, gamma = NULL, alpha = 1, ...)
{
    if (!is.null(gamma))
        warning("'gamma' is deprecated and has no effect")
    if (n < 1L)
        return(character(0L))
    h <- rep(h, length.out = 2L)
    c <- rep(c., length.out = 2L)
    l <- rep(l, length.out = 2L)
    power <- rep(power, length.out = 2L)
    rval <- seq(1, 0, length = n)
    rval <- hex(polarLUV(L = l[2L] - diff(l) * rval^power[2L],
        C = c[2L] - diff(c) * rval^power[1L], H = h[2L] - diff(h) *
            rval), fixup = fixup, ...)
    if (!missing(alpha)) {
        alpha <- pmax(pmin(alpha, 1), 0)
        alpha <- format(as.hexmode(round(alpha * 255 + 1e-04)),
            width = 2L, upper.case = TRUE)
        rval <- paste(rval, alpha, sep = "")
    }
    return(rval)
}

# default bdm colour-palette
pltt.dflt <- function(bdm, w = NULL, l = 32, i = 12, layer = 1){
  if (!is.null(bdm$wtt[[layer]]))
  {
    if (is.null(w)) w <- min(bdm$wtt[[layer]]$s, 256)
    h.seq <- floor(seq(1, 256, length.out = w))
    pltt <- sapply(h.seq, function(k) cs.palette(l, h=k)[i])
    pltt.rows <- ceiling(sqrt(length(pltt)))
    pltt <- c(pltt, pltt[1:(pltt.rows**2 -length(pltt))])
    clr.pltt <- t(matrix(pltt, nrow = pltt.rows))
  }
  else clr.pltt <- pltt.get(s = 4)[2]
  return(clr.pltt)
}

# Base Palette
pltt.base <- function(w, l){
	h.seq <- floor(seq(1, 256, length.out=w))
	return(sapply(h.seq, function(h) cs.palette(l, h=h)))
}

# Get Palette
pltt.get <- function(s = 16, l = 16){
	if (s <= 16){
		w <- c(1,2,2,2,3,3,4,4,3,5,6,6,7,7,8,8)[s]
		i <- c(1,1,2,2,2,2,2,2,3,2,2,2,2,2,2,2)[s]
	}
	else{
		w <- min(floor(sqrt(s)), 16)
		i <- min(s%/%w, l)
		while (w < 16 && w*i < s){
			w <- w +1
			i <- s%/%w
			if (i<l && w*i<s) i <- i +1
		}
	}
	pltt <- pltt.base(w, l)[seq(1, l, length.out=(i+1))[1:i],]
	pltt <- rep(pltt, ceiling(s /length(pltt)))[1:s]
	return(pltt)
}

# Palette plot
pltt.plot <- function(pltt, ttl=''){
	l <- length(pltt)
	plot(1, 1, xlab='', ylab='', xlim=c(0,1), ylim=c(1,(l+1)), xaxt='n', yaxt="n", bty="n", type="n", main=ttl, cex.main=0.8)
	rect(0, seq(1,l), 1, seq(2,l+1), col=pltt, border=NA)
}

# heatmap palette (this is a copy of the color.jet() function)
pltt.heat <- function(n, alpha = 1)
{
	if (length(n) > 1 || !is.finite(n))
		stop("'n' must be an integer positive value.", call. = FALSE)
	if (n < 1)
		stop("'n' must be an integer positive value.", call. = FALSE)
	if (length(alpha) > 1 || !is.finite(alpha))
		stop("'alpha' must be an numeric value in the range [0, 1].", call. = FALSE)
	if (alpha < 0 || alpha > 1)
		stop("'alpha' must be an numeric value in the range [0, 1].", call. = FALSE)
	alpha = round(255 * alpha)
	ramp = colorRamp(c("#00008F", "#00009F", "#0000AF", "#0000BF",
		"#0000CF", "#0000DF", "#0000EF", "#0000FF", "#0010FF",
		"#0020FF", "#0030FF", "#0040FF", "#0050FF", "#0060FF",
		"#0070FF", "#0080FF", "#008FFF", "#009FFF", "#00AFFF",
		"#00BFFF", "#00CFFF", "#00DFFF", "#00EFFF", "#00FFFF",
		"#10FFEF", "#20FFDF", "#30FFCF", "#40FFBF", "#50FFAF",
		"#60FF9F", "#70FF8F", "#80FF80", "#8FFF70", "#9FFF60",
		"#AFFF50", "#BFFF40", "#CFFF30", "#DFFF20", "#EFFF10",
		"#FFFF00", "#FFEF00", "#FFDF00", "#FFCF00", "#FFBF00",
		"#FFAF00", "#FF9F00", "#FF8F00", "#FF8000", "#FF7000",
		"#FF6000", "#FF5000", "#FF4000", "#FF3000", "#FF2000",
		"#FF1000", "#FF0000", "#EF0000", "#DF0000", "#CF0000",
		"#BF0000", "#AF0000", "#9F0000", "#8F0000", "#800000"),
		space = "Lab")
	rgb(ramp(seq(0, 1, length = n)), alpha = alpha, maxColorValue = 255)
}

# density colour palette
pltt.pakde <- function(lvls)
{
	pltt1 <- cs.palette(lvls, h = c(0, -100), c. = c(80, 40), l = c(40, 75), power = c(1, 1))[lvls:2]
	pltt2 <- cs.palette(lvls, h = c(0, 90), c. = c(100, 30), l = c(50, 90), power = c(0.2,	1))
	c("#999999", pltt1, pltt2)
}

# class density palette
pltt.class <- function(lvls)
{
	pltt1 <- cs.palette(lvls+1, h = c(-100, 100), c. = c(60, 100), l = c(15, 95), power = c(2, 0.9))
	pltt1[1:lvls]
}

# ------------------------------------------------------------------------------
# +++ optimal multiplot layout (internal)
# ------------------------------------------------------------------------------

layout.get <- function(m)
{
	m <- min(m, 25)
	layout.rows <- 1
	layout.rows <- ceiling(sqrt(m))
	layout.cols <- ceiling(m / layout.rows)
	layout.seq <- c(seq(m), rep(m+1, (layout.rows *layout.cols -m)))
	layout.mtx <- matrix(layout.seq, nrow = layout.rows, byrow = T)
	return(layout.mtx)
}

# ------------------------------------------------------------------------------
# +++ nulL plot (internal)
# ------------------------------------------------------------------------------

plot.null <- function(){
	plot(0, 1, type='n', bty='n', xaxt='n', yaxt='n', xlab='', ylab='')
}
