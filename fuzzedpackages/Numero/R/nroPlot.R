nroPlot <- function(
    topology,
    colors,
    labels=NULL,
    subplot=NULL,
    interactive=FALSE,
    clear=NULL) {

    # Convert inputs to matrices.
    if(is.factor(colors)) colors <- as.character(colors)
    if(is.factor(labels)) labels <- as.character(labels)
    if(is.vector(colors)) colors <- as.matrix(colors)
    if(is.vector(labels)) labels <- as.matrix(labels)
    if(nrow(colors)*ncol(colors) < 1) {
        warning("Empty input.")
        return(NULL)
    }

    # Default labels.
    if(is.null(labels)) labels <- matrix(nrow=0, ncol=0)
    if(nrow(labels) < 1) {
        labels <- matrix("", nrow=nrow(colors), ncol=ncol(colors))
        colnames(labels) <- colnames(colors)
    }

    # Default column names.
    if(length(colnames(colors)) < 1)
        colnames(colors) <- paste(1:ncol(colors))
    if(length(colnames(labels)) < 1)
        colnames(labels) <- colnames(colors)
   
    # Check if input is a list.
    if(is.list(topology) && !is.data.frame(topology))
        topology <- topology$topology

    # Check inputs.
    if(nrow(topology) != nrow(colors)) stop("Incompatible colors.")
    if(nrow(topology) != nrow(labels)) stop("Incompatible labels.")
    if(ncol(colors) != ncol(labels)) stop("Incompatible inputs.")
    subplot <- as.integer(subplot)

    # Make sure all data are available.
    topology <- data.frame(topology, stringsAsFactors=FALSE)
    if(is.null(topology$REGION)) topology$REGION <- ""
    if(is.null(topology$REGION.label)) topology$REGION.label <- ""
    if(is.null(topology$REGION.color)) topology$REGION.color <- ""

    # Check that all labels are single characters.
    flags <- sapply(topology$REGION.label, nchar)
    if(sum(flags > 1) > 0) {
        labels <- matrix("", nrow=nrow(labels), ncol=ncol(labels))
        warning("Multi-character region labels not supported.")
    }

    # Populate values for compatibility with existing code.
    values <- labels

    # Apply label visibility for compatibility.
    visible <- attr(labels, "visible")
    for(j in 1:ncol(labels)) {
        mask <- which(visible[,j] == FALSE)
        labels[mask,j] <- ""
    }

    # Collect parameters.
    param <- list()
    param$colors <- colors
    param$labels <- labels
    param$values <- values
    param$interactive <- as.logical(interactive[[1]])
    param$clear <- param$interactive
    if(!is.null(clear)) param$clear <- clear

    # Make variable names of uniform length.
    vars <- paste(colnames(colors), "                ")
    vars <- sapply(vars, substr, start=1, stop=16)
    param$variables <- as.character(vars)

    # Determine subplot geometry.
    if(length(subplot) < 2) {
        subplot <- sqrt(ncol(colors) + 1)
        subplot <- floor(c(subplot, subplot))
	while((subplot[1])*(subplot[2]) < ncol(colors))
	    subplot[2] <- (subplot[2] + 1)
    }
    subplot <- pmax(subplot, 1, na.rm=TRUE)
    param$nrows <- subplot[[1]]
    param$ncols <- subplot[[2]]

    # Set font size adjustment.
    param$plotcap <- (param$nrows)*(param$ncols)
    param$cexcoeff <- 1.3/log(0.8*sqrt(param$plotcap) + 1.72)

    # Set plot size parameters.
    param$xgap <- 0.5
    param$ygap <- 0.8
    param$rmax <- max(topology$RADIUS2, na.rm=TRUE)
    param$wplot <- 2*(param$rmax + param$xgap)
    param$hplot <- 2*(param$rmax + param$ygap)
    param$xbounds <- c(0, (param$ncols)*(param$wplot))
    param$ybounds <- c(0, (param$nrows)*(param$hplot))
    if(param$interactive == FALSE)
	param$ybounds[2] <- (param$ybounds[2] - 0.7)

    # Set colormap for highlights.
    param$colormap <- nroPlot.colormap()

    # Trigger a new plot.
    param$trigger <- 0
    param$nplots <- 0

    # Create base figure.
    param <- nroPlot.multi(topology, param)

    # Launch interactive mode.
    if(param$interactive) {
        topology <- nroPlot.interface(topology, param)
        if(param$clear) nroPlot.multi(topology, param) # RStudio remnants
    }
    return(topology)
}

#----------------------------------------------------------------------------

nroPlot.colormap <- function() {

    # Organize color palette (for 26 letters).
    hues <- matrix(grDevices::rainbow(25), nrow=5, ncol=5)
    hues <- grDevices::col2rgb(as.vector(t(hues)))/255.0
    red <- c(hues["red",], 0.3)
    green <- c(hues["green",], 0.2)
    blue <- c(hues["blue",], 0.8)

    # Estimate perceived brightness.
    lums <- (0.3*red + 0.5*green + 0.1*blue)/3
    lums <- lums/max(lums)

    # Subdue brightest colors.
    f <- sqrt(1.4 - lums)
    red <- pmin(f*red, 1)
    green <- pmin(f*green, 1)
    blue <- pmin(f*blue, 1)

    # Convert to hexa format.
    colrs <- grDevices::rgb(red, green, blue)

    # Add transparency.
    redmask <- which((red >= green) & (red >= blue))
    greenmask <- which((green >= red) & (green >= blue))
    bluemask <- which((blue >= red) & (blue >= green))
    colrs[redmask] <- paste(colrs[redmask], "50", sep="")
    colrs[greenmask] <- paste(colrs[greenmask], "70", sep="")
    colrs[bluemask] <- paste(colrs[bluemask], "50", sep="")
    return(colrs)
}

#----------------------------------------------------------------------------

nroPlot.exit <- function(elements, subgrp, name, param=NULL) {

   # Check if anything to do.
   if(length(subgrp) < 1) {
      cat("Subgroup '", name, "' is empty.\n", sep="")
      return(NULL)
   }

   # Check if name already exists.
   subnames <- elements$REGION
   if(sum(subnames == name) > 0) {
       cat("Subgroup '", name, "' already exists.", sep="")
       return(NULL)
   }

   # Check if name is reserved.
   if(name == "not_selected") {
       cat("Cannot use '", name, "' for a selected subgroup.", sep="")
       return(NULL)
   }

   # Use a subdued color palette.
   if(is.null(param)) param <- list(colormap=nroPlot.colormap())
   cmap <- param$colormap
   nmax <- length(cmap)

   # Find the next available subgroup label.
   label <- setdiff(LETTERS[1:nmax], elements$REGION.label)
   if(length(label) < 1) stop("Too many subgroups.")   
   if(length(label) > 1) label <- label[1]

   # Set the subgroup attributes.
   elements$REGION[subgrp] <- name
   elements$REGION.label[subgrp] <- label
   elements$REGION.color[subgrp] <- cmap[which(LETTERS == label)] 
   return(elements)
}

#----------------------------------------------------------------------------

nroPlot.identify <- function(topology, param) {
    pos <- graphics::locator(n=1)
    if(is.null(pos)) return(NA)

    # Check if within active area.
    xbounds <- param$xbounds
    ybounds <- param$ybounds
    if(pos$x < xbounds[1]) return(NA)
    if(pos$x > xbounds[2]) return(NA)
    if(pos$y < ybounds[1]) return(NA)
    if(pos$y > ybounds[2]) return(NA)

    # Determine the distance from lower left corner.
    xp <- (pos$x - xbounds[1])
    yp <- (pos$y - ybounds[1])

    # Check if over a coloring
    wplot <- param$wplot
    hplot <- param$hplot
    noccrows <- floor((param$nplot)/(param$ncols))
    nocctail <- (param$nplot)%%(param$ncols)
    if(((ybounds[2] - yp) > noccrows*hplot) & (xp > nocctail*wplot))
        return(NA)

    # Remove multiples of plot sizes.
    xp <- (xp - floor(xp/wplot)*wplot)
    yp <- (yp - floor(yp/hplot)*hplot)

    # Adjust for plot center.
    xp <- (xp - 0.5*wplot)
    yp <- (yp - 0.5*hplot)

    # Find the ring that contains the point.
    radius <- sqrt(xp*xp + yp*yp);
    ring <- which((radius >= topology$RADIUS1) &
                  (radius <= topology$RADIUS2))
    if(length(ring) < 1) return(NA)

    # Find the closest point on the ring (flipped y-coordinates).
    dx <- (topology[ring, "X"] - xp)
    dy <- (topology[ring, "Y"] + yp)
    delta <- sqrt(dx*dx + dy*dy)
    ind <- ring[which.min(delta)]
    return(ind)
}

#----------------------------------------------------------------------------

nroPlot.interface <- function(elements, param) {
    colrs <- param$colors
    labls <- param$labels

    # Original highlights.
    halocolrs <- elements$REGION.color
    halolabls <- elements$REGION.label

    # Greeting text.
    greet1 <- "\nPlease click on subplots to select subgroups."
    greet2 <- "\nTo update or exit, click outside the map colorings.\n"
    greet <- paste(greet1, greet2, sep="")

    # Respond to mouse clicks.
    cat(greet); subgrp <- c()
    while(TRUE) {
        ind <- nroPlot.identify(elements, param)
	if(is.na(ind)) {

            # Exit interactive mode.
            if(length(subgrp) < 1) {
	       ans <- ""
	       while(ans != "n") {
	           ans <- readline("End session (y/n): ")
	           if(ans == "y") {
		       skipped <- which(elements$REGION == "")
                       elements$REGION[skipped] <- "not_selected"
		       return(elements)
		   }
	       }
	       cat(greet); next
	    }

            # Ask for subgroup identifier.
	    res <- NULL; name <- ""
	    while(is.null(res)) {
	        cat("\nSubgroup name (empty to cancel, 'q' to exit):\n") 
                name <- readline()
		if(nchar(name) < 1) break
		if((name == "q") || (name == "'q'")) {
		    mask <- which(elements$REGION.label == "?")
		    vars <- c("REGION","REGION.label","REGION.color")
		    elements[mask,vars] <- ""
		    return(elements)
                }
                res <- nroPlot.exit(elements, subgrp, name, param)
            }
            if(nchar(name) < 1) {cat(greet); next}

            # Print subgroup sizes.
	    t <- table(res$REGION)
            t <- data.frame(NAME=names(t), N.district=as.integer(t))
            print(t, row.names=FALSE)

            # Confirm subgroup assignment.
            while(TRUE) {
	       ans <- readline("Confirm subgroup (y/n): ")
	       if(ans == "n") break
	       if(ans == "y") {
	           param <- nroPlot.multi(res, param)
                   halocolrs <- res$REGION.color
                   halolabls <- res$REGION.label
	           elements <- res
		   subgrp <- c()
		   break
	       }
            }
            cat(greet); next
        }

        # Check if already selected.
        if(sum(subgrp == ind) > 0) {
	    subgrp <- setdiff(subgrp, ind)
            elements$REGION.color[ind] <- halocolrs[ind]
            elements$REGION.label[ind] <- halolabls[ind]   
	}
        else {
	    subgrp <- c(subgrp, ind)
            elements$REGION.color[ind] <- "#00000060"
            elements$REGION.label[ind] <- "?"
	}

        # Print values from the currently selected district.
        if(!is.null(param$values)) {
	    vals <- as.double(param$values[ind,])
	    vals <- formatC(vals, format="g", digits=3, flag="#")
	    for(j in 1:length(vals))
	      cat("\n", param$variables[j], " ", vals[j], sep="")
            cat("\n")
	}

        # Redraw plots.
	param <- nroPlot.multi(elements, param, targets=ind)
    }
    return(NULL)
}

#----------------------------------------------------------------------------

nroPlot.multi <- function(elements, param, targets=c()) {
    colrs <- param$colors
    labls <- param$labels

    # Shorthands for parameters.
    colrs <- param$colors
    labls <- param$labels
    nrows <- param$nrows
    ncols <- param$ncols
    wplot <- param$wplot
    hplot <- param$hplot
    rmax <- param$rmax
    xgap <- param$xgap
    ygap <- param$ygap
    xbounds <- param$xbounds
    ybounds <- param$ybounds
    cexcoeff <- param$cexcoeff

    # Check if full update is needed.
    resetflag <- ((length(targets) < 1) | (param$trigger < 1))

    # Create a new plot.
    if(resetflag) {
        if(param$clear) grDevices::graphics.off()
        graphics::par(pty="m", mar=c(0,0,0,0))
        graphics::plot(x=NA, y=NA, asp=1, xlab=NA, ylab=NA,
             axes=FALSE, xlim=xbounds, ylim=ybounds)
	param$trigger <- 20
	targets <- c()
    }

    # Draw map colorings.
    nplots <- 0
    titles <- colnames(colrs)
    for(j in 1:ncol(colrs)) { 
        if(j > nrows*ncols) {
            warning("Too many subplots to fit in figure.")
            param$warn.subplot = TRUE
	    break
        }
	
        # Determine position offsets.
        dx <- ((nplots%%ncols)*wplot + rmax + xgap)
        dy <- (floor(nplots/ncols)*hplot + rmax + ygap)
        dy <- (ybounds[2] - dy)

        # Set colrs and labls.
	elements$COLOR <- colrs[,j]
	elements$LABEL <- labls[,j]

        # Set label colors.
	elements$LABEL.color <- "#00000060"
        #elements$LABEL.color <- "black"
        #lum <- grDevices::col2rgb(elements$COLOR)
        #lum <- apply(lum, 2, stats::median, na.rm=TRUE)
        #elements$LABEL.color[which(lum < 90)] <- "white"

        # Create subplot.
        nroPlot.single(elements, c(dx, dy), scale=cexcoeff, targets)
	nplots <- (nplots + 1)

        # Write title text if a full refresh.
	if(!resetflag) next
        graphics::text(x=dx, y=(dy + rmax + 0.5),
	               labels=titles[j], cex=cexcoeff)
    }

    # Update the event counter.
    param$nplots <- nplots
    param$trigger <- (param$trigger - 1)
    return(param)
}

#----------------------------------------------------------------------------

nroPlot.single <- function(elements, offsets=c(0,0), scale=1, targets=NULL) {

    # Draw background circle if untargeted refresh.
    rmax <- max(elements$RADIUS2, na.rm=TRUE)
    angC <- pi*(1:10)/5
    flexC <- rep(-0.77, length(angC))
    if(is.null(targets)) {
        graphics::xspline(x=(rmax*cos(angC) + offsets[1]),
            y=(rmax*sin(angC) + offsets[2]), shape=flexC,
	    lwd=3*scale, border="grey70", col="grey90", open=FALSE)
    }

    # Default selection mask.
    if(length(targets) < 1)
        targets <- 1:nrow(elements)

    # Find neighboring districts if a single target with a label.
    if(length(targets) == 1) {
        if(elements$LABEL[targets] != "") {
            xdist <- (elements$X - elements$X[targets])
            ydist <- (elements$Y - elements$Y[targets])
            delta <- sqrt(xdist^2 + ydist^2)
            targets <- which(delta < 1.5)
	}
    }

    # Select targets.
    elem <- elements[targets,]

    # Prepare paths.
    rad1 <- 0.999*(elem$RADIUS1)
    rad2 <- 1.001*(elem$RADIUS2)
    ang1 <- (pi/180)*(elem$ANGLE1 - 0.2)
    ang2 <- (pi/180)*(elem$ANGLE2 + 0.2)
    ang12 <- 0.5*(ang1 + ang2)
    radii <-  cbind(rad1,  rad1, rad1, rad2,  rad2, rad2, rad1)
    angles <- cbind(ang1, ang12, ang2, ang2, ang12, ang1, ang1)
    flex <-   c(       0,    -1,    0,    0,    -1,    0,    0)

    # Apply offsets and flip y-coordinates.
    elem$X <- (elem$X + offsets[1])
    elem$Y <- (offsets[2] - elem$Y)
    xcoord <- (radii*cos(angles) + offsets[1])
    ycoord <- (offsets[2] - radii*sin(angles))

    # Draw elements.
    colrs <- elem$COLOR
    for(i in 1:nrow(elem)) {

        # Circle.
        if(rad1[i] <= 0.0) {
	    graphics::xspline(x=((rad2[i])*cos(angC) + offsets[1]),
	        y=((rad2[i])*sin(angC) + offsets[2]),
		shape=flexC, col=colrs[i], border=NA, open=FALSE)
	    next
        }
        
        # Sector slice.
        graphics::xspline(x=xcoord[i,], y=ycoord[i,], shape=flex,
	    col=colrs[i], border=NA, open=FALSE)
    }

    # Draw highlights.
    mask <- which(elem$REGION.color != "")
    if(length(mask) > 0) {
        xp <- elem$X[mask]
        yp <- elem$Y[mask]
	cl <- elem$REGION.color[mask]
	lb <- elem$REGION.label[mask]
        graphics::points(x=xp, y=yp, pch=16, col=cl, cex=2.2*scale)
        graphics::points(x=xp, y=yp, col="white", cex=2.2*scale)
        graphics::text(x=xp, y=yp, cex=0.8*scale, font=2,
	    labels=lb, col="white")
    }
    
    # Write labels not conflicting with highlights.
    mask <- which((elem$LABEL != "") & (elem$REGION.label == ""))
    if(length(mask) > 0) {
        graphics::text(x=elem$X[mask], y=elem$Y[mask],
            labels=elem$LABEL[mask], cex=scale,
            col=elem$LABEL.color[mask])
    }
    return(NULL)
}
