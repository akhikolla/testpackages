nroPlot.save <- function(
    file,
    topology,
    colors,
    labels=NULL,
    subplot=NULL,
    font=1.0) {
  
    # Check if input is a list.
    if(is.list(topology) && !is.data.frame(topology))
        topology <- topology$topology
    if(nrow(topology) < 2) stop("Unusable input.")

    # Prepare topology.
    topo <- topology[,c("X","Y","RADIUS1","RADIUS2","ANGLE1","ANGLE2")]
    topo <- nroRcppMatrix(topo, trim=FALSE)

    # Copy attributes before touching variables.
    contrast <- nroRcppMatrix(attr(colors, "contrast"), trim=FALSE)
    visible <- nroRcppMatrix(attr(labels, "visible"), trim=FALSE)
    
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
    if(is.null(labels)) {
         labels <- matrix("", nrow=nrow(colors), ncol=ncol(colors))
         rownames(labels) <- rownames(colors)
         colnames(labels) <- colnames(colors)
    }
    
    # Default contrast.
    if(nrow(contrast) < 1) {
      contrast <- matrix(0, nrow=nrow(colors), ncol=ncol(colors))
      rownames(contrast) <- rownames(colors)
      colnames(contrast) <- colnames(colors)
    }
    
    # Default visibility.
    if(nrow(visible) < 1) {
      visible <- matrix(1, nrow=nrow(colors), ncol=ncol(colors))
      rownames(visible) <- rownames(colors)
      colnames(visible) <- colnames(colors)
    }

    # Prepare highlights.
    hlights <- data.frame(REGION=rep("(empty)", nrow(topo)),
        REGION.label=rep("", nrow(topo)),
        REGION.color=rep("#000000", nrow(topo)),
        stringsAsFactors=FALSE)
    if(is.finite(match("REGION", colnames(topology))))
        hlights$REGION <- as.character(topology[,"REGION"])
    if(is.finite(match("REGION.label", colnames(topology))))
        hlights$REGION.label <- as.character(topology[,"REGION.label"])
    if(is.finite(match("REGION.color", colnames(topology))))
        hlights$REGION.color <- as.character(topology[,"REGION.color"])

    # Check region labels and colors.
    js.regs <- nroPlotSave.colormap(hlights)

    # Determine subplot geometry.
    if(length(subplot) < 2) {
        subplot <- sqrt(ncol(colors) + 1)
        subplot <- floor(c(subplot, subplot))
        while((subplot[1])*(subplot[2]) < ncol(colors))
            subplot[2] <- (subplot[2] + 1)
    }

    # Number of subplot columns in figure.
    subplot <- as.integer(subplot[1:2])
    if(subplot[1] < 1) stop("Unusable number of plot rows.")
    if(subplot[2] < 1) stop("Unusable number of plot columns.")

    # Check font size.
    if(is.null(font)) font <- 1.0
    font <- as.double(font[[1]])
    if(!is.finite(font)) font <- 1.0
    if((font < 0.1) | (font > 100.0)) 
        stop("Unusable font size.")

    # Set file path.
    fname <- path.expand(file)
    if(nchar(fname) < 1) stop("Empty file name.")
    
    # Generate SVG code.
    svgdoc <- nroPlotSave.svg(topology=topo, colors=colors,
        labels=labels, visible=visible, contrast=contrast,
        hlights=hlights, subplot=subplot, font=font)
    codes <- svgdoc$codes
    boxes <- svgdoc$boxes

    # Import base javascript that enables interactive features.
    jsfile <- system.file("extcode", "circus.js",
        package="Numero", mustWork=TRUE)
    jscode <- readChar(con=jsfile, nchars=1e5)

    # Add visualization data.
    keys <- names(codes)
    map <- data.frame(topo, hlights, stringsAsFactors=FALSE)
    jscode <- c(nroPlotSave.json(map, colnames(map), "TOPOLOGY"), jscode)
    jscode <- c(sprintf("const NDISTRICTS = %d;\n", nrow(topo)), jscode)
    jscode <- c(nroPlotSave.json(labels, keys, "LABELS"), jscode)
    jscode <- c(js.regs, jscode)

    # Add plot identifiers.
    s <- paste0("\"", keys, sep="\"", collapse=",")
    s <- paste0("const SUBPLOTS = [", s, "];\n", collapse="")
    jscode <- c(s, jscode)

    # Encapsulate within a script section.
    jscode <- paste0(jscode, sep="", collapse="")
    jscode <- paste0("\n<script type=\"application/javascript\">\n",
                     jscode, "</script>\n", sep="", collapse="")

    # Create figure file.
    res <- nroPlotSave.write(file, codes, boxes, jscode)
    return(as.numeric(res$nbytes))
}

#---------------------------------------------------------------------------

nroPlotSave.svg <- function(topology, colors, labels,
    visible, contrast, hlights, subplot, font) {

    # Default titles.
    titles <- colnames(colors)
    if(length(titles) != ncol(colors))
        titles <- paste("Column", 1:ncol(colors))

    # Set plot identifiers.
    slots <- 0
    keys <- character()
    for(k in 1:length(titles)) {
	keys[k] <- intToUtf8((65 + slots), multiple=FALSE)

        # Increment counter.
	pos <- length(slots)
	slots[pos] <- (slots[pos] + 1)

        # Propagate increment.
        while(slots[pos] > 25) {
	    if(pos < 2) break
            slots[pos] <- 0
            pos <- (pos - 1)
	    slots[pos] <- (slots[pos] + 1)
        }

        # Check if propagation succeeded.
	if(slots[pos] <= 25) next

        # New slot needed.
	slots <- c(0, 0*slots)
    }
 
    # Maximum number of subplots.
    capacity <- (subplot[1])*(subplot[2])
    nsubcol <- subplot[2]

    # Generate SVG code.
    codes <- c()
    boxes <- c()
    offset <- c(0, 0)
    for(j in 1:length(keys)) {
        if(j > capacity) {
	    warning("Plot capacity exceeded.")
            break
        }

        # Paint map districts.
        res.p <- .Call("nro_circus_paint",
	     as.double(offset),
             as.matrix(topology),
	     as.character(colors[,j]),
	     as.character(keys[j]),
             as.character(titles[j]),
             PACKAGE="Numero")
        if(is.character(res.p)) stop(res.p)

        # Make sure key is the same across components.
        key <- res.p$key

        # Write labels.
        res.w <- .Call("nro_circus_write",
	     as.double(offset),
             as.matrix(topology),
             as.character(labels[,j]),
             as.logical(visible[,j]),
             as.logical(contrast[,j]),
             as.character(key),
	     as.double(font),
             PACKAGE="Numero")
        if(is.character(res.w)) stop(res.w)

        # Add highlights.
        res.s <- .Call("nro_circus_show",
	     as.double(offset),
             as.matrix(topology),
             as.character(hlights$REGION.color),
 	     as.character(hlights$REGION.label),
             as.character(key),
             PACKAGE="Numero")
        if(is.character(res.s)) stop(res.s)

        # Combine code segments.
	subcode <- c("<g>", res.w$code.shadow, res.p$code,
	    res.w$code.label, res.s$code, "</g>")
	subcode <- paste(subcode, collapse="\n")
	boxes <- rbind(boxes, res.p$bbox)
        codes <- c(codes, subcode)

        # Determine the size of the map in device coordinates.
	wplot <- (res.p$bbox[3] - res.p$bbox[1])
        hplot <- (res.p$bbox[4] - res.p$bbox[2])

        # Set the center point of the next plot.
        offset[1] <- (offset[1] + wplot)
	if(j%%nsubcol == 0)
            offset <- c(0.0, (offset[2] + hplot))
    }

    # Attach plot keys to code.
    rownames(boxes) <- keys
    names(codes) <- keys

    # Return results.
    output <- list()
    output$boxes <- boxes
    output$codes <- codes
    return(output)
}

#---------------------------------------------------------------------------

nroPlotSave.write <- function(fname, codes, boxes, jscode) { 

    # Set file path.
    fname <- path.expand(fname)
    if(nchar(fname) < 1) stop("Empty file name.")
    
    # Save document.
    len <- nchar(fname)
    filefmt <- substr(fname, (len - 4), len)
    res <- list(text="0", nbytes=0)
    if(filefmt == ".html") {

        # Import HTML wrapper.
        headfile <- system.file("extcode", "circus.head.html",
	                        package="Numero", mustWork=TRUE)
        tailfile <- system.file("extcode", "circus.tail.html",
	                        package="Numero", mustWork=TRUE)
        htmlhead <- readChar(con=headfile, nchars=1e5)
        htmltail <- readChar(con=tailfile, nchars=1e5)

        # Set up SVG parent elements for subplots.
        codes <- nroPlotSave.figure(codes, boxes)

        # Assemble final document code.
        codes <- c(htmlhead, jscode, codes, htmltail)

        # Save as an HTML document.
        res <- .Call("nro_webpage",
                     as.character(fname),
                     as.character(codes),
                     PACKAGE="Numero")
        if(is.character(res)) stop(res)
    }
    else {

       # Overall bounding box.
       bbox <- rep(NA, 4)
       bbox[1] <- min(boxes[,1], na.rm=TRUE)
       bbox[2] <- min(boxes[,2], na.rm=TRUE)
       bbox[3] <- max(boxes[,3], na.rm=TRUE)
       bbox[4] <- max(boxes[,4], na.rm=TRUE)

       # Save as an SVG document.
       res <- .Call("nro_figure",
                     as.character(fname),
                     as.character(codes),
                     as.numeric(bbox),
                     as.character(jscode),
                     PACKAGE="Numero")
       if(is.character(res)) stop(res)
    }

    # Return results.
    res$file <- fname
    return(res)
}

#----------------------------------------------------------------------------

nroPlotSave.figure <- function(codes, boxes) {
    output <- c()
    for(k in names(codes)) {
        bbox <- boxes[k,]
        dx <- -1*round(bbox[1])
        dy <- -1*round(bbox[2])
        w <- round(bbox[3] - bbox[1])
        h <- round(bbox[4] - bbox[2])

        s <- sprintf("\n<svg id=\"%s\"", k)
	if(length(output) < 1)
            s <- c(s, sprintf("onload=\"initPage('%s', true)\"", k))
        s <- c(s, "draggable=\"false\"")
        s <- c(s, "xmlns=\"http://www.w3.org/2000/svg\"")
        s <- c(s, "style=\"user-select: none;\"")
        s <- c(s, sprintf("x=\"0\" width=\"%d\"", w))
        s <- c(s, sprintf("y=\"0\" height=\"%d\">", h))

        s <- c(s, "\n<polygon points=\"")
        s <- c(s, "0,0")
        s <- c(s, sprintf("%d,0", w))
        s <- c(s, sprintf("%d,%d", w,h))
        s <- c(s, sprintf("0,%d\"", h))

        s <- c(s, "style=\"")
        s <- c(s, "fill: #ffffff;")
        s <- c(s, "pointer-events: none;\"")
        s <- c(s, sprintf("id=\"%s_background\"/>\n", k))

        tf <- sprintf("<g transform=\"translate(%d,%d)\"", dx, dy)
        s <- c(s, tf, sprintf("tfx=\"%d\" tfy=\"%d\"", dx, dy));
        s <- c(s, sprintf("id=\"%s_contents\">", k))
        s <- paste0(s, collapse="\n")

        output <- c(output, s, codes[[k]], "</svg>\n")
    }
    return(output)
}

#----------------------------------------------------------------------------

nroPlotSave.json <- function(data, keys, jsname) {
    words <- character()
    if(is.null(keys)) keys <- colnames(data)
    for(j in 1:ncol(data)) {
        x <- data[,j]
	if(is.logical(x)) x <- as.integer(x)
	if(is.character(x))
           s <- paste0("\"", x, "\"", sep="", collapse=",")
        else
           s <- paste0(x, sep="", collapse=",")
        words[j] <- paste0("\"", keys[j], "\":[", s, "]", sep="")
    }
    js <- paste0(words, sep="", collapse=",\n")
    js <- sprintf("const %s = {\n%s\n};\n", jsname, js)
    return(js);
}

#----------------------------------------------------------------------------

nroPlotSave.colormap <- function(hlights=NULL) {

    # Set color model.
    red <- c(255, 255, 200, 185, 070, 010, 033)/255.0;
    gre <- c(071, 065, 140, 215, 200, 185, 130)/255.0;
    blu <- c(189, 050, 045, 000, 025, 213, 255)/255.0;

    # Set interpolation sequence for 26 letters.
    samples <- seq(from=0, to=1, length.out=32)
    samples <- matrix(samples, nrow=4, ncol=8)
    samples <- as.vector(t(samples))
    samples <- samples[1:26]/max(samples[1:26])

    # Interpolate for 26 letters.
    pivots <- seq(from=0, to=1, length.out=length(red))
    red <- stats::approx(x=pivots, y=red, xout=samples)$y
    gre <- stats::approx(x=pivots, y=gre, xout=samples)$y
    blu <- stats::approx(x=pivots, y=blu, xout=samples)$y

    # Convert to hexa format.
    colrs <- grDevices::rgb(red, gre, blu)
    colrs <- paste0(colrs, "A0") # opacity

    # Automatic text.
    labls <- intToUtf8((65 + 0:25), multiple=TRUE)
    texts <- paste("Subgroup", labls)

    # Finish default regions.
    info <- data.frame(REGION=texts, REGION.label=labls,
        REGION.color=colrs, stringsAsFactors=FALSE)
    if(is.null(hlights)) hlights <- info

    # Check highlight labels.
    hlights <- unique(hlights)
    if(anyDuplicated(hlights$REGION.label) > 0) {
        warning("Duplicated region details.")
	rows <- which(!duplicated(hlights$REGION.label))
	hlights <- hlights[rows,]
    }
    if(nrow(hlights) > nrow(info)) {
        warning("Too many region definitions.")
        hlights <- hlights[1:nrow(info),]
    }

    # Remove redundant highlights.
    hlights <- rbind(hlights, info)
    rows <- which(!duplicated(hlights$REGION.label))
    hlights <- hlights[rows,]

    # Remove empty elements.
    rows <- which(hlights$REGION.label != "")
    hlights <- hlights[rows,]

    # Re-check number of regions.
    hlights <- hlights[1:nrow(info),]
    hlights <- hlights[order(hlights$REGION.label),]

    # Convert names to JSON.
    sA <- paste0("\"", hlights$REGION.label, sep="\"")
    sB <- paste0("\"", hlights$REGION, sep="\"")
    js1 <- paste(sA, sB, sep=":", collapse=",")
    js1 <- paste0("const REGIONS = {", js1, "};\n", sep="", collapse="")
    
    # Convert colors to JSON.
    sA <- paste0("\"", hlights$REGION.label, sep="\"")
    sB <- paste0("\"", hlights$REGION.color, sep="\"")
    js2 <- paste(sA, sB, sep=":", collapse=",")
    js2 <- paste0("const REGIONCOLORS = {", js2, "};\n", sep="", collapse="")
    js <- paste0(js1, js2)
    return(js)
}
