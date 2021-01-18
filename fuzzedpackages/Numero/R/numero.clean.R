numero.clean <- function(
    ...,
    identity=NULL,
    na.freq=0.9,
    num.only=TRUE,
    select="") {

    # Start processing.
    stamp <- date()
    cat("\n*** numero.clean ***\n", stamp, "\n", sep="")

    # Check if anything to do.
    datasets <- list(...)
    if(length(datasets) < 1) {
        warning("No data.")
        return(NULL)
    }

    # If a single dataset, name is not used.
    if(length(datasets) == 1) names(datasets) <- "dataset"

    # Require names to be set for all inputs.
    if(sum(names(datasets) != "") < length(datasets)) {
        cat("\nPlease name all datasets in the function call.\n")
        cat("For example, to pass datasets 'x' and 'y', ")
        cat("you can use 'x=x' and 'y=y'.\n\n")
        stop("Unnamed arguments detected.")
    }

    # Check that names do not clash with inputs.
    reserv <- c("identity", "na.freq", "select")
    if(anyDuplicated(c(names(datasets), reserv)) > 0)
        stop("Dataset name overlaps with input argument.")

    # Check NA frequency.
    na.freq <- as.double(na.freq[[1]])
    if(!is.finite(na.freq)) stop("Unusable NA frequency threshold.")
    if(na.freq < 0.0) stop("Unusable NA frequency threshold.")
    if(na.freq >= 1.0) stop("Unusable NA frequency threshold.")

    # Check numeric flag.
    if(!is.logical(num.only)) stop("Unusable 'num.only'.")
    num.only <- (num.only[[1]] == TRUE)

    # Check selection mode.
    select <- as.character(select[[1]])
    if(anyDuplicated(c(select, "", "union",
        "intersection", "exclusion")) < 1)
        stop("Unknown selection mode.")

    # Set row names.
    rnames <- character()
    for(dn in names(datasets)) {
        if(length(datasets) == 1) cat("\nIdentities:\n", sep="")
	else cat("\nIdentities in '", dn, "':\n", sep="")
	datasets[[dn]] <- numero.clean.names(datasets[[dn]], identity)
        rnames <- c(rnames, rownames(datasets[[dn]]))
    }

    # Select preliminary rownames.
    if(length(rnames) > 0) {
        h <- table(rnames)
        rnames <- names(h)
        if(select == "exclusion") rnames <- rnames[which(h == min(h))]
        if(select == "intersection") rnames <- rnames[which(h == max(h))]
    }

    # Process data.
    processed <- list()
    for(k in 1:100) {
        shared <- rnames
        for(dn in names(datasets)) {
            if(length(datasets) > 1)
	        cat("\nValues in '", dn, "', scan ", k, ":\n", sep="")
	    else
	        cat("\nValues, scan ", k, ":\n", sep="")

            # Clean unusable entries.
	    ds <- datasets[[dn]]
            ds <- numero.clean.filter(ds, rnames, na.freq, num.only)
            if(!is.null(ds)) shared <- intersect(shared, rownames(ds))

            # Add data frame to processed.
	    processed[[dn]] <- data.frame(ds, stringsAsFactors=FALSE)
        }

        # One scan is enough if identities not shared. 
	if(select != "intersection") break

        # Enforce shared identities.
        for(dn in names(processed)) {
	    ds <- processed[[dn]]
	    if(is.null(ds)) next
	    processed[[dn]] <- ds[shared,]
        }

        # Check if another scan is needed.
	if(length(shared) == length(rnames)) break
        rnames <- shared
    }

    # Collect names of usable rows and columns.
    output <- list()
    nrows <- 0; ncols <- 0
    rnames <- c(); cnames <- c()
    for(dn in names(processed)) {
        ds <- processed[[dn]]
	if(is.null(ds)) next
	nrows <- (nrows + nrow(ds))
	ncols <- (ncols + ncol(ds))
        rnames <- c(rnames, rownames(ds))
        cnames <- c(cnames, colnames(ds))
        output[[dn]] <- ds
    }

    # Find distinct names.
    rnames <- unique(rnames)
    cnames <- unique(cnames)

    # Expand all datasets to the same rows.
    if(select == "union") {
	nkeys <- length(rnames)
        for(dn in names(processed)) {
            ds <- processed[[dn]]
	    if(is.null(ds)) next

            # Copy data.
            dsnew <- list()
	    pos <- match(rownames(ds), rnames)
            for(vn in colnames(ds)) {
	        x <- rep(NA, nkeys)
                x[pos] <- ds[,vn]
		dsnew[[vn]] <- x
            }

            # Convert to data frame or matrix.
            dsnew <- as.data.frame(dsnew, stringsAsFactors=FALSE)
            if(is.matrix(ds)) dsnew <- as.matrix(dsnew)

            # Set row names.
            rownames(dsnew) <- rnames
            output[[dn]] <- dsnew
	}
    }

    # Final report.
    cat("\nSummary:\n", sep="")
    if(length(processed) > 1) {
        cat(length(output), " usable datasets\n", sep="")
        cat(length(rnames), " unique row names\n", sep="")
        cat(length(cnames), " unique column names\n", sep="")
    }
    cat(nrows, " usable rows\n", sep="")
    cat(ncols, " usable columns\n", sep="")

    # Return results.
    if(length(processed) > 1) return(output)
    if(length(output) == 1) return(output[[1]])
    return(NULL)
}

#--------------------------------------------------------------------------

numero.clean.names <- function(ds, keys) {
    if(length(ds) < 1) {
        cat("empty dataset\n")
        return(NULL)
    }
    if(is.vector(ds)) {
        cat("vector input excluded\n")
        return(NULL)
    }
    
    # Set row names according to identity columns.
    if(length(keys) > 0) {

        # Check that keys are available.
        dskeys <- intersect(keys, colnames(ds))
        if(length(dskeys) < length(keys)) {
            cat("one or more identity columns not found\n")
            return(NULL)
        }

        # Construct identities.
        rnames <- ds[,keys[[1]]]
        if(length(keys) > 1) {
            for(k in 2:length(keys)) {
	        rk <- ds[,keys[[k]]]
                rnames <- paste(rnames, rk, sep="_")
            }
        }

        # Remove missing values.
	mask <- which(!is.na(rnames))
	nmiss <- (length(rnames) - length(mask))
        if(nmiss > 0) {
	    cat(nmiss, " unusable identities removed\n", sep="")
	    rnames <- rnames[mask]
            ds <- ds[mask,]
        }

        # Remove duplicated identities.
	mask <- which(!duplicated(rnames))
	ndupl <- (length(rnames) - length(mask))
        if(ndupl > 0) {
	    cat(ndupl, " duplicates removed\n", sep="")
	    rnames <- rnames[mask]
            ds <- ds[mask,]
        }

        # Remove identity columns from dataset.
	rownames(ds) <- rnames
	ds <- ds[,setdiff(colnames(ds), keys)]
	cat(length(keys), " identity columns consumed\n", sep="")
    }

    # Set default row names.
    if(length(rownames(ds)) != nrow(ds)) {
        rownames(ds) <- (nrows + 1:nrow(ds))
        nrows <- (nrows + nrow(ds))
        cat("row names set to consecutive numbers\n")
    }

    # Set default column names.
    if(length(colnames(ds)) != ncol(ds)) {
        colnames(ds) <- paste("V.", 1:ncol(ds), sep="")
        cat("column names set based on consecutive numbering\n")
    }

    # Check column names.
    vars <- stats::na.omit(colnames(ds))
    vars <- vars[which(vars != "")]
    vars <- vars[which(vars != "_")] # numero.prepare.adjust
    if(length(vars) < ncol(ds)) {
        nskip <- (ncol(ds) - length(vars))
        cat(nskip, " columns with unusable name\n", sep="")
        ds <- ds[,vars]
    }

    # Remove duplicated rows.
    dupl <- which(duplicated(rownames(ds)))
    if(length(dupl) > 0) {
        u <- unique(rownames(ds))
        ds <- ds[u,]
    }
    cat(nrow(ds), " unique row names\n", sep="")

    # Check if any data.
    if(nrow(ds) < 1) {
        cat("no usable rows\n")
        return(NULL)
    }
    return(ds)
}

#--------------------------------------------------------------------------

numero.clean.filter <- function(ds, rnames, na.freq, num.only) {
    if(length(ds) < 1) return(NULL)
    if(is.vector(ds)) {
        cat("vector input excluded\n")
        return(NULL)
    }

    # Select data points.
    rnames <- intersect(rnames, rownames(ds))
    cat(length(rnames), " / ", nrow(ds), " rows selected\n", sep="")
    ds <- ds[rnames,]
    if(nrow(ds) < 2) {
        cat("less than two usable rows")
        return(NULL)
    }

    # Remove empty columns.
    flags <- apply(ds, 2, numero.clean.check, flimit=na.freq)
    if(sum(flags) > 0) {
        ds <- ds[,which(!flags)]
        cat(sum(flags), " unusable columns\n", sep="")
	if(ncol(ds) < 2) {
	    cat("less than two usable columns")
	    return(NULL)
	}
    }

    # Remove empty rows.
    flags <- apply(ds, 1, numero.clean.check, flimit=na.freq)
    if(sum(flags) > 0) {
        ds <- ds[which(!flags),]
        cat(sum(flags), " unusable rowss\n", sep="")
	if(nrow(ds) < 2) {
	    cat("less than two usable rows")
	    return(NULL)
	}
    }

    # Detect numeric and binary columns.
    binary <- character()
    factors <- character()
    numerics <- character()
    for(vn in colnames(ds)) {
      x <- ds[,vn]
      if(is.factor(x)) {
          factors <- c(factors, vn)
          next
      }
      bits.na <- is.na(x)
      suppressWarnings(x <- as.numeric(x))
      bits.fin <- is.finite(x)
      nfinite <- sum(is.finite(x))
      if(sum(bits.fin) < 0.5*sum(!bits.na)) next
      n0 <- sum((x == 0), na.rm=TRUE)
      n1 <- sum((x == 1), na.rm=TRUE)
      if((n0 + n1) == sum(bits.fin)) binary <- c(binary, vn)
      numerics <- c(numerics, vn)
      ds[,vn] <- x # enforce correct class
    }

    # Print report.
    nother <- (ncol(ds) - length(factors) - length(numerics))
    cat(length(binary), " binary columns\n", sep="")
    cat(length(factors), " factor columns\n", sep="")
    cat(length(numerics), " numeric columns\n", sep="")
    cat(nother, " other columns\n", sep="")

    # Select only numeric columns.
    if(num.only) {
        nremov <- (ncol(ds) - length(numerics))
	if(nremov > 0) {
            ds <- ds[,numerics]
            factors <- character()
	    cat(nremov, " non-numeric columns removed\n", sep="")
	}
    }

    # Return results.
    attr(ds, "numerics") <- numerics
    attr(ds, "factors") <- factors
    attr(ds, "binary") <- binary
    return(ds)
}

#---------------------------------------------------------------------------

numero.clean.check <- function(x, flimit) {
    if(is.numeric(x)) {
        f <- sum(!is.finite(x))/length(x)
        return(f > flimit)
    }
    f <- sum(is.na(x))/length(x)
    return(f > flimit)
}