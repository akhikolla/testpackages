numero.prepare <- function(
    data,
    variables=NULL,
    confounders=NULL,
    batch=NULL,
    method="standard",
    pipeline=NULL) {

    # Start processing.
    stamp <- date()
    cat("\n*** numero.prepare ***\n", stamp, "\n", sep="")

    # Use all data columns.
    if(length(variables) < 1) variables <- colnames(data)
    if(length(batch) > 1) batch <- batch[[1]]

    # Create a new pipeline.
    if(is.null(pipeline)) {
        cat("\nSetting up:\n")
        pipeline <- list(batch=batch)

        # Make sure variable groups are distinct.
        confounders <- setdiff(confounders, batch)
        variables <- setdiff(variables, c(confounders, batch))

        # Select variables.
        convars <- intersect(confounders, colnames(data))
        batvars <- intersect(batch, colnames(data))
        datvars <- intersect(variables, colnames(data))
        convars <- setdiff(convars, batch)

        # Convert to numeric and trim unusable rows and columns.
        dsT <- nroRcppMatrix(data[,c(convars, datvars)], trim=TRUE)
        convars <- intersect(convars, colnames(dsT))
        datvars <- intersect(datvars, colnames(dsT))

        # Check that enough columns.
        cat(length(convars), " / ", length(confounders),
            " confounder columns\n", sep="")
        cat(length(batvars), " / ", length(batch),
            " batch columns\n", sep="")
        cat(length(datvars), " / ", length(variables),
            " data columns\n", sep="")
        if(length(datvars) < 3) {
            cat("too few data columns\n")
            return(NULL)
        }

        # First round of preprocessing.
        suppressWarnings(dsT <- nroPreprocess(data=dsT,
	    method=method, trim=TRUE))
        pipeline$mapping1 <- attr(dsT, "mapping")

        # Regression model of confounding.
        if(length(convars) > 0) {
            dsT <- numero.prepare.regress(dsT, convars)
            pipeline$adjustment <- attr(dsT, "adjustment")
        }

        # Correction model for batch differences.
        if(length(batvars) > 0) {
            dsT <- numero.prepare.flatten(dsT, data, batvars)
            pipeline$correction <- attr(dsT, "correction")
        }

        # Second round of preprocessing.
        if((length(convars) + length(batvars)) > 0) {
            suppressWarnings(dsT <- nroPreprocess(data=dsT,
	        method=method, trim=TRUE))
            pipeline$mapping2 <- attr(dsT, "mapping")
        }
    }

    # Check if pipeline is in the attribute.
    if(is.data.frame(pipeline) || is.matrix(pipeline)) 
        pipeline <- attr(pipeline, "pipeline")

    # First round of preprocessing, keep empty rows.
    cat("\nProcessing:\n"); ds <- data
    if(!is.null(pipeline$mapping1)) {
        suppressWarnings(ds <- nroPostprocess(data=ds,
            mapping=pipeline$mapping1, trim=FALSE))
        cat(ncol(ds), " columns standardized\n", sep="")
    }

    # Adjust for confounding.
    if(!is.null(pipeline$adjustment))
        ds <- numero.prepare.adjust(ds, pipeline$adjustment)

    # Adjust for batches.
    subsets <- NULL
    if(!is.null(pipeline$correction)) {
        ds <- numero.prepare.correct(ds, data, pipeline$correction)
        subsets <- attr(ds, "subsets")
    }

    # Second round of preprocessing with removal of empty rows.
    if(!is.null(pipeline$mapping2) && (length(ds) > 0)) {
        suppressWarnings(ds <- nroPostprocess(data=ds,
            mapping=pipeline$mapping2, trim=TRUE))
        cat(ncol(ds), " columns re-standardized\n", sep="")
    }

    # Convert to matrix.
    if(is.null(ds)) ds <- matrix(nrow=0, ncol=0)
    else ds <- as.matrix(ds)

    # Final report.
    cat("\nSummary:\n", sep="")
    cat(nrow(ds), " / ", nrow(data), " usable rows\n", sep="")
    cat(ncol(ds), " / ", ncol(data), " usable columns\n", sep="")

    # Return results.
    attr(ds, "pipeline") <- pipeline
    attr(ds, "subsets") <- subsets
    return(ds)
}

#--------------------------------------------------------------------------

numero.prepare.regress <- function(ds, convars) {
    if(length(ds) < 1) return(NULL)
    if(length(convars) < 1) return(ds)

    # Impute missing values.
    suppressWarnings(confs <- nroImpute(data=ds[,convars]))

    # Add intercept to confounder matrix.
    confs <- cbind(rep(1, nrow(ds)), confs)
    confs <- as.matrix(confs) # prevent stats::lm.fit() fail
    colnames(confs) <- c("_", convars)

    # Find target columns and usable rows.
    vars <- setdiff(colnames(ds), convars)
    rows <- which(is.finite(rowMeans(confs)))

    # Prepare coefficient matrix.
    coeff <- matrix(NA, nrow=length(vars), ncol=ncol(confs))
    rownames(coeff) <- vars
    colnames(coeff) <- colnames(confs)

    # Fit models.
    for(vn in vars) {
        y <- ds[,vn]; ds[,vn] <- NA
        mask <- intersect(which(is.finite(y)), rows)
        if(length(mask) < 10) next

        # Regress confounding variance.
        m <- stats::lm.fit(x=confs[mask,], y=y[mask])
        n <- sum(is.finite(m$residuals))
        if(n < 10) next

        # Update results.
        coeff[vn,] <- as.double(m$coefficients)
        ds[mask,vn] <- m$residuals
    }

    # Collect results.
    ds <- ds[,rownames(coeff)]
    attr(ds, "adjustment") <- coeff
    return(ds)
}

#--------------------------------------------------------------------------

numero.prepare.adjust <- function(ds, coeff) {
    if(length(ds) < 1) return(NULL)
    if(is.null(coeff)) return(ds)

    # Prepare confounder matrix.
    convars <- setdiff(colnames(coeff), "_")
    convars <- intersect(convars, colnames(ds))
    confs <- cbind(rep(1, nrow(ds)), ds[,convars])
    confs <- as.matrix(confs)
    colnames(confs) <- c("_", convars)

    # Check that confounders match parameters.
    if(ncol(coeff) != ncol(confs)) {
        cat("incompatible confounders\n")
        return(NULL)
    }
    if(sum(colnames(coeff) != colnames(confs)) > 0) {
        cat("incompatible confounders\n")
        return(NULL)
    }

    # Apply confounder adjustments.
    failed <- character()
    vars <- intersect(colnames(ds), rownames(coeff))
    for(vn in vars) {
        x <- as.double(confs %*% coeff[vn,])
        ds[,vn] <- (ds[,vn] - x)
        if(sum(is.finite(ds[,vn])) > 0) next
	failed <- c(failed, vn)
    }

    # Number of failed adjustments.
    nfail <- length(failed)
    if(nfail > 0) cat(nfail, " adjustments failed\n", sep="")
    if(nfail == length(vars)) return(NULL)

    # Show report.
    cat(length(vars), " columns adjusted\n", sep="")

    # Return results.
    return(ds[,vars,drop=FALSE])
}

#--------------------------------------------------------------------------

numero.prepare.flatten <- function(ds, ds.orig, batch) {
    if(length(ds) < 1) return(NULL)
    if(length(batch) < 1) return(ds)

    # Check batch info.
    batch <- intersect(colnames(ds.orig), batch[[1]])
    if(length(batch) < 1) {
        cat("batch column missing\n")
        return(NULL)
    }

    # Assign batch labels.
    labels <- rep(NA, nrow(ds))
    pos <- match(rownames(ds), rownames(ds.orig))
    rows <- which(pos > 0)
    labels[rows] <- ds.orig[pos[rows],batch]

    # De-stratify values.
    suppressWarnings(ds.new <- nroDestratify(ds, labels))

    # Exclude incomplete corrections.
    incomp <- attr(ds.new, "incomplete")
    vars <- setdiff(colnames(ds.new), incomp)
    vars <- intersect(vars, colnames(ds))
    ds <- ds[,vars]; ds.new <- ds.new[,vars]
    cat(length(incomp), " incomplete corrections\n", sep="")

    # Check if anything to do.
    if(length(vars) < 2) {
        cat("batch correction failed\n")
        return(NULL)
    }

    # Create mappings by subgroup.
    mappings <- list()
    q <- seq(from=0, to=1, length.out=100)
    subsets <- split(1:length(labels), labels)
    for(sn in names(subsets)) {
        rows <- subsets[[sn]]
        model.in <- ds[rows,]
        model.out <- ds.new[rows,]
        model.in <- apply(model.in, 2, stats::quantile, probs=q, na.rm=T)
        model.out <- apply(model.out, 2, stats::quantile, probs=q, na.rm=T)
        mappings[[sn]] <- list(input=model.in, output=model.out)
    }

    # Return results.
    results <- list()
    results$mappings <- mappings
    results$variables <- colnames(ds.new)
    results$batch <- batch
    attr(ds.new, "correction") <- results
    return(ds.new)
}

#--------------------------------------------------------------------------

numero.prepare.correct <- function(ds, ds.orig, param) {
    if(length(ds) < 1) return(NULL)

    if(nrow(ds) != nrow(ds.orig)) stop("Incompatible inputs.")
    if(is.null(param)) return(ds)

    # Check batch info.
    batch <- intersect(colnames(ds.orig), param$batch)
    if(length(batch) < 1) {
        cat("batch column missing\n")
        return(NULL)
    }

    # Separate batches.
    labels <- ds.orig[,param$batch]
    subsets <- split(1:length(labels), labels)

    # Check compatibility of mappings.
    mappings <- param$mappings
    pos <- match(names(subsets), names(mappings))
    if(sum(is.na(pos)) > 0) {
        cat("incompatible batch models\n")
        return(NULL)
    }

    # Apply correction models.
    vars <- intersect(colnames(ds), param$variables)
    for(sn in names(subsets)) {
        rows <- subsets[[sn]]
        map <- mappings[[sn]]
        suppressWarnings(
        ds[rows,vars] <- nroPostprocess(ds[rows,vars], map))
    }

    # Show report.
    cat(length(vars), " columns corrected\n", sep="")

    # Convert row indices to names.
    rnames <- rownames(ds)
    for(sn in names(subsets))
        subsets[[sn]] <- rnames[subsets[[sn]]]

    # Finish results.
    ds <- ds[,vars]
    attr(ds, "subsets") <- subsets
    return(ds)
}
