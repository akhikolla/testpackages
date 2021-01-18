nroSummary <- function(
    data,
    districts,
    regions=NULL,
    categlim=8,
    capacity=10) {

    # Convert vector to data matrix.
    if(is.vector(data)) {
        keys <- names(data)
        data <- data.frame(X=data, stringsAsFactors=FALSE)
        if(!is.null(keys)) rownames(data) <- keys
    }
    if(nrow(data)*ncol(data) < 1) {
        warning("Empty input.")
        return(NULL)
    }

    # Check input sizes.
    if(length(districts) != nrow(data))
        stop("Incompatible inputs.")

    # Check if districts data are named.
    if(is.null(names(districts)))
        names(districts) <- rownames(data)

    # Check data consistency.
    keys <- intersect(names(districts), rownames(data))
    if(length(keys) < 1) stop("Incompatible inputs.")
    districts <- districts[keys]
    data <- data[keys,]

    # Check if regions are available.
    if(is.null(regions)) regions <- sort(unique(districts))
    if(is.vector(regions)) {
        labls <- as.integer(as.factor(regions))
        regions <- data.frame(REGION=regions, REGION.label=labls,
            stringsAsFactors=FALSE)
	rownames(regions) <- regions$REGION
    }

    # Check other inputs.
    categlim <- as.integer(categlim[[1]])
    capacity <- as.integer(capacity[[1]])

    # Replace districts with regions.
    pos <- match(districts, rownames(regions))
    mask <- which(pos > 0)
    if(length(mask) < length(pos))
        warning("Unusable districts or regions excluded.")
    g <- regions[pos[mask],"REGION"]
    g.label <- regions[pos[mask],"REGION.label"]
    data <- data[mask,]

    # Check subgroups.
    t <- table(g)
    if(length(t) < 2) {
        warning("Less than two subgroups.")
        return(NULL)
    }
    if(length(t) > capacity) {
        warning("Subgroup capacity exceeded.")
        return(NULL)
    }

    # Process data columns.
    output <- NULL
    for(vn in colnames(data)) {
        x <- data[,vn]
        nlev <- nlevels(as.factor(x))

        # Unusable values.
        if(nlev < 2) {
	    warning(paste(vn, ": Unusable values.", sep=""))
            next
        }

        # Categorical data.
        if((nlev <= categlim) || (is.numeric(x) == FALSE)) {
	    stats <- nroSummary.categ(x, g)
	    if(is.character(stats)) {
	        warning(paste(vn, ": ", stats, sep=""))
            }
	    else {
	        stats <- data.frame(VARIABLE=vn, stats,
		    stringsAsFactors=FALSE)
	        output <- rbind(output, stats)
	    }
	    next
        }

        # Numerical data.
	stats <- nroSummary.real(x, g)
	if(is.character(stats)) {
	    warning(paste(vn, ": ", stats, sep=""))
        }
	else {
	    stats <- data.frame(VARIABLE=vn, stats,
	        stringsAsFactors=FALSE)
	    output <- rbind(output, stats)
	}
    }

    # Add region labels.
    pos <- match(output$SUBGROUP, regions$REGION)
    rows <- which(pos > 0)
    output$LABEL[rows] <- regions[pos[rows],"REGION.label"]

    # Finish results.
    attr(output, "regions") <- g
    attr(output, "labels") <- g.label
    attr(output, "subgroups") <- split(1:nrow(data), g)
    rownames(output) <- NULL
    return(output)
}

#---------------------------------------------------------------------------

nroSummary.categ <- function(x, g) {

    # Check if enough data.
    mask <- which((is.na(x) == FALSE) & (is.na(g) == FALSE))
    if(length(mask) < 10) return("Too few usable data.")
    x <- x[mask]
    g <- g[mask]

    # Split into subgroups.
    xsets <- split(x, g)
    if(length(xsets) < 2) return("Only one subgroup.")

    # Recode if binary.
    xfact <- as.factor(x)
    xlevs <- levels(xfact)
    nlevs <- length(xlevs)
    if(nlevs == 2) {
        xsets <- lapply(xsets, function(v, l) {
	    v <- (v == l)
        }, xlevs[2])
	xlevs <- c(TRUE, FALSE)
	xfact <- (as.integer(xfact) - 1)
    }

    # Subgroup sizes.
    stats <- list()
    stats$N <- lapply(xsets, function(v) {
        return(sum(is.na(v) == FALSE))
    })

    # Add extra class labels to every subgroup to ensure tests
    # work (this will slightly dilute the results).
    for(j in 1:length(xsets))
        xsets[[j]] <- c(xsets[[j]], xlevs)

    # Estimate basic subgroup stats.
    stats$MEAN <- rep(NA, length(xsets))
    stats$SD <- rep(NA, length(xsets))
    stats$MEDIAN <- rep(NA, length(xsets))
    stats$MAD <- rep(NA, length(xsets))
    stats$Q025 <- rep(NA, length(xsets))
    stats$Q250 <- rep(NA, length(xsets))
    stats$Q750 <- rep(NA, length(xsets))
    stats$Q975 <- rep(NA, length(xsets))
    if(nlevs == 2) stats$MEAN <- lapply(xsets, mean, na.rm=TRUE)

    # Convert to data frame.
    stats <- lapply(stats, as.double)
    stats <- as.data.frame(stats, stringsAsFactors=FALSE)
    stats <- data.frame(SUBGROUP=names(xsets), LABEL=NA, stats,
        stringsAsFactors=FALSE)

    # Add P-value columns.
    stats$TYPE <- "categ"
    if(nlevs == 2) stats$TYPE <- "binary"
    stats$P.chisq <- NA
    stats$P.t <- NA
    stats$P.anova <- NA

    # Find the biggest subgroup.
    ind <- which.max(stats$N)
    if(length(ind) < 1) return("Unusable subgroups.")

    # Chi-squared test by subgroup.
    suppressWarnings({
    stats$P.chisq <- lapply(xsets, function(v, v0) {
        bits <- c(rep(0, length(v)), rep(1, length(v0)))
        st <- stats::chisq.test(c(v, v0), bits)
	return(st$p.value)
    }, xsets[[ind]])})
    stats$P.chisq <- as.double(stats$P.chisq)

    # Chi-squared test for whole data.
    suppressWarnings(st <- stats::chisq.test(x, g))

    # Add extra row for full test.
    stats <- rbind(stats[1,], stats)
    stats[1,c("SUBGROUP", "MEAN","SD",
        "MEDIAN","MAD","P.t","P.anova")] <- NA
    if(nlevs == 2) stats$MEAN[1] <- mean(xfact, na.rm=TRUE)
    stats$N[1] <- length(x)
    stats$P.chisq[1] <- st$p.value

    # Return results.
    return(stats)
}

#---------------------------------------------------------------------------

nroSummary.real <- function(x, g) {
    suppressWarnings(x <- as.numeric(x))
  
    # Check if enough data.
    mask <- which(0*x == 0)
    if(length(mask) < 10) return("Too few data.")
    sigma <- stats::sd(x[mask], na.rm=TRUE)
    if(!is.finite(sigma)) return("Unusable data.")
    if(sigma < 1e-9) return("Too low variance.")
    x <- x[mask]
    g <- g[mask]

    # Convert to tapered ranks.
    z <- (rank(x, na.last="keep") - 1)
    z <- (2*z/max(z, na.rm=TRUE) - 1)
    z <- 0.5*(z + z^3)

    # Split into subgroups.
    xsets <- split(x, g)
    zsets <- split(z, g)
    if(length(xsets) < 2) return("Only one subgroup.")

    # Subgroup sizes.
    stats <- list()
    stats$N <- lapply(xsets, function(v) {
        return(sum(0*v == 0))
    })

    # Estimate basic subgroup stats.
    stats$MEAN <- lapply(xsets, mean, na.rm=TRUE)
    stats$SD <- lapply(xsets, stats::sd, na.rm=TRUE)
    stats$MEDIAN <- lapply(xsets, stats::median, na.rm=TRUE)
    stats$MAD <- rep(NA, length(xsets))
    stats$Q025 <- rep(NA, length(xsets))
    stats$Q975 <- rep(NA, length(xsets))
    for(j in 1:length(xsets)) {
        xj <- as.double(xsets[[j]])
	muj <- as.double(stats$MEDIAN[j])
        stats$MAD[j] <- stats::median(abs(xj - muj), na.rm=TRUE)
        qi <- stats::quantile(xj, c(0.025, 0.25, 0.75, 0.975), na.rm=TRUE)
        stats$Q025[j] <- qi[1]
        stats$Q250[j] <- qi[2]
        stats$Q750[j] <- qi[3]
        stats$Q975[j] <- qi[4]
    }

    # Convert to data frame.
    stats <- lapply(stats, as.double)
    stats <- as.data.frame(stats, stringsAsFactors=FALSE)
    stats <- data.frame(SUBGROUP=names(xsets), LABEL=NA, stats,
        stringsAsFactors=FALSE)

    # Add P-value columns.
    stats$TYPE <- "real"
    stats$P.chisq <- NA
    stats$P.t <- NA
    stats$P.anova <- NA

    # Find the "most average" subgroup.
    mu <- mean(x, na.rm=TRUE)
    ind <- which.min(abs(stats$MEAN - mu))
    if(length(ind) < 1) return("Unusable subgroups.")

    # Rank-regulated T-tests against the reference group.
    stats$P.t <- lapply(zsets, function(v, v0) {
        xsigma <- stats::sd(v, na.rm=TRUE)
	ysigma <- stats::sd(v0, na.rm=TRUE)
        if(!is.finite(xsigma)) return(1.0)
	if(!is.finite(ysigma)) return(1.0)
        if(xsigma < 1e-20) return(1.0)
	if(ysigma < 1e-20) return(1.0)
        return(stats::t.test(x=v, y=v0)$p.value)
    }, zsets[[ind]])
    stats$P.t <- as.double(stats$P.t)

    # Analysis of variance with randomized design.
    p.anova <- 1.0
    try({
        tmp <- data.frame(Z=as.double(z), G=as.factor(g))
        fit <- stats::aov(formula=Z~G, data=tmp)
        p.anova <- unlist(summary(fit))["Pr(>F)1"]
        p.anova <- as.double(p.anova)
    }, silent=TRUE)

    # Add extra row for ANOVA.
    stats <- rbind(stats[1,], stats)
    stats[1,c("SUBGROUP", "MEAN","SD",
        "MEDIAN","MAD","P.t","P.chisq")] <- NA
    stats$N[1] <- length(x)
    stats$MEAN[1] <- mean(x, na.rm=TRUE)
    stats$SD[1] <- stats::sd(x, na.rm=TRUE)
    stats$MEDIAN[1] <- stats::median(x, na.rm=TRUE)
    stats$MAD[1] <- stats::median(abs(x - stats$MEDIAN[1]), na.rm=TRUE)
    stats$Q025[1] <- stats::quantile(x, 0.025, na.rm=TRUE)
    stats$Q250[1] <- stats::quantile(x, 0.250, na.rm=TRUE)
    stats$Q750[1] <- stats::quantile(x, 0.750, na.rm=TRUE)
    stats$Q975[1] <- stats::quantile(x, 0.975, na.rm=TRUE)
    stats$P.anova[1] <- p.anova

    # Return results.
    return(stats)
}
