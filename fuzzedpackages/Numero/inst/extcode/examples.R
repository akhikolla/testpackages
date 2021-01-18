library("Numero")
options(digits=3)
if(readline("Clear workspace? [y/N] ") != "y")
  stop("Examples not run.")

cat("\nnroAggregate.Rd\n")
rm(list=ls())

    # Import data.
    fname <- system.file("extdata", "finndiane.txt", package = "Numero")
    dataset <- read.delim(file = fname)
    
    # Prepare training data.
    trvars <- c("CHOL", "HDL2C", "TG", "CREAT", "uALB")
    trdata <- scale.default(dataset[,trvars]) 
    
    # K-means clustering.
    km <- nroKmeans(data = trdata)
    
    # Self-organizing map.
    sm <- nroKohonen(seeds = km)
    sm <- nroTrain(map = sm, data = trdata)
    
    # Assign data points into districts.
    matches <- nroMatch(centroids = sm, data = trdata)
    
    # District averages for one variable.
    chol <- nroAggregate(topology = sm, districts = matches,
                         data = dataset$CHOL)
    print(chol)
    
    # District averages for all variables.
    planes <- nroAggregate(topology = sm, districts = matches, data = dataset)
    print(head(planes))


cat("\nnroColorize.Rd\n")
rm(list=ls())

    # Import data.
    fname <- system.file("extdata", "finndiane.txt", package = "Numero")
    dataset <- read.delim(file = fname)
    
    # Prepare training data.
    trvars <- c("CHOL", "HDL2C", "TG", "CREAT", "uALB")
    trdata <- scale.default(dataset[,trvars]) 
    
    # K-means clustering.
    km <- nroKmeans(data = trdata)
    
    # Self-organizing map.
    sm <- nroKohonen(seeds = km)
    sm <- nroTrain(map = sm, data = trdata)
    
    # Assign data points into districts.
    matches <- nroMatch(centroids = sm, data = trdata)
    
    # District averages for all variables.
    planes <- nroAggregate(topology = sm, districts = matches, data = dataset)
    
    # District colors for cholesterol.
    chol <- nroColorize(values = planes[,"CHOL"])
    print(head(chol))
    
    # District colors for all variables.
    colrs <- nroColorize(values = planes)
    print(head(colrs))


cat("\nnroDestratify.Rd\n")
rm(list=ls())

    # Import data.
    fname <- system.file("extdata", "finndiane.txt", package = "Numero")
    dataset <- read.delim(file = fname)
    
    # Remove sex differences for creatinine.
    creat <- nroDestratify(dataset$CREAT, dataset$MALE)
    
    # Compare creatinine distributions.
    men <- which(dataset$MALE == 1)
    women <- which(dataset$MALE == 0)
    print(summary(dataset[men,"CREAT"]))
    print(summary(dataset[women,"CREAT"]))
    print(summary(creat[men]))
    print(summary(creat[women]))
    
    # Remove sex differences (produces warnings for binary traits).
    ds <- nroDestratify(dataset, dataset$MALE)
    
    # Compare HDL2C distributions.
    print(summary(dataset[men,"HDL2C"]))
    print(summary(dataset[women,"HDL2C"]))
    print(summary(ds[men,"HDL2C"]))
    print(summary(ds[women,"HDL2C"]))


cat("\nnroImpute.Rd\n")
rm(list=ls())

    # Import data.
    fname <- system.file("extdata", "finndiane.txt", package = "Numero")
    dataset <- read.delim(file = fname)
    
    # Convert identities to strings (produces a warning later).
    ds <- dataset
    ds[,"INDEX"] <- paste("K", ds[,"INDEX"], sep=".")
    
    # Introduce missing values to cholesterol.
    missing <- seq(from = 1, to = nrow(ds), length.out = 40)
    missing <- unique(round(missing))
    ds[missing,"CHOL"] <- NA
    
    # Impute missing values with and without standardization.
    ds.std <- nroImpute(data = ds, standard = TRUE)
    ds.orig <- nroImpute(data = ds, standard = FALSE)
    
    # Compare against "true" cholesterol values.
    rho.std <- cor(ds.std[missing,"CHOL"], dataset[missing,"CHOL"])
    rho.orig <- cor(ds.orig[missing,"CHOL"], dataset[missing,"CHOL"])
    cat("Correlation, standard = TRUE:  ", rho.std, "\n", sep="")
    cat("Correlation, standard = FALSE: ", rho.orig, "\n", sep="")


cat("\nnroKmeans.Rd\n")
rm(list=ls())

    # Import data.
    fname <- system.file("extdata", "finndiane.txt", package = "Numero")
    dataset <- read.delim(file = fname)
    
    # Prepare training data.
    trvars <- c("CHOL", "HDL2C", "TG", "CREAT", "uALB")
    trdata <- scale.default(dataset[,trvars]) 
    
    # Unbalanced K-means clustering.
    km0 <- nroKmeans(data = trdata, k = 5, balance = 0.0)
    print(table(km0$layout$BMC))
    print(km0$centroids)
    
    # Balanced K-means clustering.
    km1 <- nroKmeans(data = trdata, k = 5, balance = 1.0)
    print(table(km1$layout$BMC))
    print(km1$centroids)


cat("\nnroKohonen.Rd\n")
rm(list=ls())

    # Import data.
    fname <- system.file("extdata", "finndiane.txt", package = "Numero")
    dataset <- read.delim(file = fname)
    
    # Prepare training data.
    trvars <- c("CHOL", "HDL2C", "TG", "CREAT", "uALB")
    trdata <- scale.default(dataset[,trvars]) 
    
    # K-means clustering.
    km <- nroKmeans(data = trdata)
    
    # Self-organizing map.
    sm <- nroKohonen(seeds = km)
    print(head(sm$centroids))
    print(head(sm$topology))


cat("\nnroLabel.Rd\n")
rm(list=ls())

    # Import data.
    fname <- system.file("extdata", "finndiane.txt", package = "Numero")
    dataset <- read.delim(file = fname)
    
    # Prepare training data.
    trvars <- c("CHOL", "HDL2C", "TG", "CREAT", "uALB")
    trdata <- scale.default(dataset[,trvars]) 
    
    # K-means clustering.
    km <- nroKmeans(data = trdata)
    
    # Self-organizing map.
    sm <- nroKohonen(seeds = km)
    sm <- nroTrain(map = sm, data = trdata)
    
    # Assign data points into districts.
    matches <- nroMatch(centroids = sm, data = trdata)
    
    # District averages for all variables.
    planes <- nroAggregate(topology = sm, districts = matches, data = dataset)
    
    # District labels for cholesterol.
    chol <- nroLabel(topology = sm, values = planes[,"CHOL"])
    print(head(attr(chol, "visible")))
    print(head(chol))
    
    # District labels for all variables.
    colrs <- nroLabel(topology = sm, values = planes)
    print(head(attr(colrs, "visible")))
    print(head(colrs))


cat("\nnroMatch.Rd\n")
rm(list=ls())

    # Import data.
    fname <- system.file("extdata", "finndiane.txt", package = "Numero")
    dataset <- read.delim(file = fname)
    
    # Prepare training data.
    trvars <- c("CHOL", "HDL2C", "TG", "CREAT", "uALB")
    trdata <- scale.default(dataset[,trvars]) 
    
    # K-means clustering.
    km <- nroKmeans(data = trdata, k = 10)
    
    # Assign data points into districts.
    matches <- nroMatch(centroids = km, data = trdata)
    print(head(attr(matches,"quality")))
    print(table(matches))


cat("\nnroPermute.Rd\n")
rm(list=ls())

    # Import data.
    fname <- system.file("extdata", "finndiane.txt", package = "Numero")
    dataset <- read.delim(file = fname)
    
    # Set row names.
    rownames(dataset) <- paste("r", 1:nrow(dataset), sep="")
    
    # Prepare training data.
    trvars <- c("CHOL", "HDL2C", "TG", "CREAT", "uALB")
    trdata <- scale.default(dataset[,trvars])
    
    # K-means clustering.
    km <- nroKmeans(data = trdata)
    
    # Self-organizing map.
    sm <- nroKohonen(seeds = km)
    sm <- nroTrain(map = sm, data = trdata)
    
    # Assign data points into districts.
    matches <- nroMatch(centroids = sm, data = trdata)
    
    # Estimate statistics for cholesterol
    chol <- nroPermute(map = sm, districts = matches, data = dataset$CHOL)
    print(chol[,c("TRAINING", "Z", "P.z", "P.freq")])
    
    # Estimate statistics.
    stats <- nroPermute(map = sm, districts = matches, data = dataset)
    print(stats[,c("TRAINING", "Z", "P.z", "P.freq")])


cat("\nnroPlot.Rd\n")
rm(list=ls())

    # Import data.
    fname <- system.file("extdata", "finndiane.txt", package = "Numero")
    dataset <- read.delim(file = fname)
    
    # Detect binary columns.
    dataset <- nroPreprocess(dataset, method = "")
    
    # Prepare training data.
    trvars <- c("CHOL", "HDL2C", "TG", "CREAT", "uALB")
    trdata <- scale.default(dataset[,trvars])
    
    # K-means clustering.
    km <- nroKmeans(data = trdata)
    
    # Self-organizing map.
    sm <- nroKohonen(seeds = km)
    sm <- nroTrain(map = sm, data = trdata)
    
    # Assign data points into districts.
    matches <- nroMatch(centroids = sm, data = trdata)
    
    # Select a subset of variables and detect binary data.
    vars <- c("AGE", "MALE", "uALB", "CHOL", "DIAB_KIDNEY", "DECEASED")
    selected <- nroPreprocess(dataset[,vars], method = "")
    
    # Calculate district averages for seleted variables.
    vars <- c("AGE", "MALE", "uALB", "CHOL", "DIAB_KIDNEY", "DECEASED")
    planes <- nroAggregate(topology = sm, districts = matches, data = selected)
    
    # Estimate statistics.
    stats <- nroPermute(map = sm, districts = matches, data = selected)
    
    # Set visuals.
    colrs <- nroColorize(values = planes, amplitudes = stats)
    labls <- nroLabel(topology = sm, values = planes)
    
    # Add subgrouping information.
    topo <- sm$topology
    topo$REGION <- ""
    topo$REGION[1:8] <- "Center"
    topo$REGION[9:21] <- "Perimeter"
    
    # Add subgroup labels.
    topo$REGION.label <- ""
    topo$REGION.label[1:8] <- "C"
    topo$REGION.label[9:21] <- "P"
    
    # Add subgroup colors.
    topo$REGION.color <- ""
    topo$REGION.color[1:8] <- "#00f00060"
    topo$REGION.color[9:21] <- "#f000f060"
    
    # Plot colorings on screen.
    nroPlot(topology = topo, colors = colrs, labels = labls)
    
    # Save colorings in file.
    #fn <- "colorings.html"
    #n <- nroPlot.save(file = fn, topology = topo,
    #    colors = colrs, labels = labls)
    #cat(n, " bytes saved in '", fn, "'\n", sep="")


cat("\nnroPostprocess.Rd\n")
rm(list=ls())

    # Import data.
    fname <- system.file("extdata", "finndiane.txt", package = "Numero")
    dataset <- read.delim(file = fname)
    
    # Show original data characteristics.
    print(summary(dataset))
    
    # Preprocess a subset of data.
    ds.pre <- nroPreprocess(dataset[1:100,])
    print(summary(ds.pre))
    
    # Repeat preprocessing for the whole dataset (approximation).
    ds.post <- nroPostprocess(dataset, ds.pre)
    print(summary(ds.post))


cat("\nnroPreprocess.Rd\n")
rm(list=ls())

    # Import data.
    fname <- system.file("extdata", "finndiane.txt", package = "Numero")
    dataset <- read.delim(file = fname)
    
    # Show original data characteristics.
    print(summary(dataset))
    
    # Detect binary columns.
    ds <- nroPreprocess(dataset, method = "")
    print(attr(ds,"binary"))
    
    # Centering and scaling cholesterol.
    ds <- nroPreprocess(dataset$CHOL)
    print(summary(ds))
    
    # Centering and scaling.
    ds <- nroPreprocess(dataset)
    print(summary(ds))
    
    # Tapered ranks.
    ds <- nroPreprocess(dataset, method = "tapered")
    print(summary(ds))


cat("\nnroRcppMatrix.Rd\n")
rm(list=ls())

    # Fully numeric data frame.
    x <- data.frame(A=c(1,2,3,4), B=c(0,1,0,NA), C=c(2,3,4,5))
    print(nroRcppMatrix(data=x, trim=TRUE))
    
    # Matrix of characters, some of which can be converted to numbers.
    x <- matrix(c("1","2","b","4","","6","7","8"), nrow=4, ncol=2)
    print(nroRcppMatrix(data=x, trim=TRUE))
    
    # Object that can be converted to numbers.
    x <- list(text="abc", value="123")
    print(nroRcppMatrix(data=x, trim=TRUE))
    
    # Unusable object.
    x <- list(text="abc", value="123", multiple=c("a","b","c"))
    print(nroRcppMatrix(data=x, trim=TRUE))


cat("\nnroRcppVector.Rd\n")
rm(list=ls())

    # Empty input reverts to default.
    x <- nroRcppVector(data=NULL, default=NA, empty=TRUE)
    print(x)
    
    # Empty input reverts to default, then to specified type.
    x <- nroRcppVector(data=NULL, default=123, empty=TRUE, numeric=FALSE)
    print(x)
    
    # Convert a logical vector to numbers.
    x <- c(TRUE, TRUE, FALSE, TRUE)
    names(x) <- c("a","b","c","d")
    y <- nroRcppVector(data=x, numeric=TRUE)
    print(y)


cat("\nnroSummary.Rd\n")
rm(list=ls())

    # Import data.
    fname <- system.file("extdata", "finndiane.txt", package = "Numero")
    dataset <- read.delim(file = fname)
    
    # Prepare training data.
    trvars <- c("CHOL", "HDL2C", "TG", "CREAT", "uALB")
    trdata <- scale.default(dataset[,trvars])
    
    # K-means clustering.
    km <- nroKmeans(data = trdata)
    
    # Self-organizing map.
    sm <- nroKohonen(seeds = km)
    sm <- nroTrain(map = sm, data = trdata)
    
    # Assign data points into districts.
    matches <- nroMatch(centroids = sm, data = trdata)
    
    # Calculate district averages for urinary albumin.
    plane <- nroAggregate(topology = sm, districts = matches,
                          data = dataset$uALB)
    plane <- as.vector(plane)
    
    # Assign subgroups based on urinary albumin.
    regns <- rep("HighAlb", length.out=length(plane))
    regns[which(plane < quantile(plane, 0.67))] <- "MiddleAlb"
    regns[which(plane < quantile(plane, 0.33))] <- "LowAlb"
    
    # Add label info and make a data frame.
    regns <- data.frame(REGION=regns, REGION.label="",
        stringsAsFactors=FALSE)
    regns[which(regns$REGION == "HighAlb"),"REGION.label"] <- "H"
    regns[which(regns$REGION == "MiddleAlb"),"REGION.label"] <- "M"
    regns[which(regns$REGION == "LowAlb"),"REGION.label"] <- "L"
    
    # Calculate summary statistics.
    st <- nroSummary(data = dataset, districts = matches, regions = regns)
    print(st[,c("VARIABLE","SUBGROUP","MEAN","P.chisq","P.t","P.anova")])


cat("\nnroTrain.Rd\n")
rm(list=ls())

    # Import data.
    fname <- system.file("extdata", "finndiane.txt", package = "Numero")
    dataset <- read.delim(file = fname)
    
    # Prepare training data.
    trvars <- c("CHOL", "HDL2C", "TG", "CREAT", "uALB")
    trdata <- scale.default(dataset[,trvars]) 
    
    # K-means clustering.
    km <- nroKmeans(data = trdata)
    
    # Train with full data.
    sm <- nroKohonen(seeds = km)
    sm <- nroTrain(map = sm, data = trdata, subsample = nrow(trdata))
    print(sm$history)
    
    # Train with subsampling.
    sm <- nroKohonen(seeds = km)
    sm <- nroTrain(map = sm, data = trdata, subsample = 200)
    print(sm$history)


cat("\nnumero.clean.Rd\n")
rm(list=ls())

    # Import data.
    fname <- system.file("extdata", "finndiane.txt", package = "Numero")
    dataset <- read.delim(file = fname)
    
    # Create new versions for testing.
    dsA <- dataset[1:250, c("INDEX","AGE","MALE","uALB")]
    dsB <- dataset[151:300, c("INDEX","AGE","MALE","uALB","CHOL")]
    dsC <- dataset[201:500, c("INDEX","AGE","MALE","DIAB_RETINO")]
    
    # Select all rows.
    results <- numero.clean(a = dsA, b = dsB, c = dsC, identity = "INDEX")
    cat("\n\nNo selection:\n")
    print(nrow(results$a))
    print(nrow(results$b))
    print(nrow(results$c))
    
    # Select all rows and expanded for all identities.
    results <- numero.clean(a = dsA, b = dsB, c = dsC, identity = "INDEX",
                            select = "union")
    cat("\n\nUnion:\n")
    print(nrow(results$a))
    print(nrow(results$b))
    print(nrow(results$c))
    
    # Select only rows that are shared between all datasets.
    results <- numero.clean(a = dsA, b = dsB, c = dsC, identity = "INDEX",
                            select = "intersection")
    cat("\n\nIntersection:\n")
    print(nrow(results$a))
    print(nrow(results$b))
    print(nrow(results$c))
    
    # Select only rows with a unique INDEX ('dsB' has none).
    results <- numero.clean(a = dsA, b = dsB, c = dsC, identity = "INDEX",
                            select = "exclusion")
    cat("\n\nExclusion:\n")
    print(nrow(results$a))
    print(nrow(results$b))
    print(nrow(results$c))
    
    # Add extra identification information.
    dsA$GROUP <- "A"
    dsB$GROUP <- "B"
    dsC$GROUP <- "C"
    
    # Select rows with a unique identifier.
    results <- numero.clean(a = dsA, b = dsB, c = dsC,
                            identity = c("GROUP","INDEX"),
                            select = "exclusion")
    cat("\n\nMulti-identities:\n")
    print(nrow(results$a))
    print(nrow(results$b))
    print(nrow(results$c))


cat("\nnumero.create.Rd\n")
rm(list=ls())

    # Import data.
    fname <- system.file("extdata", "finndiane.txt", package = "Numero")
    dataset <- read.delim(file = fname)
    
    # Set identities and manage missing data.
    dataset <- numero.clean(dataset, identity = "INDEX")
    
    # Prepare training set.
    trvars <- c("CHOL", "HDL2C", "TG", "CREAT", "uALB")
    trdata <- numero.prepare(data = dataset, variables = trvars)
    
    # Create a self-organizing map.
    modl <- numero.create(data = trdata)


cat("\nnumero.evaluate.Rd\n")
rm(list=ls())

    # Import data.
    fname <- system.file("extdata", "finndiane.txt", package = "Numero")
    dataset <- read.delim(file = fname)
    
    # Set identities and manage missing data.
    dataset <- numero.clean(dataset, identity = "INDEX")
    
    # Prepare training variables.
    trvars <- c("CHOL", "HDL2C", "TG", "CREAT", "uALB")
    trdata <- numero.prepare(data = dataset, variables = trvars)
    
    # Create a self-organizing map.
    modl <- numero.create(data = trdata)
    
    # Evaluate map statistics.
    results <- numero.evaluate(model = modl, data = dataset)
    print(results$statistics[,c("TRAINING", "Z", "P.z", "P.freq")])


cat("\nnumero.plot.Rd\n")
rm(list=ls())

    # Import data.
    fname <- system.file("extdata", "finndiane.txt", package = "Numero")
    dataset <- read.delim(file = fname)
    
    # Set identities and manage missing data.
    dataset <- numero.clean(dataset, identity = "INDEX")
    
    # Prepare training variables.
    trvars <- c("CHOL", "HDL2C", "TG", "CREAT", "uALB")
    trdata <- numero.prepare(data = dataset, variables = trvars,
        batch = "MALE", confounders = c("AGE", "T1D_DURAT"))
    
    # Create a self-organizing map.
    modl <- numero.create(data = trdata)
    
    # Evaluate map statistics for all variables.
    stats <- numero.evaluate(model = modl, data = dataset)
    
    # Plot map colorings.
    numero.plot(results = stats)


cat("\nnumero.prepare.Rd\n")
rm(list=ls())

    # Import data.
    fname <- system.file("extdata", "finndiane.txt", package = "Numero")
    dataset <- read.delim(file = fname)
    
    # Set identities and manage missing data.
    dataset <- numero.clean(dataset, identity = "INDEX")
    
    # Prepare training variables using default standardization.
    trvars <- c("CHOL", "HDL2C", "TG", "CREAT", "uALB")
    trdata <- numero.prepare(data = dataset, variables = trvars)
    print(summary(trdata))
    
    # Prepare training values adjusted for age and sex and
    # standardized by rank-based method.
    trdata <- numero.prepare(data = dataset, variables = trvars,
                             batch = "MALE", confounders = "AGE",
    			 method = "tapered")
    print(summary(trdata))


cat("\nnumero.quality.Rd\n")
rm(list=ls())

    # Import data.
    fname <- system.file("extdata", "finndiane.txt", package = "Numero")
    dataset <- read.delim(file = fname)
    
    # Set identities and manage missing data.
    dataset <- numero.clean(dataset, identity = "INDEX")
    
    # Prepare training variables.
    trvars <- c("CHOL", "HDL2C", "TG", "CREAT", "uALB")
    trdata <- numero.prepare(data = dataset, variables = trvars)
    
    # Create a self-organizing map.
    modl <- numero.create(data = trdata)
    
    # Analyze map quality.
    qc <- numero.quality(model = modl)


cat("\nnumero.subrgoup.Rd\n")
rm(list=ls())

    # Import data.
    fname <- system.file("extdata", "finndiane.txt", package = "Numero")
    dataset <- read.delim(file = fname)
    
    # Set identities and manage missing data.
    dataset <- numero.clean(dataset, identity = "INDEX")
    
    # Prepare training variables.
    trvars <- c("CHOL", "HDL2C", "TG", "CREAT", "uALB")
    trdata <- numero.prepare(data = dataset, variables = trvars)
    
    # Create a self-organizing map.
    modl <- numero.create(data = trdata)
    
    # Evaluate map statistics for all variables.
    stats <- numero.evaluate(model = modl, data = dataset)
    
    # Define subgroups, uncomment to launch interactive window.
    #elem <- numero.subgroup(results = stats, variables = trvars)


cat("\nnumero.summary.Rd\n")
rm(list=ls())

    # Import data.
    fname <- system.file("extdata", "finndiane.txt", package = "Numero")
    dataset <- read.delim(file = fname)
    
    # Set identities and manage missing data.
    dataset <- numero.clean(dataset, identity = "INDEX")
    
    # Prepare training variables.
    trvars <- c("CHOL", "HDL2C", "TG", "CREAT", "uALB")
    trdata <- numero.prepare(data = dataset, variables = trvars)
    
    # Create a self-organizing map.
    modl <- numero.create(data = trdata)
    
    # Evaluate map statistics for all variables.
    stats <- numero.evaluate(model = modl, data = dataset)
    
    # Define subgroups.
    x <- stats$planes[,"uALB"]
    tops <- which(x >= quantile(x, 0.75, na.rm=TRUE))
    bottoms <- which(x <= quantile(x, 0.25, na.rm=TRUE))
    elem <- data.frame(stats$map$topology, stringsAsFactors = FALSE)
    elem$REGION <- "MiddleAlb"
    elem$REGION[tops] <- "HighAlb"
    elem$REGION[bottoms] <- "LowAlb"
    elem$REGION.label <- "M"
    elem$REGION.label[tops] <- "H"
    elem$REGION.label[bottoms] <- "L"
    
    # Compare subgroups.
    cmp <- numero.summary(results = stats, topology = elem, data = dataset)


cat("\nAll examples completed.\n\n")
