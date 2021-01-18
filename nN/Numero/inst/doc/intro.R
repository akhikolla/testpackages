## ----echo=FALSE, results="hide"------------------------------------------
# Save space on screen printouts.
options(digits = 3)

## ----eval=FALSE----------------------------------------------------------
#  # Install the package from a remote repository.
#  install.packages("Numero")

## ------------------------------------------------------------------------
# Activate the library.
library("Numero")
packageVersion("Numero")
ls("package:Numero")

## ----eval=FALSE----------------------------------------------------------
#  # Access function documentation (not shown in vignette).
#  ? numero.create

## ----eval=FALSE----------------------------------------------------------
#  # Run all code examples (not shown in vignette).
#  fn <- system.file("extcode", "examples.R", package = "Numero")
#  source(fn)

## ----eval=FALSE----------------------------------------------------------
#  # Show readme file on screen (not shown in vignette).
#  fn <- system.file("extdata", "finndiane.readme.txt", package = "Numero")
#  cat(readChar(fn, 1e5))

## ------------------------------------------------------------------------
# Import data.
fname <- system.file("extdata", "finndiane.txt", package = "Numero")
dataset <- read.delim(file = fname, stringsAsFactors = FALSE)
nrow(dataset)
colnames(dataset)

## ------------------------------------------------------------------------
# Manage unusable entries and data identification.
dataset <- numero.clean(data = dataset, identity = "INDEX")

## ------------------------------------------------------------------------
# Select training variables.
trvars <- c("CHOL", "HDL2C", "TG", "CREAT", "uALB")

## ----results="hide"------------------------------------------------------
# Center and scale the training data.
trdat.basic <- scale.default(dataset[,trvars])

## ------------------------------------------------------------------------
# Calculate standard deviations.
apply(trdat.basic, 2, sd, na.rm = TRUE)

## ------------------------------------------------------------------------
# Create a new self-organizing map based on the training set.
modl.basic <- numero.create(data = trdat.basic)
summary(modl.basic)

## ------------------------------------------------------------------------
# Calculate map quality measures for the training data.
qc.basic <- numero.quality(model = modl.basic)

## ----dev="svg", results="hide", fig.width=7, fig.height=3, fig.align="center", fig.cap="Figure: Distribution of model residuals."----
# Plot frequencies of data points at different quality levels.
par(mar = c(5,4,1,0), mfrow = c(1,2))
hist(x = qc.basic$layout$RESIDUAL, breaks = 50,
     main = NULL, xlab = "RESIDUAL", ylab = "Number of data points",
     col = "#FFEFA0", cex = 0.8)
hist(x = qc.basic$layout$RESIDUAL.z, breaks = 50,
     main = NULL, xlab = "RESIDUAL.z", ylab = "Number of data points",
     col = "#FFEFA0", cex = 0.8)

## ----dev="svg", results="hide", fig.width=9, fig.height=3, fig.align="center", fig.cap="Figure: Visualization of quality measures across map districts. "----
# Plot map quality measures.
numero.plot(results = qc.basic, subplot = c(1,4))

## ------------------------------------------------------------------------
# Map statistics for the whole dataset.
stats.basic <- numero.evaluate(model = qc.basic, data = dataset)
summary(stats.basic)

## ----dev="svg", results="hide", fig.width=6, fig.height=4, fig.align="center", fig.cap="Figure: Statistically normalized colorings of the training variables in the kidney disease dataset. The color intensity depends on how likely the observed regional variation would arise by chance; intense reds and intense blues indicate that these extremes would be very unlikely if the data point layout was random. The numbers show the average values in original units for selected districts."----
# Plot map colorings of training variables.
numero.plot(results = stats.basic, variables = trvars, subplot = c(2,3))

## ------------------------------------------------------------------------
stats.basic$statistics[c("CHOL","MALE","AGE","T1D_DURAT"),
                       c("TRAINING","Z","P.z")]

## ----eval=FALSE----------------------------------------------------------
#  # Interactive subgrouping based on the training set.
#  subgr.basic <- numero.subgroup(results = stats.basic, variables = trvars)

## ----echo=FALSE, results=FALSE-------------------------------------------
# Workaround for the vignette document, DO NOT USE IN REAL STUDIES!
x <- (stats.basic$planes[,"DIAB_KIDNEY"])*(stats.basic$planes[,"uALB"])
tops <- which(x >= quantile(x, 0.75, na.rm=TRUE))
bottoms <- which(x <= quantile(x, 0.25, na.rm=TRUE))
workaround <- as.data.frame(stats.basic$map$topology)
workaround$REGION.label <- " "
workaround$REGION.color <- ""
workaround$REGION <- "not_selected"
workaround$REGION.label[tops] <- "A"
workaround$REGION.color[tops] <- "#FF000050"
workaround$REGION[tops] <- "HighALB"
workaround$REGION.label[bottoms] <- "B"
workaround$REGION.color[bottoms] <- "#89AB0070"
workaround$REGION[bottoms] <- "LowALB"

## ------------------------------------------------------------------------
# Interactive selection not available in vignette.
if(!exists("subgr.basic")) subgr.basic <- workaround

## ----origgrp1, dev="svg", results="hide", fig.width=6, fig.height=4, fig.align="center", fig.cap="Figure: Statistically normalized colorings of the training variables in the kidney disease dataset. The labels show the results from the subgrouping procedure."----
# Plot results from the subgrouping procedure.
numero.plot(results = stats.basic, variables = trvars,
            topology = subgr.basic, subplot = c(2,3))

## ----dev="svg", results="hide", fig.width=6, fig.height=4, fig.align="center", fig.cap="Figure: Statistically normalized colorings of selected variables in the kidney disease dataset."----
# Plot results from subgrouping procedure for non-biochemical variables.
numero.plot(results = stats.basic,
            variables = c("AGE",
                          "MALE",
                          "DIAB_KIDNEY",
                          "DIAB_RETINO",
                          "MACROVASC",
                          "DECEASED"),
            topology = stats.basic$map$topology, subplot = c(2,3))

## ----dev="svg", results="hide", fig.width=6, fig.height=4, fig.align="center", fig.cap="Figure: Colorings of selected variables with subgroup labels."----
# Plot results from subgrouping procedure for non-biochemical variables.
numero.plot(results = stats.basic,
            variables = c("AGE",
                          "MALE",
                          "DIAB_KIDNEY",
                          "DIAB_RETINO",
                          "MACROVASC",
                          "DECEASED"),
            topology = subgr.basic, subplot = c(2,3))

## ------------------------------------------------------------------------
# Compare subgroups.
report.basic <- numero.summary(results = stats.basic, topology = subgr.basic)
colnames(report.basic)

## ------------------------------------------------------------------------
# Show results for mortality rate.
rows <- which(report.basic$VARIABLE == "DECEASED")
report.basic[rows,c("SUBGROUP","N","MEAN","P.chisq")]

## ------------------------------------------------------------------------
# Mitigate stratification and confounding factors.
trdata.adj <- numero.prepare(data = dataset, variables = trvars, batch = "MALE",
                             confounders = c("AGE", "T1D_DURAT"))
colnames(trdata.adj)

## ------------------------------------------------------------------------
subsets <- attr(trdata.adj, "subsets")
women <- subsets[["0"]]
men <- subsets[["1"]]
c(length(men), length(women))

## ------------------------------------------------------------------------
# Create a new self-organizing map based on sex-adjusted data.
modl.adj <- numero.create(data = trdata.adj)
summary(modl.adj)

## ------------------------------------------------------------------------
# Calculate map quality measures for sex-adjusted data.
qc.adj <- numero.quality(model = modl.adj)

## ----dev="svg", results="hide", fig.width=7, fig.height=3, fig.align="center", fig.cap="Figure: Distribution of model residuals. The training data were adjusted for age and sex."----
# Plot frequencies of data points at different quality levels.
par(mar = c(5,4,1,0), mfrow = c(1,2))
hist(x = qc.adj$layout$RESIDUAL, breaks = 20,
     main = NULL, xlab = "RESIDUAL", ylab = "Number of data points",
     col = "#FFEFA0", cex = 0.8)
hist(x = qc.adj$layout$RESIDUAL.z, breaks = 20,
     main = NULL, xlab = "RESIDUAL.z", ylab = "Number of data points",
     col = "#FFEFA0", cex = 0.8)

## ------------------------------------------------------------------------
# Maximum residuals from unadjusted and adjusted analyses.
c(max(qc.basic$layout$RESIDUAL.z, na.rm=TRUE),
  max(qc.adj$layout$RESIDUAL.z, na.rm=TRUE))

## ------------------------------------------------------------------------
# Variation in point density in unadjusted and adjusted analyses.
c(sd(qc.basic$planes[,"HISTOGRAM"], na.rm=TRUE),
  sd(qc.adj$planes[,"HISTOGRAM"], na.rm=TRUE))

## ----dev="svg", results="hide", fig.width=9, fig.height=3, fig.align="center", fig.cap="Figure: Visualization of data point quality measures across map districts from age and sex adjusted analysis."----
# Plot map quality measures.
numero.plot(results = qc.adj, subplot = c(1,4))

## ------------------------------------------------------------------------
# Map statistics for the whole dataset.
stats.adj <- numero.evaluate(model = qc.adj, data = dataset)

## ----results="hide"------------------------------------------------------
# Map statistics for women.
stats.adjW <- numero.evaluate(model = qc.adj, data = dataset[women,])

## ----results="hide"------------------------------------------------------
# Map statistics for men.
stats.adjM <- numero.evaluate(model = qc.adj, data = dataset[men,])

## ------------------------------------------------------------------------
stats.adj$statistics[c("MALE","AGE","T1D_DURAT"), c("Z","P.z")]

## ----dev="svg", results="hide", fig.width=6, fig.height=4, fig.align="center", fig.cap="Figure: Colorings of training variables using the data from female participants."----
numero.plot(results = stats.adjW, variables = trvars, subplot = c(2,3))

## ----dev="svg", results="hide", fig.width=6, fig.height=4, fig.align="center", fig.cap="Figure: Colorings of training variables using the data from male participants."----
numero.plot(results = stats.adjM, variables = trvars, subplot = c(2,3))

## ----dev="svg", results="hide", fig.width=6, fig.height=4, fig.align="center", fig.cap="Figure: Colorings of training variables using the data from female participants. Color scales were derived from the full dataset."----
numero.plot(results = stats.adjW, variables = trvars,
            subplot = c(2,3), reference = stats.adj)

## ----dev="svg", results="hide", fig.width=6, fig.height=4, fig.align="center", fig.cap="Figure: Colorings of training variables using the data from male participants. Color scales were derived from the full dataset."----
numero.plot(results = stats.adjM, variables = trvars,
            subplot = c(2,3), reference = stats.adj)

## ----eval=FALSE----------------------------------------------------------
#  # Interactive subgrouping based on the training set.
#  subgr.adj <- numero.subgroup(results = stats.adj, variables = trvars)

## ----echo=FALSE, results=FALSE-------------------------------------------
# Workaround for the vignette document, DO NOT USE IN REAL STUDIES!
x <- (stats.adj$planes[,"DIAB_KIDNEY"])*(stats.adj$planes[,"uALB"])
tops <- which(x >= quantile(x, 0.75, na.rm=TRUE))
bottoms <- which(x <= quantile(x, 0.25, na.rm=TRUE))
workaround <- as.data.frame(stats.adj$map$topology)
workaround$REGION.label <- " "
workaround$REGION.color <- ""
workaround$REGION <- "not_selected"
workaround$REGION.label[tops] <- "A"
workaround$REGION.color[tops] <- "#FF000050"
workaround$REGION[tops] <- "HighALB"
workaround$REGION.label[bottoms] <- "B"
workaround$REGION.color[bottoms] <- "#89AB0070"
workaround$REGION[bottoms] <- "LowALB"

## ------------------------------------------------------------------------
# Interactive selection not available in vignette.
if(!exists("subgr.adj")) subgr.adj <- workaround

## ----dev="svg", results="hide", fig.width=6, fig.height=4, fig.align="center", fig.cap="Figure: Statistically normalized colorings of the training variables in the kidney disease dataset. The labels show the results from the subgrouping procedure. The map was created from sex-adjusted data."----
# Plot results from subgrouping procedure.
numero.plot(results = stats.adj, variables = trvars,
            topology = subgr.adj, subplot = c(2,3))

## ----dev="svg", results="hide", fig.width=6, fig.height=4, fig.align="center", fig.cap="Figure: Colorings of selected variables with subgroup labels. The map was created from age and sex adjusted data."----
# Plot results from subgrouping procedure for non-biochemical variables.
numero.plot(results = stats.adj,
            variables = c("AGE",
                          "MALE",
                          "DIAB_KIDNEY",
                          "DIAB_RETINO",
                          "MACROVASC",
                          "DECEASED"),
            topology = subgr.adj, subplot = c(2,3))

## ------------------------------------------------------------------------
# District averages from unadjusted analysis.
summary(stats.basic$planes[,"MALE"])

## ------------------------------------------------------------------------
# District averages from adjusted analysis.
summary(stats.adj$planes[,"MALE"])

## ------------------------------------------------------------------------
# Compare subgroups.
report.adj <- numero.summary(results = stats.adj, topology = subgr.adj)

## ------------------------------------------------------------------------
# Show results for mortality rate.
rows <- which(report.adj$VARIABLE == "DECEASED")
report.adj[rows,c("SUBGROUP","N","MEAN","P.chisq")]

## ------------------------------------------------------------------------
ds.discov <- dataset[women,]
ds.replic <- dataset[men,]
ds.mets <- ds.replic[which(ds.replic[,"METAB_SYNDR"] == 1),]

## ----results="hide"------------------------------------------------------
# Discovery cohort and training set.
trdata.discov <- numero.prepare(data = ds.discov, variables = trvars,
                                confounders = c("AGE", "T1D_DURAT"),
                                method = "tapered")

## ------------------------------------------------------------------------
summary(trdata.discov[,"uALB"])

## ----results="hide"------------------------------------------------------
# Replication cohort, version A.
param <- attr(trdata.discov, "pipeline")
trdata.replicA <- numero.prepare(data = ds.replic, pipeline = param)

## ------------------------------------------------------------------------
summary(trdata.replicA[,"uALB"])

## ----results="hide"------------------------------------------------------
# Replication cohort, version B.
trdata.replicB <- numero.prepare(data = ds.replic, variables = trvars,
                             confounders = c("AGE", "T1D_DURAT"),
                             method = "tapered")

## ------------------------------------------------------------------------
summary(trdata.replicB[,"uALB"])

## ----results="hide"------------------------------------------------------
# Replication cohort, MetS version.
trdata.mets <- numero.prepare(data = ds.mets, variables = trvars,
                            confounders = c("AGE", "T1D_DURAT"),
                            method = "tapered")

## ------------------------------------------------------------------------
summary(trdata.mets[,"uALB"])

## ------------------------------------------------------------------------
# Create a new self-organizing map based on sex-adjusted data.
radius.basic <- attr(modl.basic$map$topology, "radius")
modl.discov <- numero.create(data = trdata.discov, radius = radius.basic)
summary(modl.discov)

## ----results="hide"------------------------------------------------------
# Calculate map quality measures.
qc.discov <- numero.quality(model = modl.discov)
qc.replicA <- numero.quality(model = modl.discov, data = trdata.replicA)
qc.replicB <- numero.quality(model = modl.discov, data = trdata.replicB)
qc.mets <- numero.quality(model = modl.discov, data = trdata.mets)

## ------------------------------------------------------------------------
# Define comparable histogram bins.
rz <- c(qc.adj$layout[,"RESIDUAL.z"], qc.discov$layout[,"RESIDUAL.z"])
rz.breaks <- seq(min(rz, na.rm=TRUE), max(rz, na.rm=TRUE), length.out=20)

## ----dev="svg", results="hide", fig.width=7, fig.height=3, fig.align="center", fig.cap="Figure: Distribution of model residuals when the training data were preprocessed by scaling & centering or by tapered ranking."----
# Plot frequencies of data points at different quality levels.
par(mar = c(5,4,1,0), mfrow = c(1,2))
hist(x = qc.adj$layout[,"RESIDUAL.z"], breaks = rz.breaks,
     main = NULL, xlab = "RESIDUAL.z (scale & center)",
     ylab = "Number of data points", col = "#FFEFA0", cex = 0.8)
hist(x = qc.discov$layout[,"RESIDUAL.z"], breaks = rz.breaks,
     main = NULL, xlab = "RESIDUAL.z (rank)",
     ylab = "Number of data points", col = "#FFEFA0", cex = 0.8)

## ------------------------------------------------------------------------
# Comparison of maximum residuals across examples.
r <- c(max(qc.basic$layout[,"RESIDUAL.z"], na.rm=TRUE),
       max(qc.adj$layout[,"RESIDUAL.z"], na.rm=TRUE),
       max(qc.discov$layout[,"RESIDUAL.z"], na.rm=TRUE),
       max(qc.replicA$layout[,"RESIDUAL.z"], na.rm=TRUE),
       max(qc.replicB$layout[,"RESIDUAL.z"], na.rm=TRUE),
       max(qc.mets$layout[,"RESIDUAL.z"], na.rm=TRUE))
names(r) <- c("basic", "adj", "discov", "replicA", "replicB", "mets")
print(r)

## ----dev="svg", results="hide", fig.width=10, fig.height=3, fig.align="center", fig.cap="Figure: Visualization of quality measures across map districts for the discovery cohort. "----
# Plot map quality measures.
numero.plot(results = qc.discov, subplot = c(1,4))

## ----dev="svg", results="hide", fig.width=10, fig.height=3, fig.align="center", fig.cap="Figure: Visualization of quality measures across map districts for the replication dataset A. "----
# Plot map quality measures.
numero.plot(results = qc.replicA, subplot = c(1,4))

## ----dev="svg", results="hide", fig.width=10, fig.height=3, fig.align="center", fig.cap="Figure: Visualization of quality measures across map districts for the replication dataset B. "----
# Plot map quality measures.
numero.plot(results = qc.replicB, subplot = c(1,4))

## ----dev="svg", results="hide", fig.width=10, fig.height=3, fig.align="center", fig.cap="Figure: Visualization of quality measures across map districts for the replication dataset with the metabolic syndrome. "----
# Plot map quality measures.
numero.plot(results = qc.mets, subplot = c(1,4))

## ----results="hide"------------------------------------------------------
# Map statistics for discovery and replication datasets.
stats.discov <- numero.evaluate(model = modl.discov, data = ds.discov)
stats.replicA <- numero.evaluate(model = qc.replicA, data = ds.replic)
stats.replicB <- numero.evaluate(model = qc.replicB, data = ds.replic)
stats.mets <- numero.evaluate(model = qc.mets, data = ds.mets)

## ----dev="svg", results="hide", fig.width=6, fig.height=4, fig.align="center", fig.cap="Figure: Colorings of training variables for the discovery dataset."----
numero.plot(results = stats.discov, variables = trvars,
            gain = 0.8, subplot = c(2,3))

## ----dev="svg", results="hide", fig.width=6, fig.height=4, fig.align="center", fig.cap="Figure: Colorings of selected variables for replication A."----
numero.plot(results = stats.replicA, variables = trvars,
            gain = 0.8, subplot = c(2,3), reference = stats.discov)

## ----dev="svg", results="hide", fig.width=6, fig.height=4, fig.align="center", fig.cap="Figure: Colorings of selected variables for replication B."----
numero.plot(results = stats.replicB, variables = trvars,
            gain = 0.8, subplot = c(2,3), reference = stats.discov)

## ----dev="svg", results="hide", fig.width=6, fig.height=4, fig.align="center", fig.cap="Figure: Colorings of selected variables for MetS dataset."----
numero.plot(results = stats.mets, variables = trvars,
            gain = 0.8, subplot = c(2,3), reference = stats.discov)

## ------------------------------------------------------------------------
# Selection of clinically interesting variables.
clinvars <- c("uALB", "AGE", "DIAB_KIDNEY", "DIAB_RETINO",
              "MACROVASC", "DECEASED")

## ----eval=FALSE----------------------------------------------------------
#  # Interactive subgrouping based on the training set.
#  subgr.discov <- numero.subgroup(results = stats.discov, variables = clinvars)

## ----echo=FALSE, results=FALSE-------------------------------------------
# Workaround for the vignette document, DO NOT USE IN REAL STUDIES!
x <- (stats.discov$planes[,"DIAB_KIDNEY"])*(stats.discov$planes[,"uALB"])
tops <- which(x >= quantile(x, 0.75, na.rm=TRUE))
bottoms <- which(x <= quantile(x, 0.25, na.rm=TRUE))
workaround <- as.data.frame(stats.discov$map$topology)
workaround$REGION.label <- " "
workaround$REGION.color <- ""
workaround$REGION <- "not_selected"
workaround$REGION.label[tops] <- "A"
workaround$REGION.color[tops] <- "#FF000050"
workaround$REGION[tops] <- "HighDiabKD"
workaround$REGION.label[bottoms] <- "B"
workaround$REGION.color[bottoms] <- "#89AB0070"
workaround$REGION[bottoms] <- "LowDiabKD"

## ------------------------------------------------------------------------
# Interactive selection not available in vignette.
if(!exists("subgr.discov")) subgr.discov <- workaround

## ----dev="svg", results="hide", fig.width=6, fig.height=4, fig.align="center", fig.cap="Figure: Statistically normalized colorings of the training variables in the kidney disease dataset. The labels show the results from the subgrouping procedure."----
# Plot results from subgrouping procedure.
numero.plot(results = stats.discov, variables = clinvars,
            topology = subgr.discov, subplot = c(2,3))

## ----eval=FALSE----------------------------------------------------------
#  # Compare subgroups.
#  report.discov <- numero.summary(results = stats.discov,
#                                  topology = subgr.discov)
#  report.replicA <- numero.summary(results = stats.replicA,
#                                   topology = subgr.discov)
#  report.replicB <- numero.summary(results = stats.replicB,
#                                   topology = subgr.discov)
#  report.mets <- numero.summary(results = stats.mets,
#                                topology = subgr.discov)

## ----echo=FALSE, results="hide"------------------------------------------
suppressWarnings({
report.discov <- numero.summary(results = stats.discov,
                                topology = subgr.discov)
report.replicA <- numero.summary(results = stats.replicA,
                                 topology = subgr.discov)
report.replicB <- numero.summary(results = stats.replicB,
                                 topology = subgr.discov)
report.mets <- numero.summary(results = stats.mets,
                              topology = subgr.discov)})

## ------------------------------------------------------------------------
# Show results for mortality rate in the discovery set.
rows <- which(report.discov$VARIABLE == "DECEASED")
report.discov[rows,c("SUBGROUP","N","MEAN","P.chisq")]

## ------------------------------------------------------------------------
# Show results for mortality rate in replication A.
rows <- which(report.replicA$VARIABLE == "DECEASED")
report.replicA[rows,c("SUBGROUP","N","MEAN","P.chisq")]

## ------------------------------------------------------------------------
# Show results for mortality rate in replication B.
rows <- which(report.replicB$VARIABLE == "DECEASED")
report.replicB[rows,c("SUBGROUP","N","MEAN","P.chisq")]

## ------------------------------------------------------------------------
# Show results for mortality rate in metabolic syndrome subset.
rows <- which(report.mets$VARIABLE == "DECEASED")
report.mets[rows,c("SUBGROUP","N","MEAN","P.chisq")]

## ----eval=FALSE----------------------------------------------------------
#  # Save map colorings of training variables.
#  numero.plot(results = stats.basic, variables = trvars,
#              folder = "/tmp/Results")

## ----echo=FALSE----------------------------------------------------------
s <- paste("\n", 
    "*** numero.plot ***\n", 
    "Thu Sep  5 15:44:02 2019\n", 
    "\n", 
    "Resources:\n", 
    "5 column(s) included\n", 
    "destination folder '/tmp/Results'\n", 
    "\n", 
    "Figure 1:\n", 
    "5 subplot(s)\n", 
    "file name '/tmp/Results/figure01.svg'\n", 
    "97622 bytes saved in '/tmp/Results/figure01.svg'\n", 
    "99982 bytes saved in '/tmp/Results/figure01.html'\n", 
    "\n", 
    "Summary:\n", 
    "1 figure(s) -> '/tmp/Results'\n", sep="")
cat(s)

## ----eval=FALSE----------------------------------------------------------
#  # Import topology and region assignments.
#  subgr <- read.delim(file = "Downloads/regions.txt", stringsAsFactors=FALSE)

## ----eval=FALSE----------------------------------------------------------
#  # Calculate subgroup statistics.
#  report <- numero.summary(results = stats.basic, topology = subgr)

## ----echo=FALSE----------------------------------------------------------
sessionInfo()
Sys.time()

