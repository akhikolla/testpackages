## ---- eval=TRUE, warning=FALSE, message=FALSE---------------------------------
# Imports
library(scoper)
library(dplyr)

# Clonal assignment using identical nucleotide sequences
results <- identicalClones(ExampleDb, method="nt")

## ---- eval=TRUE, warning=FALSE, message=FALSE---------------------------------
# Get results data.frame
results_db <- as.data.frame(results)
glimpse(results_db)

## ---- eval=TRUE, warning=FALSE, message=FALSE---------------------------------
# Plot a histogram of inter clonal distances
plot(results, binwidth=0.02)

# Get summary data.frame
glimpse(summary(results))

## ---- eval=TRUE, warning=FALSE, message=FALSE---------------------------------
# Clonal assignment using hierarchical clustering
results <- hierarchicalClones(ExampleDb, threshold=0.15)

## ---- eval=TRUE, warning=FALSE, message=FALSE---------------------------------
# Get results data.frame
results_db <- as.data.frame(results)
glimpse(results_db)

# Plot a histogram of inter and intra clonal distances
plot(results, binwidth=0.02)

# Get summary data.frame
glimpse(summary(results))

## ---- eval=TRUE, warning=FALSE, message=FALSE---------------------------------
# Clonal assignment using the spectral clustering method novj
results <- spectralClones(ExampleDb, method="novj")
# Plot a histogram of inter and intra clonal distances
plot(results, binwidth=0.02)

## ---- eval=TRUE, warning=FALSE, message=FALSE---------------------------------
# Clonal assignment using the spectral clustering method novj with threshold
results <- spectralClones(ExampleDb, method="novj",
                          threshold=0.15)
# Plot a histogram of inter and intra clonal distances
plot(results, binwidth=0.02)

## ---- eval=TRUE, warning=FALSE, message=FALSE---------------------------------
# Clonal assignment using the spectral clustering method vj with threshold
results <- spectralClones(ExampleDb, method="vj",
                          threshold=0.15,
                          germline="germline_alignment_d_mask")

## ---- eval=TRUE, warning=FALSE, message=FALSE---------------------------------
# Get results data.frame
results_db <- as.data.frame(results)
glimpse(results_db)

# Plot a histogram of inter and intra clonal distances
plot(results, binwidth=0.02)

# Get summary data.frame
glimpse(summary(results))

