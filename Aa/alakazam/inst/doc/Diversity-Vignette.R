## ---- eval=TRUE, warning=FALSE, message=FALSE---------------------------------
# Load required packages
library(alakazam)

# Load example data
data(ExampleDb)

## ---- eval=TRUE, warning=FALSE------------------------------------------------
# Partitions the data based on the sample column
clones <- countClones(ExampleDb, group="sample_id")
head(clones, 5)

## ---- eval=TRUE, warning=FALSE------------------------------------------------
# Partitions the data based on both the sample_id and c_call columns
# Weights the clone sizes by the duplicate_count column
clones <- countClones(ExampleDb, group=c("sample_id", "c_call"), copy="duplicate_count", clone="clone_id")
head(clones, 5)

## ---- eval=TRUE, results='hide', warning=FALSE, fig.width=6, fig.height=4-----
# Partitions the data on the sample column
# Calculates a 95% confidence interval via 200 bootstrap realizations
curve <- estimateAbundance(ExampleDb, group="sample_id", ci=0.95, nboot=200, clone="clone_id")

## ---- eval=TRUE, warning=FALSE, fig.width=6, fig.height=4---------------------
# Plots a rank abundance curve of the relative clonal abundances
sample_colors <- c("-1h"="seagreen", "+7d"="steelblue")
plot(curve, colors = sample_colors, legend_title="Sample")

## ---- eval=TRUE, results='hide'-----------------------------------------------
# Compare diversity curve across values in the "sample" column
# q ranges from 0 (min_q=0) to 4 (max_q=4) in 0.05 increments (step_q=0.05)
# A 95% confidence interval will be calculated (ci=0.95)
# 200 resampling realizations are performed (nboot=200)
sample_curve <- alphaDiversity(ExampleDb, group="sample_id", clone="clone_id",
                               min_q=0, max_q=4, step_q=0.1,
                               ci=0.95, nboot=200)

# Compare diversity curve across values in the c_call column
# Analyse is restricted to c_call values with at least 30 sequences by min_n=30
# Excluded groups are indicated by a warning message
isotype_curve <- alphaDiversity(ExampleDb, group="c_call", clone="clone_id",
                                min_q=0, max_q=4, step_q=0.1,
                                ci=0.95, nboot=200)

## ---- eval=TRUE, fig.width=6, fig.height=4------------------------------------
# Plot a log-log (log_q=TRUE, log_d=TRUE) plot of sample diversity
# Indicate number of sequences resampled from each group in the title
sample_main <- paste0("Sample diversity")
sample_colors <- c("-1h"="seagreen", "+7d"="steelblue")
plot(sample_curve, colors=sample_colors, main_title=sample_main, 
     legend_title="Sample")

# Plot isotype diversity using default set of Ig isotype colors
isotype_main <- paste0("Isotype diversity")
plot(isotype_curve, colors=IG_COLORS, main_title=isotype_main, 
     legend_title="Isotype")

## ---- eval=TRUE, fig.width=6, fig.height=3------------------------------------
# Test diversity at q=0, q=1 and q=2 (equivalent to species richness, Shannon entropy, 
# Simpson's index) across values in the sample_id column
# 200 bootstrap realizations are performed (nboot=200)
isotype_test <- alphaDiversity(ExampleDb, group="c_call", min_q=0, max_q=2, step_q=1, nboot=200, clone="clone_id")

# Print P-value table
print(isotype_test@tests)

# Plot results at q=0 and q=2
# Plot the mean and standard deviations at q=0 and q=2
plot(isotype_test, 0, colors=IG_COLORS, main_title=isotype_main, 
     legend_title="Isotype")
plot(isotype_test, 2, colors=IG_COLORS, main_title=isotype_main, 
     legend_title="Isotype")

