# create example data and save as RData file
library(lineup2)

# download the data file (913 MB)
url <- "https://phenomedoc.jax.org/QTL_Archive/attie_2015/Attie_2015_eqtl_raw.zip"
tmpdir <- tempdir()
file <- file.path(tmpdir, basename(url))
if(!file.exists(file)) {
    download.file(url, file)
}

# extract two data files
csvfiles <- file.path("Raw", c("gastroc_mlratio_raw.csv", "islet_mlratio_raw.csv"))
unzip(file, csvfiles, exdir=tmpdir)

# read in the data files
lineup2ex <- lapply(csvfiles, function(file) {
    x <- read.csv(file.path(tmpdir, file), header=TRUE, stringsAsFactors=FALSE)
    rownames(x) <- x[,1]
    x[,-1] })

# reduce to a selected subset: 100 highly correlated genes + 100 random others
# check that columns are aligned
stopifnot(all(colnames(lineup2ex[[1]]) == colnames(lineup2ex[[2]])))
# first align the rows
aligned <- align_matrix_rows(lineup2ex[[1]], lineup2ex[[2]])
# calculate correlations between column pairs
rho <- corr_betw_matrices(aligned[[1]], aligned[[2]], "paired")

# pick out the most-correlated pairs plus others
best <- order(rho, decreasing=TRUE)[1:100]
set.seed(20201015)
other <- sample(seq_along(rho)[-best], 100, replace=TRUE)
keep <- sort(c(other, best))

# subset the columns and add names (taken from file names)
lineup2ex <- lapply(lineup2ex, function(a) a[,keep])
names(lineup2ex) <- sapply(strsplit(csvfiles, "[/_]"), "[", 2)

# save to file
save(lineup2ex, file="../../data/lineup2ex.RData", compress="xz")

# cleanup
unlink(file)
unlink(file.path(tmpdir, "Raw"), recursive=TRUE)
