## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- cache = FALSE, include = FALSE------------------------------------------
## copied from examples from the "simmer" R package
## after: https://www.enchufa2.es/archives/suggests-and-vignettes.html
## by Iñaki Úcar
required <- c("lfa", "BEDMatrix", "snpStats", "pryr") # first is not a CRAN package, only suggested since popkin doesn't need it to run... might as well do the same for other dependencies

if (!all(sapply(required, requireNamespace, quietly = TRUE)))
  knitr::opts_chunk$set(eval = FALSE)

## -----------------------------------------------------------------------------
# Data dimensions.
# Choose non-multiples of 4 to test edge cases of BED parsers.
# Number of loci.
m_loci <- 10001
# Number of individuals
n_ind <- 1001
# Overall allele frequency
# (we don't do any actual genetics in this vignette,
# so just set to a reasonable value)
p <- 0.5
# Missingness rate
miss <- 0.1

# Total number of genotypes
n_data <- n_ind * m_loci
# Draw random genotypes from Binomial
X <- rbinom( n_data, 2, p)
# Add missing values
X[ sample(n_data, n_data * miss) ] <- NA
# Turn into matrix
X <- matrix(X, nrow = m_loci, ncol = n_ind)

# Inspect the first 10 individuals at the first 10 loci
X[1:10, 1:10]

## -----------------------------------------------------------------------------
library(genio)

## -----------------------------------------------------------------------------
# We have to specify the number of loci
bim <- make_bim( n = m_loci )

# Inspect the default values
bim

# Let's add the "chr" prefix to the chromosome values,
# so we recognize them when we see them later.
bim$chr <- paste0('chr', bim$chr)
# Make SNP IDs look like "rs" IDs
bim$id <- paste0('rs', bim$id)
# Make positions 1000 bigger
bim$pos <- bim$pos * 1000
# Select randomly between Cs and Gs for the reference alleles
bim$ref <- sample(c('C', 'G'), m_loci, replace = TRUE)
# Select randomly between As and Ts for the alternative alleles
bim$alt <- sample(c('A', 'T'), m_loci, replace = TRUE)

# Inspect the table with our changes
bim

## -----------------------------------------------------------------------------
# Specify the number of individuals
fam <- make_fam( n = n_ind )

# Inspect the default values
fam

# Add prefixes to families and IDs to recognize them later
fam$fam <- paste0('fam', fam$fam)
fam$id <- paste0('id', fam$id)
# Sex values are usually 1 and 2
fam$sex <- sample(1:2, n_ind, replace = TRUE)
# Let's make phenotypes continuous.
# Draw independently from Standard Normal.
fam$pheno <- rnorm(n_ind)
# Let's leave maternal and paternal IDs as missing

# Inspect again
fam

## -----------------------------------------------------------------------------
# Add column and row names from bim and fam tables we just created.
rownames(X) <- bim$id
colnames(X) <- fam$id
# Inspect again the first 10 individuals and loci
X[1:10, 1:10]

## -----------------------------------------------------------------------------
# Will delete at the end of the vignette
file_plink <- tempfile('vignette-random-data')

# Write genotypes, along with the BIM and FAM files we created.
# Omiting them would result in writing the original dummy version of these tables,
# before we edited them.
time_write_genio <- system.time(
    write_plink(file_plink, X, bim, fam)
)
time_write_genio

## -----------------------------------------------------------------------------
# Read the data back in memory.
# Time this step
time_read_genio <- system.time(
    data_genio <- read_plink(file_plink)
)
time_read_genio

# Inspect the data we just read back

# The same random genotypes (first 10 individuals and loci, now with row and column names):
data_genio$X[1:10, 1:10]

# The locus annotations
data_genio$bim

# The individual annotations
data_genio$fam

# Quickly test that the inputs and outputs are identical.
# Genotypes have NAs, so we have to compare this way.
stopifnot( all( X == data_genio$X, na.rm = TRUE) )
stopifnot( bim == data_genio$bim )
# FAM has mixed content (chars, ints, and doubles).
# First 5 columns should be exact match:
stopifnot( fam[,1:5] == data_genio$fam[,1:5] )
# Exact equality may fail for pheno due to precisions, so test this way instead:
stopifnot( max(abs(fam$pheno - data_genio$fam$pheno)) < 1e-4 )

## -----------------------------------------------------------------------------
# Constants
bytes_per_genotype <- 4
bytes_per_gb <- 1024 ^ 3
# Example data dimensions
num_ind <- 1000
num_loci <- 500000
# Gigabytes per 1000 individuals for a typical genotyping array
bytes_per_genotype * num_ind * num_loci / bytes_per_gb

## -----------------------------------------------------------------------------
library(BEDMatrix)
# Time it too.
# Although the BIM and FAM tables are not returned,
# they are partially parsed and kept in memory,
# which can take time for extremely large files
time_read_bedmatrix_1 <- system.time(
    X_BEDMatrix <- BEDMatrix(file_plink)
)

# Inspect the first 10 loci and individuals as usual.
# Note it is transposed compared to our X.
# Also note the column and row names are different from genio's.
X_BEDMatrix[1:10, 1:10]

## -----------------------------------------------------------------------------
# This turns it into a regular R matrix.
# Since most of the reading is actually happening now,
# we time this step now.
time_read_bedmatrix_2 <- system.time(
    X_BEDMatrix_Rmat <- as.matrix(X_BEDMatrix)
)
time_read_bedmatrix_2
# Now we can test that the BEDMatrix agrees with the original matrix we simulated.
# Need to transpose first.
stopifnot( all( X == t(X_BEDMatrix_Rmat), na.rm = TRUE) )

## -----------------------------------------------------------------------------
library(snpStats)

# Read data, time it.
time_read_snpStats_1 <- system.time(
    data_snpStats <- read.plink(file_plink)
)
time_read_snpStats_1

# Inspect the data

# Genotypes.
# Note data is not visualized this way.
# This matrix is also transposed compared to the genio matrix.
data_snpStats$genotypes

# Locus annotations
head( data_snpStats$map )

# Individual annotations
head (data_snpStats$fam )

## -----------------------------------------------------------------------------
# Transpose, then convert to a regular R matrix.
# Let's time this step too.
time_read_snpStats_2 <- system.time(
    X_snpStats <- as( t(data_snpStats$genotypes), 'numeric')
)
time_read_snpStats_2

# Now we can visualize the matrix.
# First 10 loci and individuals, as before.
# Note that, compared to (genio, BEDMatrix, lfa), alleles are encoded in reverse,
# so 0s and 2s are flipped in this matrix.
X_snpStats[1:10, 1:10]

# Again verify that the matrices are identical.
# (Here 2-X flips back 0s and 2s)
stopifnot( all( X == 2 - X_snpStats, na.rm = TRUE) )

## -----------------------------------------------------------------------------
# Let's write this to another file
file_plink_copy <- tempfile('vignette-random-data-copy')

# Copy objects to not change originals
X_snpStats <- X
bim_snpStats <- as.data.frame(bim) # to use rownames
fam_snpStats <- as.data.frame(fam) # ditto

# All data requires matching row and/or column names.
# These first two were already done above.
#rownames(X_snpStats) <- bim$id
#colnames(X_snpStats) <- fam$id
# Row names here are redundant but required.
rownames(bim_snpStats) <- bim$id
rownames(fam_snpStats) <- fam$id

# We shall time several required steps in order to write genotypes in a standard R matrix,
# and the related annotation tables, to BED.
time_write_snpStats <- system.time({
    # Transpose and convert our genotypes to SnpMatrix object.
    # We flip 0s and 2s before converting
    X_snpStats <- as(2 - t(X_snpStats), 'SnpMatrix')
    
    # This complicated command is required to write the data.
    # Although X, fam, and bim are passed as for genio's write_plink,
    # here the name of every column must be specified (there are no reasonable defaults).
    # Interestingly, the parameter names of snpStats' write.plink
    # do not agree with read.plink from the same package.
    write.plink(
        file_plink_copy,
        snps = X_snpStats,
        subject.data = fam_snpStats,
        pedigree = fam,
        id = id,
        father = pat,
        mother = mat,
        sex = sex,
        phenotype = pheno,
        snp.data = bim_snpStats,
        chromosome = chr,
        genetic.distance = posg,
        position = pos,
        allele.1 = ref,
        allele.2 = alt
    )
})

# remove the new file, no longer need it
delete_files_plink(file_plink_copy)

## -----------------------------------------------------------------------------
library(lfa)
# Parse the genotype matrix
time_read_lfa <- system.time(
    X_lfa <- read.bed(file_plink)
)
time_read_lfa

# This genotype matrix does not have column or row names:
X_lfa[1:10, 1:10]

# Again verify that the matrices are identical
stopifnot( all( X == X_lfa, na.rm = TRUE) )

## ---- fig.width = 6, fig.height = 4, fig.align = 'center'---------------------
# Extract component 3 of each time object,
# which is is total time elapsed.
# Sum the two steps it takes for each of BEDMatrix and snpStats to obtain a native R matrix.
times_read <- c(
    time_read_genio[3],
    time_read_bedmatrix_1[3] + time_read_bedmatrix_2[3],
    time_read_snpStats_1[3] + time_read_snpStats_2[3],
    time_read_lfa[3]
)
names_read <- c(
    'genio',
    'BEDMatrix',
    'snpStats',
    'lfa'
)
# Create barplot
barplot(
    times_read,
    names.arg = names_read,
    main = 'BED reader runtimes',
    xlab = 'packages',
    ylab = 'runtime (s)'
)

## ---- fig.width = 4, fig.height = 4, fig.align = 'center'---------------------
times_write <- c(
    time_write_genio[3],
    time_write_snpStats[3]
)
names_write <- c(
    'genio',
    'snpStats'
)
# Create barplot
barplot(
    times_write,
    names.arg = names_write,
    main = 'BED writer runtimes',
    xlab = 'packages',
    ylab = 'runtime (s)'
)

## ---- fig.width = 6, fig.height = 4, fig.align = 'center'---------------------
library(pryr)
# Store directly into a vector
sizes <- c(
    object_size( X ),
    object_size( data_genio$X ),
    object_size( X_BEDMatrix ),
    object_size( data_snpStats$genotypes ),
    object_size( X_lfa )
)
names_sizes <- c(
    'original',
    'genio',
    'BEDMatrix',
    'snpStats',
    'lfa'
)
# Create barplot
barplot(
    sizes,
    names.arg = names_sizes,
    main = 'Native genotype object sizes',
    xlab = 'packages',
    ylab = 'memory (bytes)'
)

## -----------------------------------------------------------------------------
delete_files_plink(file_plink)

## -----------------------------------------------------------------------------
sessionInfo()

