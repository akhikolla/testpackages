## ---- eval=TRUE, warning=FALSE, message=FALSE---------------------------------
# Load required packages
library(alakazam)
library(dplyr)

# Subset example data
data(ExampleDb)
db <- ExampleDb[ExampleDb$sample_id == "+7d", ]

## ---- eval=TRUE, warning=FALSE, fig.width=7.5, fig.height=6-------------------
db_props <- aminoAcidProperties(db, seq="junction", nt=TRUE, trim=TRUE, 
                                label="cdr3")

# The full set of properties are calculated by default
dplyr::select(db_props[1:3, ], starts_with("cdr3"))

# Define a ggplot theme for all plots
tmp_theme <- theme_bw() + theme(legend.position="bottom")

# Generate plots for all four of the properties
g1 <- ggplot(db_props, aes(x=c_call, y=cdr3_aa_length)) + tmp_theme +
    ggtitle("CDR3 length") + 
    xlab("Isotype") + ylab("Amino acids") +
    scale_fill_manual(name="Isotype", values=IG_COLORS) +
    geom_boxplot(aes(fill=c_call))
g2 <- ggplot(db_props, aes(x=c_call, y=cdr3_aa_gravy)) + tmp_theme + 
    ggtitle("CDR3 hydrophobicity") + 
    xlab("Isotype") + ylab("GRAVY") +
    scale_fill_manual(name="Isotype", values=IG_COLORS) +
    geom_boxplot(aes(fill=c_call))
g3 <- ggplot(db_props, aes(x=c_call, y=cdr3_aa_basic)) + tmp_theme +
    ggtitle("CDR3 basic residues") + 
    xlab("Isotype") + ylab("Basic residues") +
    scale_y_continuous(labels=scales::percent) +
    scale_fill_manual(name="Isotype", values=IG_COLORS) +
    geom_boxplot(aes(fill=c_call))
g4 <- ggplot(db_props, aes(x=c_call, y=cdr3_aa_acidic)) + tmp_theme +
    ggtitle("CDR3 acidic residues") + 
    xlab("Isotype") + ylab("Acidic residues") +
    scale_y_continuous(labels=scales::percent) +
    scale_fill_manual(name="Isotype", values=IG_COLORS) +
    geom_boxplot(aes(fill=c_call))

# Plot in a 2x2 grid
gridPlot(g1, g2, g3, g4, ncol=2)

## ---- eval=TRUE, warning=FALSE------------------------------------------------
db_props <- aminoAcidProperties(db, seq="junction", property=c("gravy", "charge"),
                                nt=TRUE, trim=TRUE, label="cdr3")
dplyr::select(db_props[1:3, ], starts_with("cdr3"))

## ---- eval=TRUE, warning=FALSE, message=FALSE---------------------------------
# Load the relevant data objects from the seqinr package
library(seqinr)
data(aaindex)
data(pK)
h <- aaindex[["KIDA850101"]]$I
p <- setNames(pK[["Murray"]], rownames(pK))
# Rename the hydrophobicity vector to use single-letter codes
names(h) <- translateStrings(names(h), ABBREV_AA)
db_props <- aminoAcidProperties(db, seq="junction", property=c("gravy", "charge"), 
                                nt=TRUE, trim=TRUE, label="cdr3", 
                                hydropathy=h, pK=p)
dplyr::select(db_props[1:3, ], starts_with("cdr3"))

## ---- eval=TRUE, warning=FALSE, message=FALSE---------------------------------
# Translate junction DNA sequences to amino acids and trim first and last codons
cdr3 <- translateDNA(db$junction[1:3], trim=TRUE)

# Grand average of hydrophobicity
gravy(cdr3)

# Average bulkiness
bulk(cdr3)

# Average polarity
polar(cdr3)

# Normalized aliphatic index
aliphatic(cdr3)
# Unnormalized aliphatic index
aliphatic(cdr3, normalize=FALSE)

# Normalized net charge
charge(cdr3)
# Unnormalized net charge
charge(cdr3, normalize=FALSE)

# Count of acidic amino acids
# Takes a named list of regular expressions
countPatterns(cdr3, c(ACIDIC="[DE]"), label="cdr3")

