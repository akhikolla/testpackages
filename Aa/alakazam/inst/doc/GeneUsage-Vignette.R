## ---- eval=TRUE, warning=FALSE, message=FALSE---------------------------------
# Load required packages
library(alakazam)
library(dplyr)
library(scales)

# Subset example data
data(ExampleDb)

## ---- eval=TRUE, warning=FALSE------------------------------------------------
# Quantify usage at the gene level
gene <- countGenes(ExampleDb, gene="v_call", groups="sample_id", mode="gene")
head(gene, n=4)

## ---- eval=TRUE, warning=FALSE------------------------------------------------
# Assign sorted levels and subset to IGHV1
ighv1 <- gene %>%
    mutate(gene=factor(gene, levels=sortGenes(unique(gene), method="name"))) %>%
    filter(getFamily(gene) == "IGHV1")

# Plot V gene usage in the IGHV1 family by sample
g1 <- ggplot(ighv1, aes(x=gene, y=seq_freq)) +
    theme_bw() +
    ggtitle("IGHV1 Usage") +
    theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
    ylab("Percent of repertoire") +
    xlab("") +
    scale_y_continuous(labels=percent) +
    scale_color_brewer(palette="Set1") +
    geom_point(aes(color=sample_id), size=5, alpha=0.8)
plot(g1)

## ---- eval=TRUE, warning=FALSE------------------------------------------------
# Quantify V family usage by sample
family <- countGenes(ExampleDb, gene="v_call", groups="sample_id", mode="family")

# Plot V family usage by sample
g2 <- ggplot(family, aes(x=gene, y=seq_freq)) +
    theme_bw() +
    ggtitle("Family Usage") +
    theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
    ylab("Percent of repertoire") +
    xlab("") +
    scale_y_continuous(labels=percent) +
    scale_color_brewer(palette="Set1") +
    geom_point(aes(color=sample_id), size=5, alpha=0.8)
plot(g2)

## ---- eval=TRUE, warning=FALSE------------------------------------------------
# Quantify V family clonal usage by sample and isotype
family <- countGenes(ExampleDb, gene="v_call", groups=c("sample_id", "c_call"), 
                     clone="clone_id", mode="family")
head(family, n=4)

## ---- eval=TRUE, warning=FALSE------------------------------------------------
# Subset to IGHM and IGHG for plotting
family <- filter(family, c_call %in% c("IGHM", "IGHG"))
# Plot V family clonal usage by sample and isotype
g3 <- ggplot(family, aes(x=gene, y=clone_freq)) +
    theme_bw() +
    ggtitle("Clonal Usage") +
    theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
    ylab("Percent of repertoire") +
    xlab("") +
    scale_y_continuous(labels=percent) +
    scale_color_brewer(palette="Set1") +
    geom_point(aes(color=sample_id), size=5, alpha=0.8) +
    facet_grid(. ~ c_call)
plot(g3)

## ---- eval=TRUE, warning=FALSE------------------------------------------------
# Calculate V family copy numbers by sample and isotype
family <- countGenes(ExampleDb, gene="v_call", groups=c("sample_id", "c_call"), 
                     mode="family", copy="duplicate_count")
head(family, n=4)

## ---- eval=TRUE, warning=FALSE------------------------------------------------
# Subset to IGHM and IGHG for plotting
family <- filter(family, c_call %in% c("IGHM", "IGHG"))
# Plot V family copy abundance by sample and isotype
g4 <- ggplot(family, aes(x=gene, y=copy_freq)) +
    theme_bw() +
    ggtitle("Copy Number") +
    theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
    ylab("Percent of repertoire") +
    xlab("") +
    scale_y_continuous(labels=percent) +
    scale_color_brewer(palette="Set1") +
    geom_point(aes(color=sample_id), size=5, alpha=0.8) +
    facet_grid(. ~ c_call)
plot(g4)

