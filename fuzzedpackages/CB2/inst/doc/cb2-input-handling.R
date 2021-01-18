## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## -----------------------------------------------------------------------------
library(CB2)
library(dplyr)
library(readr)

## -----------------------------------------------------------------------------
data("Evers_CRISPRn_RT112")
head(Evers_CRISPRn_RT112$count)

## -----------------------------------------------------------------------------
Evers_CRISPRn_RT112$design

## -----------------------------------------------------------------------------
sgrna_stats <- measure_sgrna_stats(Evers_CRISPRn_RT112$count, Evers_CRISPRn_RT112$design, "before", "after")
gene_stats <- measure_gene_stats(sgrna_stats)
head(gene_stats)

## -----------------------------------------------------------------------------
df <- read_csv("https://raw.githubusercontent.com/hyunhwaj/CB2-Experiments/master/01_gene-level-analysis/data/Evers/CRISPRn-RT112.csv")
df

## -----------------------------------------------------------------------------
head(measure_sgrna_stats(df, Evers_CRISPRn_RT112$design, "before", "after", ge_id = 'gene', sg_id = 'sgRNA'))

