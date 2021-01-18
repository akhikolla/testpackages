## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(ampir)

## ---- warning=FALSE, message=FALSE--------------------------------------------
my_protein_df <- read_faa(system.file("extdata/little_test.fasta", package = "ampir"))

## ---- echo=FALSE--------------------------------------------------------------
display_df <- my_protein_df
display_df$seq_aa <- paste(substring(display_df$seq_aa,1,45),"...",sep="")
knitr::kable(display_df)

## -----------------------------------------------------------------------------
my_prediction <- predict_amps(my_protein_df, model = "precursor")

## ---- echo=FALSE--------------------------------------------------------------
my_prediction$seq_aa <- paste(substring(my_prediction$seq_aa,1,45),"...",sep="")
knitr::kable(my_prediction, digits = 3)

## -----------------------------------------------------------------------------
my_predicted_amps <- my_protein_df[my_prediction$prob_AMP > 0.8,]

## ---- echo=FALSE--------------------------------------------------------------
my_predicted_amps$seq_aa <- paste(substring(my_predicted_amps$seq_aa,1,45),"...",sep="")
knitr::kable(my_predicted_amps)

## -----------------------------------------------------------------------------
df_to_faa(my_predicted_amps, tempfile("my_predicted_amps.fasta", tempdir()))

