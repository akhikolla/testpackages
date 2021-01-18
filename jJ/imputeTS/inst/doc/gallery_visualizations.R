## ----knitr-setup, include = FALSE---------------------------------------------
knitr::opts_chunk$set(fig.align = "center",
                      fig.width = 6,
                      fig.height = 4,
                      dpi = 100)

## ----ggplot-na-distribution, message=FALSE------------------------------------
library("imputeTS")
ggplot_na_distribution(tsAirgap)

## ----create-df, results=F, echo = F, warning=FALSE, fig.show='hide'-----------
df <- structure(list(date = structure(c(-21185, -20819, -20454, -20089, 
-19724, -19358, -18993, -18628, -18263, -17897, -17532, -17167, 
-16802, -16436, -16071, -15706, -15341, -14975, -14610, -14245, 
-13880, -13514, -13149, -12784, -12419, -12053, -11688, -11323, 
-10958, -10592, -10227, -9862, -9497, -9131, -8766, -8401, -8036, 
-7670, -7305, -6940, -6575, -6209, -5844, -5479, -5114, -4748, 
-4383, -4018, -3653, -3287, -2922, -2557, -2192, -1826, -1461, 
-1096, -731, -365, 0, 365), class = "Date"), value = structure(c(48.2, 
50.5, 49.4, 51.1, 49.4, 47.9, 49.8, 50.9, 49.3, 51.9, 50.8, 49.6, 
49.3, 50.6, 48.4, 50.7, 50.9, 50.6, 51.5, 52.8, 51.8, 51.1, 49.8, 
50.2, 50.4, NA, NA, NA, 48.8, 51.7, 51, 50.6, 51.7, 51.5, 
52.1, 51.3, 51, 54, 51.4, 52.7, 53.1, 54.6, NA, 52, 50.9, 52.6, 
50.2, 52.6, 51.6, 51.9, 50.5, 50.9, 51.7, 51.4, 51.7, 50.8, 51.9, 
51.8, 50.0, 49.1 ), .Tsp = c(1912, 1971, 1), class = "ts")), class = "data.frame", row.names = c(NA, 
-60L))

## ----ggplot-na-distribution2--------------------------------------------------
ggplot_na_distribution(x = df$value, x_axis_labels = df$date)

## ----ggplot-na-intervals------------------------------------------------------
ggplot_na_intervals(tsNH4)

## ----ggplot-na-intervals2-----------------------------------------------------
ggplot_na_intervals(tsNH4, measure = "count", interval_size = 144, color_missing = "gold3")

## ----ggplot-na-gapsize--------------------------------------------------------
library(imputeTS)
ggplot_na_gapsize(tsNH4)

## ----ggplot-na-gapsize2-------------------------------------------------------
library(imputeTS)
ggplot_na_gapsize(tsNH4, include_total = F, limit = 15)

## ----ggplot-na-imputations1---------------------------------------------------
library(imputeTS)
imp <- na_interpolation(tsAirgap)
ggplot_na_imputations(tsAirgap, imp)

## ----ggplot-na-imputations2---------------------------------------------------
library(imputeTS)
imp <- na_mean(tsAirgap)
ggplot_na_imputations(x_with_na = tsAirgap, x_with_imputations = imp, x_with_truth = tsAirgapComplete )

