## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ---- echo=FALSE, fig.cap="\\label{fig:wallthickness}**Left ventricular wall thickness and contributing factors.** Left ventricular wall thickness was measured at more than 27k positions in a cohort of 1,185 healthy volunteers. The average wall thickness at each position is depicted in saturated colours (the mesh depicts the right ventricle for reference). Heart wall thickness is determined by an interplay of genetic (green arrows) and non-genetic factors (grey arrows). Additional correlation between positional thickness can be observed for measures in close spatial proximity.", fig.align='center'----
knitr::include_graphics(system.file("extdata/vignettes/wallthickness-min.png",
                                    package = "PhenotypeSimulator"))

## ----load real data------------------------------------------------------
library(PhenotypeSimulator)
library(ggplot2)

phenofile <- system.file("extdata/vignettes/pheno_small.rds", 
                            package = "PhenotypeSimulator")
covsfile <- system.file("extdata/vignettes/covs.rds", 
                            package = "PhenotypeSimulator")

pheno <- readRDS(phenofile)
covs <- readRDS(covsfile)

dim(pheno)
dim(covs)
colnames(covs)

## ----distribution, fig.cap="\\label{fig:distribution_covs}**Distribution of covariates.** Histograms of the reduced cohort data (100 individuals) for sex, age, height, weight and number of 2D cardiac MRI slices.", fig.align='center'----
covs_melt <- reshape2::melt(covs, value.name="meassure", 
                            variable.name="covariate")
covs_melt$covariate <- factor(covs_melt$covariate,
                              levels=c("Age", "Height", "Weight", "Sex", 
                                       "Slices"),
                              labels= c("Age [years]", "Height [cm]", 
                                        "Weight [kg]", "Sex", "Slices"))
p_cont <- ggplot(dplyr::filter(covs_melt, covariate %in% c("Age [years]", 
                                                           "Height [cm]", 
                                                           "Weight [kg]")), 
                   aes(x=meassure))
p_cont <- p_cont + geom_histogram(binwidth=1) +
    facet_wrap(~covariate, scales = "free") +
    theme_bw() +
    theme(axis.title.x=element_blank(),
          strip.background = element_rect(fill="white", color="black"))

p_cat <- ggplot(dplyr::filter(covs_melt, covariate %in% c("Sex", "Slices")), 
                   aes(x=meassure))
p_cat <- p_cat + geom_bar(aes(x=meassure)) +
    facet_wrap(~covariate, scales = "free") +
    theme_bw() +
    theme(axis.title.x=element_blank(),
          strip.background = element_rect(fill="white", color="black"))

combine <- cowplot::plot_grid(p_cont, p_cat, nrow=2) 
print(combine)

## ----correlation phenotypes, eval=FALSE----------------------------------
#  cor_pheno <- cor(pheno)
#  gplots::heatmap.2(cor_pheno, trace="none", dendrogram="none", keysize=1,
#                    col=colorRampPalette(c("white", "#5BBCD6")),
#                    labRow=FALSE, labCol=FALSE,
#                    breaks=seq(0,1,0.001),
#                    key.title="",
#                    key.xlab="Pearson's correlation coefficient",
#                    key.ylab="",
#                    density.info="none")

## ---- echo=FALSE, out.width='\\textwidth', fig.align='center',fig.cap="\\label{fig:corpheno}\\textbf{Correlation of left ventricular wall thickness measurements.}"----
cor_pheno <- cor(pheno)
knitr::include_graphics(system.file("extdata/vignettes/cor_pheno-min.png",
            package = "PhenotypeSimulator"))

## ----data parameters-----------------------------------------------------
age_range <- range(covs$Age)

weight_mean <- mean(covs$Weight)
weight_sd <- sd(covs$Weight)

height_mean <- mean(covs$Height)
height_sd <- sd(covs$Height)

slices_categories <- length(unique(covs$Slices))

sex_proportion_male <- length(which(covs$Sex==1))/nrow(covs)

correlation <- quantile(unlist(cor_pheno[lower.tri(cor_pheno)]), 0.999)

## ----simulate data-------------------------------------------------------
set.seed(25)
covs_simulated <- noiseFixedEffects(N=100, P=1000, NrFixedEffects = 5, 
                                    NrConfounders = 1,
                  distConfounders = c("unif", "norm", "norm", "cat_unif", "bin"),
                  mConfounders = c(mean(age_range), weight_mean, height_mean),
                  sdConfounders = c(age_range[2] -mean(age_range), weight_sd, 
                                    height_sd),
                  probConfounders = sex_proportion_male,
                  catConfounders = slices_categories)

background_simulated <- correlatedBgEffects(N=100, P=1000, pcorr=correlation)


## ----simulation, fig.cap="\\label{fig:distribution_simulated}**Distribution of simulated covariates.**", fig.align='center'----
covs_simulated_melt <- reshape2::melt(covs_simulated$cov, value.name="meassure", 
                            variable.name="covariate")
covs_simulated_melt$covariate <- factor(covs_simulated_melt$covariate,
                              levels=c("sharedConfounder1_unif1", 
                                       "sharedConfounder2_norm1", 
                                       "sharedConfounder3_norm1", 
                                       "sharedConfounder5_bin1", 
                                       "sharedConfounder4_cat_unif1"),
                              labels= c("Age simulated", "Height simulated", 
                                        "Weight simulated", "Sex simulated", 
                                        "Slices simulated"))
p_cont <- ggplot(dplyr::filter(covs_simulated_melt, covariate %in% 
                                   c("Age simulated", "Height simulated", 
                                                           "Weight simulated")), 
                   aes(x=meassure))
p_cont <- p_cont + geom_histogram(binwidth=1) +
    facet_wrap(~covariate, scales = "free") +
    theme_bw() +
    theme(axis.title.x=element_blank(),
          strip.background = element_rect(fill="white", color="black"))

p_cat <- ggplot(dplyr::filter(covs_simulated_melt, covariate %in% 
                                  c("Sex simulated", "Slices simulated")), 
                   aes(x=meassure))
p_cat <- p_cat + geom_bar(aes(x=meassure)) +
    facet_wrap(~covariate, scales = "free") +
    theme_bw() +
    theme(axis.title.x=element_blank(),
          strip.background = element_rect(fill="white", color="black"))

combine <- cowplot::plot_grid(p_cont, p_cat, nrow=2) 
print(combine)

## ----correlation simulated background, eval=FALSE------------------------
#  gplots::heatmap.2(cor(background_simulated$correlatedBg), trace="none",
#                    dendrogram="none",
#                    keysize=1, col=colorRampPalette(c("white", "#5BBCD6")),
#                    labRow=FALSE, labCol=FALSE,
#                    breaks=seq(0,1,0.001),
#                    key.title="",
#                    key.xlab="Pearson's correlation coefficient",
#                    key.ylab="",
#                    density.info="none")

## ----load simulated background, echo=FALSE, out.width='\\textwidth', fig.align='center', fig.cap="\\label{fig:corsimulated}\\textbf{Correlation of simulated spatial proximity.}"----
knitr::include_graphics(system.file("extdata/vignettes/cor_simulated-min.png",
            package = "PhenotypeSimulator"))

