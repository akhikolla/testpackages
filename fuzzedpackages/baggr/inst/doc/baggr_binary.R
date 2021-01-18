## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE, fig.width=7, 
  fig.height = 5, comment = "#>")
library(baggr)
library(ggplot2)
library(gridExtra)

## -----------------------------------------------------------------------------
df_yusuf <- read.table(text="
       trial  a n1i  c n2i
      Balcon 14  56 15  58
     Clausen 18  66 19  64
 Multicentre 15 100 12  95
      Barber 10  52 12  47
      Norris 21 226 24 228
      Kahler  3  38  6  31
     Ledwich  2  20  3  20
", header=TRUE)

## -----------------------------------------------------------------------------
df <- df_yusuf
df$b <- df$n1i-df$a
df$d <- df$n2i-df$c
df$tau <- log((df$a*df$d)/(df$b*df$c))
df$se <- sqrt(1/df$a + 1/df$b + 1/df$c + 1/df$d)

df

## ---- include = F-------------------------------------------------------------
bg_model_agg <- baggr(df, 
                    effect = "logarithm of odds ratio",
                    group = "trial")

## ---- eval = F, echo = T------------------------------------------------------
#  bg_model_agg <- baggr(df,
#                      effect = "logarithm of odds ratio",
#                      group = "trial")

## -----------------------------------------------------------------------------
a <- 9; b <- 1; c <- 99; d <- 1
cat("Risk ratio is", (a/(a+b))/(c/(c+d)), "\n" )
cat("Odds ratio is", a*d/(b*c), "\n")

## -----------------------------------------------------------------------------
a <- 10; b <- 20; c <- 100; d <- 100
cat("Risk ratio is", (a/(a+b))/(c/(c+d)), "\n" )
cat("Odds ratio is", a*d/(b*c), "\n")

## -----------------------------------------------------------------------------
bg_model_agg

## -----------------------------------------------------------------------------
forest_plot(bg_model_agg, show = "both", print = "inputs")

## -----------------------------------------------------------------------------
effect_plot(bg_model_agg)

## ---- warning = FALSE---------------------------------------------------------
gridExtra::grid.arrange(
  plot(bg_model_agg, transform = exp) + xlab("Effect on OR"),
  effect_plot(bg_model_agg, transform = exp) + xlim(0, 3) + xlab("Effect on OR"),
  ncol = 2)

## ---- echo=T, eval = F--------------------------------------------------------
#  # Instead of writing...
#  # bg1 <- baggr(df, pooling = "none")
#  # bg2 <- baggr(df, pooling = "partial")
#  # bg3 <- baggr(df, pooling = "full")
#  
#  # ...we use this one-liner
#  bg_c <- baggr_compare(df, effect = "logarithm of odds ratio")

## ---- include = FALSE---------------------------------------------------------
bg_c <- baggr_compare(df, what = "pooling", effect = "logarithm of odds ratio")

## -----------------------------------------------------------------------------
plot(bg_c)

## -----------------------------------------------------------------------------
effect_plot(
  "Partial pooling, default prior" = bg_c$models$partial,
  "Full pooling, default prior" = bg_c$models$full) +
  theme(legend.position = "bottom")

## ---- include = F, warning = F------------------------------------------------
a <- loocv(df, pooling = "partial")
b <- loocv(df, pooling = "full")
#a; b; #you can print out individual loocv() calculations
loo_compare(a,b) #...but typically we compare them to each other

## ---- echo = T, eval = F------------------------------------------------------
#  a <- loocv(df, pooling = "partial")
#  b <- loocv(df, pooling = "full")
#  #a; b; #you can print out individual loocv() calculations
#  loo_compare(a,b) #...but typically we compare them to each other

## -----------------------------------------------------------------------------
df_ind <- data.frame()
for(i in 1:nrow(df_yusuf)) {
  df_ind <- rbind(df_ind, data.frame(group = df_yusuf$trial[i],
             treatment = c(rep(1, df_yusuf$n1i[i]), rep(0, df_yusuf$n2i[i])),
             outcome = c(rep(1, df_yusuf$a[i]), rep(0, df_yusuf$n1i[i] - df_yusuf$a[i]),
                         rep(1, df_yusuf$c[i]), rep(0, df_yusuf$n2i[i] - df_yusuf$c[i]))))
}
head(df_ind)

## ----logit model, include = F-------------------------------------------------
bg_model_ind <- baggr(df_ind, model = "logit", effect = "logarithm of odds ratio")

## ---- echo = T, eval = F------------------------------------------------------
#  bg_model_ind <- baggr(df_ind, model = "logit", effect = "logarithm of odds ratio")

## -----------------------------------------------------------------------------
baggr_compare(bg_model_agg, bg_model_ind)

## -----------------------------------------------------------------------------
prepare_ma(df_ind, effect = "logOR")
prepare_ma(df_ind, effect = "logRR")

## -----------------------------------------------------------------------------
df_rare <- data.frame(group = paste("Study", LETTERS[1:4]),
                      a = c(0, 2, 1, 3), c = c(2, 2, 3, 3),
                      n1i = c(120, 300, 110, 250),
                      n2i = c(120, 300, 110, 250))

## -----------------------------------------------------------------------------
df_rare_ind <- data.frame()
                      
for(i in 1:nrow(df_rare)) {
  df_rare_ind <- rbind(df_rare_ind, data.frame(group = df_rare$group[i],
             treatment = c(rep(1, df_rare$n1i[i]), rep(0, df_rare$n2i[i])),
             outcome = c(rep(1, df_rare$a[i]), rep(0, df_rare$n1i[i] - df_rare$a[i]),
                         rep(1, df_rare$c[i]), rep(0, df_rare$n2i[i] - df_rare$c[i]))))
}

## -----------------------------------------------------------------------------
df_rare_logor <- prepare_ma(df_rare_ind, effect = "logOR")
df_rare_logor

## -----------------------------------------------------------------------------
pma01 <- prepare_ma(df_rare_ind, effect = "logOR", 
                            rare_event_correction = 0.1)
pma1 <- prepare_ma(df_rare_ind, effect = "logOR", 
                            rare_event_correction = 1)

## ----rare event comparison, include=F-----------------------------------------
bg_correction01 <- baggr(pma01, effect = "logOR")
bg_correction025 <- baggr(df_rare_logor, effect = "logOR")
bg_correction1 <- baggr(pma1, effect = "logOR")
bg_rare_ind <- baggr(df_rare_ind, effect = "logOR")

## ---- echo=T, eval=F----------------------------------------------------------
#  bg_correction01 <- baggr(pma01, effect = "logOR")
#  bg_correction025 <- baggr(df_rare_logor, effect = "logOR")
#  bg_correction1 <- baggr(pma1, effect = "logOR")
#  bg_rare_ind <- baggr(df_rare_ind, effect = "logOR")

## -----------------------------------------------------------------------------
bgc <- baggr_compare(
  "Correct by .10" = bg_correction01,
  "Correct by .25" = bg_correction025,
  "Correct by 1.0" = bg_correction1,
  "Individual data" = bg_rare_ind
)
bgc
plot(bgc)

## -----------------------------------------------------------------------------
df_rare <- data.frame(group = paste("Study", LETTERS[1:4]),
                      a = c(1, 2, 1, 3), c = c(2, 2, 3, 3),
                      n1i = c(120, 300, 110, 250),
                      n2i = c(120, 300, 110, 250))

df_rare_ind <- data.frame()
                      
for(i in 1:nrow(df_rare)) {
  df_rare_ind <- rbind(df_rare_ind, data.frame(group = df_rare$group[i],
             treatment = c(rep(1, df_rare$n1i[i]), rep(0, df_rare$n2i[i])),
             outcome = c(rep(1, df_rare$a[i]), rep(0, df_rare$n1i[i] - df_rare$a[i]),
                         rep(1, df_rare$c[i]), rep(0, df_rare$n2i[i] - df_rare$c[i]))))
}
df_rare_logor <- prepare_ma(df_rare_ind, effect = "logOR")

## ---- include=F---------------------------------------------------------------
bg_rare_agg <- baggr(df_rare_logor, effect = "logOR")
bg_rare_ind <- baggr(df_rare_ind, effect = "logOR")

## ---- eval=F, echo=T----------------------------------------------------------
#  bg_rare_agg <- baggr(df_rare_logor, effect = "logOR")
#  bg_rare_ind <- baggr(df_rare_ind, effect = "logOR")

## -----------------------------------------------------------------------------
bgc <- baggr_compare(
  "Summary-level (Rubin model on logOR)" = bg_rare_agg,
  "Individual-level (logistic model)" = bg_rare_ind
)
bgc
plot(bgc)

## ---- include = FALSE, echo = FALSE-------------------------------------------
#let's use the data.frame we created from Yusuf et al earlier
df$study_grouping      <- c(1,1,1,0,0,0,0)
df$different_contrasts <- c(1,1,1,0,0,0,0) - .5
bg_cov1 <- baggr(df, covariates = c("study_grouping"), effect = "logarithm of odds ratio")
bg_cov2 <- baggr(df, covariates = c("different_contrasts"), effect = "logarithm of odds ratio")

## ---- echo = TRUE, eval = FALSE-----------------------------------------------
#  #let's use the data.frame we created from Yusuf et al earlier
#  df$study_grouping      <- c(1,1,1,0,0,0,0)
#  df$different_contrasts <- c(1,1,1,0,0,0,0) - .5
#  bg_cov1 <- baggr(df, covariates = c("study_grouping"), effect = "logarithm of odds ratio")
#  bg_cov2 <- baggr(df, covariates = c("different_contrasts"), effect = "logarithm of odds ratio")

## -----------------------------------------------------------------------------
baggr_compare("No covariate" = bg_model_agg, 
              "With covariates, 0-1 coding" = bg_cov1,
              "With covariates, +-1/2 coding" = bg_cov2)

## -----------------------------------------------------------------------------
bg_cov1

