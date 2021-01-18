## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(knitr)
library(ggplot2)
library(gridExtra)

## ----setup--------------------------------------------------------------------
# install.packages('WRI')
library(WRI)

## -----------------------------------------------------------------------------
data(strokeCTdensity)
?strokeCTdensity

predictor = strokeCTdensity$predictors
dSup = strokeCTdensity$densitySupport
densityCurves = strokeCTdensity$densityCurve

## -----------------------------------------------------------------------------
equal_t = den2Q_qd(densityCurves, dSup, t_vec = seq(0, 1, length.out = 120))
nonequal_t = den2Q_qd(densityCurves, dSup, t_vec = unique(c(seq(0, 0.05, 0.001), seq(0.05, 0.95, 0.05), seq(0.95, 1, 0.001))))

## ---- fig.align='center', fig.height=3, fig.width= 7, echo=F------------------
i = 1
df_t = data.frame(densityFun = c(densityCurves[i, ], equal_t$fobs[i, ], nonequal_t$fobs[i, ]),
                  dSup = c(dSup,equal_t$Qobs[i, ], nonequal_t$Qobs[i, ]),
                  density = c(rep('original f(x)', 101), rep('equal t: 1/q(t)', 120), rep('nonequal t: 1/q(t)', 119)))
   
ggplot(data = df_t, aes(x = dSup, y = densityFun, color = density)) +
         geom_line() +
         theme_linedraw()+
         geom_point(size = 0.3) +
         ggtitle('Comparison of the effect of equally and nonequally spaced t grid vector')


## -----------------------------------------------------------------------------
res = wass_regress(rightside_formula = ~., Xfit_df = predictor, Ytype = 'density', Ymat = densityCurves, Sup = dSup)

## -----------------------------------------------------------------------------
summary(res)

## -----------------------------------------------------------------------------
globalF_res =  globalFtest(res, alpha = 0.05, permutation = TRUE, numPermu = 200)
kable(globalF_res$summary_df, digits = 3)
sprintf('The wasserstein F-statistic is %.3f on %.3f degrees of freedom', globalF_res$wasserstein.F_stat, globalF_res$chisq_df)

## ---- eval=T------------------------------------------------------------------
# the reduced model only has four radiological variables
reduced_res = wass_regress(~ log_b_vol + b_shapInd + midline_shift + B_TimeCT, Xfit_df = predictor, Ymat = densityCurves, Ytype = 'density', Sup = dSup)
full_res = wass_regress(rightside_formula = ~., Xfit_df = predictor, Ymat = densityCurves, Ytype = 'density', Sup = dSup)

partialFtable = partialFtest(reduced_res, full_res, alpha = 0.05)
kable(partialFtable, digits = 3)

## ---- fig.width=7.2, fig.height=2.5-------------------------------------------
xpred = colMeans(predictor)
confidence_Band = confidenceBands(res, Xpred_df = xpred, type = 'both')

## -----------------------------------------------------------------------------
mean_Mode <- function(vec) {
  return(ifelse(length(unique(vec)) < 3, modeest::mfv(vec), mean(vec)))
}
mean_mode_vec = apply(predictor, 2, mean_Mode)
predictorDF = rbind(mean_mode_vec, mean_mode_vec)
predictorDF[ , 1] = quantile(predictor$log_b_vol, probs = c(1/4, 3/4))

## ---- fig.width=7.2, fig.height=3---------------------------------------------
res_cb = confidenceBands(res, predictorDF, level = 0.95, delta = 0.01, type = 'both', figure = F)
m = ncol(res_cb$quan_list$Q_lx)
na.mat = matrix(NA, nrow = 2, ncol = m - ncol(res_cb$den_list$f_lx))

cb_plot_df = with(res_cb, data.frame(
  fun = rep(c('quantile function', 'density function'), each = 2*m), 
  Q1Q3 = rep(rep(c('Q1', 'Q3'), each = m), 2),
  value_m = c(as.vector(t(quan_list$Qpred)), as.vector(t(cbind(cdf_list$fpred)))),
  value_u = c(as.vector(t(quan_list$Q_ux)), as.vector(t(cbind(den_list$f_ux, na.mat)))),
  value_l = c(as.vector(t(quan_list$Q_lx)), as.vector(t(cbind(den_list$f_lx, na.mat)))),
  support_full = c(rep(quan_list$t, 2), as.vector(t(cbind(cdf_list$Fsup)))),
  support_short = c(rep(quan_list$t, 2), as.vector(t(cbind(den_list$Qpred, na.mat))))
  ))
  
ggplot(data = cb_plot_df, aes(color = Q1Q3)) +
  theme_linedraw()+
  geom_line(aes(x = support_full, y = value_m)) +
  geom_ribbon(aes(x = support_short, ymin = value_l, ymax = value_u, fill = Q1Q3), alpha = 0.25) +
  facet_wrap( ~ fun, scales = "free_y") +
  ylab('Confidence band') +
  xlab('Support')

