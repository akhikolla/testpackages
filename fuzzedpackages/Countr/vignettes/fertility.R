library(Countr)

library(dplyr) 
library(xtable)

data(fertility)

dataCaption <- "First few rows of fertility data."
print(xtable(head(fertility), caption = dataCaption, label = "tbl:data"), 
      rotate.colnames = TRUE)

freqtable <- 
    count_table(count = fertility$children, breaks = 0:9, formatChar = TRUE)
print(xtable(freqtable, caption = "Fertility data: Frequency distribution of column \\texttt{children}.",
             label = "tbl:freq"))

nam_fac <- sapply(fertility, function(x) !is.numeric(x))
fert_factor <- summary(fertility[ , nam_fac])
fert_num <- t(sapply(fertility[ , !nam_fac], summary))

print(xtable(fert_factor, caption = "Summary of the factor variables", label = "tbl:frecfac"))

print(xtable(fert_num, caption = "Summary of the numeric explanatory variables", 
                       label = "tbl:frecnum"))

form <- children ~ german + years_school + voc_train + university + religion +
    year_birth + rural + age_marriage
pois <- glm(formula = form, data = fertility, family = poisson())
summary(pois)

form <- children ~ german + years_school + voc_train + university + religion +
    year_birth + rural + age_marriage
wei <- renewalCount(formula = form, data = fertility, dist = "weibull",
                    weiMethod = "conv_dePril",
                    control = renewal.control(trace = 0, method = "nlminb")
                    )
summary(wei)

residuals_plot(wei, type = "pearson")

par(mfrow = c(1, 2))
res_wei <- residuals(wei, type = "pearson")
qqnorm(res_wei, ylim = range(res_wei), main = "Weibull Renewal Model")
qqline(res_wei, ylim = range(res_wei))
grid()
pois <- glm(formula = form, data = fertility, family = poisson())
res_pois <- residuals(pois, type = "pearson")
qqnorm(res_pois, ylim = range(res_wei), main = "GLM Poisson")
qqline(res_pois, ylim = range(res_wei))
grid()

form <- children ~ german + years_school + voc_train + university + religion +
                   year_birth + rural + age_marriage
wei <- renewalCount(formula = form, data = fertility, dist = "weibull",
                    control = renewal.control(trace = 0, method = "nlminb")
                    )
wei_summary <- summary(wei)$coef

t_shape <- wei_summary[rownames(wei_summary) == "shape_",, drop = FALSE]
print(xtable(t_shape), floating = FALSE)

library(lmtest)
lrtest(pois, wei)

form <- children ~ german + years_school + voc_train + university + religion +
    year_birth + rural + age_marriage
pois <- glm(formula = form, data = fertility, family = poisson())
wei <- renewalCount(formula = form, data = fertility, dist = "weibull",
                    control = renewal.control(trace = 0, method = "nlminb")
                    )
tab <- compareToGLM(poisson_model = pois, breaks = c(0:5, 7, 9), weibull = wei)
tab_tbl <- xtable(tab, caption = "Comparison between some models.")
print(tab_tbl)

colSums(dplyr::select(tab, contains("_pearson")))

frequency_plot(tab$Counts, tab$Actual,
               dplyr::select(tab, contains("_predicted"))
)

  form <- children ~ german + years_school + voc_train + university + religion +
      year_birth + rural + age_marriage
  wei <- renewalCount(formula = form, data = fertility, dist = "weibull",
                      control = renewal.control(trace = 0, method = "nlminb")
                      )
  gof <- chiSq_gof(wei, breaks = c(0:5, 7, 9))
  gof
