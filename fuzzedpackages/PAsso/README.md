# PAsso: an R Package for Assessing Partial Association between Ordinal Variables

<!-- badges: start -->

[![CRAN checks](https://cranchecks.info/badges/summary/PAsso)](https://cran.r-project.org/web/checks/check_results_PAsso.html)
[![](http://cranlogs.r-pkg.org/badges/grand-total/PAsso?color=blue)](https://cran.r-project.org/package=PAsso)
[![](http://cranlogs.r-pkg.org/badges/last-month/PAsso?color=green)](https://cran.r-project.org/package=PAsso)
[![](http://cranlogs.r-pkg.org/badges/last-week/PAsso?color=yellow)](https://cran.r-project.org/package=PAsso)
[![](https://travis-ci.com/XiaoruiZhu/PAsso.svg?branch=master)](https://travis-ci.com/XiaoruiZhu/PAsso.svg)

<!-- badges: end -->

Overview
--------

An R package of a unified framework for assessing **P**arrtial **Asso**ciation between Ordinal variables. It includes quantification, visualization, and hypothesis testing. All the products are based on the paper by Dungang Liu, 
Shaobo Li, Yan Yu and Irini Moustaki (2020) and the approach described in [Dungang and Zhang
(2017)](http://www.tandfonline.com/doi/abs/10.1080/01621459.2017.1292915?journalCode=uasa20).

## Installation

The `PAsso` package is currently not available on [PAsso CRAN]() and wait for future updates.

### Install the development version from GitHub

``` r
if (!requireNamespace("devtools")) install.packages("devtools")
devtools::install_github("XiaoruiZhu/PAsso")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(PAsso)

PAsso_1 <- PAsso(responses = c("PreVote.num", "PID"),
                 adjustments = c("income.num", "age", "edu.year"),
                 data = ANES2016,
                 uni.model = "probit",
                 method = c("kendall")
                )
                
# Print the partial association matrix only
print(PAsso_1, 5)

# Provide partial association matrix, marginal association matrix, and summary of models' coefficients
summary(PAsso_1, 4)

# Plot partial association regression plot: residuals
plot(PAsso_1)

# Retrieve residuals
test_resids <- residuals(PAsso_1, draw = 1)
head(test_resids)
dim(test_resids)

# Association analysis between three ordinal responses ----------

PAsso_2 <- PAsso(responses = c("PreVote.num", "PID", "selfLR"),
                adjustments = c("income.num", "age", "edu.year"),
                data = ANES2016,
                uni.model <- "probit",
                method = c("kendall"),
                resids.type = "surrogate")
                
# Compare marginal correlation and partial correlation.
summary(PAsso_2, digits=4)
plot(PAsso_2)

# test function: Conduct inference based on object of "PAsso.test" class ----------------------------
library(progress); #library(doParallel)

system.time(Pcor_SR_test1 <- test(object = PAsso_2, boot_SE=100, H0=0, parallel=F))
print(Pcor_SR_test1, digits=6)
print(PAsso_2, 6)

# diagnostic.plot function -----------------------------------------------------
check_qq <- diagnostic.plot(object = PAsso_2, output = "qq")

check_fitted <- diagnostic.plot(object = PAsso_2, output = "fitted")

check_covar <- diagnostic.plot(object = PAsso_2, output = "covariate")

# general association measure and 3-D plot for VOTE and PID ------------------
library("copula")
library("plotly")

# Draw all pairs
testPlots <- plot3D(PAsso_1)
testPlots$`PreVote.num v.s. PID`

# Draw just one pair
testPlots2 <- plot3D(object = PAsso_2, y1 = "selfLR", y2 = "PID")
testPlots2

# "PAsso" advanced using of the function: Input a few models directly ------------------------------

fit_vote <- glm(PreVote.num ~ income.num + age + edu.year, data = nes2016,
               family = binomial(link = "probit"))
summary(fit_vote)

fit_PID <- polr(as.factor(PID) ~ income.num + age + edu.year, data = nes2016,
               method = "probit", Hess = TRUE)

summary(fit_PID)

system.time(PAsso_adv1 <- PAsso(fitted.models=list(fit_vote, fit_PID),
                                association = c("partial"),
                                method = c("kendall"),
                                resids.type = "surrogate")
)

# Partial association coefficients (Parts of Table 7 in paper)
print(PAsso_adv1, digits = 3)
summary(PAsso_adv1, digits = 3)
```

References
----------

Dungang Liu, Shaobo Li, Yan Yu and Irini Moustaki (2019). "Assessing partial association between ordinal variables: quantification, visualization, and hypothesis testing", accepted by JASA.

Dungang Liu & Heping Zhang (2018). "Residuals and Diagnostics for Ordinal Regression Models: A Surrogate Approach", Journal of the American Statistical Association, 113:522, 845-854, DOI: 10.1080/01621459.2017.1292915, URL
<http://www.tandfonline.com/doi/abs/10.1080/01621459.2017.1292915?journalCode=uasa20>

Greenwell, B.M., McCarthy, A.J., Boehmke, B.C. & Dungang, L. (2018)
"Residuals and diagnostics for binary and ordinal regression models: An
introduction to the sure package." The R Journal. URL
<https://journal.r-project.org/archive/2018/RJ-2018-004/index.html>

