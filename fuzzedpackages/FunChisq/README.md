The 'FunChisq' R package
===============================

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/FunChisq)](https://cran.r-project.org/package=FunChisq)
[![CRAN_latest_release_date](https://www.r-pkg.org/badges/last-release/FunChisq)](https://cran.r-project.org/package=FunChisq)
[![metacran downloads](https://cranlogs.r-pkg.org/badges/FunChisq)](https://cran.r-project.org/package=FunChisq)
[![metacran downloads](https://cranlogs.r-pkg.org/badges/grand-total/FunChisq)](https://cran.r-project.org/package=FunChisq)
[![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)



### Overview

The package provides statistical hypothesis testing methods for inferring model-free functional dependency. Functional test statistics are asymmetric and functionally optimal, unique from other related statistics. The test significance is based on either asymptotic chi-squared or exact distributions.

The tests include asymptotic *functional chi-squared tests* (Zhang & Song, 2013) [<arXiv:1311.2707>](https://arxiv.org/pdf/1311.2707v3.pdf) and an *exact functional test* (Zhong & Song, 2019) [<10.1109/TCBB.2018.2809743>](https://doi.org/10.1109/TCBB.2018.2809743). The *normalized* functional chi-squared test was used by Best Performer NMSUSongLab in HPN-DREAM (DREAM8) Breast Cancer Network Inference Challenges (Hill et al., 2016) [<10.1038/nmeth.3773>](https://doi.org/10.1038/nmeth.3773). 

To measure the effect size, one can use the asymmetric *function index* (Zhong & Song, 2019) [<10.1186/s12920-019-0565-9>](https://doi.org/10.1186/s12920-019-0565-9) (Kumar et al., 2018) [<10.1109/BIBM.2018.8621502>](https://doi.org/10.1109/BIBM.2018.8621502). Its value is minimized to 0 by perfectly independent patterns and maximized to 1 by perfect non-constant functions.

A simulator (Sharma et al., 2017) [<10.32614/RJ-2017-053>](https://doi.org/10.32614/RJ-2017-053) can generate functional, non-functional, and independent patterns as contingency tables. The simulator provides options to control row and column marginal distributions and the noise level.

### When to use the package

Tests in this package can be used to reveal evidence for causality based on the causality-by-functionality principle. They target model-free inference without assuming a parametric model. For continuous data, these tests offer an advantage over regression analysis when a parametric functional form cannot be assumed. Data can be first discretized, e.g., by R packages ['Ckmeans.1d.dp'](https://cran.r-project.org/package=Ckmeans.1d.dp) or ['GridOnClusters'](https://cran.r-project.org/package=GridOnClusters). For categorical data, they provide a novel means to assess directional dependency not possible with symmetrical Pearson's chi-squared or Fisher's exact tests. They are a better alternative to conditional entropy in many aspects.

### To download and install the package

```{r}
install.packages("FunChisq")
```
