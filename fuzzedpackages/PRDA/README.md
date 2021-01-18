
<!-- README.md is generated from README.Rmd. Please edit that file -->

# PRDA: Prospective and Retrospective Design Analysis

<!-- badges: start -->

[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/ClaudioZandonella/PRDA?branch=master&svg=true)](https://ci.appveyor.com/project/ClaudioZandonella/PRDA/branch/master)
[![Travis build
status](https://travis-ci.org/ClaudioZandonella/PRDA.svg?branch=master)](https://travis-ci.org/ClaudioZandonella/PRDA)
[![Codecov test
coverage](https://codecov.io/gh/ClaudioZandonella/PRDA/branch/master/graph/badge.svg)](https://codecov.io/gh/ClaudioZandonella/PRDA/branch/master)
[![DOI](https://zenodo.org/badge/212573857.svg)](https://zenodo.org/badge/latestdoi/212573857)

<hr>

<!-- badges: end -->

{PRDA} allows performing a prospective or retrospective design analysis
to evaluate inferential risks (i.e., power, Type M error, and Type S
error) in a study considering Pearson’s correlation between two
variables or mean comparisons (one-sample, paired, two-sample, and
Welch’s *t*-test).

For an introduction to design analysis and a general overview of the
package see `vignette("PRDA")`. Examples for retrospective design
analysis and prospective design analysis are provided in
`vignette("retrospective")` and `vignette("prospective")` respectively.

All the documentation is available at
<https://claudiozandonella.github.io/PRDA/>.

## Installation

<!-- You can install the released version of PRDA from [CRAN](https://CRAN.R-project.org) with: -->

<!-- ``` r -->

<!-- install.packages("PRDA") -->

<!-- ``` -->

<!-- And the development version from [GitHub](https://github.com/) with: -->

You can install the development version from
[GitHub](https://github.com/ClaudioZandonella/PRDA/tree/master) with:

``` r
# install.packages("devtools")
devtools::install_github("ClaudioZandonella/PRDA",
                         build_vignettes = TRUE)
```

## The Package

{PRDA} package can be used for Pearson’s correlation between two
variables or mean comparisons (i.e., one-sample, paired, two-sample, and
Welch’s t-test) considering an hypothetical value of *ρ* or Cohen’s *d*
respectively. See `vignette("retrospective")` and
`vignette("prospective")` to know how to set function arguments for the
different effect types.

### Functions

In {PRDA} there are two main functions `retrospective()` and
`prospective()`.

#### • `retrospective()`

Given the hypothetical population effect size and the study sample size,
the function `retrospective()` performs a retrospective design analysis.
According to the defined alternative hypothesis and the significance
level, the inferential risks (i.e., Power level, Type M error, and Type
S error) are computed together with the critical effect value (i.e., the
minimum absolute effect size value that would result significant).

Consider a study that evaluated the correlation between two variables
with a sample of 30 subjects. Suppose that according to the literature
the hypothesized effect is *ρ* = .25. To evaluate the inferential risks
related to the study we use the function `retrospective()`.

``` r
set.seed(2020) # set seed to make results reproducible

retrospective(effect_size = .25, sample_n1 = 30, 
              test_method = "pearson")
#> 
#>  Design Analysis
#> 
#> Hypothesized effect:  rho = 0.25 
#> 
#> Study characteristics:
#>    test_method   sample_n1   sample_n2   alternative   sig_level   df
#>    pearson       30          NULL        two_sided     0.05        28
#> 
#> Inferential risks:
#>    power   typeM   typeS
#>    0.27    1.826   0.003
#> 
#> Critical value(s): rho  =  ± 0.361
```

In this case, the statistical power is almost 30% and the associated
Type M error and Type S error are respectively around 1.80 and 0.003.
That means, statistical significant results are on average an
overestimation of 80% of the hypothesized population effect and there is
a .3% probability of obtaining a statistically significant result in the
opposite direction.

To know more about function arguments and further examples see the
function documentation `?retrospective` and `vignette("retrospective")`.

#### • `prospective()`

Given the hypothetical population effect size and the required power
level, the function `prospective()` performs a prospective design
analysis. According to the defined alternative hypothesis and the
significance level, the required sample size is computed together with
the associated Type M error, Type S error, and the critical effect value
(i.e., the minimum absolute effect size value that would result
significant).

Consider a study that will evaluate the correlation between two
variables. Knowing from the literature that we expect an effect size of
*ρ* = .25, the function `prospective()` can be used to compute the
required sample size to obtain a power of 80%.

``` r
prospective(effect_size = .25, power = .80, test_method = "pearson",
            display_message = FALSE)
#> 
#>  Design Analysis
#> 
#> Hypothesized effect:  rho = 0.25 
#> 
#> Study characteristics:
#>    test_method   sample_n1   sample_n2   alternative   sig_level   df 
#>    pearson       122         NULL        two_sided     0.05        120
#> 
#> Inferential risks:
#>    power   typeM   typeS
#>    0.797   1.119   0    
#> 
#> Critical value(s): rho  =  ± 0.178
```

The required sample size is \(n=122\), the associated Type M error is
around 1.10 and the Type S error is approximately 0.

To know more about function arguments and further examples see the
function documentation `?prospective` and `vignette("prospective")`.

### Hypothetical effect size

The hypothetical population effect size can be defined as a single value
according to previous results in the literature or experts indications.
Alternatively, {PRDA} allows users to specify a distribution of
plausible values to account for their uncertainty about the hypothetical
population effect size. To know how to specify the hypothetical effect
size according to a distribution and an example of application see
`vignette("retrospective")`.

## Contributing to PRDA

The PRDA package is still in the early stages of its life. Thus, surely
there are many bugs to fix and features to propose. Anyone is welcome to
contribute to the PRDA package.

Please note that this project is released under a [Contributor Code of
Conduct](https://www.contributor-covenant.org/). By contributing to this
project, you agree to abide by its terms.

#### Bugs and New Features

To propose a new feature or to report a bug, please open an issue on
[GitHub](https://github.com/ClaudioZandonella/PRDA/issues). See
[Community
guidelines](https://github.com/ClaudioZandonella/PRDA/blob/master/CONTRIBUTING.md).

#### Future Plans

  - Improve compute time by parallelizing the code
  - Implement design analysis in the case of linear regression models

## Citation

To cite {PRDA} in publications use:

Zandonella Callegher, C., Pastore, M., Andreella, A., Vesely, A.,
Toffalini, E., Bertoldo, G., & Altoè G. (2020). PRDA: Prospective and
Retrospective Design Analysis (Version 1.0.0). Zenodo.
<https://doi.org/10.5281/zenodo.4044214>

A BibTeX entry for LaTeX users is

    @Misc{,
        author       = {Zandonella Callegher, Claudio and Pastore, Massimiliano and Andreella, Angela and 
                        Vesely, Anna and Toffalini, Enrico and Bertoldo, Giulia and Altoè, Gianmarco},
        title        = {PRDA: Prospective and Retrospective Design 
                       Analysis},
        year         = 2020,
        publisher    = {Zenodo},
        version      = {1.0.0},
        doi          = {10.5281/zenodo.4044214},
        url          = {https://doi.org/10.5281/zenodo.4044214}
      }

## References

<div id="refs" class="references">

<div id="ref-altoeEnhancingStatisticalInference2020">

Altoè, Gianmarco, Giulia Bertoldo, Claudio Zandonella Callegher, Enrico
Toffalini, Antonio Calcagnì, Livio Finos, and Massimiliano Pastore.
2020. “Enhancing Statistical Inference in Psychological Research via
Prospective and Retrospective Design Analysis.” *Frontiers in
Psychology* 10. <https://doi.org/10.3389/fpsyg.2019.02893>.

</div>

<div id="ref-bertoldoDesigningStudiesEvaluating2020">

Bertoldo, Giulia, Claudio Zandonella Callegher, and Gianmarco Altoè.
2020. “Designing Studies and Evaluating Research Results: Type M and
Type S Errors for Pearson Correlation Coefficient.” Preprint. PsyArXiv.
<https://doi.org/10.31234/osf.io/q9f86>.

</div>

<div id="ref-gelmanPowerCalculationsAssessing2014">

Gelman, Andrew, and John Carlin. 2014. “Beyond Power Calculations:
Assessing Type S (Sign) and Type M (Magnitude) Errors.” *Perspectives on
Psychological Science* 9 (6): 641–51.
<https://doi.org/10.1177/1745691614551642>.

</div>

</div>
