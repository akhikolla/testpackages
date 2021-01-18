SCORNET
================

<!-- badges: start -->

[![R build
status](https://github.com/celehs/SCORNET/workflows/R-CMD-check/badge.svg)](https://github.com/celehs/SCORNET/actions)
<!-- badges: end -->

## Overview

The Semi-supervised Calibration of Risk with Noisy Event Times (SCORNET)
is a consistent, semi-supervised, non-parametric survival curve
estimator optimized for efficient use of Electronic Health Record (EHR)
data with a limited number of current status labels. Derived from van
der Laan and Robins’ Inverse Probability of Censoring Weighted (IPCW)
estimator, it achieves locally efficient survival curve estimation using
current status labels – binary indicators of phenotype status at
censoring time – rather than more expensive event time labels. SCORNET
boosts efficiency over IPCW in the typical EHR setting by (1) utilizing
unlabeled patients in a semi-supervised fashion, and (2) leveraging
information-dense engineered EHR features to maximize imputation
precision in the unlabeled set.

See Ahuja et al. (2020) for details.

Click [HERE](https://celehs.github.io/SCORNET/) to view a simulated example.

## References

  - Ahuja Y, Liang L, Huang S, Cai T (2020). Semi-supervised Calibration
    of Risk with Noisy Event Times (SCORNET) Using Electronic Health
    Record Data. BioArxiv.

  - Mark J. van der Laan & James M. Robins (1998) Locally Efficient
    Estimation with Current Status Data and Time-Dependent Covariates,
    Journal of the American Statistical Association, 93:442, 693-701,
    DOI: 10.1080/01621459.1998.10473721
