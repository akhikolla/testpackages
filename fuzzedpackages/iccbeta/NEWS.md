# iccbeta 1.2.0

## Feature

- Added automatic calculation of icc on `lmer` model objects.

## Changes

- Added improved error messaging related to presence of missing values and
  dimension mis-match.
- Added ORCiDs to `DESCRIPTION`
- Increased version dependencies.
- Enabled use of OpenMP and C++11.
- Added `CITATION` information for the _R_ package.

## Documentation

- Switched documentation to use Markdown.
- Improved contents of documentation.

## Deployment

- Enable the default TMSALab Travis-CI configuration.
- Enable code coverage checks.

# iccbeta 1.1.0

- Adds instructions on how to create `simICCdata2` to relevant helpdocs
- Switched over to using Rcpp to handle function registration instead of `src/init.c`.
- Added namespace requirement check for `lme4` and `RLRsim` to documentation.\
- Modified Travis-CI to use a later gcc

# iccbeta 1.0.1

- Added `src/init.c` to comply with R 3.4 requirements
- Added GitHub repository location information.
- Improved documentation flow.

# iccbeta 1.0.0

- Released initial version of iccbeta



