# cIRT 1.3.1

## Changes

- Added `LazyData: true` to `DESCRIPTION` to match with how the data 
  documentation is called. 
  ([#3](https://github.com/tmsalab/cIRT/pull/3))

## Documentation

- Added a `pkgdown` website that deploys to <https://tmsalab.github.io/cIRT>.
  ([#3](https://github.com/tmsalab/cIRT/pull/3))

## Deployment

- Switched from Travis-CI to GitHub Actions for R.
  ([#5](https://github.com/tmsalab/cIRT/pull/5))

# cIRT 1.3.0

## Changes

- Updated package dependencies
- Enabled C++11 and OpenMP.
- Switched to allowing Rcpp to handle native registration

## Bugfix

- Addressed issues in the choice generation procedure.

## Documentation

- Improved in-line documentation.
- Added Authors' ORCIDs to `DESCRIPTION`.

## Testing

- Enabled TMSA Lab's configuration for Travis-CI.

# cIRT 1.2.1

- Added `src/init.c` for R 3.4 compatibility
- Added GitHub project page link

# cIRT 1.2.0

- Added two vignettes that cover the model estimation and simulation results in the package. 
- Added a `NEWS.md` file to track changes to the package.

# cIRT 1.1.0

- Adds two columns to choice matrix: `hard_q_id` and `easy_q_id`.

# cIRT 1.0.0

## Modeling Framework
- Implementation of the hierarchical framework described in "A Hierarchical Model for Accuracy and Choice on Standardized Tests"
- Specifically, a choice inclusive Probit HLM and a Two Parameter Ogive Model.

## C++ Functions
- Random Number Generation for the following distributions: Wishart, Inverse Wishart, and Multivariate Normal
- Matrix Centering
- Direct Sum calculation

## Data
- Student Performance on Revised Purdue Spatial Visualization Test (Revised PSVT:R) by Yoon, 2011 in `trial_matrix`
- The choices students made among items presented to them in `choice_matrix`
- The end payout results for students based on their choices in `payout_matrix`
- One additional data set exists containing the student's sex response in `survey_data`
