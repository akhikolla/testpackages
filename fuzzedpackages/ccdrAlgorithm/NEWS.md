# ccdrAlgorithm 0.0.5

## Notes

* Minor patch to fix failing CRAN check

# ccdrAlgorithm 0.0.4

## Features

* `ccdr.run` now supports black lists and white lists (sparsebn #6)
* Users can now use prior knowledge of variances via the new `sigmas` argument to `ccdr.run` (sparsebn #12)

## Bug fixes

* `ivn.rand = FALSE` not handled correctly (#4)

# ccdrAlgorithm 0.0.3

## Features

* New `generate_mvn_data()` method to generate multivariate normal data from a DAG.

## Notes

* Added warning when dataset contains more than 10,000 columns: This requires building from source. The CCDr algorithm has been safely tested on datasets with up to 8,000 variables.
* By default, `ccdr.run()` includes the node names in the `sparsebnPath` output.

# ccdrAlgorithm 0.0.2

## Features

* `ccdr.run()` is now compatible with interventional data

# ccdrAlgorithm 0.0.1

* Initial stable release

