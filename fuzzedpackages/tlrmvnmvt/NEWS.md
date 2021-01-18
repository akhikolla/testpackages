# tlrmvnmvt 1.0.0

The initial release of the package.

For a complete description of the package functions, please refer to the help file. To install a highly optimized version of this package, please refer to the README.md file. 

# tlrmvnmvt 1.0.1

This version fixed four bugs:
	1. The memory allocation size error in functions `mvn_internal`, `mvn_internal2`, `mvt_internal`, and `mvt_internal2`
	2. Change the position of the memory release command
	3. Correct the functioning of the mean parameter. It was assumed to be "shifted" but should be "Kshirsagar"; Refer to `mvtnorm::pmvt` function's documentation for the meaning of the two names
	4. Implement the internal Matern function so that it produces the same result as the `geoR::matern` function

# tlrmvnmvt 1.1.0

This version makes changes to the R interfaces of the package
	1. The previous `pmvn.genz` and `pmvn.tlr` functions are combined into the new `pmvn` function
	2. The previous `pmvt.genz` and `pmvt.tlr` functions are combined into the new `pmvt` function
	3. The interfaces of the `pmvn` and the `pmvt` functions resembles those of the `pmvnorm` and the `pmvt` functions from the `mvtnorm` package at the best effort. Differences include:
		1. Input correlation matrices are treated the same as input covariance matrices
		2. Users can specify whether the results should be returned in terms of logarithm
		3. The sample size, block size (for tile-low-rank only), and the truncation level can be specified through the `algorithm` parameter
		4. Covariance matrices can be constructed by the package with the optional parameters but these parameters are not specified to align the function interfaces with those from the `mvtnorm` package