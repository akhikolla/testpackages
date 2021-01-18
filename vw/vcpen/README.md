[![CRAM_Status_Badge](http://www.r-pkg.org/badges/version/vcpen)](https://CRAN.R-project.org/package=vcpen)
[![Downloads](http://cranlogs.r-pkg.org/badges/vcpen)](https://CRAN.R-project.org/package=vcpen)
[![Total-Downloads](https://cranlogs.r-pkg.org/badges/grand-total/vcpen)](https://CRAN.R-project.org/package=vcpen)

# The `vcpen` Package
Penalized Variance Component analysis is performed by the vcpen function. The calculations do..Dan please edit...


# The `vcpen()` Function

`vcpen()` is a function that performs penalized model fitting of variance components. The variance components are pre-calculated matrices of similarity between subjects, including residual error matrix.

# The `kernel_linear()` Function

`kernel_linear()` calculates a kernel matrix of genomic similarity from SNP dosage, where the rows are subjects, and columns are the minor allele dosage (0/1/2), assuming bi-allelic SNPs.
