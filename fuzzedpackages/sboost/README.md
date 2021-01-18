# sboost
**Machine learning package used to build and test classifiers using AdaBoost on decision stumps.**

Creates classifier for binary outcomes using Adaptive Boosting (AdaBoost)
on decision stumps with a fast C++ implementation. Feature vectors may be a 
combination of continuous (numeric) and categorical (string, factor) elements. 
Methods for classifier assessment, predictions, and cross-validation also included. 
The advantage of this type of classifier is that it is non-linear but it is more 
interpretable than random forests, neural-nets, and other non-linear classifiers. 

See [jadonwagstaff.github.io/sboost](https://jadonwagstaff.github.io/sboost.html) for a description
of how the classifier functions, and what makes this classifier more interpretable than others.

For original paper describing AdaBoost see:

Freund, Y., Schapire, R.E.: A decision-theoretic generalization of on-line learning and an application to boosting. Journal of Computer and System Sciences 55(1), 119-139 (1997)

## Installation
Install this package from the CRAN repository.

```
install.packages("sboost")
```

Alternatively, use devtools to install the development version of this package.

To install devtools on R run:

```
install.packages("devtools")
```

After devtools is installed, to install the sboost package on R run:

```
devtools::install_github("jadonwagstaff/sboost")
```

## Functions

*sboost* - Main machine learning algorithm, uses categorical or continuous features to build a classifier that predicts a binary outcome.  Run ```?sboost::sboost``` to see documentation in R.

*validate* - Uses k-fold cross validation on a training set to validate the classifier.

*assess* - Shows performance of a classifier on a set of feature vectors and outcomes.

*predict* - Outputs predictions of a classifier on a set of feature vectors.

## Author
Jadon Wagstaff

## Licence
MIT
