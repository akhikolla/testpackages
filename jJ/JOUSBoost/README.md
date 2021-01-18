<!-- README.md is generated from README.Rmd. Please edit that file -->
JOUSBoost
=========

The JOUSBoost package implements under/oversampling with jittering for probability estimation. Its intent is to be used to improve probability estimates that come from boosting algorithms (such as AdaBoost), but is modular enough to be used with virtually any classification algorithm from machine learning. See Mease (2007) for more information.

You can install:

-   the latest released version from CRAN with

``` r
install.packages("JOUSBoost")
```

-   the latest development version from github with

``` r
devtools::install_github("molson2/JOUSBoost")
```

Illustration
------------

The following example gives a useage case for JOUSBoost. This example illustrates the improvement in probability estimates on gets from applying the JOUS procedure to AdaBoost on a simulated data set. First, we'll train AdaBoost applied to depth three decision trees, and then we'll get the estimated probabilities.

``` r
# Generate data from Friedman model #
library(JOUSBoost)
set.seed(111)
dat = friedman_data(n = 1000, gamma = 0.5)
train_index = sample(1:1000, 800)

# Train AdaBoost classifier using depth 3 decision tree
ada = adaboost(dat$X[train_index,], dat$y[train_index], tree_depth = 3, n_rounds = 400)

# get probability estimate on test data
phat_ada = predict(ada, dat$X[-train_index, ], type="prob")
```

Next, we'll compute probabilities by using the JOUS procedure.

``` r
# Apply jous to adaboost classifier
class_func = function(X, y) adaboost(X, y, tree_depth = 3, n_rounds = 400)
pred_func = function(fit_obj, X_test) predict(fit_obj, X_test)

jous_fit = jous(dat$X[train_index,], dat$y[train_index], class_func,
                pred_func, type="under", delta=10, keep_models=TRUE)

# get probability estimate on test data
phat_jous = predict(jous_fit, dat$X[-train_index, ], type="prob")
```

Finally, we can see the benefit of using JOUSBoost!

``` r
# compare MSE of probability estimates
p_true = dat$p[-train_index]
mean((p_true - phat_jous)^2)
#> [1] 0.05455999
mean((p_true - phat_ada)^2)
#> [1] 0.1277416
```

Mease, D., Wyner, A., and Buja, A. 2007. “Costweighted Boosting with Jittering and over/under-Sampling. JOUS-Boost.” Journal of Machine Learning Research 8 409-439.
