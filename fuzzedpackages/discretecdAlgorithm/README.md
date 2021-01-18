discretecdAlgorithm
===================

[![Build Status](https://travis-ci.org/gujyjean/discretecdAlgorithm.svg?branch=master)](https://travis-ci.org/gujyjean/discretecdAlgorithm)

An algorithm to learn structure of discrete Bayesian network, this package can deal with observational data, interventional data, or a mixture of both.

algorithm
---------

-   `cd.run` is the main function to run coordinate descent algorithm. With the `adaptive` option, users may choose to use regular group lasso penalty, or adaptive group lasso penalty.
-   `max_lambda` is a function to calculate the maximum value of lambda that will penalized all edges to zero.
