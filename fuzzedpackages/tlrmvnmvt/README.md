## About the package

This package offers functionalities for computing multivariate normal and Student-t probabilities. The meanings of this package:
1. Faster implementation of the classic Genz algorithm (compared with `pmvnorm` and `pmvt` functions)
2. Able to return the results in the log2 form, which is useful when the true probability is smaller than machine precision
3. Accepts both a matrix and a geometry for building the covariance matrix
4. Implement the tile-low-rank method with block reordering, which can compute problems in tens of thousands dimensions
5. Provide the interface for users to adjust the number of Monte Carlo sample size to make a balance between accuracy and computation time specific to their applications.

## Speed improvement

For better performance, the package should be compiled with a proper optimization flag. To achieve this, you can do either of the following in the Makevars file under ~/.R directory:

1. Set CXX14FLAGS to empty, CXX14FLAGS=, which will leave the configure file to find the optimization option based on the CXX14 compiler your R uses.
2. Set the CXX14 compiler and flags by yourself, e.g. 
    CXX14 = g++
    CXX14FLAGS = -O3

The goal here is to override the default compiler options used by R. The default options are fine, just not the fastest.

## On Mac OS

After some installations on Mac OS, two issues may happen:

1. The binary gfortran compiler is not available
2. The gfortran library is not available

The explanation can be sought from: https://cran.r-project.org/bin/macosx/tools/

And the solution is:

1. Choose and install gfortran from the above-mentioned website. For issue 2, this should be sufficient
2. For issue 1, a soft link should be created:
    ln -s /usr/local/gfortran/bin/gfortran /usr/local/bin/

## On Windows

If you are encountering the issue "C++14 standard requested but CXX14 is not defined", one option is to append the Makevars.win file under the .R directory with:

CXX11FLAGS=-O3 -Wno-unused-variable -Wno-unused-function
CXX14FLAGS=-O3 -Wno-unused-variable -Wno-unused-function
CXX14 = $(BINPREF)g++ -m$(WIN) -std=c++1y

Thanks to Robert Aue for providing this solution.

Then hopefully the installation works!
