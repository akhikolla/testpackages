# MrSGUIDE

**M**ultiple **R**esponse **S**ubgroup Identification using GUIDE algorithm.

The current package is still under heavily development.

For subgroup identification use [GUIDE](http://pages.stat.wisc.edu/~loh/guide.html) algorithm.

Here are the links and paper for reference:

Loh, W.-Y. and Zhou, P. (2020), [The GUIDE approach to subgroup identification](http://www.stat.wisc.edu/~loh/treeprogs/guide/LZ20.pdf). In Design and Analysis of Subgroups with Biopharmaceutical Applications, N. Ting, J. C. Cappelleri, S. Ho, and D.-G. Chen (Eds.) Springer, in press.

Loh, W.-Y., Cao, L. and Zhou, P. (2019), [Subgroup identification for precision medicine: a comparative review of thirteen methods](http://www.stat.wisc.edu/~loh/treeprogs/guide/wires19.pdf), Wiley Interdisciplinary Reviews: Data Mining and Knowledge Discovery, vol. 9, 5, e1326. [DOI](http://dx.doi.org/10.1002/widm.1326)

Loh, W.-Y., Man, M. and Wang, S. (2019), [Subgroups from regression trees with adjustment for prognostic effects and po3st-selection inference](http://pages.stat.wisc.edu/~loh/treeprogs/guide/sm19.pdf), Statistics in Medicine, vol. 38, 545-557. [DOI](https://onlinelibrary.wiley.com/doi/10.1002/sim.7677)

Loh, W.-Y., Fu, H., Man, M., Champion, V. and Yu, M. (2016), [Identification of subgroups with differential treatment effects for longitudinal and multiresponse variables](http://www.stat.wisc.edu/~loh/treeprogs/guide/LFMCY16.pdf), Statistics in Medicine, vol. 35, 4837-4855. [DOI](https://onlinelibrary.wiley.com/doi/full/10.1002/sim.7020)

## Dependencies

R packages:

- Rcpp
- RcppArmardillo
- BH

## Package install

```
library(devtools)
install_github('baconzhou/MrSGUIDE')
```

### MacOS

Please refer to the following link if you have problem with install the package from source.

1. [R Compiler Tools for Rcpp on macOS](https://thecoatlessprofessor.com/programming/r-compiler-tools-for-rcpp-on-macos/) 
2. [Rcpp FAQ section 2.10](https://CRAN.R-project.org/package=Rcpp)

Here are my personal set up for macOS building `MrSGUIDE` used homebrew for libomp

```
brew update
brew install llvm libomp gcc
```

In the `~/.R/Makevars` file

```
CC=/usr/local/opt/llvm/bin/clang
CXX=/usr/local/opt/llvm/bin/clang++
CXX1X=/usr/local/opt/llvm/bin/clang++
CXX11=/usr/local/opt/llvm/bin/clang++

FLIBS=-L/usr/local/Cellar/gcc/8.2.0 -L/usr/local/Cellar/gcc/8.2.0/lib/gcc/8 -lgfortran -lquadmath -lm
```
