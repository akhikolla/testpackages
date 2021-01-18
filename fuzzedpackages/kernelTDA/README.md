
<!-- README.md is generated from README.Rmd. Please edit that file -->

# kernelTDA

<!-- badges: start -->

<!-- badges: end -->

This `R`-package provides an implementation of the most famous kernels
found in the framework of Topological Data Analysis (TDA), more
specifically:

  - [Persistence Scale Space
    Kernel](http://openaccess.thecvf.com/content_cvpr_2015/papers/Reininghaus_A_Stable_Multi-Scale_2015_CVPR_paper.pdf)
  - [Sliced Wasserstein
    Kernel](https://dl.acm.org/citation.cfm?id=3305450)
  - [Persistence Fisher
    Kernel](http://papers.nips.cc/paper/8205-persistence-fisher-kernel-a-riemannian-manifold-kernel-for-persistence-diagrams.pdf)
  - [Geodesic Wasserstein Kernel(s)](https://arxiv.org/abs/1709.07100)
  - [Persistence
    Images](http://www.jmlr.org/papers/volume18/16-337/16-337.pdf)


Here you can also find an `R` interface to the C++ library
[HERA](https://bitbucket.org/grey_narn/hera/src/master/), which contains
an efficient implementation of the L_p q-Wasserstein distance
between persistence diagrams.


Finally, this package contains a solver for kernelized Support Vector Machine problems with indefinite kernels, based on the algorithm proposed by [Loosli et al.](https://hal.archives-ouvertes.fr/hal-01593553/document). The implementation is largely based on the C++ library [LIBSVM](https://www.csie.ntu.edu.tw/~cjlin/libsvm/), and on its R interface in the package [e1071](https://CRAN.R-project.org/package=e1071).


This package is now on CRAN, you can install it with:

``` r
install.packages("kernelTDA")
```
