afCEC
===

Active function cross-entropy clustering partitions the n-dimensional data into the clusters by finding the parameters of the mixed generalized multivariate normal distribution, that optimally approximates the scattering of the data in the n-dimensional space, whose density. The above-mentioned generalization is performed by introducing so called "f-adapted Gaussian densities" (i.e. the ordinary Gaussian densities adapted by the "active function"). Additionally, the active function cross-entropy clustering performs the automatic reduction of the unnecessary clusters. For more information please refer to P. Spurek, J. Tabor, K.Byrski, "Active function Cross-Entropy Clustering" (2017) .
The afCEC package is a part of CRAN repository and it can be installed by the following command:

```R
install.packages("afCEC")
library("afCEC")
```

The basic usage comes down to the function ` afCEC ` with two required arguments: input data (`points`) and the initial number of centers (`maxClusters `):

```R
afCEC (points= , maxClusters= )
```
Below, a simple session with **R** is presented, where the component
(waiting) of the Old Faithful dataset is split into two clusters:

```R
library(afCEC)
data(fire)
plot(fire, asp=1, pch=20)

result <- afCEC(fire, 5,  numberOfStarts=10);
print(result)
plot(result)
```

As the main result, afCEC returns data cluster membership `cec$cluster`. The following parameters of 
clusters can be obtained as well:

- means (`result$means`)
- covariances (`result$covariances`)
- cardinalities (`result$cardinalities`)

