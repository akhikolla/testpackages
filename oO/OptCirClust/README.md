The 'OptCirClust' R package
=================================

### Overview

The package provides fast, optimal, and reproducible clustering for circular, periodic, or framed data. The performance is guaranteed by embedding divide-and-conquer, dynamic programming, and another divide-and-conquer. The runtime of the algorithms is $O(K N \log^2 N)$, where $K$ is the number of clusters and $N$ is the number of circular data points.

The core algorithm solves the framed clustering problem by minimizing the within-cluster sum of squared distances. In contrast to heuristic clustering, the efficiency and accuracy of fast optimal framed clustering is evident both theoretically and practically. The fast algorithm use divide-and-conquer inside dynamic programming inside another divide-and-conquer, guaranteeing cluster optimality---the total of within-cluster sums of squared distances is always the minimum given the number of clusters. Two plot functions visualize circular and framed clusters to reveal patterns in data. On a desktop computer using a single processor core, millions of data points can be clustered within seconds. The algorithms offer a general high-performance solution to circular, periodic, or framed data clustering.

### Main methods

The function to perform fast optimal circular data clustering is `CirClust()`. The function to perform the fast framed clustering is `FramedClust()`. Both functions with default `method` have a runtime of $O(K N \log^2N)$, linear polylogarithmic in sample size $N$ and linear in the number of clusters $K$. The space complexity of both functions is $O(KN)$. Their implementations scale large input data to improve numerical precision.

### When to use the package

One can use the algorithms to cluster circular, periodic, or framed data. Here circular data broadly refer to data points on any non-self-intersecting loop. Periodic data clustering can be formulated as circular clustering. Framed clustering finds one frame of fixed size that contains best clusters. One can use them to characterize events along circular DNA molecules, circular RNA molecules, and circular genomes of bacteria, chloroplast, and mitochondria. One can also cluster climate data along any given longitude or latitude.

### To download and install the package
```{r}
install.packages("OptCirClust")
```
