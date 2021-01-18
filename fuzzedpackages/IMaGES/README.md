# IMaGES

This is the repository for the R implementation of the IMaGES algorithm. This project was initiated by SJ Hanson (RUBIC, Rutgers University). The repository started as a fork of [pcalg](https://cran.r-project.org/package=pcalg) and is now a standalone product. The additional code and changes were written/made by Noah Frazier-Logue.

IMaGES is based on the paper 

Ramsey JD, Hanson SJ, Hanson C, Halchenko YO, Poldrack RA, Glymour C (2010). Six problems for causal inference from fMRI. Neuroimage, 49, 1545-1558.

This algorithm elaborates on the GES algorithm by using a global score across the supplied datasets and operating over the datasets concurrently to determine the representative graph(s) with the best goodness of fit.


**NOTE**: This software is in beta! If you come across any issues while using this package or have any suggestions for improvement, submit a pull request.

### Installation

To install from this repository, simply run these commands in an R shell:

```
> library(devtools)
> install_github("noahfl/IMaGES")
```

TODO: Add stuff about CRAN when that becomes relevant.

### Usage

```R
#matrices should be a list of >= 1 datasets with an optional header
matrices <- list(matrix1, matrix2,...)

#load supplied sample data
data(IMdata)

im.results <- IMaGES(matrices=IMData, penalty=3, num.markovs=5)

#plot individual graph, in this case the global graph
plotIMGraph(im.results$.global)

#plot Markov Equivalence Class (size specified by num.markovs)
plotMarkovs(im.results)

#plot global graph with SEM data, and all individual datasets' SEM data
#imposed on the global graph
plotAll(im.results)

#compare IMaGES result against individual graphs
data(IMTrue)

for (i in 1:length(IMTrue)) {
    plot(IMTrue[[i]])
}

```
