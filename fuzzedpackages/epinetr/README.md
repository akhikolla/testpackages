# epinetr
An R Package for Epistastic Network Modelling with Forward-Time Simulation

## About this project
*epinetr* is a package for the R statistical computing environment designed to facilitate the modelling of epistasis and epistatic networks of arbitrary complexity in populations across generations. Our hope is that this software will aid researchers in uncovering the genetic architecture of complex traits and bridging the conceptual divide between quantitative and molecular genetics.

Using *epinetr*, you can test  the impacts of various mixes of additive and epistatic effects against different population structures and selection criteria on populations. Our primary goal is to investigate the relationship between biological epistasis in individuals and additive models in populations.

## Installation
Installation is straightforward, provided you already have the `devtools` package installed. Simply run the command

```r
install_github("diondetterer/epinetr")
```

and the *epinetr* package will be installed into your R library.

## Usage
There is a vignette in the package which provides a fairly comprehensive tutorial, and we encourage all users to read it. However, here are some minimal commands to get you started:

```r
pop <- Population(
  popSize = 500, map = map100snp, QTL = 20,
  alleleFrequencies = runif(100), broadH2 = 0.9,
  narrowh2 = 0.75, traitVar = 40
)
```

This will create a `Population` object called `pop` with 500 individuals, a chromosome map given by `map100snp`, 20 randomly selected QTLs, randomly-generated allele frequencies, broad-sense heritability at 0.9, narrow-sense heritability at 0.75 and trait variance at 40.

```r
pop <- addEffects(pop)
pop <- attachEpiNet(pop)
```

These commands will attach additive and epistatic effects to the population.

```r
plot(getEpiNet(pop))
```

This will provide a visualisation of the epistatic network.

```r
pop <- runSim(pop, generations = 150)
```

This will run the simulator for 150 generations.

Finally, plot the run:

```r
plot(pop)
```

## Authors and support
Dion Detterer, Paul Kwan and Cedric Gondro wrote the *epinetr* package, with Dion as the maintainer.

Issues can be reported via the [issues tab](https://github.com/diondetterer/epinetr/issues), or you can email Dion at ddettere@myune.edu.au for assistance.

## Contributing
We welcome contributions to the project; please see the [project wiki](https://github.com/diondetterer/epinetr/wiki) for details on the codebase.

For advice on setting up an appropriate R development environment, see Hadley Wickham's advice on system setup at https://r-pkgs.org/setup.html

## License
*epinetr* is released under the GPLv3 license. See the file `LICENSE` for more details.
