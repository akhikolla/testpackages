# cbq: An R Package for Conditional Binary Quantile Models


The package `cbq` provides basic functionalities of conditional binary quantile models using Markov chain Monte Carlo methods. The estimation is conducted through pre-compiled stan codes. The conditional binary quantile (CBQ) models extend the simple version of binary quantile models for analyzing discrete choices (including but beyond binary choices) with varying choice alternatives. Each quantile estimation represents a local inspection of effects at the specified quantile of the choice probabilities. In the simple binary setting, the features of the choice alternatives within each choice set are assumed to be the same. However, in reality we oftentimes observe individuals who are facing different sets and numbers of choice alterantives. The CBQ models are developed to solve this problem by introducing a conditional multinomial structure for modeling varying choice alternatives. Even though the CBQ models are called "binary", they actually belong to a more general family of dicrete choice models. I refer the readers to [Lu (2020)](https://doi.org/10.1017/pan.2019.29)
for the details of the estimation.


## Installation

```r
# Make sure that the following packages have been installed in your local R environment
if(!require(rstan)) install.packages("rstan")

# Install cirque from github
if(!require(devtools)) install.packages("devtools")
devtools::install_github("xiao-lu-research/cbq")
```


## Usage

```r

# Load the package
library(cbq)

# Get help
?cbq

# Simulate the data
x <- rnorm(50)
y <- ifelse(x > 0, 1, 0)
dat <- as.data.frame(cbind(y, x))

# Estimate the CBQ model
model <- cbq(y ~ x, dat, 0.5)

# Show the results
print(model)
coef(model)
plot(model)

```

## References

Lu, Xiao. (2020). Discrete Choice Data with Unobserved Heterogeneity: A Conditional Binary Quantile Model. Political Analysis, 28(2), 147-167.
