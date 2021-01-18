## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

library(armspp)

## ----normal_example------------------------------------------------------
output <- arms(
  5000, function(x) -x ^ 2 / 2,
  -1000, 1000,
  metropolis = FALSE, include_n_evaluations = TRUE
)
print(str(output))
hist(output$samples, br = 50, freq = FALSE, main = 'Normal samples')
shapiro.test(output$samples)

## ----mixture_of_normals_density------------------------------------------
dnormmixture <- function(x) {
  parts <- log(c(0.4, 0.6)) + dnorm(x, mean = c(-1, 4), log = TRUE)
  log(sum(exp(parts - max(parts)))) + max(parts)
}

curve(exp(Vectorize(dnormmixture)(x)), -7, 10)

## ----mixture_of_normals_example------------------------------------------
samples <- arms(1000, dnormmixture, -1000, 1000)
hist(samples, freq = FALSE, br = 50)
curve(exp(Vectorize(dnormmixture)(x)), -7, 10, col = 'red', add = TRUE)

