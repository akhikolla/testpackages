## ---- SETTINGS-knitr, include=FALSE-------------------------------------------
stopifnot(require(knitr))
opts_chunk$set(
  comment=NA, 
  message = FALSE, 
  warning = FALSE, 
  eval = identical(Sys.getenv("NOT_CRAN"), "true"),
  dev = "png",
  dpi = 150,
  fig.asp = 0.618,
  fig.width = 5,
  out.width = "60%",
  fig.align = "center"
)

## ---- SETTINGS-gg, include=TRUE-----------------------------------------------
library(ggplot2)
library(bayesplot)
theme_set(bayesplot::theme_default())

## ---- default-prior-1, results="hide", eval = TRUE----------------------------
library("rstanarm")
default_prior_test <- stan_glm(mpg ~ wt + am, data = mtcars, chains = 1)

## ---- default-prior-summary, eval = TRUE--------------------------------------
prior_summary(default_prior_test)

## ---- echo=FALSE, eval = TRUE-------------------------------------------------
priors <- prior_summary(default_prior_test)
fr2 <- function(x) format(round(x, 2), nsmall = 2)

## ---- no-autoscale, results="hide", eval = TRUE-------------------------------
test_no_autoscale <-
  update(
    default_prior_test,
    prior = normal(0, 5),
    prior_intercept = student_t(4, 0, 10),
    prior_aux = cauchy(0, 3)
  )

## ---- no-autoscale-prior-summary, eval = TRUE---------------------------------
prior_summary(test_no_autoscale)

## ---- eval = TRUE-------------------------------------------------------------
p <- 1 - 2 * pnorm(-250, mean = 0, sd = 500)
print(paste("Pr(-250 < theta < 250) =", round(p, 2)))

## ---- fig.cap="_There is much more probability mass outside the interval (-250, 250)._", eval = TRUE----
theta <- rnorm(1e5, mean = 0, sd = 500)
p_approx <- mean(abs(theta) < 250)
print(paste("Pr(-250 < theta < 250) =", round(p_approx, 2)))

d <- data.frame(theta, clr = abs(theta) > 250)
library(ggplot2)
ggplot(d, aes(x = theta, fill = clr)) + 
  geom_histogram(binwidth = 5, show.legend = FALSE) + 
  scale_y_continuous(name = "", labels = NULL, expand = c(0,0)) + 
  scale_x_continuous(name = expression(theta), breaks = c(-1000, -250, 250, 1000))

## ---- flat-prior-1, echo=FALSE, results="hide", eval = TRUE-------------------
flat_prior_test <- stan_glm(mpg ~ wt, data = mtcars, prior = NULL, iter = 10, chains = 1)

## ---- flat-prior-2, eval=FALSE, eval = TRUE-----------------------------------
flat_prior_test <- stan_glm(mpg ~ wt, data = mtcars, prior = NULL)

## ---- flat-prior-summary, eval = TRUE-----------------------------------------
prior_summary(flat_prior_test)

## ---- eval=FALSE--------------------------------------------------------------
#  my_prior <- normal(location = c(-10, 0), scale = c(5, 2))
#  stan_glm(y ~ x1 + x2, data = dat, prior = my_prior)

