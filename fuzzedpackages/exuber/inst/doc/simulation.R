## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE, warning = FALSE, message = FALSE,
  comment = "#>"
)

## ----load, echo=FALSE, message=FALSE------------------------------------------
library(exuber)

## ----generate simulated data--------------------------------------------------
set.seed(125)
# The fundamental value from the Lucas pricing model
pf <- sim_div(400)

# The Evans bubble term
pb <- sim_evans(400)

# the scaling factor for the bubble
kappa <- 20
# 
# The simulated price
p <- pf + kappa*pb

## ----plot simulation, eval = TRUE, fig.width = 9, echo=FALSE------------------
library(ggplot2)
library(gridExtra)
library(dplyr)

tibble(
  index = 1:NROW(p),
  "Price" = p
  ) %>% 
  ggplot(aes(x = index, y = Price)) + 
  geom_line() +
  theme_bw() +
  labs(y = "", x = "")

## ----simulate-----------------------------------------------------------------
library(ggplot2)
library(purrr)

sims <- tibble(
  sim_psy1 = sim_psy1(100),
  sim_psy2 = sim_psy2(100),
  sim_evans = sim_evans(100),
  sim_blan = sim_blan(100)
)

## ----autoplot-sim, fig.width = 9, fig.height=7--------------------------------
sims %>%
  map(autoplot) %>%
  do.call(grid.arrange, .)

