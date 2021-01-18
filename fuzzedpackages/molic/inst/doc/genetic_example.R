## ----global_options, include=FALSE--------------------------------------------
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE, include = TRUE)

## ----e------------------------------------------------------------------------
library(molic)
library(ess)   # for fit_components
library(dplyr)
set.seed(7)    # for reproducibility

## Data

# The components - here microhaplotypes
haps <- tgp_haps

# All the Europeans
eur <- tgp_dat %>%
  as_tibble() %>%
  filter(pop_meta == "EUR")

# Extracting the two databases for each copy of the chromosomes
eur_a <- eur %>%
  filter(grepl("a$", sample_name)) %>%
  select(-c(1:2))

eur_b <- eur %>%
  filter(grepl("b$", sample_name)) %>%
  select(-c(1:2))

## -----------------------------------------------------------------------------
print(eur, n_extra = 0)

## ----fig.align='center', fig.width = 6, fig.height = 6,fig.show='hold'--------
ga <- fit_components(eur_a, comp = haps, trace = FALSE)
gb <- fit_components(eur_b, comp = haps, trace = FALSE)
print(ga)
plot(ga, vertex.size = 3, vertex.label = NA)

## -----------------------------------------------------------------------------
# Only 500 simulations is used here to exeplify
# The default number of simulations is 10,000
m1 <- fit_outlier(eur_a, ga, nsim = 500) # consider using more cores (ncores argument)
m2 <- fit_outlier(eur_b, gb, nsim = 500) # consider using more cores (ncores argument)
m  <- fit_mixed_outlier(m1, m2)

## -----------------------------------------------------------------------------
print(m)

## ----fig.align='center', fig.width = 6, fig.height = 6,fig.show='hold'--------
plot(m)

## -----------------------------------------------------------------------------
outs <- outliers(m)
eur_a_outs <- eur_a[which(outs), ]
eur_b_outs <- eur_b[which(outs), ]
print(eur_a_outs, n_extra = 0)

## -----------------------------------------------------------------------------
x1 <- rbind(eur_a_outs[1, ], eur_b_outs[1, ])
x2 <- rbind(eur_a[1, ], eur_b[1, ])
dev1 <- deviance(m, x1) # falls within the critical region in the plot (the red area)
dev2 <- deviance(m, x2) # falls within the acceptable region in the plot

dev1
dev2

# Retrieving the pvalues
pval(m, dev1)
pval(m, dev2)

## ----fig.align='center', fig.width = 6, fig.height = 6,fig.show='hold'--------
amr <- tgp_dat %>%
  as_tibble() %>%
  filter(pop_meta == "AMR")

z1  <- amr %>%
  filter(grepl("a$", sample_name)) %>% 
  select(unname(unlist(haps))) %>%
  slice(1) %>%
  unlist()

z2  <- amr %>%
  filter(grepl("b$", sample_name)) %>% 
  select(unname(unlist(haps))) %>%
  slice(1) %>%
  unlist()

# Only 500 simulations is used here to exemplify
# The default number of simulations is 10,000
m3 <- fit_outlier(eur_a, ga, z1, nsim = 500) # consider using more cores (ncores argument)
m4 <- fit_outlier(eur_b, gb, z2, nsim = 500) # consider using more cores (ncores argument)
m5 <- fit_mixed_outlier(m3, m4)
print(m5)
plot(m5)

