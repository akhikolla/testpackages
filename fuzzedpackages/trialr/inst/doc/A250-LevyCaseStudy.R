## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- message=FALSE-----------------------------------------------------------
library(dplyr)

tibble(
  `Dose-level` = 1:5,
  `Dose (mg/m2 per day)` = c(0.5, 1, 3, 5, 6),
  `Prior Pr(DLT)` = c(0.05, 0.1, 0.15, 0.33, 0.5)
) %>% knitr::kable()

## ---- results='hide', warning=FALSE-------------------------------------------
library(trialr)

target <- 0.33
skeleton <- c(0.05, 0.1, 0.15, 0.33, 0.5)

fit0 <- crm_prior_beliefs(skeleton, target, model = 'logistic_gamma',
                          a0 = 1, beta_shape = 1, beta_inverse_scale = 1)

## -----------------------------------------------------------------------------
fit0

## ---- fig.width=7, fig.height=5, results='hide', warning=FALSE----------------
library(tidyr)
library(purrr)
library(ggplot2)

get_prior_fit <- function(a0) {
  crm_prior_beliefs(skeleton, target, 
                    model = 'logistic_gamma', a0 = a0, 
                    beta_shape = 1, beta_inverse_scale = 1)
}

tibble(a0 = c(-1, 0, 1, 2, 3, 4)) %>% 
  mutate(Mod = map(a0, get_prior_fit)) %>% 
  mutate(
    Dose = Mod %>% map("dose_indices"),
    ProbTox = Mod %>% map("prob_tox"),
  ) %>% 
  select(-Mod) %>% 
  unnest(cols = c(Dose, ProbTox)) %>% 
  mutate(a0 = factor(a0)) %>% 
  ggplot(aes(x = Dose, y = ProbTox, group = a0, col = a0)) + 
  geom_line() + 
  ylim(0, 1) + 
  labs(title = 'Prior Prob(DLT) location is affected by the fixed intercept, a0')

## ---- eval=FALSE--------------------------------------------------------------
#  tibble(a0 = c(-1, 0, 1, 2, 3, 4)) %>%
#    mutate(Mod = map(a0, get_prior_fit)) %>%
#    mutate(
#      Dose = Mod %>% map("dose_indices"),
#      ProbTox = Mod %>% map("prob_tox"),
#      ) %>%
#    select(-Mod) %>%
#    unnest(cols = c(Dose, ProbTox)) %>%
#    mutate(a0 = factor(a0)) %>%
#    ggplot(aes(x = Dose, y = ProbTox, group = a0, col = a0)) +
#    geom_line() +
#    ylim(0, 1) +
#    labs(title = 'Prior Prob(DLT) location is affected by the fixed intercept, a0')

## ---- results='hide', warning=FALSE-------------------------------------------
fit1 <- stan_crm(outcome_str = '1NNN', 
                 skeleton = skeleton, target = target, model = 'logistic_gamma', 
                 a0 = 1, beta_shape = 1, beta_inverse_scale = 1,
                 seed = 123, control = list(adapt_delta = 0.99))

## -----------------------------------------------------------------------------
fit1

## ---- results='hide', warning=FALSE-------------------------------------------
fit2 <- stan_crm(outcome_str = '1NNN 3NNT', 
                 skeleton = skeleton, target = target, model = 'logistic_gamma', 
                 a0 = 1, beta_shape = 1, beta_inverse_scale = 1,
                 seed = 123, control = list(adapt_delta = 0.99))

## -----------------------------------------------------------------------------
fit2

## ---- results='hide', warning=FALSE-------------------------------------------
fit3 <- stan_crm(outcome_str = '1NNN 3NNT 4NNT', 
                 skeleton = skeleton, target = target, model = 'logistic_gamma', 
                 a0 = 1, beta_shape = 1, beta_inverse_scale = 1,
                 seed = 123, control = list(adapt_delta = 0.99))

## -----------------------------------------------------------------------------
fit3

## ---- results='hide', warning=FALSE-------------------------------------------
fit4 <- stan_crm(outcome_str = '1NNN 3NNT 4NNT 4NNN', 
                 skeleton = skeleton, target = target, model = 'logistic_gamma', 
                 a0 = 1, beta_shape = 1, beta_inverse_scale = 1,
                 seed = 123, control = list(adapt_delta = 0.99))

## -----------------------------------------------------------------------------
fit4

## ---- results='hide', warning=FALSE-------------------------------------------
fit5 <- stan_crm(outcome_str = '1NNN 3NNT 4NNT 4NNN 4NTN', 
                 skeleton = skeleton, target = target, model = 'logistic_gamma', 
                 a0 = 1, beta_shape = 1, beta_inverse_scale = 1,
                 seed = 123, control = list(adapt_delta = 0.99))

## -----------------------------------------------------------------------------
fit5

## ---- results='hide', warning=FALSE-------------------------------------------
fit6 <- stan_crm(outcome_str = '1NNN 3NNT 4NNT 4NNN 4NTN 4TNT', 
                 skeleton = skeleton, target = target, model = 'logistic_gamma', 
                 a0 = 1, beta_shape = 1, beta_inverse_scale = 1,
                 seed = 123, control = list(adapt_delta = 0.99))

## -----------------------------------------------------------------------------
fit6

## -----------------------------------------------------------------------------
apply(as.data.frame(fit6, pars = 'prob_tox'), 2, quantile, 
      probs = c(0.025, 0.975))

## -----------------------------------------------------------------------------
prob_tox_mtd <- fit6$prob_tox[fit6$recommended_dose]
prob_tox_mtd

## ---- results = 'hide', warning=FALSE-----------------------------------------
paths <- crm_dtps(skeleton = skeleton, target = target, 
                  model = 'logistic_gamma', cohort_sizes = c(3), 
                  previous_outcomes = '1NNN 3NNT 4NNT 4NNN 4NTN 4TNT',
                  a0 = 1, beta_shape = 1, beta_inverse_scale = 1,
                  seed = 123, control = list(adapt_delta = 0.99), refresh = 0)

library(tibble)
df <- as_tibble(paths)

## -----------------------------------------------------------------------------
library(dplyr)
library(purrr)
library(tidyr)

df %>% 
  mutate(prob_tox = map(fit, 'prob_tox')) %>% 
  select(-fit, -parent_fit) %>% 
  unnest(cols = c(dose_index, prob_tox)) %>% 
  filter(dose_index == 4)

## -----------------------------------------------------------------------------
df %>% 
  filter(.depth > 0) %>% 
  mutate(prob_tox = map(fit, 'prob_tox')) %>% 
  select(-fit, -parent_fit) %>% 
  unnest(cols = c(dose_index, prob_tox)) %>% 
  filter(dose_index == 4) %>% 
  select(outcomes, prob_tox) %>% 
  bind_cols(
    lik = dbinom(x = 0:3, size = 3, prob = prob_tox_mtd)) -> future_scenario

future_scenario

## -----------------------------------------------------------------------------
future_scenario %>% 
  mutate(prob_tox_change = abs(prob_tox - prob_tox_mtd)) %>% 
  summarise(expected_change = sum(lik * prob_tox_change))

## -----------------------------------------------------------------------------
levy_reported <- tibble(
  Dose = 1:5,
  ProbTox = c(0.06, 0.12, 0.17, 0.36, 0.53),
)

## -----------------------------------------------------------------------------
fit_levy_crm <- function(outcomes, a0, beta_inverse_scale) {
  stan_crm(outcome_str = outcomes, 
           skeleton = skeleton, target = target, 
           model = 'logistic_gamma', 
           a0 = a0, beta_shape = 1, 
           beta_inverse_scale = beta_inverse_scale,
           control = list(adapt_delta = 0.99), 
           seed = 123, refresh = 0)
}

## ---- results='hide', warning=FALSE-------------------------------------------
expand.grid(a0 = 1:4, beta_inverse_scale = c(0.5, 1, 2)) %>% 
  mutate(Series = rownames(.)) %>% 
  mutate(Mod = map2(a0, beta_inverse_scale, fit_levy_crm, 
                    outcomes = '1NNN 3NNT 4NNT 4NNN 4NTN 4TNT')) %>% 
  mutate(
    Dose = Mod %>% map("dose_indices"),
    ProbTox = Mod %>% map("prob_tox"),
  ) %>% 
  select(-Mod) %>% 
  unnest(cols = c(Dose, ProbTox)) %>% 
  mutate(a0 = factor(a0), 
         beta_inverse_scale = factor(beta_inverse_scale)) -> all_fits

## ---- fig.width=7, fig.height=5-----------------------------------------------
all_fits %>% 
  ggplot(aes(x = Dose, y = ProbTox)) + 
  geom_line(aes(group = Series, col = a0)) + 
  geom_line(data = levy_reported, col = 'green', size = 1.2) + 
  facet_wrap(~ beta_inverse_scale) + 
  ylim(0, 0.7) + 
  labs(title = "None of the exponential models quite matches the investigators' inferences")

