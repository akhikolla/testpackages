## -----------------------------------------------------------------------------
outcomes <- '2NN 3NN 4TT'

## -----------------------------------------------------------------------------
skeleton <- c(0.05, 0.15, 0.25, 0.4, 0.6)
target <- 0.25

## ---- message=FALSE, warning=FALSE--------------------------------------------
library(trialr)

## ---- results = "hide", warning=FALSE, message=FALSE--------------------------
fit <- stan_crm(outcomes, skeleton = skeleton, target = target, 
                beta_sd = sqrt(1.34), seed = 123)
fit

## ---- warning=FALSE, message=FALSE--------------------------------------------
library(dplyr)
library(tidybayes)

prob_tox_samp_tall <- fit %>% 
  gather_draws(prob_tox[dose]) %>% 
  rename(prob_dlt = .value) %>% 
  ungroup

## -----------------------------------------------------------------------------
prob_tox_samp_tall %>% head(10)

## ---- fig.width=7, fig.height=7-----------------------------------------------
library(ggplot2)

prob_tox_samp_tall %>% 
  ggplot(aes(x = dose, y = prob_dlt, group = dose)) +
  geom_boxplot() + 
  ylim(0, 1) + 
  labs(title = 'Boxplot of Prob(DLT) under CRM')

## ---- fig.width=7, fig.height=7-----------------------------------------------
prob_tox_samp_tall %>% 
  ggplot(aes(x = dose, y = prob_dlt, group = dose)) +
  geom_violin(fill = 'orange') + 
  ylim(0, 1) + 
  labs(title = 'Violin plot of Prob(DLT) under CRM')

## ---- fig.width=7, fig.height=7, message=FALSE--------------------------------
library(ggridges)

prob_tox_samp_tall %>% 
  mutate(dose = factor(dose)) %>% 
  ggplot(aes(x = prob_dlt, y = dose, fill = dose)) +
  geom_density_ridges() + 
  theme(legend.position = 'none') +
  labs(title = 'Joyplot of Prob(DLT) under CRM') + 
  theme(legend.position = 'bottom')

## ---- fig.width=7, fig.height=7, message=FALSE--------------------------------
prob_tox_samp_tall  %>% 
  group_by(.draw) %>% 
  summarise(mtd = dose[which.min(abs(prob_dlt - target))]) %>% 
  mutate(mtd = factor(mtd)) -> mtd_candidates

prob_tox_samp_tall %>% 
  left_join(mtd_candidates, by = '.draw') %>% 
  filter(.draw <= 200) %>% 
  ggplot(aes(x = dose, y = prob_dlt, group = .draw)) +
  geom_line(aes(col = mtd), alpha = 0.5) + 
  geom_hline(yintercept = target, col = 'red', linetype = 'dashed') + 
  labs(title = 'The identify of the MTD is still shrouded in mystery', 
       y = 'Prob(DLT)', col = 'MTD') +
  theme(legend.position = 'bottom')

## ---- fig.width=7, fig.height=7, message=FALSE, warning=FALSE-----------------
mtd_candidates %>% 
  count(mtd) %>% 
  mutate(prob_mtd = n / sum(n)) %>% 
  ggplot(aes(x = mtd, y = prob_mtd, fill = mtd)) + 
  geom_col() +
  labs(x = 'MTD') +
  theme(legend.position = 'bottom')

## ---- fig.width=7, fig.height=7-----------------------------------------------
fit %>% 
  gather_draws(prob_tox[dose]) %>% 
  group_by(dose) %>% 
  summarise(prob_too_toxic = mean(.value > target)) %>%
  ggplot(aes(x = dose, y = prob_too_toxic, fill = dose)) + 
  geom_col() + 
  scale_fill_gradient(low="green", high="red") + 
  labs(title = 'Posterior probability that each dose is too toxic',
       y = 'Prob(DLT risk > target)', fill = 'Probability dose is too toxic') +
  theme(legend.position = 'bottom')

