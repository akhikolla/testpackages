## -----------------------------------------------------------------------------
outcome_str <- '1NNN 2NTN 2NNN'

## -----------------------------------------------------------------------------
skeleton <- c(0.05, 0.15, 0.25, 0.4, 0.6)
target <- 0.25

## ---- message=FALSE-----------------------------------------------------------
library(trialr)

path <- crm_path_analysis(
  outcome_str = outcome_str,
  skeleton = skeleton, target = target, model = 'empiric',
  beta_sd = 1, seed = 123, refresh = 0)

## -----------------------------------------------------------------------------
names(path)

## -----------------------------------------------------------------------------
library(tibble)

df <- as_tibble(path)
df

## ---- message=FALSE-----------------------------------------------------------
library(tidyr)
library(purrr)
library(dplyr)

df %>% 
  mutate(prob_tox = fit %>% map('prob_tox')) %>% 
  select(outcomes, dose_index, prob_tox) %>% 
  unnest(cols = c(dose_index, prob_tox)) %>% 
  filter(dose_index == 2)

## ---- message=FALSE, eval=FALSE-----------------------------------------------
#  library(tidyr)
#  library(purrr)
#  library(dplyr)
#  
#  df %>%
#    mutate(prob_tox = fit %>% map('prob_tox')) %>%
#    select(outcomes, dose_index, prob_tox) %>%
#    unnest %>%
#    filter(dose_index == 2)

## -----------------------------------------------------------------------------
df %>% 
  mutate(
    recommended_dose = fit %>% map_int('recommended_dose'),
    careful_dose = fit %>% map_dbl(careful_escalation, 
                                   tox_threshold = target + 0.1, 
                                   certainty_threshold = 0.7)
  ) %>% 
  select(outcomes, recommended_dose, careful_dose)

## -----------------------------------------------------------------------------
paths <- crm_path_analysis(
  outcome_str = '1NNN 2NTN 2NNN 3TTT 1TTT 1TNT',
  skeleton = skeleton, target = target, model = 'logistic',
  a0 = 3, beta_mean = 0,beta_sd = 1,
  seed = 123, refresh = 0)

df <- as_tibble(paths)

df %>% 
  mutate(
    recommended_dose = fit %>% map_int('recommended_dose'),
    careful_dose = fit %>% map_dbl(careful_escalation, 
                                   tox_threshold = target + 0.1, 
                                   certainty_threshold = 0.7)
  ) %>% 
  select(outcomes, recommended_dose, careful_dose)

## -----------------------------------------------------------------------------
paths <- crm_dtps(skeleton = skeleton,
                  target = target, 
                  model = 'empiric', 
                  cohort_sizes = c(2, 2), 
                  next_dose = 2, 
                  beta_sd = 1,
                  refresh = 0)

## -----------------------------------------------------------------------------
df <- as_tibble(paths)
df

## -----------------------------------------------------------------------------
spread_paths(df %>% select(-fit, -parent_fit, -dose_index))

## -----------------------------------------------------------------------------
paths2 <- crm_dtps(skeleton = skeleton,
                   target = target,
                   model = 'empiric',
                   cohort_sizes = c(3, 3),
                   previous_outcomes = '2NN 3TN',
                   next_dose = 2, 
                   beta_sd = 1,
                   refresh = 0)

## -----------------------------------------------------------------------------
spread_paths(as_tibble(paths2) %>% select(-fit, -parent_fit, -dose_index))

## -----------------------------------------------------------------------------
paths3 <- crm_dtps(
  skeleton = skeleton,
  target = target,
  model = 'empiric',
  cohort_sizes = c(3, 3),
  previous_outcomes = '2NN 3TN',
  next_dose = 2, 
  beta_sd = 1,
  user_dose_func = function(x) {
    careful_escalation(x, tox_threshold = target + 0.1, 
                       certainty_threshold = 0.7)
  }, 
  seed = 123, refresh = 0)

df3 <- as_tibble(paths3)
spread_paths(df3 %>% select(-fit, -parent_fit, -dose_index))

## ---- fig.width=6, fig.height=6-----------------------------------------------
# This section of code outputs to the Viewer pane in RStudio
if(Sys.getenv("RSTUDIO") == "1") {
  
  library(DiagrammeR)
  
  df3 %>%
    transmute(id = .node,
              type = NA,
              label = case_when(
                is.na(next_dose) ~ 'Stop',
                TRUE ~ next_dose %>% as.character()),
              shape = 'circle',
              fillcolor = case_when(
                next_dose == 1 ~ 'slategrey',
                next_dose == 2 ~ 'skyblue1',
                next_dose == 3 ~ 'royalblue1',
                next_dose == 4 ~ 'orchid4',
                next_dose == 5 ~ 'royalblue4',
                is.na(next_dose) ~ 'red'
              )
    ) -> ndf
  
  df3 %>% 
    filter(!is.na(.parent)) %>% 
    select(from = .parent, to = .node, label = outcomes) %>% 
    mutate(rel = "leading_to") -> edf
  
  graph <- create_graph(nodes_df = ndf, edges_df = edf)
  render_graph(graph)
}

