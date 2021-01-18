## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  message = FALSE,
  fig.width = 10
)

## ----load, echo=FALSE, message=FALSE------------------------------------------
library(exuber)
options(exuber.show_progress = FALSE)
library(dplyr)
library(ggplot2)
library(tidyr)

## ----options, echo=FALSE------------------------------------------------------
options(exuber.parallel = FALSE)

## ----tstats-cv----------------------------------------------------------------
set.seed(123)
sims <- tibble(
  sim_psy1 = sim_psy1(100),
  sim_psy2 = sim_psy2(100),
  sim_evans = sim_blan(100),
  sim_blan = sim_evans(100),
) 

# Esimation
estimation <- radf(sims, lag = 1)
  
# Critical Values
crit_values <- radf_mc_cv(nrow(sims))

## ----autoplot-basic-----------------------------------------------------------
autoplot(estimation, crit_values)

## ----autoplot-color-theme-----------------------------------------------------
autoplot(estimation, crit_values) +
  scale_color_manual(values = c("grey","black")) + 
  theme_classic()

## ----autoplot-shade-----------------------------------------------------------
autoplot(estimation, crit_values, shade_opt = shade(fill = "pink", opacity = 0.3))

## ----join-sets----------------------------------------------------------------
joined <- augment_join(estimation, crit_values)
joined

## ----facet-joined-------------------------------------------------------------
joined %>% 
  ggplot(aes(x = index)) +
  geom_line(aes(y = tstat)) +
  geom_line(aes(y = crit)) +
  facet_grid(sig + name ~  id  , scales = "free_y")

## ----facet-joined-theme-exuber, warning=FALSE---------------------------------
joined %>%
  pivot_longer(cols = c("tstat", "crit"), names_to = "nms") %>% 
  ggplot(aes(x = index, y = value, col = nms)) +
  geom_line() +
  facet_grid(sig + name ~  id  , scales = "free_y") +
  scale_exuber_manual() +
  theme_exuber()

## ----distributions------------------------------------------------------------
distr <- radf_mc_distr(n = 300)
autoplot(distr)

## ----ecdf---------------------------------------------------------------------
library(tidyr)
distr %>%
  tidy() %>%
  rename_all(~ stringr::str_to_upper(.)) %>%
  gather(Statistic, value, factor_key = TRUE) %>%
  ggplot(aes(value, color = Statistic)) +
  stat_ecdf() +
  ggtitle("Empirical Cumulative Distribution") +
  geom_hline(yintercept = 0.95, linetype = "dashed") + theme_bw()

## ----lapply-arrange-----------------------------------------------------------
library(gridExtra)

# To choose only positive series (i.e. statistically significant for 5%)
positive_series <- diagnostics(estimation, crit_values)$positive 

# Through a loop on positive series 
plot_list1 <- list()
for (as in positive_series) {
  plot_list1[[as]] <- autoplot(estimation, crit_values, select_series = as)
}

# Alternatively  with lapply
plot_list2 <- lapply(positive_series, function(x) autoplot(estimation, crit_values, select_series = x))
names(plot_list2) <- positive_series

do.call(gridExtra::grid.arrange, plot_list1)

## ----example-old--------------------------------------------------------------
plot_list1[[1]] <- plot_list1[[1]] + theme_classic()

