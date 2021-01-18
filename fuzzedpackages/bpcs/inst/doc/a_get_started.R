## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, echo=T, results='hide', warning=F, message=F----------------------
library(bpcs)
library(ggplot2)
library(dplyr)
library(tibble)
library(kableExtra)
library(bayesplot)
library(knitr)

## ----eval=FALSE, echo=T-------------------------------------------------------
#  remotes::install_github('davidissamattos/bpcs')

## -----------------------------------------------------------------------------
library(bpcs)

## -----------------------------------------------------------------------------
knitr::kable(tennis_agresti) %>% 
  kableExtra::kable_styling()

## -----------------------------------------------------------------------------
m1 <- bpc(data = tennis_agresti,
          player0 = 'player0',
          player1 = 'player1',
          result_column = 'y',
          model_type = 'bt',
          solve_ties = 'none', #there are no ties
          show_chain_messages = T)

## ----eval=F-------------------------------------------------------------------
#  launch_shinystan(m1)

## ----eval=F-------------------------------------------------------------------
#  stanfit <- get_stanfit(m1)
#  shinystan::launch_shinystan(stanfit)

## -----------------------------------------------------------------------------
knitr::kable(m1$lookup_table)

## -----------------------------------------------------------------------------
stanfit <- get_stanfit(m1)
posterior<-rstan::extract(stanfit,inc_warmup=T,permuted=F)

## ----eval=F-------------------------------------------------------------------
#  bayesplot::mcmc_trace(posterior,pars = c("lambda[1]","lambda[2]","lambda[3]","lambda[4]"), n_warmup=1000)

## -----------------------------------------------------------------------------
rstan::summary(stanfit ,pars=c('lambda'))$summary

## -----------------------------------------------------------------------------
y<-as.vector(tennis_agresti$y)
yrep<-predict(m1,tennis_agresti,n=100,return_matrix = T)
yrep<-yrep[,1:46] #from column 47 we have if it was a tie or not. We just need to remove this

## -----------------------------------------------------------------------------
bayesplot::ppc_bars(y=y, yrep=yrep) +
  labs(title = 'Bar plot with medians and uncertainty\n intervals superimposed')

## -----------------------------------------------------------------------------
summary(m1)

## -----------------------------------------------------------------------------
knitr::kable(get_hpdi_parameters(m1), caption = 'Parameter distribution and the High Posterior Density intervals', digits = 2) %>% 
  kable_styling()

## -----------------------------------------------------------------------------
hpdi <- get_hpdi_parameters(m1) %>%
  dplyr::filter(startsWith(Parameter, "lambda"))

ggplot2::ggplot(hpdi, aes(x = Parameter)) +
  ggplot2::geom_pointrange(aes(ymin = HPD_lower,
                               ymax = HPD_higher,
                               y = Mean)) +
  ggplot2::labs(y = "Estimate", x = "Player", title = "HPDI interval of the strength of the players") +
  ggplot2::coord_flip()


## -----------------------------------------------------------------------------
prob_table<-get_probabilities(m1)$Table
knitr::kable(prob_table, caption = 'Probabilities of one player beating the other', digits = 2) %>% 
  kableExtra::kable_styling()

## -----------------------------------------------------------------------------
ranking <- get_rank_of_players(m1)

## -----------------------------------------------------------------------------
t <- ranking %>% dplyr::select(Parameter, MedianRank, StdRank)
knitr::kable(t, caption = 'Rank of the players') %>%
  kable_styling()

## -----------------------------------------------------------------------------
ggplot2::ggplot()+
  ggplot2::geom_histogram(aes(x=ranking$PosteriorRank[1]$rank),bins = 5)+
  ggplot2::labs(title = 'Posterior distribution of the rank for Graf', x='Rank')

## -----------------------------------------------------------------------------
tennis_new_games<- tibble::tribble(~player0, ~player1,
                                  'Seles', 'Graf',
                                  'Seles', 'Sabatini',
                                  'Seles', 'Navratilova',
                                  'Seles', 'Sanchez')
y_seles<-predict(m1,tennis_new_games,n=100)
#Now let's summarize the posterior
y_seles <- dplyr::mutate(y_seles, avg_win_player1 = rowMeans(select(y_seles, starts_with('y_pred')))) 
y_seles %>% 
  dplyr::select(player0, player1,avg_win_player1) %>%
  knitr::kable()

