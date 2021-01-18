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

## ----packages-1, message=FALSE------------------------------------------------
library(rstanarm)
library(ggplot2)
library(bayesplot)
theme_set(bayesplot::theme_default())
# options(mc.cores = 4) 

## ----packages-2, eval=FALSE, message=FALSE------------------------------------
#  library(dplyr)
#  library(tidyr)

## ---- include=FALSE, collapse=TRUE--------------------------------------------
simulate_mrp_data <- function(n) {
  J <- c(2, 3, 7, 3, 50) # male or not, eth, age, income level, state
  poststrat <- as.data.frame(array(NA, c(prod(J), length(J)+1))) # Columns of post-strat matrix, plus one for size
  colnames(poststrat) <- c("male", "eth", "age","income", "state",'N')
  count <- 0
  for (i1 in 1:J[1]){
    for (i2 in 1:J[2]){
      for (i3 in 1:J[3]){
        for (i4 in 1:J[4]){
          for (i5 in 1:J[5]){
              count <- count + 1
              # Fill them in so we know what category we are referring to
              poststrat[count, 1:5] <- c(i1-1, i2, i3,i4,i5) 
          }
        }
      }
    }
  }
  # Proportion in each sample in the population
  p_male <- c(0.52, 0.48)
  p_eth <- c(0.5, 0.2, 0.3)
  p_age <- c(0.2,.1,0.2,0.2, 0.10, 0.1, 0.1)
  p_income<-c(.50,.35,.15)
  p_state_tmp<-runif(50,10,20)
  p_state<-p_state_tmp/sum(p_state_tmp)
  poststrat$N<-0
  for (j in 1:prod(J)){
    poststrat$N[j] <- round(250e6 * p_male[poststrat[j,1]+1] * p_eth[poststrat[j,2]] *
      p_age[poststrat[j,3]]*p_income[poststrat[j,4]]*p_state[poststrat[j,5]]) #Adjust the N to be the number observed in each category in each group
  }
  
  # Now let's adjust for the probability of response
  p_response_baseline <- 0.01
  p_response_male <- c(2, 0.8) / 2.8
  p_response_eth <- c(1, 1.2, 2.5) / 4.7
  p_response_age <- c(1, 0.4, 1, 1.5,  3, 5, 7) / 18.9
  p_response_inc <- c(1, 0.9, 0.8) / 2.7
  p_response_state <- rbeta(50, 1, 1)
  p_response_state <- p_response_state / sum(p_response_state)
  p_response <- rep(NA, prod(J))
  for (j in 1:prod(J)) {
    p_response[j] <-
      p_response_baseline * p_response_male[poststrat[j, 1] + 1] *
      p_response_eth[poststrat[j, 2]] * p_response_age[poststrat[j, 3]] *
      p_response_inc[poststrat[j, 4]] * p_response_state[poststrat[j, 5]]
  }
  people <- sample(prod(J), n, replace = TRUE, prob = poststrat$N * p_response)
  
  ## For respondent i, people[i] is that person's poststrat cell,
  ## some number between 1 and 32
  n_cell <- rep(NA, prod(J))
  for (j in 1:prod(J)) {
    n_cell[j] <- sum(people == j)
  }
  
  coef_male <- c(0,-0.3)
  coef_eth <- c(0, 0.6, 0.9)
  coef_age <- c(0,-0.2,-0.3, 0.4, 0.5, 0.7, 0.8, 0.9)
  coef_income <- c(0,-0.2, 0.6)
  coef_state <- c(0, round(rnorm(49, 0, 1), 1))
  coef_age_male <- t(cbind(c(0, .1, .23, .3, .43, .5, .6),
                           c(0, -.1, -.23, -.5, -.43, -.5, -.6)))
  true_popn <- data.frame(poststrat[, 1:5], cat_pref = rep(NA, prod(J)))
  for (j in 1:prod(J)) {
    true_popn$cat_pref[j] <- plogis(
      coef_male[poststrat[j, 1] + 1] +
        coef_eth[poststrat[j, 2]] + coef_age[poststrat[j, 3]] +
        coef_income[poststrat[j, 4]] + coef_state[poststrat[j, 5]] +
        coef_age_male[poststrat[j, 1] + 1, poststrat[j, 3]]
      )
  }
  
  #male or not, eth, age, income level, state, city
  y <- rbinom(n, 1, true_popn$cat_pref[people])
  male <- poststrat[people, 1]
  eth <- poststrat[people, 2]
  age <- poststrat[people, 3]
  income <- poststrat[people, 4]
  state <- poststrat[people, 5]
  
  sample <- data.frame(cat_pref = y, 
                       male, age, eth, income, state, 
                       id = 1:length(people))
  
  #Make all numeric:
  for (i in 1:ncol(poststrat)) {
    poststrat[, i] <- as.numeric(poststrat[, i])
  }
  for (i in 1:ncol(true_popn)) {
    true_popn[, i] <- as.numeric(true_popn[, i])
  }
  for (i in 1:ncol(sample)) {
    sample[, i] <- as.numeric(sample[, i])
  }
  list(
    sample = sample,
    poststrat = poststrat,
    true_popn = true_popn
  )
}

## ----include=FALSE, eval=FALSE------------------------------------------------
#  mrp_sim <- simulate_mrp_data(n=1200)
#  save(mrp_sim, file = "mrp-files/mrp_sim.rda", version = 2)

## ----eval=FALSE---------------------------------------------------------------
#  mrp_sim <- simulate_mrp_data(n=1200)
#  str(mrp_sim)

## ---- echo=FALSE--------------------------------------------------------------
load("mrp-files/mrp_sim.rda")
str(mrp_sim)

## ---- message=FALSE-----------------------------------------------------------
sample <- mrp_sim[["sample"]]
rbind(head(sample), tail(sample))

## ----message=FALSE------------------------------------------------------------
poststrat <- mrp_sim[["poststrat"]]
rbind(head(poststrat), tail(poststrat))

## ----message=FALSE------------------------------------------------------------
true_popn <- mrp_sim[["true_popn"]]
rbind(head(true_popn), tail(true_popn))

## ----order-states-------------------------------------------------------------
sample$state <- factor(sample$state, levels=1:50)
sample$state <- with(sample, factor(state, levels=order(table(state))))
true_popn$state <- factor(true_popn$state,levels = levels(sample$state))
poststrat$state <- factor(poststrat$state,levels = levels(sample$state))

## ----state-and-pop-data-for-plots, eval=FALSE, include=FALSE------------------
#  # not evaluated to avoid tidyverse dependency
#  income_popn <- poststrat %>%
#    group_by(income) %>%
#    summarize(Num=sum(N)) %>%
#    mutate(PROP=Num/sum(Num),TYPE='Popn',VAR='Income',CAT=income) %>%
#    ungroup()
#  income_data <- sample %>%
#    group_by(income) %>%
#    summarise(Num=n()) %>%
#    mutate(PROP=Num/sum(Num),TYPE='Sample',VAR='Income',CAT=income) %>%
#    ungroup()
#  income<-rbind(income_data[,2:6],income_popn[,2:6])
#  
#  age_popn <- poststrat%>%
#    group_by(age)%>%
#    summarize(Num=sum(N))%>%
#    mutate(PROP=Num/sum(Num),TYPE='Popn',VAR='Age',CAT=age)%>%
#    ungroup()
#  age_data <- sample%>%
#    group_by(age)%>%
#    summarise(Num=n())%>%
#    mutate(PROP=Num/sum(Num),TYPE='Sample',VAR='Age',CAT=age)%>%
#    ungroup()
#  age <- rbind(age_data[,2:6],age_popn[,2:6] )
#  
#  eth_popn <- poststrat%>%
#    group_by(eth)%>%
#    summarize(Num=sum(N))%>%
#    mutate(PROP=Num/sum(Num),TYPE='Popn',VAR='Ethnicity',CAT=eth)%>%
#    ungroup()
#  eth_data <- sample%>%
#    group_by(eth)%>%
#    summarise(Num=n())%>%
#    mutate(PROP=Num/sum(Num),TYPE='Sample',VAR='Ethnicity',CAT=eth)%>%
#    ungroup()
#  eth<-rbind(eth_data[,2:6],eth_popn[,2:6])
#  
#  male_popn <- poststrat%>%
#    group_by(male)%>%
#    summarize(Num=sum(N))%>%
#    mutate(PROP=Num/sum(Num),TYPE='Popn',VAR='Male',CAT=male)%>%
#    ungroup()
#  male_data <- sample%>%
#    group_by(male)%>%
#    summarise(Num=n())%>%
#    mutate(PROP=Num/sum(Num),TYPE='Sample',VAR='Male',CAT=male)%>%
#    ungroup()
#  male <- rbind(male_data[,2:6],male_popn[,2:6])
#  
#  state_popn <- poststrat%>%
#    group_by(state)%>%
#    summarize(Num=sum(N))%>%
#    mutate(PROP=Num/sum(poststrat$N),TYPE='Popn',VAR='State',CAT=state)%>%
#    ungroup()
#  
#  state_plot_data <- sample%>%
#    group_by(state)%>%
#    summarise(Num=n())%>%
#    mutate(PROP=Num/nrow(sample),TYPE='Sample',VAR='State',CAT=state)%>%
#    ungroup()
#  
#  state_plot_data <- rbind(state_plot_data[,2:6],state_popn[,2:6])
#  state_plot_data$TYPE <- factor(state_plot_data$TYPE, levels = c("Sample","Popn"))
#  
#  plot_data <- rbind(male,eth,age,income)
#  plot_data$TYPE <- factor(plot_data$TYPE, levels = c("Sample","Popn"))
#  
#  save(state_plot_data, file = "mrp-files/state_plot_data.rda", version = 2)
#  save(plot_data, file = "mrp-files/plot_data.rda", version = 2)

## ----plot-data, echo=FALSE, fig.height = 4, fig.width = 7, fig.align = "center"----
load("mrp-files/plot_data.rda") # created in previous chunk
ggplot(data=plot_data, aes(x=as.factor(CAT), y=PROP, group=as.factor(TYPE), linetype=as.factor(TYPE))) +
  geom_point(stat="identity",colour='black')+
  geom_line()+
  facet_wrap( ~ VAR, scales = "free",nrow=1,ncol=5)+
  theme_bw()+
  scale_fill_manual(values=c('#1f78b4','#33a02c',
                             '#e31a1c','#ff7f00','#8856a7'),guide=FALSE)+
  scale_y_continuous(breaks=c(0,.25,.5,.75,1), labels=c('0%','25%',"50%","75%","100%"))+
  scale_alpha_manual(values=c(1, .3))+
  ylab('Proportion')+
  labs(alpha='')+
  theme(legend.position="bottom",
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        legend.title=element_blank(),
        legend.text=element_text(size=10),
        axis.text=element_text(size=10),
        strip.text=element_text(size=10),
        strip.background = element_rect(fill='grey92'))

load("mrp-files/state_plot_data.rda") # created in previous chunk
ggplot(data=state_plot_data, aes(x=as.factor(CAT), y=PROP, group=as.factor(TYPE),    linetype=as.factor(TYPE))) +
  geom_point(stat="identity",colour='black')+
  geom_line()+
  facet_wrap( ~ VAR)+
  theme_bw()+
  scale_fill_manual(values=c('#1f78b4','#33a02c',
                             '#e31a1c','#ff7f00','#8856a7'),guide=FALSE)+
  scale_y_continuous(breaks=c(0,.025,.05,1), labels=c('0%','2.5%',"5%","100%"),expand=c(0,0),limits=c(0,.06))+
  scale_alpha_manual(values=c(1, .3))+
  ylab('Proportion')+
  labs(alpha='')+
  theme(legend.position="bottom",
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        legend.title=element_blank(),
        legend.text=element_text(size=10),
        axis.text.y=element_text(size=10),
        axis.text.x=element_text(size=8,angle=90),
        strip.text=element_text(size=10),
        strip.background = element_rect(fill='grey92'))

## ---- eval=FALSE, echo=FALSE--------------------------------------------------
#  # not evaluated to avoid dependency on tidyverse
#  
#  #Summarise
#  summary_by_poststrat_var <- sample %>%
#    gather(variable,category,c("income","eth","age","male")) %>%
#    group_by(variable,category) %>%
#    #Wald confidence interval
#    summarise(y_mean=mean(cat_pref),y_sd=sqrt(mean(cat_pref)*(1-mean(cat_pref))/n())) %>%
#    ungroup()
#  summary_by_poststrat_var$variable <- as.factor(summary_by_poststrat_var$variable)
#  levels(summary_by_poststrat_var$variable) <- list('Age'='age','Ethnicity'='eth','Income'='income','Male'='male')
#  
#  save(summary_by_poststrat_var, file = "mrp-files/summary_by_poststrat_var.rda",
#       version = 2)

## ----plot-summary-by-poststrat-var, echo=FALSE, fig.height = 4, fig.width = 7, fig.align = "center"----
load("mrp-files/summary_by_poststrat_var.rda") # created in previous chunk
ggplot(data=summary_by_poststrat_var, aes(x=as.factor(category), y=y_mean,group=1)) +
  geom_errorbar(aes(ymin=y_mean-y_sd, ymax=y_mean+y_sd), width=0)+
  geom_line()+
  geom_point()+
  scale_colour_manual(values=c('#1f78b4','#33a02c','#e31a1c','#ff7f00',
                             '#8856a7'))+theme_bw()+
facet_wrap(~variable,scales = "free_x",nrow=1,ncol=5)+
    scale_y_continuous(breaks=c(.5,.75,1), labels=c("50%","75%",
                                        "100%"), limits=c(0.4-.4*.05,.9),expand = c(0,0))+
  labs(x="",y="Cat preference")+
  theme(legend.position="none",
        axis.title.y=element_text(size=10),
        axis.title.x=element_blank(),
        axis.text=element_text(size=10),
        strip.text=element_text(size=10),
        strip.background = element_rect(fill='grey92'))

## ----interaction-summary, eval=FALSE, echo=FALSE------------------------------
#  # not evaluated to avoid dependency on tidyverse
#  
#  #Summarise
#  interaction <- sample %>%
#    gather(variable, category, c("age", "eth")) %>%
#    group_by(variable, category, male) %>%
#    summarise(y_mean = mean(cat_pref),
#              y_sd = sqrt(mean(cat_pref) * (1 - mean(cat_pref)) / n())) %>%
#    ungroup()
#  
#  #Tidy for nice facet labels
#  interaction$variable <- as.factor(interaction$variable)
#  levels(interaction$variable) <- list('Ethnicity' = 'eth', 'Age' = 'age')
#  save(interaction, file = "mrp-files/interaction.rda", version = 2)

## ----plot-interaction, echo=FALSE, fig.height = 4, fig.width = 7, fig.align = "center"----
load("mrp-files/interaction.rda") # created in previous chunk
ggplot(data=interaction, aes(x=as.factor(category), y=y_mean, colour=as.factor(male),group=as.factor(male))) +
  geom_errorbar(aes(ymin=y_mean-y_sd, ymax=y_mean+y_sd),width=0 )+
  geom_line(aes(x=as.factor(category), y=y_mean,colour=as.factor(male)))+
  geom_point()+
  facet_wrap(~variable,scales = "free_x",nrow=1,ncol=2)+
  labs(x="",y="Cat preference",colour='Gender')+
  scale_y_continuous(breaks=c(0,.25,.5,.75,1), labels=c("0%",'25%',"50%","75%",
                                        "100%"), limits=c(0,1),expand=c(0,0))+
  scale_colour_manual(values=c('#4575b4','#d73027'))+theme_bw()+
  theme(axis.title=element_text(size=10),
        axis.text=element_text(size=10),
        legend.position='none',
        strip.text=element_text(size=10),
        strip.background = element_rect(fill='grey92'))


## ---- eval=FALSE, echo=FALSE--------------------------------------------------
#  # not evaluated to avoid dependency on tidyverse
#  
#  #Summarise by state
#  preference_by_state <- sample %>%
#    group_by(state) %>%
#    summarise(y_mean = mean(cat_pref),
#              y_sd = sqrt(mean(cat_pref) * (1 - mean(cat_pref)) / n())) %>%
#    ungroup()
#  
#  save(preference_by_state, file = "mrp-files/preference_by_state.rda", version = 2)

## ---- echo=FALSE, fig.height = 4, fig.width = 8, fig.align = "center"---------
load("mrp-files/preference_by_state.rda")
compare <- ggplot(data=preference_by_state, aes(x=state, y=y_mean,group=1)) +
  geom_ribbon(aes(ymin=y_mean-y_sd,ymax=y_mean+y_sd,x=state),fill='lightgrey',alpha=.7)+
  geom_line(aes(x=state, y=y_mean))+
  geom_point()+
  scale_y_continuous(breaks=c(0,.25,.5,.75,1), 
                     labels=c("0%","25%","50%","75%","100%"), 
                     limits=c(0,1), expand=c(0,0))+
  scale_x_discrete(drop=FALSE)+
  scale_colour_manual(values=c('#1f78b4','#33a02c','#e31a1c','#ff7f00',
                               '#8856a7'))+
  theme_bw()+
  labs(x="States",y="Cat preference")+
  theme(legend.position="none",
        axis.title=element_text(size=10),
        axis.text.y=element_text(size=10),
        axis.text.x=element_text(angle=90,size=8),
        legend.title=element_text(size=10),
        legend.text=element_text(size=10))

compare2 <- ggplot()+
  geom_hline(yintercept = mean(sample$cat_pref),size=.8)+
  geom_text(aes(x = 5.2, y = mean(sample$cat_pref)+.025, label = "Sample"))+
  scale_y_continuous(breaks=c(0,.25,.5,.75,1), 
                     labels=c("0%","25%","50%","75%","100%"),
                     limits=c(-0.25,1.25),expand=c(0,0))+
  theme_bw()+
  labs(x="Popn",y="")+
   theme(legend.position="none",
        axis.title.y=element_blank(),
        axis.title.x=element_text(size=10),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        legend.title=element_text(size=10),
        legend.text=element_text(size=10))

bayesplot_grid(compare,compare2, 
               grid_args = list(nrow=1, widths = c(8,1)))

## ---- message=FALSE, warning=FALSE, results='hide'----------------------------
fit <- stan_glmer(
  cat_pref ~ factor(male) + factor(male) * factor(age) + 
    (1 | state) + (1 | age) + (1 | eth) + (1 | income),
  family = binomial(link = "logit"),
  data = sample
)

## -----------------------------------------------------------------------------
print(fit)

## ---- message=FALSE-----------------------------------------------------------
posterior_prob <- posterior_linpred(fit, transform = TRUE, newdata = poststrat)
poststrat_prob <- posterior_prob %*% poststrat$N / sum(poststrat$N)
model_popn_pref <- c(mean = mean(poststrat_prob), sd = sd(poststrat_prob))
round(model_popn_pref, 3)

## ---- message=FALSE-----------------------------------------------------------
sample_popn_pref <- mean(sample$cat_pref)
round(sample_popn_pref, 3)

## ---- message=FALSE,fig.height = 4, fig.width = 8, fig.align = "center"-------
compare2 <- compare2 +
  geom_hline(yintercept = model_popn_pref[1], colour = '#2ca25f', size = 1) +
  geom_text(aes(x = 5.2, y = model_popn_pref[1] + .025), label = "MRP", colour = '#2ca25f')
bayesplot_grid(compare, compare2, 
               grid_args = list(nrow = 1, widths = c(8, 1)))

## ---- message=FALSE-----------------------------------------------------------
true_popn_pref <- sum(true_popn$cat_pref * poststrat$N) / sum(poststrat$N)
round(true_popn_pref, 3)

## ---- echo=FALSE, message=FALSE,fig.height = 4, fig.width = 8, fig.align = "center"----
compare2 <- compare2 +
  geom_hline(yintercept = mean(true_popn_pref), linetype = 'dashed', size = .8) +
  geom_text(aes(x = 5.2, y = mean(true_popn_pref) - .025), label = "True")
bayesplot_grid(compare, compare2, 
               grid_args = list(nrow = 1, widths = c(8, 1)))

## ---- message=FALSE-----------------------------------------------------------
state_df <- data.frame(
  State = 1:50,
  model_state_sd = rep(-1, 50),
  model_state_pref = rep(-1, 50),
  sample_state_pref = rep(-1, 50),
  true_state_pref = rep(-1, 50),
  N = rep(-1, 50)
)

for(i in 1:length(levels(as.factor(poststrat$state)))) {
  poststrat_state <- poststrat[poststrat$state == i, ]
    posterior_prob_state <- posterior_linpred(
    fit,
    transform = TRUE,
    draws = 1000,
    newdata = as.data.frame(poststrat_state)
  )
  poststrat_prob_state <- (posterior_prob_state %*% poststrat_state$N) / sum(poststrat_state$N)
  #This is the estimate for popn in state:
  state_df$model_state_pref[i] <- round(mean(poststrat_prob_state), 4)
  state_df$model_state_sd[i] <- round(sd(poststrat_prob_state), 4)
  #This is the estimate for sample
  state_df$sample_state_pref[i] <- round(mean(sample$cat_pref[sample$state == i]), 4)
  #And what is the actual popn?
  state_df$true_state_pref[i] <-
    round(sum(true_popn$cat_pref[true_popn$state == i] * poststrat_state$N) /
            sum(poststrat_state$N), digits = 4)
  state_df$N[i] <- length(sample$cat_pref[sample$state == i])
}

state_df[c(1,3:6)]
state_df$State <- factor(state_df$State, levels = levels(sample$state))

## -----------------------------------------------------------------------------
round(100 * c(
  mean = mean(abs(state_df$sample_state_pref-state_df$true_state_pref), na.rm = TRUE),
  max = max(abs(state_df$sample_state_pref-state_df$true_state_pref), na.rm = TRUE)
))

## -----------------------------------------------------------------------------
round(100 * c(
  mean = mean(abs(state_df$model_state_pref-state_df$true_state_pref)),
  max = max(abs(state_df$model_state_pref-state_df$true_state_pref))
))

## ---- message=FALSE, echo=FALSE, fig.height = 4, fig.width = 8, fig.align = "center",warning=FALSE, fig.align = "center"----
#Summarise by state
compare <- compare +
  geom_point(data=state_df, mapping=aes(x=State, y=model_state_pref),
             inherit.aes=TRUE,colour='#238b45')+
  geom_line(data=state_df, mapping=aes(x=State, y=model_state_pref,group=1),
            inherit.aes=TRUE,colour='#238b45')+
  geom_ribbon(data=state_df,mapping=aes(x=State,ymin=model_state_pref-model_state_sd,
                                        ymax=model_state_pref+model_state_sd,group=1), 
              inherit.aes=FALSE,fill='#2ca25f',alpha=.3)+
  geom_point(data=state_df, mapping=aes(x=State, y=true_state_pref),
             alpha=.5,inherit.aes=TRUE)+
  geom_line(data=state_df, mapping=aes(x=State, y=true_state_pref),
            inherit.aes = TRUE,linetype='dashed')

bayesplot_grid(compare, compare2, 
               grid_args = list(nrow = 1, widths = c(8, 1)))

## ---- eval=FALSE--------------------------------------------------------------
#  # not evaluated to avoid dependency on tidyverse
#  sample_alt <- sample %>%
#    group_by(male, age, income, state, eth) %>%
#    summarise(N_cat_pref = sum(cat_pref), N = n()) %>%
#    ungroup()

## ---- include=FALSE-----------------------------------------------------------
load("mrp-files/sample_alt.rda")

## ---- message=FALSE, warning=FALSE, results='hide'----------------------------
fit2 <- stan_glmer(
  cbind(N_cat_pref, N - N_cat_pref) ~ factor(male) + factor(male) * factor(age) + 
    (1 | state) + (1 | age) + (1 | eth) + (1 | income),
  family = binomial("logit"),
  data = sample_alt,
  refresh = 0
)

## -----------------------------------------------------------------------------
print(fit2)

## ---- message=FALSE-----------------------------------------------------------
posterior_prob_alt <- posterior_linpred(fit2, transform = TRUE, newdata = poststrat)
poststrat_prob_alt <- posterior_prob_alt %*% poststrat$N / sum(poststrat$N)
model_popn_pref_alt <- c(mean = mean(poststrat_prob_alt), sd = sd(poststrat_prob_alt))
round(model_popn_pref_alt, 3)

## -----------------------------------------------------------------------------
print(simulate_mrp_data)

