library(Countr)
library(dplyr)

data(football)
table(football$awayTeamGoals)

away_poiss <- glm(formula = awayTeamGoals ~ 1, family = poisson, data = football)
away_wei <- renewalCount(formula = awayTeamGoals ~ 1, data = football,
                         dist = "weibull", weiMethod = "conv_dePril",
                         computeHessian = FALSE, 
                         control = renewal.control(trace = 0,
                                                   method = "nlminb")
                         )

breaks_ <- 0:5
pears <- compareToGLM(poisson_model = away_poiss,
                      breaks = breaks_, weibull = away_wei)

frequency_plot(pears$Counts, pears$Actual,
               dplyr::select(pears, contains("_predicted")),
               colours = c("grey", "blue", "green", "black")
               )

library(lmtest)
lr <- lrtest(away_poiss, away_wei)

lr

gof_wei <- chiSq_gof(away_wei, breaks = breaks_)
gof_pois <- chiSq_gof(away_poiss, breaks = breaks_)
print(gof_wei)

save.image()
