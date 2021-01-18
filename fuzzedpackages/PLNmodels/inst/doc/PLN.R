## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  screenshot.force = FALSE, 
  echo = TRUE,
  rows.print = 5,
  message = FALSE, 
  warning = FALSE)

## ----requirement--------------------------------------------------------------
library(PLNmodels)
library(ggplot2)
library(corrplot)

## ----data_load----------------------------------------------------------------
data(trichoptera)
trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)

## ----geometricalInsight, echo = FALSE, message = FALSE, warning = FALSE, fig.cap = "PLN: geometrical view", fig.width=7, fig.height=7, fig.align='center'----

library(grid)
library(gridExtra)
library(dplyr)

set.seed(20171110)
x <- rnorm(100)
y <- rnorm(100)
b <- data.frame(x = x + y, y = y) / 1
mu <- 0
##
data.perfect <- as.data.frame((b + matrix(rep(mu, each = length(x)), ncol = 2)))
p.latent <- ggplot(data.perfect, aes(x, y)) + geom_point() + ggtitle(expression(Latent~Space~(Z)))
.rpois <- function(lambda) {
  unlist(lapply(exp(lambda), function(x) {rpois(1, x)}))
}
observation <- as.data.frame(lapply(data.perfect, .rpois))
mapped.parameter <- as.data.frame(lapply(data.perfect, exp))
## segment between mapped and observed data
segment.data <- cbind(mapped.parameter, observation)
names(segment.data) <- c("x", "y", "xend", "yend")
## Mapped parameters
p.mapped <- ggplot(mapped.parameter, aes(x, y)) + geom_point(col = "red") + ggtitle(expression(Observation~Space~(exp(Z))))
## Observations only
obs <- group_by(observation, x, y)
obs <- dplyr::summarize(obs, count = n())
p.observation.only <- ggplot(obs, aes(x, y)) +
  geom_point(aes(size = count)) +
  ggtitle(Observation~Space~(Y)~+'noise') +
  theme(legend.position = c(.95, .95), legend.justification = c(1, 1),
        legend.background = element_rect(color = "transparent"),
        legend.box.background = element_blank())
## Observations and latent parameters
p.observation.mixed <- p.observation.only +
  geom_point(data = mapped.parameter, color = "red", alpha = 0.5) +
  geom_segment(data = segment.data, aes(xend = xend, yend = yend), color = "black", alpha = 0.2) +
  ggtitle(Observation~Space~(Y==P(exp(Z)))~+'noise')

grid.arrange(p.latent + labs(x = "species 1", y = "species 2"),
             p.mapped  + labs(x = "species 1", y = "species 2"),
             p.observation.mixed + labs(x = "species 1", y = "species 2"),
             p.observation.only + labs(x = "species 1", y = "species 2"),
             ncol = 2)

## ----simple PLN---------------------------------------------------------------
myPLN <- PLN(Abundance ~ 1, trichoptera)

## ----show-method--------------------------------------------------------------
myPLN

## ----fields-access------------------------------------------------------------
c(myPLN$loglik, myPLN$BIC, myPLN$ICL, myPLN$R_squared)
myPLN$criteria

## ----fitted, fig.cap = "fitted value vs. observation", fig.dim=c(7,5)---------
data.frame(
  fitted   = as.vector(fitted(myPLN)),
  observed = as.vector(trichoptera$Abundance)
) %>% 
  ggplot(aes(x = observed, y = fitted)) + 
    geom_point(size = .5, alpha =.25 ) + 
    scale_x_log10() + 
    scale_y_log10() + 
    theme_bw() + annotation_logticks()

## ----coef---------------------------------------------------------------------
data.frame(
  rbind(t(coef(myPLN)), t(standard_error(myPLN))), 
  row.names = c("effect", "stderr")
 ) %>% select(1:5) %>% knitr::kable()

## ----plot covariance, fig.width=7, fig.height=5-------------------------------
corrplot(sigma(myPLN), is.corr = FALSE)

## ----weighted, fig.width=7, fig.height=5--------------------------------------
myPLN_weighted <-
  PLN(
    Abundance ~ 1,
    data    = trichoptera,
    weights = runif(nrow(trichoptera)),
    control = list(trace = 0)
  )
data.frame(
  unweighted = as.vector(fitted(myPLN)),
  weighted   = as.vector(fitted(myPLN_weighted))
) %>%
  ggplot(aes(x = unweighted, y = weighted)) +
    geom_point(size = .5, alpha =.25 ) +
    scale_x_log10() +
    scale_y_log10() +
    theme_bw() + annotation_logticks()

## ----PLN offset---------------------------------------------------------------
myPLN_offsets <- 
  PLN(Abundance ~ 1 + offset(log(Offset)), 
      data = trichoptera, control = list(trace = 0))

## ----compare w/wo offset------------------------------------------------------
rbind(
  myPLN$criteria,
  myPLN_offsets$criteria
) %>% knitr::kable()

## ----PLN wind-----------------------------------------------------------------
myPLN_wind <- PLN(Abundance ~ 1 + Wind + offset(log(Offset)), data = trichoptera)

## ----compare models-----------------------------------------------------------
rbind(
  myPLN_offsets$criteria,
  myPLN_wind$criteria
) %>% knitr::kable()

## ----covariances models spherical---------------------------------------------
myPLN_spherical <-
  PLN(
    Abundance ~ 1 + offset(log(Offset)),
    data = trichoptera, control = list(covariance = "spherical", trace = 0)
  )

## ----covariances model diagonal-----------------------------------------------
myPLN_diagonal <-
  PLN(
    Abundance ~ 1 + offset(log(Offset)),
    data = trichoptera, control = list(covariance = "diagonal", trace = 0)
  )

## ----PLN covariance full, evaluate = FALSE------------------------------------
myPLN_default <-
  PLN(Abundance ~ 1, data = trichoptera, )
myPLN_full <-
  PLN(Abundance ~ 1, data = trichoptera, control = list(covariance = "full"))

## ----compare covariances------------------------------------------------------
rbind(
  myPLN_offsets$criteria,
  myPLN_diagonal$criteria,
  myPLN_spherical$criteria
) %>%
  as.data.frame(row.names = c("full", "diagonal", "spherical")) %>%
  knitr::kable()

## ----final--------------------------------------------------------------------
myPLN_final <-
  PLN(
    Abundance ~ 1 + Wind + offset(log(Offset)),
    data    = trichoptera, control = list(covariance = "diagonal", trace = 0)
  )
rbind(
  myPLN_wind$criteria,
  myPLN_diagonal$criteria,
  myPLN_final$criteria
) %>% knitr::kable()

