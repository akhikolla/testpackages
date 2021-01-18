## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  screenshot.force = FALSE, 
  echo = TRUE,
  rows.print = 5,
  message = FALSE, 
  warning = FALSE)

## ----requirement, cache = FALSE-----------------------------------------------
library(PLNmodels)

## ----data_load----------------------------------------------------------------
data(trichoptera)
trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)

## ----description, echo = FALSE------------------------------------------------
data.frame(Label = 1:12, 
       `Number of Consecutive Nights` = c(12, 5, 5, 4, 4, 1, 3, 4, 5, 4, 1, 1), 
       Date  = paste(rep(c("June", "July", "June", "July"), times = c(4, 1, 6, 1), sep = " "), 
                     rep(c(59, 60), times = c(6, 6)),
                     sep = " ")) %>% 
  knitr::kable(align = "c")

## ----LDA-nocov----------------------------------------------------------------
myLDA_nocov <- PLNLDA(Abundance ~ 0 + offset(log(Offset)),
                      grouping = Group, 
                      data = trichoptera)

## ----show nocov---------------------------------------------------------------
myLDA_nocov

## ----vcov---------------------------------------------------------------------
sigma(myLDA_nocov) %>% corrplot::corrplot(is.corr = FALSE)

## ----coef---------------------------------------------------------------------
coef(myLDA_nocov)

## ----group-means--------------------------------------------------------------
myLDA_nocov$group_means %>% head() %>% knitr::kable(digits = 2)

## ----plot_model, fig.width=7, fig.height=5------------------------------------
plot(myLDA_nocov)

## ----extract-scores-----------------------------------------------------------
myLDA_nocov$scores %>% head %>% knitr::kable(digits = 2)

## ----extract-corr-------------------------------------------------------------
myLDA_nocov$corr_map %>% head %>% knitr::kable(digits = 2)

## ----predict_class_posterior--------------------------------------------------
predicted.class <- predict(myLDA_nocov, newdata = trichoptera)
## equivalent to 
## predicted.class <- predict(myLDA_nocov, newdata = trichoptera,  type = "posterior")
predicted.class %>% head() %>% knitr::kable(digits = 2)

## ----predict_class_posterior_prob---------------------------------------------
predicted.class <- predict(myLDA_nocov, newdata = trichoptera, scale = "prob")
predicted.class %>% head() %>% knitr::kable(digits = 3)

## ----predict_class------------------------------------------------------------
predicted.class <- predict(myLDA_nocov, newdata = trichoptera,  type = "response")
predicted.class

## ----check_predicted_class----------------------------------------------------
table(predicted.class, trichoptera$Group, dnn = c("predicted", "true"))

## ----predicted_scores, fig.width=7, fig.height=5------------------------------
library(ggplot2)
predicted.scores <- predict(myLDA_nocov, newdata = trichoptera,  type = "scores")
colnames(predicted.scores) <- paste0("Axis.", 1:ncol(predicted.scores))
predicted.scores <- as.data.frame(predicted.scores)
predicted.scores$group <- trichoptera$Group
plot(myLDA_nocov, map = "individual", nb_axes = 2, plot = FALSE) + 
  geom_point(data = predicted.scores, 
             aes(x = Axis.1, y = Axis.2, color = group, label = NULL))

## ---- warning=FALSE-----------------------------------------------------------
myLDA_cov <- PLNLDA(Abundance ~ Wind + 0 + offset(log(Offset)), 
                    grouping = Group, 
                    data = trichoptera)

## ----coef-cov-----------------------------------------------------------------
coef(myLDA_cov) %>% head %>% knitr::kable()

## ----group-means-cov----------------------------------------------------------
myLDA_cov$group_means %>% head %>% knitr::kable(digits = 2)

## ----plot_model_cov, fig.width=7, fig.height=5--------------------------------
plot(myLDA_cov)

## ----predict_class_cov--------------------------------------------------------
predicted.class_cov <- predict(myLDA_cov, newdata = trichoptera, type = "response")

## ----check_predicted_class_cov------------------------------------------------
table(predicted.class_cov, trichoptera$Group, dnn = c("predicted", "true"))

