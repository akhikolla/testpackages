## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval=FALSE---------------------------------------------------------------
#  
#  library(SSLR)
#  library(tidymodels)
#  library(caret)

## ----include=FALSE------------------------------------------------------------
knitr::opts_chunk$set(
  digits = 3,
  collapse = TRUE,
  comment = "#>"
)
options(digits = 3)

library(SSLR)
library(tidymodels)
library(caret)

## ----wine, results="hide"-----------------------------------------------------
data(wine)

set.seed(1)

#Train and test data
train.index <- createDataPartition(wine$Wine, p = .7, list = FALSE)
train <- wine[ train.index,]
test  <- wine[-train.index,]

cls <- which(colnames(wine) == "Wine")

# 20 % LABELED
labeled.index <- createDataPartition(wine$Wine, p = .2, list = FALSE)
train[-labeled.index,cls] <- NA

## ----fitformula, results="hide", eval=FALSE-----------------------------------
#  m <- SSLRDecisionTree() %>% fit(Wine ~ ., data = train)
#  

## ----fitxy, results="hide", eval=FALSE----------------------------------------
#  m <- SSLRDecisionTree() %>% fit_xy(x = train[,-cls], y = train$Wine)
#  

## ----fitxyu, results="hide", eval=FALSE---------------------------------------
#  m <- SSLRDecisionTree() %>% fit_x_u(x = train[labeled.index,-cls], y = train[labeled.index,cls],
#                                       x_U = train[-labeled.index,-cls])
#  

