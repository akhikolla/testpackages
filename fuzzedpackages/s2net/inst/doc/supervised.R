## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = FALSE,
  results = "hold",
  comment = " ",
  warning = FALSE,
  message = FALSE
)

## -----------------------------------------------------------------------------
library(s2net)
data("auto_mpg")

# Preprocess the data using the s2Data function
train = s2Data(auto_mpg$P1$xL, auto_mpg$P1$yL, preprocess = TRUE)

## -----------------------------------------------------------------------------
lm.fit = lm( y~ 0 + ., data = data.frame(train$xL, y = train$yL))

## -----------------------------------------------------------------------------
obj = s2netR(train, s2Params(0))
# We set all the hyper-parameters to 0

## -----------------------------------------------------------------------------
library(Metrics)
# Training error
ypred = predict.lm(lm.fit, data.frame(train$xL))
print("OLS error:")
mse(ypred, train$yL)

ypred = predict(obj, train$xL)
print("s2net error:")
mse(ypred, train$yL)

#Estimations
data.frame(mle = lm.fit$coefficients, s2net = obj$beta)

## -----------------------------------------------------------------------------
library(glmnet)

lasso.fit = glmnet(train$xL, train$yL, family = "gaussian", 
                           alpha = 1, lambda = 0.01, intercept = F)
ypred = predict(lasso.fit, train$xL)
print("Lasso error:")
mse(ypred, train$yL)

obj = s2netR(train, s2Params(lambda1 = 0.01))
ypred = predict(obj, train$xL)
print("s2net error")
mse(ypred, train$yL)

print("Coefficients")
data.frame(lasso = as.numeric(lasso.fit$beta), s2net = obj$beta)

## -----------------------------------------------------------------------------
enet.fit = glmnet(train$xL, train$yL, family = "gaussian", 
                          alpha = 0.3333, lambda = 0.03, intercept = F)
ypred = predict(enet.fit, train$xL)
print("glmnet error")
mse(ypred, train$yL)

obj = s2netR(train, s2Params(lambda1 = 0.01, lambda2 = 0.01))
ypred = predict(obj, train$xL)
print("s2net error")
mse(ypred, train$yL)

print("Coefficients")
data.frame(enet = as.matrix(enet.fit$beta), s2net = obj$beta)

