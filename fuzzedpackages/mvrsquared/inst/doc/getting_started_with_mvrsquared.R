## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>", warning = FALSE
)

## -----------------------------------------------------------------------------

library(mvrsquared)

data(mtcars)

# fit a linear model
f <- lm(mpg ~ cyl + disp + hp + wt, data = mtcars)

# extract r-squared 
f_summary <- summary(f)

r2_lm <- f_summary$r.squared

r2_lm

# calculate univariate r-squared using mvrsquared
r2_mv <- calc_rsquared(y = mtcars$mpg, yhat = f$fitted.values)

r2_mv

# just to be 100% sure...
r2_lm == r2_mv


## -----------------------------------------------------------------------------

x <- cbind(1, f$model[, -1]) # note, you have to add 1's for the intercept and
                             # I'm removing the first column of f$model as it
                             # is the outcome we are predicting

x <- as.matrix(x) # x needs to be a matrix, not a data.frame or tibble for now

w <- matrix(f$coefficients, ncol = 1) # w also has to be a matrix

# this calculates yhat as the dot product x %*% w
r2_mv2 <- calc_rsquared(y = mtcars$mpg, 
                        yhat = list(x = x,
                                    w = w))

r2_mv2


## -----------------------------------------------------------------------------
r2_mv2 == r2_lm

## -----------------------------------------------------------------------------
round(r2_mv2, 14) == round(r2_lm, 14)

## -----------------------------------------------------------------------------
calc_rsquared(y = cbind(mtcars$mpg, mtcars$mpg),
              yhat = cbind(f$fitted.values, f$fitted.values))

## -----------------------------------------------------------------------------
library(nnet)

# let's generate some synthetic data
set.seed(666)

# Some continuous variables
x1 <- rnorm(n = 10000, mean = 1, sd = 2)
x2 <- rnorm(n = 10000, mean = 1.5, sd = 2.5)
x3 <- rnorm(n = 10000, mean = .6, sd = 1.5)

# linear combinations used to generate outcomes with logit functions
z1 <- 1 - 1.2 * x1 + 4.5 * x2 - 2.8 * x3 

z2 <- -2 + 2.8 * x1 - 1.2 * x2 + 4.5 * x3 

y1 <- rbinom(10000, 1, prob = 1 / (1 + exp(-z1)))

y2 <- rbinom(10000, 1, prob = 1 / (1 + exp(-z2)))

# fit a multinomial model using a 1-layer neural net with 10 nodes
f_mv <- nnet(cbind(y1, y2) ~ x1 + x2 + x3, 
             size = 10)

yhat <- predict(f_mv, data.frame(x1 = x1, x2 = x2, x3 = x3), type = "raw")

# and now calculate r-squared
calc_rsquared(y = cbind(y1, y2), yhat = yhat)


## -----------------------------------------------------------------------------
library(tidytext)
library(textmineR)
library(dplyr)
library(stringr)

# load documents in a data frame
docs <- nih_sample 

# tokenize using tidytext's unnest_tokens
tidy_docs <- docs %>% 
  select(APPLICATION_ID, ABSTRACT_TEXT) %>% 
  unnest_tokens(output = word, 
                input = ABSTRACT_TEXT,
                stopwords = stop_words$word,
                token = "ngrams",
                n_min = 1, n = 2) %>% 
  count(APPLICATION_ID, word) %>% 
  filter(n>1) #Filtering for words/bigrams per document, rather than per corpus

tidy_docs <- tidy_docs %>% # filter words that are just numbers
  filter(! str_detect(tidy_docs$word, "^[0-9]+$"))

# turn a tidy tbl into a sparse dgCMatrix for use in textmineR
dtm <- tidy_docs %>% 
  cast_sparse(APPLICATION_ID, word, n)


# create a topic model
lda <- FitLdaModel(dtm = dtm, 
                   k = 20,
                   iterations = 200,
                   burnin = 175)


## -----------------------------------------------------------------------------
r2_lda <- calc_rsquared(y = dtm, 
                        yhat = list(x = rowSums(dtm) * lda$theta, w = lda$phi))

r2_lda


## -----------------------------------------------------------------------------

lsa <- FitLsaModel(dtm = dtm, k = 20)

r2_lsa <- calc_rsquared(y = dtm,
                        yhat = list(x = lsa$theta %*% diag(lsa$sv), w = lsa$phi))

r2_lsa


## ----eval = FALSE-------------------------------------------------------------
#  library(parallel)
#  
#  batch_size <- 10
#  
#  batches <- mclapply(X = seq(1, nrow(dtm), by = batch_size),
#                      FUN = function(b){
#  
#                        # rows to select on
#                        rows <- b:min(b + batch_size - 1, nrow(dtm))
#  
#                        # rows of the dtm
#                        y_batch <- dtm[rows, ]
#  
#                        # rows of theta multiplied by document length
#                        x_batch <- rowSums(y_batch) * lda$theta[rows, ]
#  
#                        list(y = y_batch,
#                             x = x_batch)
#                      }, mc.cores = 2)
#  
#  
#  # calculate ybar for the data
#  # in this case, lazily doing colMeans, but you could divide this problem up too
#  ybar <- colMeans(dtm)
#  
#  # MAP: calculate sums of squares
#  ss <- mclapply(X = batches,
#                 FUN = function(batch){
#                   calc_rsquared(y = batch$y,
#                                 yhat = list(x = batch$x, w = lda$phi),
#                                 ybar = ybar,
#                                 return_ss_only = TRUE)
#                 }, mc.cores = 2)
#  
#  
#  # REDUCE: get SST and SSE by summation
#  ss <- do.call(rbind, ss) %>% colSums()
#  
#  r2_mapreduce <- 1 - ss["sse"] / ss["sst"]
#  
#  # should be the same as above
#  r2_mapreduce
#  

