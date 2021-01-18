## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(sbo)

## -----------------------------------------------------------------------------
head(sbo::twitter_train, 3)

## -----------------------------------------------------------------------------
p <- sbo_predictor(object = sbo::twitter_train, # preloaded example dataset
                   N = 3, # Train a 3-gram model
                   dict = target ~ 0.75, # cover 75% of training corpus
                   .preprocess = sbo::preprocess, # Preprocessing transformation 
                   EOS = ".?!:;", # End-Of-Sentence tokens
                   lambda = 0.4, # Back-off penalization in SBO algorithm
                   L = 3L, # Number of predictions for input
                   filtered = "<UNK>" # Exclude the <UNK> token from predictions
                   )

## -----------------------------------------------------------------------------
predict(p, "i love")

## -----------------------------------------------------------------------------
set.seed(840)
babble(p)
babble(p)
babble(p)

## -----------------------------------------------------------------------------
t <- sbo_predtable(object = sbo::twitter_train, # preloaded example dataset
                   N = 3, # Train a 3-gram model
                   dict = target ~ 0.75, # cover 75% of training corpus
                   .preprocess = sbo::preprocess, # Preprocessing transformation 
                   EOS = ".?!:;", # End-Of-Sentence tokens
                   lambda = 0.4, # Back-off penalization in SBO algorithm
                   L = 3L, # Number of predictions for input
                   filtered = "<UNK>" # Exclude the <UNK> token from predictions
                   )

## -----------------------------------------------------------------------------
p <- sbo_predictor(t) # This is the same as 'p' created above

## ---- eval=FALSE--------------------------------------------------------------
#  save(t)
#  # ... and, in another session:
#  load("t.rda")

## -----------------------------------------------------------------------------
summary(p)

## -----------------------------------------------------------------------------
head(t[[3]])

## -----------------------------------------------------------------------------
t[[1]]

## ----message=FALSE, warning=FALSE---------------------------------------------
library(dplyr) # installed with `sbo`

## -----------------------------------------------------------------------------
set.seed(840)
(evaluation <- eval_sbo_predictor(p, test = sbo::twitter_test))

## -----------------------------------------------------------------------------
evaluation %>% summarise(accuracy = sum(correct)/n(), 
                   uncertainty = sqrt(accuracy * (1 - accuracy) / n())
                   )

## -----------------------------------------------------------------------------
evaluation %>% # Accuracy for in-sentence predictions
        filter(true != "<EOS>") %>%
        summarise(accuracy = sum(correct) / n(),
                  uncertainty = sqrt(accuracy * (1 - accuracy) / n())
                  )

## ---- fig.align = "center"----------------------------------------------------
if (require(ggplot2)) {
        evaluation %>%
                filter(correct, true != "<EOS>") %>%
                select(true) %>%
                transmute(rank = match(true, table = attr(p, "dict"))) %>%
                ggplot(aes(x = rank)) + geom_histogram(binwidth = 25)
}

## -----------------------------------------------------------------------------
dict <- sbo_dictionary(corpus = sbo::twitter_train, 
                       max_size = 100, 
                       target = 0.5, 
                       .preprocess = sbo::preprocess,
                       EOS = ".?!:;")

## -----------------------------------------------------------------------------
(c <- word_coverage(p, sbo::twitter_train))

## -----------------------------------------------------------------------------
summary(c)

## ---- fig.align = "center", fig.width=5---------------------------------------
plot(c)

## -----------------------------------------------------------------------------
f <- kgram_freqs(corpus = sbo::twitter_train, 
                 N = 3, 
                 dict = target ~ 0.75,
                 .preprocess = sbo::preprocess,
                 EOS = ".?!:;"
                 )

## -----------------------------------------------------------------------------
predict(f, "i love")

## -----------------------------------------------------------------------------
predict(p, "i love")

## -----------------------------------------------------------------------------
size_in_MB <- function(x) format(utils::object.size(x), units = "MB")
sapply(list(sbo_predtable = t, kgram_freqs = f), size_in_MB)

## -----------------------------------------------------------------------------
chrono_predict <- function(x) system.time(predict(x, "i love"), gcFirst = TRUE)
lapply(list(sbo_predictor = p, kgram_freqs = f), chrono_predict)

