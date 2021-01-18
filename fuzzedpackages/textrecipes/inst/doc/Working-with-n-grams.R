## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(textrecipes)
library(tokenizers)

## -----------------------------------------------------------------------------
abc <- c("The Bank is a place where you put your money;",
         "The Bee is an insect that gathers honey.")

tokenize_words(abc)

## -----------------------------------------------------------------------------
tokenize_ngrams(abc, n = 2)

## -----------------------------------------------------------------------------
tokenize_ngrams(abc, n = 3)

tokenize_ngrams(abc, n = 1)

## -----------------------------------------------------------------------------
tokenize_ngrams(abc, n = 3, ngram_delim = "_")

## -----------------------------------------------------------------------------
abc_tibble <- tibble(text = abc)

rec <- recipe(~ text, data = abc_tibble) %>%
  step_tokenize(text, token = "ngrams") %>%
  step_tokenfilter(text) %>%
  step_tf(text)

abc_ngram <- rec %>%
  prep() %>%
  juice()

abc_ngram

names(abc_ngram)

## -----------------------------------------------------------------------------
abc_tibble <- tibble(text = abc)

rec <- recipe(~ text, data = abc_tibble) %>%
  step_tokenize(text, token = "ngrams", options = list(n = 2, 
                                                       ngram_delim = "_")) %>%
  step_tokenfilter(text) %>%
  step_tf(text)

abc_ngram <- rec %>%
  prep() %>%
  juice()

abc_ngram

names(abc_ngram)

## -----------------------------------------------------------------------------
abc_tibble <- tibble(text = abc)

bigram <- function(x) {
  tokenizers::tokenize_ngrams(x, lowercase = FALSE, n = 2, ngram_delim = ".")
}

rec <- recipe(~ text, data = abc_tibble) %>%
  step_tokenize(text, custom_token = bigram) %>%
  step_tokenfilter(text) %>%
  step_tf(text)

abc_ngram <- rec %>%
  prep() %>%
  juice()

abc_ngram

names(abc_ngram)

## -----------------------------------------------------------------------------
abc_tibble <- tibble(text = abc)

rec <- recipe(~ text, data = abc_tibble) %>%
  step_tokenize(text) %>%
  step_ngram(text, num_tokens = 3) %>%
  step_tokenfilter(text) %>%
  step_tf(text)

abc_ngram <- rec %>%
  prep() %>%
  juice()

abc_ngram

names(abc_ngram)

## -----------------------------------------------------------------------------
abc_tibble <- tibble(text = abc)

rec <- recipe(~ text, data = abc_tibble) %>%
  step_tokenize(text) %>%
  step_stem(text) %>%
  step_ngram(text, num_tokens = 3) %>%
  step_tokenfilter(text) %>%
  step_tf(text)

abc_ngram <- rec %>%
  prep() %>%
  juice()

abc_ngram

names(abc_ngram)

