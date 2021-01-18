## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- message=FALSE-----------------------------------------------------------
library(dplyr)
library(recipes)
library(textrecipes)
library(modeldata)
data("okc_text")

## -----------------------------------------------------------------------------
words <- c("you", "i", "sad", "happy")

okc_rec <- recipe(~ ., data = okc_text) %>%
  step_tokenize(essay0) %>%
  step_stopwords(essay0, custom_stopword_source = words, keep = TRUE) %>% 
  step_tf(essay0)

okc_obj <- okc_rec %>%
  prep()
   
bake(okc_obj, okc_text) %>%
  select(starts_with("tf_essay0"))

## -----------------------------------------------------------------------------
stopwords_list <- c("was", "she's", "who", "had", "some", "same", "you", "most", 
                    "it's", "they", "for", "i'll", "which", "shan't", "we're", 
                    "such", "more", "with", "there's", "each")

words <- c("sad", "happy")

okc_rec <- recipe(~ ., data = okc_text) %>%
  step_tokenize(essay0) %>%
  step_stopwords(essay0, custom_stopword_source = stopwords_list) %>% 
  step_stopwords(essay0, custom_stopword_source = words) %>% 
  step_tfidf(essay0)

okc_obj <- okc_rec %>%
  prep()
   
bake(okc_obj, okc_text) %>%
  select(starts_with("tfidf_essay0"))

## -----------------------------------------------------------------------------
okc_rec <- recipe(~ ., data = okc_text) %>%
  step_tokenize(essay0, token = "characters") %>%
  step_stopwords(essay0, custom_stopword_source = letters, keep = TRUE) %>%
  step_tf(essay0)

okc_obj <- okc_rec %>%
  prep()
   
bake(okc_obj, okc_text) %>%
  select(starts_with("tf_essay0"))

## -----------------------------------------------------------------------------
okc_rec <- recipe(~ ., data = okc_text) %>%
  step_tokenize(essay0, token = "words") %>%
  step_stem(essay0) %>%
  step_untokenize(essay0) %>%
  step_tokenize(essay0, token = "ngrams") %>%
  step_tokenfilter(essay0, max_tokens = 500) %>%
  step_tfidf(essay0)

okc_obj <- okc_rec %>%
  prep()
   
bake(okc_obj, okc_text) %>%
  select(starts_with("tfidf_essay0"))

