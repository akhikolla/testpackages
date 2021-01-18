## ----setup, include=FALSE, cache=FALSE----------------------------------------
knitr::opts_chunk$set(echo = TRUE, message = FALSE, comment = NA, eval = TRUE)

## -----------------------------------------------------------------------------
library(text.alignment)

## -----------------------------------------------------------------------------
a <- "Gaspard	Tournelly cardeur à laine"
b <- "Gaspard	Bourelly cordonnier"
smith_waterman(a, b)

a <- "Gaspard	T.	cardeur à laine"
b <- "Gaspard	Tournelly cardeur à laine"
smith_waterman(a, b, type = "characters")

## -----------------------------------------------------------------------------
a <- system.file(package = "text.alignment", "extdata", "example1.txt")
a <- readLines(a)
a <- paste(a, collapse = "\n")
b <- system.file(package = "text.alignment", "extdata", "example2.txt")
b <- readLines(b)
b <- paste(b, collapse = "\n")
cat(a, sep = "\n")
cat(b, sep = "\n")

## -----------------------------------------------------------------------------
smith_waterman(a, b, type = "words")

## -----------------------------------------------------------------------------
x <- smith_waterman("Lange rei", b)
x$b$tokens[x$b$alignment$from:x$b$alignment$to]
overview <- as.data.frame(x)
overview$b_from
overview$b_to
substr(overview$b, overview$b_from, overview$b_to)

## -----------------------------------------------------------------------------
x <- smith_waterman(a, b)
x <- as.data.frame(x, alignment_id = "matching-a-to-b")
str(x)

