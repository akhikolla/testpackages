## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(MASS)
library(ranger)
library(tree.interpreter)

## ----reg----------------------------------------------------------------------
# Setup
set.seed(42L)
rfobj <- ranger(medv ~ ., Boston, keep.inbag = TRUE, importance = 'impurity')
tidy.RF <- tidyRF(rfobj, Boston[, -14], Boston[, 14])

# MDI
t(Boston.MDI <- MDI(tidy.RF, Boston[, -14], Boston[, 14]))
all.equal(as.vector(Boston.MDI),
          as.vector(importance(rfobj) /
                      sum(rfobj$inbag.counts[[1]])))

# MDI-oob
t(MDIoob(tidy.RF, Boston[, -14], Boston[, 14]))

## ----class--------------------------------------------------------------------
# Setup
set.seed(42L)
rfobj <- ranger(Species ~ ., iris, keep.inbag = TRUE, importance = 'impurity')
tidy.RF <- tidyRF(rfobj, iris[, -5], iris[, 5])

# MDI
(iris.MDI <- rowSums(MDI(tidy.RF, iris[, -5], iris[, 5])))
all.equal(as.vector(iris.MDI),
          as.vector(importance(rfobj) /
                      sum(rfobj$inbag.counts[[1]])))

# MDI-oob
rowSums(MDIoob(tidy.RF, iris[, -5], iris[, 5]))

