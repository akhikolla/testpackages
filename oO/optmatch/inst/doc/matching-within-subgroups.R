## ---- echo = FALSE-------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, prompt=TRUE)
library(optmatch)

## ------------------------------------------------------------------------
data(infert)
head(infert)

## ------------------------------------------------------------------------
table(infert$case)
table(infert$education, infert$case)

## ------------------------------------------------------------------------
f1 <- fullmatch(case ~ age, data = infert[infert$education == "0-5yrs", ])
f2 <- fullmatch(case ~ age, data = infert[infert$education == "6-11yrs", ])
f3 <- fullmatch(case ~ age, data = infert[infert$education == "12+ yrs", ])
summary(f1)
summary(f2)
summary(f3)

## ------------------------------------------------------------------------
f2 <- fullmatch(case ~ age, data = infert[infert$education == "6-11yrs", ],
                max.controls = 4)
f3 <- fullmatch(case ~ age, data = infert[infert$education == "12+ yrs", ],
                max.controls = 4)
summary(f2)
summary(f3)

## ------------------------------------------------------------------------
fcombine <- c(f1, f2, f3)
summary(fcombine)
infert$match <- fcombine

## ------------------------------------------------------------------------
fwithin <- fullmatch(case ~ age, data = infert, max.controls = 4,
                     within = exactMatch(case ~ education, data = infert))
summary(fwithin)

## ---- eval = FALSE-------------------------------------------------------
#  f1 <- fullmatch(z ~ x, data = d[d$group == 1, ], max.controls = 2)
#  f2 <- fullmatch(z ~ x, data = d[d$group == 2, ], min.controls = 1/3)
#  c(f1, f2)

