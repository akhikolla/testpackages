
require(regmed)

set.seed(1000)

data(regmed_example)

 y <- regmed_example$y
 x <- regmed_example$x
 med <- regmed_example[, -c(1,2)]
 fit.grid <- regmed.grid(x, med, y, lambda.vec= c(seq(from=1, to=0, by = -.1)), 
frac.lasso=.8)
fit.grid

## get best fit
fit.trim <- trim.best(fit.grid)
which.med <- colnames(med) %in% dimnames(fit.trim$alpha)[[1]]
med.selected <- med[, which.med]

fit.regmed <- regmed.fit(x, med.selected, y, lambda = 0.2, frac.lasso=.8)
summary(fit.regmed)


