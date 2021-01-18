library(tlrmvnmvt)
n = 49
a = rep(-10, n)
b = rnorm(n, 5, 2)
m = sqrt(n)
vec1 = 1 : m
vec2 = rep(1, m)
geom = cbind(kronecker(vec1, vec2), kronecker(vec2, vec1))
geom = geom / m
beta = 0.3
distM = as.matrix(dist(geom))
covM = exp(-distM / beta)
val1 = mvtnorm::pmvnorm(a, b, sigma = covM)
val2 = pmvn(lower = a, upper = b, mean = 0, sigma = covM, uselog2 = FALSE, 
            algorithm = GenzBretz(N = 521))
val3 = pmvn(lower = a, upper = b, mean = 0, uselog2 = TRUE, geom = geom, 
            kernelType = "matern", para = c(1.0, 0.3, 0.5, 0.0))
show("Check the (relative) difference with the mvtnorm package (pmvn)")
show(abs(val2 - val1) / val1)
show(- abs(val3 - log2(val1)) / val3)

nu = 10
val1 = mvtnorm::pmvt(lower = a, upper = b, delta = 2, df = nu, 
                     sigma = covM, type = "Kshirsagar")
val2 = pmvt(lower = a, upper = b, delta = 2, df = nu, sigma = covM, 
            uselog2 = FALSE, algorithm = GenzBretz(N = 521), 
            type = "Kshirsagar")
val3 = pmvt(lower = a, upper = b, delta = 2, df = nu, uselog2 = TRUE,
           algorithm = GenzBretz(N = 521), type = "Kshirsagar",  geom = geom, 
           kernelType = "matern", para = c(1.0, 0.3, 0.5, 0.0))
show("Check the (relative) difference with the mvtnorm package (pmvt Kshirsagar)")
show(abs(val2 - val1) / val1)
show(- abs(val3 - log2(val1)) / val3)

val1 = mvtnorm::pmvt(lower = a, upper = b, delta = 2, df = nu, 
                     sigma = covM, type = "shifted")
val2 = pmvt(lower = a, upper = b, delta = 2, df = nu, sigma = covM, 
            uselog2 = FALSE, algorithm = GenzBretz(N = 521), 
            type = "shifted")
val3 = pmvt(lower = a, upper = b, delta = 2, df = nu, uselog2 = TRUE,
           algorithm = GenzBretz(N = 521), type = "shifted",  geom = geom, 
           kernelType = "matern", para = c(1.0, 0.3, 0.5, 0.0))
show("Check the (relative) difference with the mvtnorm package (pmvt shifted)")
show(abs(val2 - val1) / val1)
show(- abs(val3 - log2(val1)) / val3)