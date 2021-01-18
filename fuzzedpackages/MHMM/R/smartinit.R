# oneEM.gaussian <- function(my_x, my_w, g){
#   prod_xw <- my_x * my_w
#   n <- length(my_x)
#   mu <- sample(my_x, g, prob = my_w/sum(my_w))
#   tmpmean <- sum(prod_xw) / sum(my_w)
#   vartotal <- sum( ((my_x - tmpmean)**2) * my_w) / sum(my_w)
#   varinter <- mean( (tmpmean - mu)**2)/sum(my_w)
#   s <- rep(sqrt(vartotal - varinter), g)
#   logprobs <- sapply(1:g, function(k) log(1/g) + dnorm(my_x, mean = mu[k], sd = s[k], log = TRUE))
#   maxprobs <- apply(logprobs, 1, max)
#   logdensity <- maxprobs + (log(rowSums(exp(sweep(logprobs, 1, maxprobs, "-")))))
#   loglike <- sum(my_w * logdensity)
#   repeat{
#     # Estep
#     tik <- exp(sweep(logprobs, 1, logdensity, "-"))
#     # Mstep
#     norm <- sapply(1:g, function(k) sum(my_w * tik[,k]))
#     mu <- as.numeric(t(prod_xw)%*%tik) / norm
#     s <- sqrt(sapply(1:g, function(k) sum(((my_x - mu[k])**2) * my_w * tik[,k])) / norm)
#     # loglike
#     logprobs <- sapply(1:g, function(k) log(1/g) + dnorm(my_x, mean = mu[k], sd = s[k], log = TRUE))
#     maxprobs <- apply(logprobs, 1, max)
#     logdensity <- maxprobs + (log(rowSums(exp(sweep(logprobs, 1, maxprobs, "-")))))
#     prec <- loglike
#     loglike <- sum(my_w * logdensity)
#     if ((loglike - prec)<0.000001) break
#   }
#   list(a = mu, b = s, loglike = loglike)
# }

oneEM.gamma <- function(my_x, my_w, g, val){
  res <- oneEMgammaCPP(my_x, my_w, g, val)
  res$a <- as.numeric(res$a)
  res$b <- as.numeric(res$b)
  res
}

smartinitparam <- function(y, M, nbcores){
  x <- na.omit(unlist(y))
  x <- x[which(x!=0)]
  res <- density(x, n = 10**3, from = min(x) , to = max(x), kernel = "r")
  tmp <- mclapply(1:50, function(it) try(oneEM.gamma(res$x, res$y, M, res$x[sample(1:length(res$x), M, prob = res$y/sum(res$y))])), mc.cores = nbcores)
  tmp <- tmp[which(sapply(tmp, function(u) class(u))=="list")]
  tmp <- tmp[[which.max(sapply(tmp, function(res) res$loglike))]]
  list(eps = rep(0, M), a=as.numeric(tmp$a), b=as.numeric(tmp$b))
}
