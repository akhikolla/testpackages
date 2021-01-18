 kombinuj1 <- function(l)
{
   stopifnot(is.list(l) && length(l) > 0)
   stopifnot(any(lapply(l,is.numeric) == 1))
   stopifnot(length(l[[1]]) > 0)
   stopifnot(all(lapply(l, length) == length(l[[1]])) )
   
   n <- length(l[[1]])
   k <- length(l)
   new <- vector("list", n)
   
   for (i in 1:n)
      for (j in 1:k)
         new[[i]][j] <- l[[j]][i]
   
   new
}
