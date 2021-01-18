# perm <- function(transpositions){
#   x <- seq_along(transpositions)
#   for(i in x){
#     x[c(i,transpositions[i])] <- x[c(transpositions[i],i)]
#   }
#   x
# }

isStrictPositiveInteger <- function(x){
  all(x > 0 & (floor(x) == x))
}