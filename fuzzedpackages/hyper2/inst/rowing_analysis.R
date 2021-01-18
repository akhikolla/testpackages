## This script creates hyper2 objects "H" and "I" which are identical
## but object "I" uses slicker idiom.  These objects are identical to
## sculls2016 [use 'data(rowing)' at the prompt to access this].

## More documentation is provided at rowing.Rd.

library("hyper2")



filename <- "rowing.txt"  # could be rowing_minimal.txt
o <- strsplit(readLines(filename)," ")

rowers <- sort(unique(unlist(o)))
H <- hyper2(list(),0,pnames=rowers)

for(v in o){
    v <- rev(v)
    for(i in seq_along(v)){
        H[v[i]] %<>% inc           # H[v[i]] %<>% `+`(1)
        H[v[seq_len(i)]] %<>% dec  # H[v[seq_len(i)]] %<>% `-`(1)
    }
}

I <- hyper2(pnames=rowers)
for(v in o){
  I <- I+order_likelihood(character_to_number(v,rowers))
}


data("rowing")
stopifnot(H == sculls2016)
stopifnot(I == sculls2016)
