dimension <- 3
signature(dimension + 1)
eplus <- basis(dimension+1)
eminus <- basis(dimension + 2)
e0 <-  (eminus - eplus)/2
einf <- eminus + eplus
E <- e0 %^% einf

point <- function(x){  # IPNS
    stopifnot(length(x)==dimension)
    as.1vector(x) + sum(x^2)*einf/2 + e0
}

sphere <- function(x,r){  # IPNS
    point(x) - r^2*einf/2
}

spherestar <- function(...){ # OPNS
    jj <- list(...)
    stopifnot(length(jj) == dimension+1)
    Reduce(`%^%`,jj)
}

plane <- function(n,d){ # IPNS
    stopifnot(length(n)==dimension)
    stopifnot(length(d) == 1)
    as.1vector(n/sqrt(sum(n^2))) + d*einf
}

planestar <- function(...){ # OPNS; A^B^C^Inf
    circlestar(list(...)) %^% einf  
}

circle     <- function(S1,S2){  # IPNS
    S1 %^% S2
}

circlestar <- function(...){  # OPNS; A^B^C
    stopifnot(length(jj) == dimension)
    Reduce(`%^%`,jj)
}

line <- function(P1,P2){  # IPNS
    P1 %^% P2
}

pointpair <- function(A,B){ # OPNS
    A %^% B
}




