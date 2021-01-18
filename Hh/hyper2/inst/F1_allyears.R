library("hyper2")
library("magrittr")

## Analysis of different years of F1 results.  File formula1.R is
## functionally identical to this but uses slightly different idiom.

files <- c(
    "formula1_2014.txt",
    "formula1_2015.txt",
    "formula1_2016.txt",
    "formula1_2017.txt"
    )

allyears <- list()
drivers <- c()
for(i in seq_along(files)){
    jj <- read.table(files[i],header=TRUE)
  allyears[[i]] <- jj
  drivers <- c(drivers, as.character(jj$driver))
}

drivers %<>% unique %>% sort


likelihood_from_finishing_order <- function(H, df){
  ## 'H' a pre-existing hyper2 object [needed because we have to know
  ## the drivers' names ahead of time], 'df' a dataframe such as
  ## read.table("formula1_2017.txt",header=T)

    drivers <- as.character(df$driver)

    for(i in seq(from=2,to=ncol(df)-1)){  # first column is driver, last is points
        d <- as.numeric(as.character(df[,i,drop=TRUE])) # coerces text to NA
        jj <- is.na(d)|(d==0)
        nonfinishers <- drivers[jj]
        finishers <- drivers[!jj]
        finishers <- finishers[order(d[!jj])]
        while(length(finishers)>1){
            H[finishers[1 ]] %<>% "+"(1)
            H[c(nonfinishers,finishers)]%<>% "-"(1)
            finishers <- finishers[-1]
        }
    }
    return(H)
}


H <- hyper2(pnames=drivers)

if(FALSE){  # too slow!
    for(p in allyears){
        H %<>% likelihood_from_finishing_order(p)
    }
} else {  ## Much faster (identical result):
  jj <- list(
    H1 = likelihood_from_finishing_order(hyper2(pnames=drivers),allyears[[1]])
    H2 = likelihood_from_finishing_order(hyper2(pnames=drivers),allyears[[2]])
    H3 = likelihood_from_finishing_order(hyper2(pnames=drivers),allyears[[3]])
    H4 = likelihood_from_finishing_order(hyper2(pnames=drivers),allyears[[4]])
    
    ## HH = H1+H2+H3+H4
    HH <- Reduce("+",jj)
}

