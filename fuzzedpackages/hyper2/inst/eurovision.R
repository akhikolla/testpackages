## This file creates hyper2 object 'euro2009'.

## The dataset wiki_matrix, defined below, is copied from "Eurovision
## Song Contest 2009," Wikipedia, accessed May 13, 2018.

## More documentation is given in euro.Rd [type help(euro2009) at the
## R prompt]

library("hyper2")

abbreviated <- TRUE
## change to FALSE to use full country names rather than two-letter abbreviations.


## First specify the matrix as appearing in the Wikipedia page:
wiki_matrix <- matrix(c(
##  ME  CZ  BE  BY  SW  AM  AD  CH  TR  IL  BG  IS  MK  RO  FI  PT  MT  BA  DE  UK
    NA, 00, 00, 03, 00, 05, 01, 02, 05, 01, 00, 00, 08, 00, 00, 01, 06, 10, 02, 00, # ME
    00, NA, 00, 00, 00, 00, 00, 00, 00, 00, 00, 00, 00, 00, 00, 00, 00, 00, 00, 00, # CZ
    00, 00, NA, 00, 00, 01, 00, 00, 00, 00, 00, 00, 00, 00, 00, 00, 00, 00, 00, 00, # BE
    02, 01, 00, NA, 01, 04, 00, 00, 00, 04, 01, 01, 06, 00, 04, 00, 01, 00, 00, 00, # BY
    00, 06, 04, 07, NA, 08, 07, 04, 04, 07, 00, 10, 03, 04, 10, 08, 08, 04, 04, 07, # SW
    04, 12, 10, 10, 05, NA, 00, 01, 10, 10, 08, 02, 02, 08, 01, 00, 00, 01, 10, 05, # AM
    00, 00, 00, 00, 00, 00, NA, 00, 01, 00, 00, 00, 00, 00, 00, 04, 03, 00, 00, 00, # AD
    00, 00, 00, 02, 02, 00, 02, NA, 00, 00, 00, 00, 00, 00, 05, 02, 00, 02, 00, 00, # CH
    08, 05, 12, 06, 07, 10, 05, 12, NA, 06, 12, 07, 12, 12, 07, 05, 10, 12, 12, 12, # TR
    05, 04, 03, 04, 06, 07, 08, 05, 03, NA, 04, 06, 01, 03, 06, 00, 04, 00, 05, 01, # IL
    00, 00, 00, 00, 00, 00, 00, 00, 02, 00, NA, 00, 05, 00, 00, 00, 00, 00, 00, 00, # BG
    07, 10, 07, 12, 12, 12, 10, 07, 08, 12, 06, NA, 04, 10, 12, 12, 12, 07, 06, 08, # IS
    10, 03, 00, 00, 00, 00, 00, 06, 06, 00, 10, 00, NA, 02, 00, 00, 00, 08, 00, 00, # MK
    06, 00, 02, 01, 00, 02, 04, 00, 07, 08, 05, 04, 07, NA, 00, 10, 02, 06, 01, 02, # RO
    03, 00, 01, 00, 10, 00, 03, 00, 00, 00, 00, 12, 00, 01, NA, 03, 05, 00, 00, 04, # FI
    00, 02, 06, 00, 03, 00, 12, 10, 00, 02, 02, 08, 00, 07, 02, NA, 00, 03, 07, 06, # PT
    01, 07, 08, 08, 04, 03, 06, 03, 00, 05, 03, 05, 00, 06, 03, 06, NA, 05, 03, 10, # MT
    12, 08, 05, 05, 08, 06, 00, 08, 12, 03, 07, 03, 10, 05, 08, 07, 07, NA, 08, 03) # BA
, nrow=18,byrow=TRUE)

## Each row corresponds to a competitor and each column to a
## judge. Note that the matrix is not square: Germany and the UK were
## judges but not competitors.

## I have removed the first column which is the total of the points,
## replaced self-voting entries with NA, and replaced blanks with
## '00'; every entry is padded to two digits.


points <- c(12,10,8,7,6,5,4,3,2,1)
## Variable 'points' gives the number of points awarded, under
## Eurovision rules, to voters' first, second, third, etc choice.
## Points for any competitor are added and the winner is the
## competitor with the most points.  However, in this model, the
## numerical values themselves do not affect the likelihood function;
## only the order of the voters' preferences matters.


preference <- wiki_matrix*0  
for(i in seq_along(points)){ preference[wiki_matrix == points[i]] <- i }
## Matrix 'preference' records voters' first, second, third, etc
## choice.  A zero entry means no points (nul punkte!) and NA means
## that voter was forbidden from voting for that player (countries
## cannot vote for themselves).


countries <- data.frame(
    fullname =
        c("Montenegro", "Czech rep", "Belgium", "Belarus", "Sweden",
          "Armenia", "Andorra", "Switzerland", "Turkey", "Israel",
          "Bulgaria", "Iceland", "Macedonia", "Romania", "Finland",
          "Portugal", "Malta", "Bosnia Herz", "Germany", "UK"),
    abbreviation = c("ME","CZ","BE", "BY", "SW", "AM", "AD", "CH",
                     "TR", "IL", "BG", "IS", "MK", "RO", "FI", "PT",
                     "MT", "BA", "DE", "UK")
)

if(abbreviated){ 
    jj <- countries$abbreviation
} else { 
    jj <- countries$fullname
}

competitors <- jj[1:18]
colnames(preference) <- jj      # voters; 20 countries (18 + DE + UK)
rownames(preference) <- competitors  

rownames(wiki_matrix) <- jj[1:18]
colnames(wiki_matrix) <- jj

## The competitors were the first 18 countries (the last two
## countries, Germany and the UK, voted but did not compete)

## We need the voters' choices to be the *rows* for consistency with
## the rest of the package (and, for that matter, the aylmer package
## and indeed the emulator package), so we take the transpose:
preference <- t(preference)

## Now, take the first row of 'preference'.  This represents the votes
## cast BY (sic) Montenegro ("ME").  Their favourite was [last column]
## Bosnia & Herzegovina who they gave rank 1 to.  Their second
## favourite was Macedonia, their third was Turkey, and so on.  They
## were not allowed to vote for themselves, which is why the first
## column is NA.  So the order was:
## bh, mc, tu, ic, ro, is, ar, fi, br, ma

## Define an empty hyper2 object:
euro2009 <- hyper2(d=18)

for(i in seq_len(nrow(preference))){   # cycle through the rows; each row is a voter
    d <- preference[i,,drop=TRUE]
    d[is.na(d)] <- -1  # kludge: make the voting country ineligible to vote.
    while(any(d>0)){
        eligible <- which(d>=0)   # This is why we set NA values to -1
        euro2009[which(d==1)] %<>% "+"(1)
                                        # The first choice among
                                        # eligible players has +1
                                        # power on the numerator
        
        euro2009[eligible] %<>% dec()   # denominator of all eligible players

        d[d==1] <- -1  # once you've won, you are ineligible to be chosen again
                       # NB, pipe idiom is ugly:  d %<>% `[<-`(. == 1, -1)

        d[d>0] %<>% dec()  # everyone moves down the list, so who
                            # *was* second choice becomes first
                            # choice, who *was* third choice becomes
                            # second, and so on.

    } # while() loop closes
} # i loop closes


## syntatic sugar:
pnames(euro2009) <- competitors

## plot points vs strength:
plot(rowSums(wiki_matrix,na.rm=TRUE),maxp(euro2009),
     pch=16,xlab="points scored",ylab="MLE strength")
