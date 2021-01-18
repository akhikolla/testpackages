library("hyper2")
library("magrittr")

## This script takes about a minute to run.  It creates hyper2 objects
## 'F1_2014' through 'F1_2017', which are the likelihood functions for
## the F1 results table for each of the years.

## Function F1_likelihood() is defined here; this creates objects such
## as F1_2017 by processing formula1_2017.txt.

## At the end of the script, hyper2 object 'F1_total' is created,
## which is a likelihood function for all four years' data.  This
## object is called 'formula1' in the package.  File F1_allyears.R is
## functionally identical but uses slightly different R idiom.


## use-case for F1_likelihood():

## R> F1_likelihood(wiki_table=read.table("formula1_2017.txt",header=TRUE))

## Files like 'formula1_2017.txt' are directly copied from Wikipedia
## (with slight whitespace changes)


`F1_likelihood` <- function(wiki_table){

  ## columns of wiki_table are assumed to be: driver, venue_1,
  ## venue_2, ..., venue_n, points

  noscore <- c("Ret", "WD", "DNS", "DSQ", "DNP", "NC")

  first_and_last <- c(1,ncol(wiki_table))
  
  racers <- wiki_table[,1]
  venues <- colnames(wiki_table)[-first_and_last]

  

  ## Now create a numeric matrix, fmat.  Two steps.  First step, strip
  ## out no-score entries;  we need to count any of noscore
  ## [Ret=retired, WD=withdrawn, DNS=did not start, etc] as a zero:
  
  f <- function(x){
    if(any(x %in% noscore)){x[x%in%noscore] <- 0}
    return(x)
  }

  jj <- apply(wiki_table,2,f)
  

  ## Second step: convert to numeric and strip out names; transpose of
  ## wiki_table (because we want each row to be a venue):
  fmat <- matrix(as.numeric(jj[,-first_and_last]),byrow=TRUE,ncol=nrow(wiki_table)) 
  colnames(fmat) <- racers
  rownames(fmat) <- venues


  ## Considering Formula1, 2017 as an example: taking the first row of
  ## fmat is AUS (Australia), in which Hamilton came second, Vettel
  ## came first, etc.  The first column of fmat is Hamilton's results.
  ## He came second in AUS, first in CHN, second in BHR, etc.


  ## Following is similar to, but slightly different from, the
  ## analysis in eurovision.R: Define an empty hyper2 object:

  F1 <- hyper2(d=ncol(fmat))

  for(i in seq_len(nrow(fmat))){   # cycle through the rows; each row is a venue [voter]
    d <- fmat[i,,drop=TRUE]
    print(d)
    while(any(d>0)){
      eligible <- which(d>=0)  
      
      ## The first choice among eligible players has +1 power on the
      ## numerator:
      F1[which(d==1)] %<>% "+"(1)

      ## denominator of all eligible players; power -1
      F1[eligible] %<>% "-"(1)

      ## once you've come first in the field, you are ineligible to be first again:
      d[d==1] <- -1  
      
      ## everyone moves down the list, so who *was* in second place
      ## becomes first place, who *was* third place becomes second,
      ## and so on:
      d[d>0] %<>% "-"(1)

    } # while() loop closes
  } # i loop closes
  

  ## syntatic sugar:
  pnames(F1) <- racers

  return(F1)
}  # function F1_likelihood() closes



F1_2014 <- "formula1_2014.txt" %>% read.table(header=TRUE) %>% F1_likelihood
F1_2015 <- "formula1_2015.txt" %>% read.table(header=TRUE) %>% F1_likelihood
F1_2016 <- "formula1_2016.txt" %>% read.table(header=TRUE) %>% F1_likelihood
F1_2017 <- "formula1_2017.txt" %>% read.table(header=TRUE) %>% F1_likelihood


## Do the 2017 season:
wiki_table <- read.table("formula1_2017.txt",header=TRUE)
points <- wiki_table$points
names(points) <- wiki_table$driver

m <- maxp(F1_2017)
dotchart(m,pch=16,main='2017 season')

dev.new()

ox <- order(points,decreasing=TRUE)
oy <- order(m,decreasing=TRUE)
png(file="formula1.png")
par(pty='s') # square plot
plot(ox,oy,asp=1,pty='s',xlim=c(0,25),ylim=c(0,25),pch=16,xlab="official order",ylab="my order",main='Formula 1, 2017 season')
par(xpd=TRUE) # allow drivers' names to appear outside plotting region
for(i in seq_along(ox)){  text(ox[i],oy[i],names(m)[ox[i]],pos=4,col='gray') }
par(xpd=FALSE) # stop diagonal line from protruding beyond plotting region
abline(0,1)
dev.off()

a <- list(F1_2014, F1_2015, F1_2016, F1_2017)
alldrivers <- all_pnames(a)

F1_total <- hyper2(pnames=alldrivers)

F1_total %<>% add(change_pnames(F1_2014,alldrivers)) 
F1_total %<>% add(change_pnames(F1_2015,alldrivers)) 
F1_total %<>% add(change_pnames(F1_2016,alldrivers)) 
F1_total %<>% add(change_pnames(F1_2017,alldrivers)) 

mallyears <- maxp(F1_total)

dev.new()
dotchart(mallyears,pch=16,main='Formula 1, 2014-7')
