## ----runtime-n, results='hide', message=FALSE, warning=FALSE, echo=FALSE------
library(ggplot2)
library(OptCirClust)
library(knitr)
library(reshape2)

opts_chunk$set(fig.width=6, fig.height=4) 






FOCC <- c()

BOCC <- c()

HEUC <- c()

number <- c()


FOCC_SSQ <- c()

BOCC_SSQ <- c()

HEUC_SSQ <- c()




First <- 1

Prev <- -1

Next <- -1

K <- 3



for(i in seq(10,100, 10 ) )
{


n <- i

set.seed(1)

x <- c(rnorm(n, sd=0.3), rnorm(n, mean=100, sd=0.3), rnorm(n, mean=200, sd=0.3))


Data_Points <- c(x,(x + length(x)))

width <- length(Data_Points)/2

Last <- length(Data_Points) - width 


time <- system.time(result1 <- FramedClust(Data_Points,  K, width, First, Last, "linear.polylog"))

FOCC <- c(FOCC,as.double(time[1]))


time <- system.time(result2 <- FramedClust(Data_Points,  K, width, First, Last, "Ckmeans.1d.dp"))

BOCC <- c(BOCC,as.double(time[1]))


time <- system.time(result3 <- FramedClust(Data_Points,  K, width, First, Last, "kmeans"))

HEUC <- c(HEUC,as.double(time[1]))




number <- c(number, n*3)



}

FOCC <- FOCC * 1000
BOCC <- BOCC * 1000
HEUC <- HEUC * 1000


df1 <- data.frame("No_Points" = number, 
                  "FOCC" = FOCC, "BOCC" = BOCC, "HEUC" = HEUC)

df <- melt(df1, id.var='No_Points')

plot1 <- ggplot(df, aes(x=No_Points, y=value, col=variable)) + geom_line() + geom_point(alpha=3) +
labs(y = "Runtime (millisecond)", x = "Number of points (N)") +
labs(colour = "Methods")

## ----runtime-k,  results='hide', message=FALSE, warning=FALSE, echo=FALSE-----

FOCC <- c()

BOCC <- c()

HEUC <- c()

clusters <- c()


FOCC_SSQ <- c()

BOCC_SSQ <- c()

HEUC_SSQ <- c()




K <- 3


set.seed(1)

x <- c(rnorm(50, sd=0.3), rnorm(50, mean=100, sd=0.3), rnorm(50, mean=200, sd=0.3))


Data_Points <- c(x,(x + length(x)))

width <- length(Data_Points)/2

Last <- length(Data_Points) - width - 1


for(K in seq(10,100,10))
{

  ptm <- proc.time()
  
  result1 <- FramedClust(Data_Points,  K, width, First, Last, "linear.polylog")
  
  time <- proc.time() - ptm
  
  FOCC <- c(FOCC,as.double(time[3]))
  
  FOCC_SSQ <- c(FOCC_SSQ,result1$tot.withinss)
  
  
  
  
  
  ptm <- proc.time()
  
  result2 <-  FramedClust(Data_Points,  K, width, First, Last, "Ckmeans.1d.dp")
  
  time <- proc.time() - ptm
  
  BOCC <- c(BOCC,as.double(time[3]))
  
  BOCC_SSQ <- c(BOCC_SSQ,result2$tot.withinss)
  
  
  
  ptm <- proc.time()
  
  result3 <- FramedClust(Data_Points,  K, width, First, Last, "kmeans")
  
  time <- proc.time() - ptm
  
  HEUC <- c(HEUC,as.double(time[3]))
  
  HEUC_SSQ <- c(HEUC_SSQ,result3$tot.withinss)
  
  
  clusters <- c(clusters, K)
  

}



FOCC <- FOCC * 1000
BOCC <- BOCC * 1000
HEUC <- HEUC * 1000

df1 <- data.frame("No_Clusters" = clusters, 
                  "FOCC" = FOCC, "BOCC" = BOCC, "HEUC" = HEUC)

df <- melt(df1, id.var='No_Clusters')

plot2 <- ggplot(df, aes(x=No_Clusters, y=value, col=variable)) + geom_line() + geom_point(alpha=3) +
labs(y = "Runtime (milliseccond)", x = "Number of clusters (K)") + labs(colour = "Methods") 

## ----SSQ,  results='hide', message=FALSE, warning=FALSE, echo=FALSE-----------


df1 <- data.frame("No_Clusters" = clusters, "FOCC" = FOCC_SSQ, "BOCC" = BOCC_SSQ, "HEUC" = HEUC_SSQ)

df <- melt(df1, id.var='No_Clusters')

df$mysize <- rep(0, nrow(df))

df$mysize[df$variable=="FOCC"] <- 1
plot3 <- ggplot(df, aes(x=No_Clusters, y=value, col=variable, size=mysize)) + geom_line() + geom_point(alpha=3) +
labs(y="Within-cluster sum of squared distances", x = "Number of clusters (K)") + labs(colour = "Methods") +
scale_size(range = c(2, 4), guide="none") + 
  scale_y_continuous(trans = 'log2')








## ----results='hide', message=FALSE, warning=FALSE, echo=FALSE-----------------
plot(plot1)

## ----  results='hide', message=FALSE, warning=FALSE, echo=FALSE---------------
plot(plot2)

## ----  results='hide', message=FALSE, warning=FALSE, echo=FALSE---------------
plot(plot3)

