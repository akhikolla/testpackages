#' Create Parameters 
#'
#' Structures the parameters into the correct format for use in  \link{community}
#'
#' @param Ds		An array of death rates for the structured population, of length \code{n}
#' @param Gs		An array of growth rates for the structured population, of length \code{n-1}
#' @param Rs		Either the seed production rate of adults in the population, or an array of seed production rates, of length \code{n}.
#' @param dispersal Dispersal distances
#' @param radius		Optional (use if there are any interactions). Either one radius of
#' interactions or an array of interaction radiuses, of length \code{n}.
#' @param stress	Optional (use to create a stress gradient). An array of
#' values of stress gradient slope. The full value will be added to death rate at the right
#' of the plot, half value at the middle of the plot, and so on, proportionally.
#' @param n         Number of stages in the population
#' @export
#' @examples
#' # create a sample parameters
#' create.parameters(n=3)
#'
#' # structure parameters from arrays
#' create.parameters(Ds=c(10,5,2),Gs=c(2,2),Rs=20,radius=2)
#'
create.parameters <- function(Ds, Gs, Rs, dispersal, radius, stress, n){
    if(missing(n)) n=length(Ds)
    else{
        if(missing(Ds)) Ds = runif(n,0,5)
        if(missing(Gs)) Gs = runif(n-1,0,5)
        if(missing(Rs)) Rs = runif(n,0,5)
    }
    if(missing(radius)) radius=rep(0,n)
    if(missing(stress)) stress=rep(0,n)
    if(missing(dispersal)) dispersal=rep(1,n)

    if(length(radius)==1){radius <- rep(radius,n)}
    if(length(dispersal)==1){dispersal <- rep(dispersal,n)}
    if(length(Rs)==1){Rs <- c(rep(0,n-1),Rs)}
    Gs[n]<-0
	data.frame(D=Ds,G=Gs,R=Rs,dispersal=dispersal,radius=radius,stress=stress)
}

#' abundance matrix
#' 
#' Returns a matrix with abundances of each life stage/species over time
#'
#' The rows in the matrix are the lifestages/species id. The times are in the row names. To
#' visualize the abundance matrix data we recomment the function \code{\link{stackplot}}.
#' @param data	result of a simulation, created by \code{\link{community}}
#' @param times	array of times at which the abundances will be calculated
#' @param by.age T/F. Use this option to get the number of individuals to reach each age, instead of
#' abundances for each time.
#' @param cap.living Logical. Use this option with by.age=T, to set the time of death of living individuals to max
#' simulation time. Otherwise, living individuals are excluded from the data. Either way, this data
#' will be more representative if only a small fraction of total individuals is living at the end of
#' simulation.
#' @examples
#' data(malthusian)
#' times <- seq(0,malthusian$maxtime,by=0.1)
#' ab <- abundance.matrix(malthusian,times)
#' 
#' ab.by.age <- abundance.matrix(malthusian,times,by.age=TRUE)
#' 
#' @export
abundance.matrix <- function(data,times=seq(0,data$maxtime,length.out=50),by.age=FALSE,cap.living=FALSE){
    ## check if array of times is appropriate to simulation
	if(max(times) > data$maxtime){ "Warning: array of times goes further than simulation maximum time" }

    ## check if by.age
    if(by.age){d <- age.data(data,cap.living)}
    else{d<-data$data}

    ## gather data from time points
	subs <- lapply(times,function(t){subset(d,d$begintime <= t & (d$endtime >= t | is.na(d$endtime)),select=c(1,2))})

    ## number of stage/species id
	n <- data$num.total
    if(n==1) { ## treat n==1 separately because R is weird
        ab <- matrix(sapply(subs,nrow),ncol=1)
    }
    else {
        abmatline <- function(x){
            l <- tapply(x$id,x$sp,length)
            # complete the rows that are missing
            if(length(l) == n){
                abl = l
            }
            else {
                abl <- rep(0,n)
                names(abl) <- 1:n
                for(i in 1:n){
                    if(i %in% names(l)){
                        c <- which(names(l)==i)
                        abl[i] <- l[c]
                    }
                }		
            }
            # if that didn't work because of NA's
            abl[is.na(abl)] <- 0

            abl
        }
        ab <- t(sapply(subs,abmatline))
    }
    rownames(ab) <- times

    ab
}

#' longevity
#' 
#' Calculates the lifespan of each individual. Returns a data.frame with the individual's id,
#' the last stage reached by that individual, the time of birth, time of death (if dead), and
#' longevity (if dead).
#'
#' @param data	result of a simulation, created by \code{\link{community}}
#' @export
#' @examples
#' data(malthusian)
#' longevity(malthusian)
longevity <- function(data){
    d <- data$data[with(data$data,order(-endtime)),] #orders by endtime, last to first
    ind.life <- function(i){c(i$sp[1],i$id[1],min(i$begintime),i$endtime[1])}
    b <- by(d,d$id,ind.life)
    m <- do.call("rbind",b)  # TURNS THE OUTPUT INTO A MATRIX
    r <- data.frame(m)
    names(r)=c("last.stage","id","birth","death")
    r$longevity <- r$death-r$birth
    r
}

# age data
# 
# lifehistory relative to birthdate of each individual (used by abundance.matrix with by.age)
# 
# @param data	result of a simulation, created by \code{\link{community}}
age.data <- function(data,cap.living=F){
    d <- data$data
    if(cap.living){
        d[is.na(d)]<-data$maxtime
    }
    else{
        d<-na.exclude(d)
    }
    relative <- function(i){
        birtht<-min(i$begintime)
        i$begintime <- i$begintime - birtht
        i$endtime <- i$endtime - birtht
        i
    }
    b <- by(d,d$id,relative)
    data$age.data <- do.call("rbind",b)
}

