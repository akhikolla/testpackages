#' An internal function to obtain column index by its name
#' @param x Dataset
#' @param col.name Column name 
#' @return column index(es) in the provided dataset
get.column.index <- function(x, col.name) {
    col.ind <- grep(paste("\\b", col.name, "\\b", sep=""), colnames(x))
    return(col.ind)
}

linear.interpolation <- function(times, points, n) { # Between known pair
    res.times <- c(times[1])
    res <- c(points[1])
    #y = dy/dt * (t - t1) + y1
    for(i in 1:(length(points)-1)) {
        y1 <- points[i]
        y2 <- points[i+1]
        dy <- y2 - y1
        t1 <- times[i]
        t2 <- times[i+1]
        if(t1 == t2) t2 <- 1.01*t1
        dt <- t2 - t1
        t <- seq(t1, t2, by=dt/n)
        y <- dy/dt * (t - t1) + y1
        res <- c(res, y[-1])
        res.times <- c(res.times, t[-1])
    }
    
    return(list(res, res.times))
}

#'Data pre-processing for analysis with stochastic process model methodology.
#'@param x A path to the file with table of follow-up oservations (longitudinal table). 
#'File formats: csv, sas7bdat
#'@param col.id A name of column containing subject ID. 
#'This ID should be the same in both x (longitudinal) and y (vital statistics) tables.
#'None: if col.id not provided, the first column of the x and 
#'first column of the y will be used by default.
#'@param col.status A name of the column containing status variable 
#'(0/1, which is an indicator of death/censoring). 
#'Note: if not provided - then the column #2 from the y (vital statistics) dataset will be used.
#'@param col.age A name of age column (also called 't1'). 
#'This column represents a time (age) of measurement.
#'If not provided then the 3rd column from the longitudinal dataset (x) will be used.
#'@param col.age.event A name of 'event' column.
#'The event column indicates a time when the even occured (e.g. system failure).
#'Note: if not provided then the 3rd column from the y (vital statistics) dataset will be used.
#'@param covariates A list of covariates (physiological variables). 
#'If covariates not provided, then all columns from longitudinal table having index > 3 will be used as covariates. 
#'@param interval A number of breaks between observations for data for discrete model. 
#'This interval must be integer and should be equal or greater than 1.
#'Default = 1 unit of time. 
#'@param verbose A verbosing output indicator. Default=FALSE.
#'@return A list of two elements: first element contains a preprocessed data for continuous model, with arbitrary intervals between observations  and 
#'second element contains a prepocessed data table for a discrete model (with constant intervals between observations).
#'@export
#'@examples \dontrun{ 
#'library(stpm) 
#'data <- prepare_data(x=system.file("extdata","longdat.csv",package="stpm"))
#'head(data[[1]])
#'head(data[[2]])
#'}
prepare_data <- function(x,
                         col.id=NA, 
                         col.status=NA,
                         col.age=NA, 
                         col.age.event=NA, 
                         covariates=NA, 
                         interval=1, 
                         verbose=FALSE) {
  
    ### Constant variables ###
    col.id.ind <- 1
    col.status.ind <- 2
    col.age.ind <- 3
    col.covar.ind <- c()
  
    if(interval < 1) 
    {
        warning("Interval must be more or equal to 1.")
    } 
    else if(interval != round(interval)) 
    {
        stop("Interval must be integer.")
    }
  
    if(class(x) == "character")
    {
        if(file_ext(x) == "csv") 
        {
            merged.data <- read.csv(x)
        } else if(file_ext(x) == "sas7bdat") 
        {
            merged.data <- read.sas7bdat(x)
        } else {
            stop(paste(x, ":", "unknown file format, it must be csv or sas7bdat."))
        }
    } 
    else if(class(x) == "data.frame")
    {
        merged.data <- x
    } else
    {
        stop("Unknown format of input x.")
    }
  
    # Parsing input parameters in order to check for errors:
    if( get.column.index(merged.data, col.id) == 0 ) {
        warning(paste("ID column",col.id, "not found in data table!"))
    } else {
        col.id.ind <- get.column.index(merged.data, col.id)
    }
    
    if( get.column.index(merged.data, col.status) == 0 ) {
        warning(paste("Status column",col.status, "not found in data table!"))
    } else {
        col.status.ind <- get.column.index(merged.data, col.status)
    }
  
    if( get.column.index(merged.data, col.age) == 0 ) {
        warning(paste("Age column",col.age, "not found in longdat table!"))
    } else {
        col.age.ind <- get.column.index(merged.data, col.age)
    }
    
    if( is.null(get.column.index(merged.data, col.age.event)) ) {
            warning(paste("Event column", col.age.event, "not found in data table!"))
            col.age.event <- col.age
            col.age.event.ind <- col.age.ind
    } else {
        col.age.event.ind <- get.column.index(merged.data, col.age.event)
    }
  
    for(c in covariates) {
        if( get.column.index(merged.data, c) == 0 ) {
            stop(paste("Covariate",c, "not found."))
        } else {
            col.covar.ind <- c(col.covar.ind, get.column.index(merged.data, c))
        }
    } 
  
    if(interval < 1)
    {
        interval <- 1
    }
  
    #-----------Done parsing imput parameters---------------------#
    
    # Remove records in which id = NA
    #merged.data <- merged.data[which(!is.na(merged.data[ , col.id.ind])),]
    #merged.data <- merged.data[which(!is.na(merged.data[ , col.status.ind])),]
    #merged.data <- merged.data[which(!is.na(merged.data[ , col.age.ind])),]
  
    # Prepare data for continuous optimisation:
    data_cont <- prepare_data_cont(merged.data, col.status.ind, col.id.ind, col.age.ind, col.age.event.ind, col.covar.ind, verbose, interval)
  
    # Prepare data for fast discrete optimization:
    #data_discr <- prepare_data_discr(merged.data, interval, col.status.ind, col.id.ind, col.age.ind, col.age.event.ind, col.covar.ind, verbose)
  
    list(model.continuous=data_cont, model.discrete=NULL)
}

#'Prepares continuouts-time dataset.
#'@param merged.data a longitudinal study dataset.
#'@param col.status.ind index of "status" column.
#'@param col.id.ind subject id column index.
#'@param col.age.ind index of the age column.
#'@param col.age.event.ind an index of the column which represents the time in which event occured.
#'@param col.covar.ind a set of column indexes which represent covariates.
#'@param verbose turns on/off verbosing output.
#'@param dt interval between observations.
prepare_data_cont <- function(merged.data, 
                              col.status.ind, 
                              col.id.ind, 
                              col.age.ind, 
                              col.age.event.ind, 
                              col.covar.ind, 
                              verbose,
                              dt) {
    #merged.data = long
    #
    #col.id.ind = 1
    #col.status.ind = 2 
    #col.age.ind = 3
    #col.age.event.ind = 3 
    #col.covar.ind = 4
    #dt <- 2
    # Split records by ID:
    #splitted <- split(merged.data, factor(merged.data[ , col.id.ind]))
    splitted <- split(merged.data, merged.data[[col.id.ind]])
    
  
    for(iii in 1:length(splitted)) {
        nrows <- length(splitted[[iii]][ , col.id.ind])
        if(nrows > 1) {
            id <- splitted[[iii]][ -1, col.id.ind]
            case <- splitted[[iii]][-1, col.status.ind]
            t1 <- splitted[[iii]][ 1:(nrows-1), col.age.ind]
            #t2 <- c(splitted[[iii]][ , col.age.ind][-1], tail(splitted[[iii]][ , col.age.event.ind],n=1))
            t2 <- splitted[[iii]][ , col.age.ind][-1]
            
        } else {
            id <- splitted[[iii]][ 1 , col.id.ind]
            case <- splitted[[iii]][1, col.status.ind]
            t1 <- splitted[[iii]][ 1:(nrows-1), col.age.ind]
            #t2 <- c(splitted[[iii]][ , col.age.ind][-1], tail(splitted[[iii]][ , col.age.event.ind],n=1))
            t2 <- t1 + 0.01*t1
            
        }
        
        tmp.frame <- cbind(id, case, t1, t2)
        # Adding covariates:
        if(nrows > 1) {
            for(ind in col.covar.ind) {
                tmp.frame <- cbind(tmp.frame, 
                               splitted[[iii]][1:(nrows-1), ind], 
                               splitted[[iii]][-1, ind])
      
            }
        } else {
            for(ind in col.covar.ind) {
                tmp.frame <- cbind(tmp.frame, 
                               splitted[[iii]][1:(nrows-1), ind], 
                               splitted[[iii]][1, ind])
            
            }
        }
        
        splitted[[iii]] <- tmp.frame
    }
    
    prep.dat <- do.call("rbind", splitted)
    rownames(prep.dat) <- c()
    #prep.dat <- prep.dat[rowSums( matrix(is.na(prep.dat[,5:dim(prep.dat)[2]]), ncol=2*length(covariates),byrow=T)) !=2*length(covariates),]
    prep.dat <- prep.dat[which(!is.na(prep.dat[,4])),]
  
    if(verbose) {
        head(prep.dat)  
    }
  
    ans_final <- prep.dat
  
    if(verbose)
        cat("Making final table...\n")
  
    # Finalizing:
    colnames(ans_final) <- c("id", "case", "t1", "t2", unlist(lapply(1:length(col.covar.ind), function(n) {c(names(merged.data)[col.covar.ind[n]], paste(names(merged.data)[col.covar.ind[n]],".next",sep=""))} )) )
    #print(ans_final)
    ans_final <- data.frame(ans_final[which(ans_final[,3] <= ans_final[,4]),]) # t1 must be less than t3
    
    # t1 must be equal t3 on previous step, if status = 0 and id is the same
    ndim <- dim(ans_final)[2] - 4
    for(i in 2:(dim(ans_final)[1]-1)) {
        if( (ans_final$case[i] == 0) & (ans_final$id[i] == ans_final$id[(i-1)]) & (ans_final$id[i] == ans_final$id[(i+1)]) ) {
            for(ii in seq(0,(ndim-1),2)) {
                ans_final[(i+1),(5+ii)] <- ans_final[i,(6+ii)]
            }
        } else if(ans_final$case[i] == 1) {
            for(ii in seq(0,(ndim-1),2)) {
                ans_final[i,(6+ii)] <- NA
            }
        } else if( ans_final$id[i] != ans_final$id[(i+1)] ) {
            if(ans_final$t1[i] >= ans_final$t2[i]) {
                ans_final$t2[i] <- ans_final$t1[i] + dt/2
            }
        }
    }
  
    return(data.frame(ans_final))
}

#'Prepares discrete-time dataset.
#'@param merged.data a longitudinal study dataset.
#'@param interval interval between observations.
#'@param col.status.ind index of "status" column.
#'@param col.id.ind subject id column index.
#'@param col.age.ind index of the age column.
#'@param col.age.event.ind an index of the column which represents the time in which event occured.
#'@param col.covar.ind a set of column indexes which represent covariates.
#'@param verbose turns on/off verbosing output.
prepare_data_discr <- function(merged.data, interval, col.status.ind, col.id.ind, col.age.ind, col.age.event.ind, col.covar.ind, verbose) {
  #---DEBUG---#
  #merged.data = data
  #col.status.ind = 2 
  #col.id.ind = 1
  #col.age.ind = 3
  #col.age.event.ind = 3 
  #col.covar.ind = 4
  #interval <- 1
  #---END DEBUG---#
  
    #'Filling the last cell
    fill_last <- function(x) {
        na_idx <- which(is.na(x))
        unique_elements <- unique(x[-na_idx])
        set_diff <- unique_elements[length(unique_elements)]
        x[na_idx] <- set_diff
        x
    }
  
  
    # Interpolation
    dt <- interval
    ans <- matrix(nrow=0, ncol=(4+2*length(col.covar.ind)))
    # Split records by ID:
    splitted <- split(merged.data, merged.data[, col.id.ind])
  
    # For each particular person's record:
    for(iii in 1:length(splitted)) {
        if( !is.na(tail(splitted[[iii]][ , col.age.event.ind],n=1)) & !is.na(tail(splitted[[iii]][ , col.status.ind],n=1)) ) {
            if(verbose) {
                print(paste(iii, "individual processed."))
            }
            # Individual ID:
            id <- splitted[[iii]][ , col.id.ind][1]
      
            #nrows <- (dim(splitted[[iii]])[1]-1)*dt
            nrows <- (dim(splitted[[iii]])[1])*dt
            
            if(nrows > 1) {
                t1.approx <- matrix(ncol=4, nrow=nrows)
                t1.approx[,1] <- id
                t1.approx[,2] <- 0
                t1.approx[nrows,2] <- tail(splitted[[iii]][ , col.status.ind],n=1) #Last value
                    
                aprx <- linear.interpolation(splitted[[iii]][, col.age.ind], splitted[[iii]][, col.covar.ind[1]], n=dt)
                
                print(splitted[[iii]])
                print(nrows)
                print(aprx)
                print(t1.approx)
                
                #t1.approx[,3] <- aprx[[2]][1:(length(aprx[[2]])-1)]
                t1.approx[,3] <- aprx[[2]][1:(length(aprx[[2]]))]
                
                if(col.age.ind == col.age.event.ind) {
                    t1.approx[,4] <- c(aprx[[2]][-1], tail(aprx[[2]],1)*0.01)
                } else {
                    t1.approx[,4] <- c(aprx[[2]][-1], tail(splitted[[iii]][, col.age.event.ind], 1))
                }
            } else {
                print(splitted[[iii]])
                t1.approx <- matrix(ncol=4, nrow=1)
                t1.approx[,1] <- id
                t1.approx[,2] <- tail(splitted[[iii]][ , col.status.ind],n=1) #Last value
                t1.approx[,3] <- splitted[[iii]][ , col.age.ind]
                t1.approx[,4] <- 1.01*t1.approx[,3]
            }
          
            
            for(ind in col.covar.ind) {
                if ( (length(splitted[[iii]][, ind]) > 1) & (length(which(!is.na(splitted[[iii]][, ind]))) > 0) ) {
                    if(length(which(!is.na(splitted[[iii]][, ind]))) == 1) {
                        splitted[[iii]][, ind] <- fill_last(splitted[[iii]][, ind])
                    }
                    
                    aprx <- linear.interpolation(splitted[[iii]][, col.age.ind], splitted[[iii]][, ind], n=dt)
                    
                    t1.approx <- cbind(t1.approx, aprx[[1]][1:(length(aprx[[1]])-1)], aprx[[1]][-1])
                } else {
                    t1.approx <- cbind(t1.approx, splitted[[iii]][, ind], splitted[[iii]][, ind])
                }
            }
            
            ans <- rbind(ans, t1.approx)
        }
    }
  
    colnames(ans) <- c("id", "case", "t1", "t2", 
                     unlist(lapply(1:length(col.covar.ind), 
                                   function(n) {c(names(merged.data)[col.covar.ind[n]], 
                                                paste(names(merged.data)[col.covar.ind[n]],".next",sep=""))})))
    rownames(ans) <- 1:dim(ans)[1]
    return(data.frame(ans))
}
