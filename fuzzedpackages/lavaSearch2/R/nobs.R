### nobs.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: maj  2 2018 (09:15) 
## Version: 
## Last-Updated: jun 14 2019 (13:41) 
##           By: Brice Ozenne
##     Update #: 3
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## not defined by lava
nobs.lvmfit <- function(object){
    return(object$data$n)
}

nobs.gls2 <- function(object){
    return(NROW(object$sCorrect$score))
}
nobs.lme2 <- nobs.gls2

######################################################################
### nobs.R ends here
