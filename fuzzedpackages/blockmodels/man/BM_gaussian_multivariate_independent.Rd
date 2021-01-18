\name{BM_gaussian_multivariate_independent}
\alias{BM_gaussian_multivariate_independent}
\alias{\S4method{BM_gaussian_multivariate_independent}{new}}
\title{Perform estimation on blockmodels for multivariate independent
    homoscedastic gaussian probability distribution}
\description{
    With the provided network and blockmodel type, estimate number of groups,
parameters and node membership}
\usage{\S4method{BM_gaussian_multivariate_independent}{new}(
    membership_type, 
    adj, 
    verbosity=6,
    autosave='',
    plotting=character(0),
    exploration_factor=1.5,
    exploration_direction=numeric(0),
    explore_min=4,
    explore_max=Inf,
    ncores=detectCores())}
\arguments{
    \item{membership_type}{The type of node membership, i.e. 'SBM', 'SBM_sym' or 'LBM'}
    \item{adj}{The list of adjacency matrices. All matrices must have the same size}
    \item{verbosity}{The verbosity level, 0 means quiet. Level 1 display the phase of reinitialization. Level 2 display the level 1 and the ascending and descending phase for the number of groups. Level 3 display the level 2 and the number current number of groups which is estimated. Level 4 display the level 3 and the steps inside the estimation. Level 5 display the level 4, the current status of parallel running jobs and the current sub-step. Level 6 display level 5 and informations about ICL criteria found. Default is level 6. This parameter can be changed by accessing to the field $verbosity of the object.}
\item{autosave}{If \var{autosave} != '', after each estimation, the model object is writed into file \var{autosave}. The model object is readable by the function \var{readRDS}. Use-it for long computation to allow restarting the estimation on system crash. You can use it to alanyze the partial results when the estimation is running. This parameter can be changed by accessing to the field $autosave of the object.}
\item{plotting}{Control plot of ICL values while the estimation is running. If plotting==character(0) (the default), plots are done on screen, if plotting=='', no plot are done, if plotting is a filename, plots are done in this filename. This parameter can be changed by accessing the field $plotting of the object.}\item{exploration_factor}{Control the exploration of the number of groups. The exploration is stop when the number of groups reach exploration factor times the current maximum. By default 1.5. This parameter can be changed by accessing the field $exploration_factor of the object.}
\item{explore_min}{Explore to the explore_min number of groups even if the exploration_factor rule is satisfied. By default 4. This parameter can be changed by accessing the field $explore_min of the object.}
\item{explore_max}{Stop exploration after explore_max number of group in any case. By default Inf. This parameter can be changed by accessing the field $explore_max of the object.}
\item{exploration_direction}{Only for LBM membership. Control the exploration direction for groups number. When provided, the exploring strategy is made to explore the provided group number. Must be a vector of two integer value representing the row group number and the column group number.}
\item{ncores}{Number of parallel jobs to launch different EM intializations. By default detectCores(). This parameter can be changed by accessing the field $ncores of the object. This parameters is used only on Linux. Parallism is disabled on other plateform. (Not working on Windows, not tested on Mac OS, not tested on *BSD.)}

}
\examples{\dontrun{

##
## SBM
##

## generation of one SBM network
npc <- 30 # nodes per class
Q <- 3 # classes
n <- npc * Q # nodes
Z<-diag(Q)\%x\%matrix(1,npc,1)
Mu1<-8*matrix(runif(Q*Q),Q,Q)
Mu2<-8*matrix(runif(Q*Q),Q,Q)
M1<-matrix(rnorm(n*n,sd=5),n,n)+Z\%*\%Mu1\%*\%t(Z) ## adjacency
M2<-matrix(rnorm(n*n,sd=10),n,n)+Z\%*\%Mu2\%*\%t(Z) ## adjacency

## estimation
my_model <- BM_gaussian_multivariate_independent("SBM",list(M1,M2) )
my_model$estimate()
which.max(my_model$ICL)

##
## SBM symmetric
##

## generation of one SBM_sym network
npc <- 30 # nodes per class
Q <- 3 # classes
n <- npc * Q # nodes
Z<-diag(Q)\%x\%matrix(1,npc,1)
Mu1<-8*matrix(runif(Q*Q),Q,Q)
Mu2<-8*matrix(runif(Q*Q),Q,Q)
Mu1[lower.tri(Mu1)]<-t(Mu1)[lower.tri(Mu1)]
Mu2[lower.tri(Mu2)]<-t(Mu2)[lower.tri(Mu2)]
M1<-matrix(rnorm(n*n,sd=5),n,n)+Z\%*\%Mu1\%*\%t(Z) ## adjacency
M2<-matrix(rnorm(n*n,sd=10),n,n)+Z\%*\%Mu2\%*\%t(Z) ## adjacency
M1[lower.tri(M1)]<-t(M1)[lower.tri(M1)]
M2[lower.tri(M2)]<-t(M2)[lower.tri(M2)]

## estimation
my_model <- BM_gaussian_multivariate_independent("SBM_sym",list(M1,M2) )
my_model$estimate()
which.max(my_model$ICL)

##
## LBM
##

## generation of one LBM network
npc <- c(50,40) # nodes per class
Q <- c(2,3) # classes
n <- npc * Q # nodes
Z1<-diag(Q[1])\%x\%matrix(1,npc[1],1)
Z2<-diag(Q[2])\%x\%matrix(1,npc[2],1)
Mu1<-8*matrix(runif(Q[1]*Q[2]),Q[1],Q[2])
Mu2<-8*matrix(runif(Q[1]*Q[2]),Q[1],Q[2])
M1<-matrix(rnorm(n[1]*n[2],sd=5),n[1],n[2])+Z1\%*\%Mu1\%*\%t(Z2) ## adjacency
M2<-matrix(rnorm(n[1]*n[2],sd=10),n[1],n[2])+Z1\%*\%Mu2\%*\%t(Z2) ## adjacency

## estimation
my_model <- BM_gaussian_multivariate_independent("LBM",list(M1,M2) )
my_model$estimate()
which.max(my_model$ICL)
}}
