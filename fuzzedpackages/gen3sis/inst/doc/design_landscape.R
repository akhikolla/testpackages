## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width=7,
  tidy=T,
  fig.align='center',
  tidy.opts = list(width.cutoff=80)
)

## ----eval=FALSE---------------------------------------------------------------
#  Ncol<-31
#  Nrow<-31
#  Ntimesteps<-140

## ----eval=FALSE---------------------------------------------------------------
#  #Create the temperature array and matrix:
#  temp_array<-array(NA,dim=c(Nrow,Ncol,Ntimesteps))
#  temp_dataframe<-matrix(NA,nrow=Nrow*Ncol, ncol=length(c(Ncol,Nrow))+Ntimesteps)
#  
#  #Create the precipitation array and matrix:
#  prec_array<-array(NA,dim=c(Nrow,Ncol,Ntimesteps))
#  prec_dataframe<-matrix(NA,nrow=Nrow*Ncol, ncol=length(c(Ncol,Nrow))+Ntimesteps)
#  
#  #Creating a vector for the names to be used to name input files (one per time step)
#  stringtimestepsnames<-vector(mode="character",length=Ntimesteps)

## ----tidy=T, out.width='100%', echo=F, fig.cap='A1: Animation of the conceptual island that we will create as an input landscape for gen3sis.', fig.margin=T----
knitr::include_graphics("../inst/extdata/CaseStudy1/landscape/case_study_landscape.gif")

## ----eval=FALSE---------------------------------------------------------------
#  for (timestep in 1:Ntimesteps){ # temporal loop
#    counting<-1
#    stringtimestepsnames[timestep]<-paste("X",timestep,"", sep="")
#    for (y in 1:Nrow){ #loop over the first spatial dimension
#      for (x in 1:Ncol){ #loop over the second spatial dimension
#        if((timestep<=10)||(timestep>120 && timestep<=140)){#time steps with only four (2x2) suitable sites
#          if((y>=15 && y<=16)&&(x>=15 && x<=16)) {  #suitable sites
#            temp_array[x,y,timestep] <- rnorm(1,20,0.5) #temperature
#            prec_array[x,y,timestep] <- rnorm(1,500,50) #precipitation
#          }
#        }
#        if((timestep>10 && timestep<=20)||(timestep>100 && timestep<=120)){#time steps with nine (3x3) suitable sites
#          if((y>=15 && y<=17)&&(x>=15 && x<=17)) {
#            temp_array[x,y,timestep] <- rnorm(1,20,0.5)
#            prec_array[x,y,timestep] <- rnorm(1,500,50)
#          }
#        }
#        if((timestep>20 && timestep<=30)||(timestep>80 && timestep<=100)){#time steps with 25 (5x5) suitable sites
#          if((y>=14 && y<=18)&&(x>=14 && x<=18)) {
#            temp_array[x,y,timestep] <- rnorm(1,20,0.5)
#            prec_array[x,y,timestep] <- rnorm(1,500,50)
#          }
#        }
#        if((timestep>30 && timestep<=40)||(timestep>60 && timestep<=80)){#time steps with 49 (7x7) suitable sites
#          if((y>=13 && y<=19)&&(x>=13 && x<=19)) {
#            temp_array[x,y,timestep] <- rnorm(1,20,0.5)
#            prec_array[x,y,timestep] <- rnorm(1,500,50)
#          }
#        }
#        if(timestep>40 && timestep<=60){#time steps with 81 (9x9) suitable sites
#          if((y>=12 && y<=20)&&(x>=12 && x<=20)) {
#            temp_array[x,y,timestep] <- rnorm(1,20,0.5)
#            prec_array[x,y,timestep] <- rnorm(1,500,50)
#          }
#        }
#        #Saving the environmental variables in a dataframe format for distance matrices
#        if(timestep==1){
#          temp_dataframe[counting,1]<-x
#          temp_dataframe[counting,2]<-y
#          prec_dataframe[counting,1]<-x
#          prec_dataframe[counting,2]<-y
#        }
#        temp_dataframe[counting,2+timestep]<-temp_dataframe[x,y,timestep]
#        prec_dataframe[counting,2+timestep]<-prec_dataframe[x,y,timestep]
#        counting<-counting+1
#      }
#    }
#  }

## ----eval=FALSE---------------------------------------------------------------
#  library(raster)
#  landscapes_list <- list()
#  for (timestep in 1:Ntimesteps){
#    temp_raster <- rasterFromXYZ(temp_dataframe[, c(1,2, timestep+2)])
#    prec_raster <- rasterFromXYZ(prec_dataframe[, c(1,2, timestep+2)])
#  
#    landscapes_list$temp <- c(landscapes_list$temp, temp_raster)
#    landscapes_list$prec <- c(landscapes_list$prec, prec_raster)
#  }
#  
#  ##saving the list of rasters into .rds format to be used as input
#  saveRDS(landscapes_list, "inputfolder/my_experiment/landscapes.rds")

