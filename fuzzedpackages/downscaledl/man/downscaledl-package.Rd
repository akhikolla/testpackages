\name{downscaledl-package}
\alias{downscaledl-package}
\alias{downscaledl}
\docType{package}
\title{
\packageTitle{downscaledl}
}
\description{
\packageDescription{downscaledl}
}
\details{
This package provides the downscaling functions for the coarse-resolution remote sensing images. Currently, we used the deep learning algorithm as the base learner in downscaling. 

The major functions defined in this package include: 

\code{\link{AutoEncoderModel}}:The residual network function used in downscaling;

\code{\link{ResautoDownscale}}:The major downscaling function that uses \code{\link{AutoEncoderModel}};

\code{\link{r2_squ}}:The function calculating R2 ;

\code{\link{rmse}}:The function calculate RMSE.

}
\author{
\packageAuthor{downscaledl}

Maintainer: \packageMaintainer{downscaledl}
}
\references{
~~ Literature or other references for background information ~~
}
\keyword{ package }

\examples{
  #Load the high-resolution raster of a covariate, 
  # elevation to be used for downscaling  
  eleFile=file.path(system.file(package = "downscaledl"), "extdata", "sample_ele.tif")
  ele=raster::raster(eleFile)

  #Load the coarse-resolution raster of the target variable to be downscaled  
  coarseFile=file.path(system.file(package = "downscaledl"), "extdata", "sample_coarse_res.tif")
  coarseRst=raster::raster(coarseFile)

  #Extract x and y to be used as two predictors in downscaling 
  xRast=ele
  yRast=ele
  pos=raster::rasterToPoints(ele,spatial=TRUE)
  cell=raster::cellFromXY(xRast,pos)
  xyc=sp::coordinates(pos)
  xRast[cell]=xyc[,"x"]
  yRast[cell]=xyc[,"y"]

  #Merge the covariates 
  covStk=raster::stack(xRast,yRast,ele)
  names(covStk)=c("x","y","ele")

  #Use the fine-resolution covariate (elevation) as the target image 
  fineTarget=ele
\donttest{
  #Set the paramneters and start to downscale ...  
  ares=ResautoDownscale(covStk,fineTarget,coarseRst,ss= 0.2, cores= 5, thresh = 0.01,ntime=3)
  
  #Show the iteration results 
  message(paste(capture.output(ares$diogRMSE), collapse = "\n",sep=""))
  
  #Show the optimal results in the final predictions 
  message(paste("test R2:",round(ares$r2,2),",test RMSE:",round(ares$rmse,4),sep=""))  
  
  #Obtain the downscaled
  downscaled_img=ares$raster
  
  #Save the current par setting 
  curpar = par(no.readonly = TRUE)
  
  #Set the new par setting 
  par(mfrow=c(1,2),mar=c(1,1,1,1))
  
  #Show the final predictions of fine resolution and 
  # original coarse-resolution image for a comparison 
  raster::plot(coarseRst)
  raster::plot(downscaled_img)
  
  #Restore the previous par setting 
  par(curpar) 
}
}
