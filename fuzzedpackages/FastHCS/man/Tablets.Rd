\name{Tablets}
\alias{Tablets}
\title{Near-infrared (NIR) spectroscopy of a sample of 310 tablets.}
\description{
  The original data set contains  near-infrared (NIR) spectroscopy data for 
310 tablets of four different dosages from pilot, laboratory and full scale 
production settings are included in the study. In this subset, we combine all 80
samples of 80mg tablets with the first 50 samples of 250mg tablets. 
}
\usage{
data(Tablets)
}
\format{
  \describe{
	\item{Labels}{The observations with label '1' correspond to the 80mg Tablets samples and the '0' to the 250mg ones.}
	\item{Column 2--405}{Near Infrared Transmittance; 404 variables; 7400 to 10507 cm-1.}
  }
}
\source{
M. Dyrby, S.B. Engelsen, L. Norgaard, M. Bruhn and L. Lundsberg Nielsen (2002).
Chemometric Quantitation of the Active Substance in a Pharmaceutical Tablet Using Near Infrared (NIR) 
Transmittance and NIR FT Raman Spectra Applied Spectroscopy 56(5): 579--585.
}
\examples{
data(Tablets)
alpha<-0.5
Q<-15
p<-ncol(Tablets[,-1])
ns<-FHCSnumStarts(q=Q,eps=(1-alpha)*4/5)
RunExample<-FALSE
if(RunExample){
	Fit<-FastHCS(x=Tablets[,-1],q=Q,nSamp=ns,seed=1,alpha=0.5)
	colvec<-rep("orange",nrow(Tablets))
	colvec[Tablets[,1]==1]<-"blue"
	plot(Fit,col=colvec,pch=16)
}

}
\keyword{datasets}
