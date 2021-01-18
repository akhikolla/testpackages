\name{DnaAlteration}
\alias{DnaAlteration}
\title{Cytosine methylation beta values for a sample of 198 non-pathological human tissue specimens.}
\description{
  A data frame with the subset of the 'Dna Alteration'
 data set corresponding to the sample of 'blood' and 'non-blood, non placenta' tissues. 
}
\usage{
data(DnaAlteration)
}
\format{
  \describe{
	\item{Labels}{Observations with label "0" correspond to the subset of 'blood' tissues.}
	\item{Column 2--1414}{Cytosine methylation beta values collected at 1413 autosomal CpG loci.}
  }
}
\source{
Christensen, B.C Houseman, E.A. Marsit, C.J. Zheng, S. Wrench, M.R. Wiemels, J.L. Nelson, H.H. Karagas, M.R.
Padbury, J.F. Bueno, R. Sugarbaker, D.J Yeh, R., Wiencke, J.K. Kelsey, K.T. (2009). Aging and Evnironemental
Exposure Alter Tissue-Specific DNA Methylation Dependent upon CpG Island Context. PLoS Genet 5(8), e1000602.
}
\examples{
data(DnaAlteration)
alpha<-0.5
Q<-15
p<-ncol(DnaAlteration[,-1])
ns<-FHCSnumStarts(q=Q,eps=(1-alpha)*4/5)
RunExample<-FALSE
if(RunExample){
	Fit<-FastHCS(x=DnaAlteration[,-1],q=Q,nSamp=ns,seed=0)
	colvec<-rep("orange",nrow(DnaAlteration))
	colvec[DnaAlteration[,1]==1]<-"blue"
	plot(Fit,col=colvec,pch=16)
}

}
\keyword{datasets}
