\name{inputData}
\alias{inputData}
\title{Input data which have been transformed}
\description{
Input all the dataset which have been transformed}
\usage{
inputData(readraw,Genformat,method,trait,PopStrType)
}
\arguments{
	\item{readraw}{genotype matrix.}
	\item{Genformat}{genotype format.}
	\item{method}{which method to analysis.}
	\item{trait}{which trait to analysis.}
      \item{PopStrType}{The type of population structure.}
      }
\author{
Zhang Ya-Wen, Li Pei, Zhang Yuan-Ming\cr
Maintainer: Yuan-Ming Zhang<soyzhang@mail.hzau.edu.cn>
}
\examples{
G1=data(Gen)
P1=data(Phe)
Readraw=ReadData(fileGen=Gen,filePhe=Phe,fileKin=NULL,filePS =NULL,
fileCov=NULL,Genformat=1)
result=inputData(readraw=Readraw,Genformat=1,method="mrMLM",trait=1,
PopStrType=NULL)
}



