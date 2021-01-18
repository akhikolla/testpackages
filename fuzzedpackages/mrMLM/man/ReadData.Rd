\name{ReadData}
\alias{ReadData}
\title{read raw data}
\description{
read raw data which have not been transformed}
\usage{
ReadData(fileGen,filePhe,fileKin,filePS,fileCov,Genformat)
}
\arguments{
	\item{fileGen}{genotype matrix.}
	\item{filePhe}{phenotype matrix.}
	\item{fileKin}{kinship matrix.}
	\item{filePS}{population structure matrix.}
      \item{fileCov}{Covariate matrix.}
      \item{Genformat}{genotype format.}
      }
\author{
Zhang Ya-Wen, Li Pei, Zhang Yuan-Ming\cr
Maintainer: Yuan-Ming Zhang<soyzhang@mail.hzau.edu.cn>
}
\examples{
G1=data(Gen)
P1=data(Phe)
result=ReadData(fileGen=Gen,filePhe=Phe,fileKin=NULL,filePS =NULL,
fileCov=NULL,Genformat=1)
}



