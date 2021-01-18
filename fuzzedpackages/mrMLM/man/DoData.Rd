\name{DoData}
\alias{DoData}
\title{process raw data}
\description{
process raw data for later use}
\usage{
DoData(genRaw,Genformat,pheRaw1q,kkRaw,psmatrixRaw,covmatrixRaw,trait,
type,PopStrType)
}
\arguments{
	\item{genRaw}{raw genotype matrix.}
	\item{Genformat}{genotype format.}
	\item{pheRaw1q}{raw phenotype matrix.}
	\item{kkRaw}{raw kinship matrix.}
      \item{psmatrixRaw}{raw population structure matrix.}
      \item{covmatrixRaw}{raw covariate matrix.}
      \item{trait}{which trait to analysis.} 
      \item{type}{which type to transform.}
      \item{PopStrType}{The type of population structure.}
      }
\author{
Zhang Ya-Wen, Li Pei, Zhang Yuan-Ming\cr
Maintainer: Yuan-Ming Zhang<soyzhang@mail.hzau.edu.cn>
}
\examples{
G1=data(Gen)
P1=data(Phe)
readraw=ReadData(fileGen=Gen,filePhe=Phe,fileKin=NULL,filePS =NULL,
fileCov=NULL,Genformat=1)
result=DoData(readraw$genRaw,Genformat=1,readraw$pheRaw1q,readraw$kkRaw,
readraw$psmatrixRaw,readraw$covmatrixRaw,trait=1,type=2,PopStrType=NULL)
}



