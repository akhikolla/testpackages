\name{pLARmEB}
\alias{pLARmEB}
\title{To perform GWAS with pLARmEB method}
\description{
polygene-background-control-based least angle regression plus Empirical Bayes}
\usage{
pLARmEB(gen,phe,outATCG,genRaw,kk,psmatrix,CriLOD,lars1,Genformat,Bootstrap,CLO)
}
\arguments{
	\item{gen}{genotype matrix.}
	\item{phe}{phenotype matrix.}
	\item{outATCG}{genotype for code 1.}
	\item{genRaw}{raw genotype.}
	\item{kk}{kinship matrix.}
	\item{psmatrix}{population structure matrix.}
      	\item{CriLOD}{Critical LOD score for significant QTN.}
      \item{lars1}{No. of potentially associated variables selected by LARS.}
      \item{Genformat}{Format for genotypic codes.}
      \item{Bootstrap}{Bootstrap=FALSE indicates the analysis of only real dataset, Bootstrap=TRUE indicates the analysis of both real dataset and four resampling datasets.}
      \item{CLO}{number of CPU.}
      }
\author{
Zhang Ya-Wen, Li Pei, Zhang Yuan-Ming\cr
Maintainer: Yuan-Ming Zhang<soyzhang@mail.hzau.edu.cn>
}
\examples{
G1=data(Gen)
P1=data(Phe)
Readraw=ReadData(fileGen=Gen,filePhe=Phe,fileKin=NULL,filePS =NULL,
Genformat=1)
InputData=inputData(readraw=Readraw,Genformat=1,method="pLARmEB",trait=1)
result=pLARmEB(InputData$doMR$gen,InputData$doMR$phe,InputData$doMR$outATCG,
InputData$doMR$genRaw,InputData$doMR$kk,InputData$doMR$psmatrix,
CriLOD=3,lars1=20,Genformat=1,Bootstrap=FALSE,CLO=1)   
}




