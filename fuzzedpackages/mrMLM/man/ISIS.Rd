\name{ISIS}
\alias{ISIS}
\title{To perform GWAS with ISIS EM-BLASSO method}
\description{
Iterative Sure Independence Screening EM-Bayesian LASSO}
\usage{
ISIS(gen,phe,outATCG,genRaw,kk,psmatrix,svpal,svmlod,Genformat,CLO)
}
\arguments{
	\item{gen}{genotype matrix.}
	\item{phe}{phenotype matrix.}
	\item{outATCG}{genotype for code 1.}
	\item{genRaw}{raw genotype.}
	\item{kk}{kinship matrix.}
	\item{psmatrix}{population structure matrix.}
	\item{svpal}{Critical P-value for selecting variable.}
	\item{svmlod}{Critical LOD score for significant QTN.}
	\item{Genformat}{Format for genotypic codes.}
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
InputData=inputData(readraw=Readraw,Genformat=1,method="ISIS EM-BLASSO",
trait=1)
result=ISIS(InputData$doMR$gen,InputData$doMR$phe,InputData$doMR$outATCG,
InputData$doMR$genRaw,InputData$doMR$kk,InputData$doMR$psmatrix,
0.01,svmlod=3,Genformat=1,CLO=1) 
}


