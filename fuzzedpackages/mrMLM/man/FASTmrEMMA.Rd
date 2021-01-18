\name{FASTmrEMMA}
\alias{FASTmrEMMA}
\title{To perform GWAS with FASTmrEMMA method}
\description{
FAST multi-locus random-SNP-effect EMMA}
\usage{
FASTmrEMMA(gen,phe,outATCG,genRaw,kk,psmatrix,svpal,svmlod,Genformat,Likelihood,CLO)
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
      \item{Likelihood}{restricted maximum likelihood (REML) and maximum likelihood (ML).}
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
InputData=inputData(readraw=Readraw,Genformat=1,method="FASTmrEMMA",trait=1)
result=FASTmrEMMA(InputData$doFME$gen,InputData$doFME$phe,
InputData$doFME$outATCG,InputData$doFME$genRaw,
InputData$doFME$kk,InputData$doFME$psmatrix,0.005,
svmlod=3,Genformat=1,Likelihood="REML",CLO=1)   
}


