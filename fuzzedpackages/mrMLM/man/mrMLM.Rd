\name{mrMLM}
\alias{mrMLM}
\alias{mrMLM}
\docType{package}
\title{Multi-Locus Random-SNP-Effect Mixed Linear Model Tools for Genome-Wide Association Study 
}
\description{
Conduct multi-locus genome-wide association study under the framework of multi-locus random-SNP-effect mixed linear model (mrMLM). First, each marker on the genome is scanned. Bonferroni correction is replaced by a less stringent selection criterion for significant test. Then, all the markers that are potentially associated with the trait are included in a multi-locus genetic model, their effects are estimated by empirical Bayes and all the nonzero effects were further identified by likelihood ratio test for ture QTL.
}
\details{
\tabular{ll}{
Package: \tab mrMLM\cr
Type: \tab Package\cr
Version: \tab 4.0.2\cr
Date: \tab 2020-10-8\cr
Depends: \tab lars\cr
Imports: \tab methods,foreach,ncvreg,coin,sampling,data.table,doParallel
\cr
License: \tab GPL version 2 or newer\cr
LazyLoad: \tab yes\cr
}
}
\usage{
mrMLM(fileGen,filePhe,fileKin,filePS,PopStrType,fileCov,Genformat,
method,Likelihood,trait,SearchRadius,CriLOD,SelectVariable,Bootstrap,
DrawPlot,Plotformat,dir)
}
\arguments{
	\item{fileGen}{File path and name in your computer of Genotype,
       i.e.,"D:/Users/Genotype_num.csv".}
	\item{filePhe}{File path and name in your computer of Phenotype,
       i.e.,"D:/Users/Phenotype.csv".}
	\item{fileKin}{File path and name in your computer of Kinship,
       i.e.,"D:/Users/Kinship.csv".}
	\item{filePS}{File path and name in your computer of Population 
      Structure,i.e.,"D:/Users/PopStr.csv".}
      \item{PopStrType}{The type of population structure,i.e.,Q (Q matrix), 
      PCA (principal components), EvolPopStr (evolutionary population 
      structure).}                    
      \item{fileCov}{File path and name in your computer of covariate,  
       i.e.,"D:/Users/Covariate.csv".} 
	\item{Genformat}{Format for genotypic codes, Num (number), Cha   
       (character) and Hmp (Hapmap).}
	\item{method}{Six multi-locus GWAS methods. Users may select one to 
       six methods, including mrMLM, FASTmrMLM, FASTmrEMMA, pLARmEB, pKWmEB 
       and ISIS EM-BLASSO.}
	\item{Likelihood}{This parameter is only for FASTmrEMMA, including 
        REML(restricted maximum likelihood) and  ML(maximum likelihood).}
      \item{trait}{Traits analyzed from number 1 to number 2,i.e.,1:2.}
	\item{SearchRadius}{This parameter is only for mrMLM and FASTmrMLM, 
       indicating Search Radius in search of potentially associated QTN,the 
       default is 20.}
	\item{CriLOD}{Critical LOD score for significant QTN.} 
      \item{SelectVariable}{This parameter is only for pLARmEB. 
       SelectVariable=50 indicates that 50 potentially associated variables 
       are selected from each chromosome. Users may change this number in 
        real data analysis in order to obtain the best results as final 
       results,the default is 50.}
      \item{Bootstrap}{This parameter is only for pLARmEB, including FASLE 
      and TRUE, Bootstrap=FALSE indicates the analysis of only real 
      dataset, Bootstrap=TRUE indicates the analysis of both real dataset 
      and four resampling datasets,the default is FALSE.}
       \item{DrawPlot}{This parameter is for all the six methods, including 
        FALSE and TRUE, DrawPlot=FALSE indicates no figure output, 
        DrawPlot=TRUE indicates the output of the Manhattan, QQ figures,the 
      default is TRUE.}
     \item{Plotformat}{This parameter is for all the figure files,
      including *.jpeg, *.png, *.tiff and *.pdf,the default is "tiff".}
     \item{dir}{This parameter is for the save path,i.e.,"D:/Users"} 
}
\author{
Zhang Ya-Wen, Li Pei, Zhang Yuan-Ming\cr
Maintainer: Yuan-Ming Zhang<soyzhang@mail.hzau.edu.cn>
}
\note{
     Once the running of the software mrMLM v4.0.2 is ended, 
     the "results" files should appear on the Directory, which was 
     set up by users before running the software. The results for 
     each trait include "*_intermediate result.csv", "*_Final result.csv", 
     Manhattan plot, and QQ plot. If only pLARmEB and ISIS EM-BLASSO
     methods are selected, there will be no intermediate results and 
     figures output. Users can decompress the mrMLM package and find 
     the User Manual file (name: Instruction.pdf) in the folder of 
     ".../mrMLM/inst/doc".
}
\references{
1. Zhang YM, Mao Y, Xie C, Smith H, Luo L, Xu S. Genetics 2005,169:2267-2275.
2. Wang SB, Feng JY, Ren WL, Huang B, Zhou L, Wen YJ, Zhang J, Dunwell JM, Xu S, Zhang YM. Sci Rep 2016,6:19444.
3. Tamba CL, Ni YL, Zhang YM. PLoS Comput Biol 2017,13(1):e1005357.
4. Zhang J, Feng JY, Ni YL, Wen YJ, Niu Y, Tamba CL, Yue C, Song Q, Zhang YM. Heredity 2018,118(6):517-524.
5. Ren WL, Wen YJ, Dunwell JM, Zhang YM. Heredity 2018,120(3): 208-218.
6. Wen YJ, Zhang H, Ni YL, Huang B, Zhang J, Feng JY, Wang SB, Dunwell JM, Zhang YM, Wu R. Brief Bioinform 2018,19(4): 700-712.
7. Tamba CL, Zhang YM. bioRxiv,preprint first posted online Jun. 7, 2018,
doi:https://doi.org/10.1101/341784.
8. Zhang YW, Tamba CL, Wen YJ, Li P, Ren WL, Ni YL, Gao J, Zhang YM. 
Genomics, Proteomics & Bioinformatics, Accept.
}
\examples{
Ge1=data(Genotype)
Ph1=data(Phenotype)
mrMLM(fileGen=Genotype,filePhe=Phenotype,Genformat="Num",
method=c("FASTmrMLM"),trait=1,CriLOD=3,DrawPlot=FALSE,
dir=tempdir())
}