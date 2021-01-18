\name{MultiManhattan}
\alias{MultiManhattan}
\title{Drawing multi-locus Manhattan plot}
\description{
Using the results of the mrMLM software to draw a multi-locus Manhattan plot }
\usage{
MultiManhattan(ResultIntermediate,ResultFinal,mar=c(2.9,2.8,0.7,2.8), 
LabDistance=1.5,ScaleDistance=0.4,LabelSize=0.8,ScaleSize=0.7,
AxisLwd=5,TckLength=-0.03,LogTimes=2,LODTimes=1.2,lodline=3, 
dirplot=getwd(), PlotFormat="tiff", 
width=28000,height=7000,pointsize = 60,res=600,
MarkGene=FALSE,Pos_x=NULL,Pos_y=NULL,GeneName=NULL,
GeneNameColour=NULL,...)
}
\arguments{
	\item{ResultIntermediate}{Intermediate results obtained by the mrMLM
       software,"D:/Users/ResultIntermediate.csv".}
      \item{ResultFinal}{Final results obtained by the mrMLM 
       software,"D:/Users/ResultFinal.csv".}
	\item{mar}{A numerical vector of the form c(bottom, left, top, 
        right) which gives the number of lines of margin to be specified on 
        the four sides of the plot, and the default is c(2.9, 2.8, 0.7, 
        2.8).}
	\item{LabDistance}{Distance between label and axis; the default is  
       1.5.}
	\item{ScaleDistance}{Distance between scale values and axis; the 
       default is 0.4.}
    \item{LabelSize}{Size of all the three labels; the default is 0.8.}
	\item{ScaleSize}{Size of scale values; the default is 0.7.}
    \item{AxisLwd}{The width of axis, a positive number; the default 
       is 5.}
	\item{TckLength}{The length of tick marks; the default is -0.03.}
	\item{LogTimes}{Magnification of {-log10(P-value)}; the default is
       2.}
    \item{LODTimes}{Magnification of {LOD score}; the default is 1.2.}
    \item{lodline}{The significant LOD score; the default is 3.}
    \item{dirplot}{Path to save plot; the default is current working
       directory}
    \item{PlotFormat}{Format of the plot.i.e., *.tiff, *.png, *.jpeg, 
      *.pdf}
    \item{width}{Figure width; the default is 28000.}
    \item{height}{Figure height; the default is 7000.}
    \item{pointsize}{Word resolution, with the unit of 1/72 inch, being 
       pixels per inch (ppi); the default is 60.}
    \item{res}{Figure resolution, with the unit of pixels per inch 
      (ppi); the default is 600.}
    \item{MarkGene}{To mark genes in plot or not; if "TRUE" is selected, a  
       file, namely "Reference information to mark gene.csv", that contains 
       the x and y axis information of all the significant QTNs will  
       generate. The default is "FALSE", indicating that no candidate or  
       known gene names are marked in Manhattan plot.}
    \item{Pos_x}{Numeric vectors of x axis where the text labels   
      should be written.}
    \item{Pos_y}{Numeric vectors of y axis where the text labels 
      should be written.}
    \item{GeneName}{A character vector or expression specifying the text 
      to be written.}
    \item{GeneNameColour}{The colour of gene names.}
    \item{...}{Arguments passed to points, axis, text.}
   }
\author{
Zhang Ya-Wen, Li Pei, and Zhang Yuan-Ming\cr
Maintainer: Yuan-Ming Zhang<soyzhang@mail.hzau.edu.cn>
}
\examples{
inter<-data(ResultIntermediate)
fin<-data(ResultFinal)
MultiManhattan(ResultIntermediate=ResultIntermediate,ResultFinal=ResultFinal,dirplot=tempdir())
}