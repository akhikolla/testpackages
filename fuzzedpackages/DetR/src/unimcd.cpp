#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <functional>
#include <fstream>
#include <iostream>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <vector>
#include <Eigen/Dense>

using namespace Eigen;
using Eigen::VectorXd;

extern int unimcd(
			VectorXd& y,
			const int& n,
			const int& h,
			const int& len,
			double& initmean,
			double& initcov
		);

extern "C"{
	void R_unimcd(
			double* yin,
			int* h,
			int* len,
			int* n,
			double* inmean,
			double* incov,
			int* inloc
		){
		VectorXd y=Map<VectorXd>(yin,*n);
		int iloc;
		const int In=*n,Ih=*h,Ilen=*len;
		double initmean,initcov;
		iloc=unimcd(y,In,Ih,Ilen,initmean,initcov);
		*incov=initcov;
 		*inmean=initmean;
		*inloc=iloc;
	}
}
/*
Sys.getenv("R_INCLUDE_DIR")	#/usr/share/R/include
Sys.getenv("R_LIBS_SITE")	#/usr/local/lib/R/site-library
###
q("no")
###
cd /home/kaveh/Desktop/work/mia/detlts/detltsR
g++ -I/usr/share/R/include -I/home/kaveh/Desktop/work/p1/geqw4/vi3/out/sp/ccode/eigen -fpic -O3 -c -DEIGEN_NO_DEBUG unimcd.cpp -o unimcd.o
g++ -shared -o unimcd.so unimcd.o -L/usr/local/lib/R/site -Wl,-rpath,/usr/local/lib/R/site -lR
###
R
library(robustbase)
setwd("/home/kaveh/Desktop/work/mia/detlts/detltsR/")
source("UNIMCD.R")
dyn.load("unimcd.so");

unimcd<-function(y,h,len){
#This code calls the UNIMCD.C routine (it's a C translation of the one found in robustbase). 
	inloc<-inmean<-incov<-0.0;
	out<-.C("R_unimcd",as.double(y),as.integer(h),as.integer(len),as.integer(length(y)),as.double(inmean),as.double(incov),as.integer(inloc))
	list(initmean=as.numeric(out[[5]]),initcov=as.numeric(out[[6]]),iloc=as.numeric(out[[7]]))
}
n<-100
set.seed(123)
y<-rnorm(n);
h<-ceiling((n+2)/2)
len<-n-h+1
aa<-UniMCD(y=y,h=h,len=len)
bb<-unimcd(y,h=h,len=len)
library(robustbase)
cc <-robustbase:::.fastmcd(as.matrix(y), as.integer(h), nsamp = 0, nmini = 300)
*/
