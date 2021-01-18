#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <functional>
#include <fstream>
#include <iostream>
#include <limits>
#include <vector>

#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <Eigen/Dense>
#include <Eigen/LU>
#include <Eigen/SVD>

using namespace std;
using namespace Eigen;
using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::RowVectorXd;

extern double whimed_i(
				VectorXd& a,
				VectorXi& w,
				int n,
				VectorXd& a_cand,
				VectorXd& a_psrt,
				VectorXi& w_cand
			);
extern double qn(VectorXd& y);
extern double Fmedian(VectorXd& x);
extern double pCalc(
				VectorXd& x,
				double (*pCalcMethod)(VectorXd&)
			);
extern double scaleTau2(VectorXd& x);
MatrixXd FFOgkBasis(
			const MatrixXd x,
			const int& calcM
		){
	//http://forum.codecall.net/topic/48576-function-pointers/#axzz2HU9vi2IB
	double (*pFo[])(VectorXd&)={&qn,&scaleTau2}; 
	const int p=x.cols();
	const int n=x.rows();
	double b1=0.0,b2=0.0;
	VectorXd dvec1=VectorXd::Ones(p);
	MatrixXd U=dvec1.asDiagonal();
	VectorXd sY(n);
	VectorXd dY(n);
	VectorXd sy(n);
	VectorXd sz(n);
	VectorXd sYi(n);
	VectorXd sYj(n);
	for(int i=0;i<p;++i){
		sYi=x.col(i);
		for(int j=0;j<i;++j){
			sYj=x.col(j);
			sY=sYi+sYj;
			dY=sYi-sYj;
			b1=pCalc(sY,pFo[calcM]);
			b1*=b1;
			b2=pCalc(dY,pFo[calcM]);
			b2*=b2;
			U(i,j)=0.25*(b1-b2);
			U(j,i)=U(i,j);	
		}		
	}
	return(U);
}
extern "C"{
	void R_FastOGK(
				int* n,
				int* p,
				double* X,
				double* Q,
				int* cMet
			){
		const int CalcMet=*cMet;
		MatrixXd xi=Map<MatrixXd>(X,*n,*p);		
		MatrixXd Um=FFOgkBasis(xi,CalcMet);
		Map<MatrixXd>(Q,Um.rows(),Um.cols())=Um;
	}
}
