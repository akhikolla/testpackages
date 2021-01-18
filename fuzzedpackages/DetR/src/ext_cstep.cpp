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
#include <Eigen/LU>
#include <Eigen/SVD>

using namespace std;
using namespace Eigen;
using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::RowVectorXd;

extern VectorXi CStep(
			const VectorXd& dPin,
			const MatrixXd& x,
			const VectorXd& y,
			const VectorXi& h,
			VectorXd& objfun,
			const int& I,
			VectorXi& citer
		);

extern "C"{
	void R_extCstep(int* rn,	//1
			int* rp,	//2
			double* X,	//3
			int* h,		//4
			double* dst,	//5
			double* Ofuns,	//6
			int* Citer,	//7
			int* P		//8
		){
		const int n=*rn,p=*rp,rh=*h;
		VectorXd dY=Map<VectorXd>(dst,n);
		MatrixXd xi=Map<MatrixXd>(X,n,p);
		VectorXi hi(2);
		hi<<rh,rh;
		VectorXd yi=xi.col(p-1);
		xi.col(p-1).setOnes();
		VectorXi hsub(rh);
		VectorXi ct(2);
		VectorXd of(2);
		hsub=CStep(dY,xi,yi,hi,of,0,ct);
		Map<VectorXi>(P,rh)=hsub.array()+1;
		Map<VectorXd>(Ofuns,2)=of.array().exp();
		Map<VectorXi>(Citer,2)=ct;
	}
}
