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

extern double Fmedian(
			Ref<VectorXd> y,
			int retMed
		);

extern "C"{
	void R_inFM(int* n,double* X,double* Q){
		VectorXd xi=Map<VectorXd>(X,*n);	
		*Q=Fmedian(xi.col(0),0);
	}
}

