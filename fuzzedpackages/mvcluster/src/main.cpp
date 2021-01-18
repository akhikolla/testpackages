#ifdef INSIDE
#define ARMA_DONT_USE_CXX11
#include <RcppArmadillo.h>
#include <RInside.h>                    // for the embedded R via RInside
#include "cluster.h"
using namespace Rcpp;
using namespace std;
int main(int argc, char **argv) {
    RInside R(argc, argv);              // create an embedded R instance
    int s = 2;//cluster(320,0,0);
    Language call("print",s);
    call.eval();
    return 0;
}
#endif


/*
 * main.cpp
 *
 *  Created on: Aug 18, 2015
 *      Author: Javon
 */

//#ifdef INSIDE
//#include <RcppArmadillo.h>
//#include <RInside.h>                    // for the embedded R via RInside
////#include "rcpp_hello_world.h"
//
//#include "MvSsvd.h"
//#include "MvLrmaL1.h"
//#include "utils.h"
//#include "MvLrmaL0.h"
//
//using namespace Rcpp;
//using namespace arma;
//using namespace std;
//
//int main(int argc, char *argv[]) {
//	mat phe;
//	phe.load("/home/jrm10008/Downloads/mfeat-pix.csv", csv_ascii); // load phenotype
//	mat gen;
//	gen.load("/home/jrm10008/Downloads/mfeat-fou.csv", csv_ascii); // load genotype
//	mat thi;
//	thi.load("/home/jrm10008/Downloads/mfeat-fac.csv", csv_ascii); // load genotype
//	vector<mat> datasets(3);
//	datasets[0] = phe;
//	datasets[1] = gen;
//	datasets[2] = thi;
//
//	uint32_t sz = 200;
//	uvec svs(3);
//	svs[0] = 200;
//	svs[1] = 60;
//	svs[2] = 22;
//	MvLrmaL0 mvlrmal0(datasets, sz, svs);
//	mvlrmal0.setDebugLevel(1);
////	mvlrmal0.setISeedFeat(0);
//
//	clock_t begin = clock();
//	mvlrmal0.clustering(); // run clustering
//	clock_t end = clock();
//	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
//	printf("time elapsed in seconds: %.2f\n", elapsed_secs);
//
//    return 0;
//}
//#endif
