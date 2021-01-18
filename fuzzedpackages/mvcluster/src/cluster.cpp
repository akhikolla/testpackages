#include "MvSsvd.h"
#include "MvLrmaL1.h"
#include "utils.h"
#include "MvLrmaL0.h"

#define ARMA_DONT_USE_CXX11
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::export]]
List clusterl0(vector<mat> data, uvec svs, int sz, int seed, int maxIter, double thres, int debug)
{
	MvLrmaL0 mvlrmal0(data, sz, svs);
	mvlrmal0.setDebugLevel(debug);
	mvlrmal0.setISeedFeat(seed);
	mvlrmal0.setMaxIter(maxIter);
	mvlrmal0.setThreshold(thres);

	clock_t begin = clock();
	mvlrmal0.clustering(); // run clustering
	clock_t end = clock();
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	if(debug>0)
	{
		Rprintf("time elapsed in seconds: %.2f\n", elapsed_secs);
	}

	return List::create(Named("Cluster") = mvlrmal0.getCluster(),
			Named("FeatClusters") = mvlrmal0.getFeatClus(),
			Named("U") = mvlrmal0.getMatU(),
			Named("V") = mvlrmal0.getPVecV(),
			Named("z") = mvlrmal0.getVecZ());
}

List clusterl1(vector<mat> data, fvec lus, fvec lvs, float lz, int maxOuter, double thresOuter, int maxInner, double thresInner, int debug)
{
	MvLrmaL1 mvlrmal1(data, lz, lus, lvs);
	mvlrmal1.setMaxOuterIter(maxOuter);
	mvlrmal1.setThresholdOuter(thresOuter);
	mvlrmal1.setMaxInnerIter(maxInner);
	mvlrmal1.setThresholdInner(thresInner);
	mvlrmal1.setDebugLevel(debug);

	clock_t begin = clock();
	mvlrmal1.clustering(); // run clustering
	clock_t end = clock();
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	if(debug>0)
	{
		Rprintf("time elapsed in seconds: %.2f\n", elapsed_secs);
	}

	return List::create(Named("Cluster") = mvlrmal1.getCluster(),
				Named("FeatClusters") = mvlrmal1.getFeatClus(),
				Named("U") = mvlrmal1.getMatU(),
				Named("V") = mvlrmal1.getV(),
				Named("z") = mvlrmal1.getVecZ());
}
List cMvSsvd(vector<mat> data, fvec lvs, float lz, int maxOuter, double thresOuter, int maxInner, double thresInner, int debug)
{
	MvSsvd mvssvd(data, lz, lvs);
	mvssvd.setMaxOuterIter(maxOuter);
	mvssvd.setThresholdOuter(thresOuter);
	mvssvd.setMaxInnerIter(maxInner);
	mvssvd.setThresholdInner(thresInner);
	mvssvd.setDebugLevel(debug);

	clock_t begin = clock();
	mvssvd.clustering(); // run clustering
	clock_t end = clock();
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	if(debug>0)
	{
		Rprintf("time elapsed in seconds: %.2f\n", elapsed_secs);
	}

	vector<vec> z;
	vector<vec*,allocator<vec*>> za = mvssvd.getV();
	for(uint32_t i=0;i<za.size();i++)
	{
		vec ze = *za.at(0);
		z.push_back(ze);
	}

	return List::create(Named("Cluster") = mvssvd.getCluster(),
				Named("FeatClusters") = mvssvd.getFeatClus(),
				Named("U") = mvssvd.getU(),
				Named("V") = z,
				Named("z") = mvssvd.getZ());
}

