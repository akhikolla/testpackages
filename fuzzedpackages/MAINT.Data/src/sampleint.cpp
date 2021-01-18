#include <cmath>
#include "sampleint.h"

#ifndef SCOPE
RNGScope scope;
#define SCOPE
#endif

using namespace std;

void sampleint(const int n,const int size,vector<int>& res,const bool C0cnv)
{
//	static vector<double> rndnmbs;
//	if (rndnmbs.size()!=n) rndnmbs.resize(n);
	vector<double> rndnmbs(n);
	rndnmbs = as<vector<double> >(runif(n,0.,1.));
	highestval<double>(n,size,rndnmbs,res,C0cnv);
	return;
}

