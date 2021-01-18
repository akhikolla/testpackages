#include "aag.h"

using namespace arma;
using namespace Rcpp;

int getsign(double x) {
	if (x > 0) {
		return 1;
	}
	else if (x < 0) {
		return -1;
	}
	else {
		return 0;
	}
}

uvec findzerocol(mat l, int p) {
	uvec zerc(p, fill::zeros);
	uvec::iterator zp = zerc.begin();
	for (int i = 0; i < p; ++i) {
		l(i, i) = 0;
		*zp++ = all(l.col(i) == 0);
	}
	return zerc;
}
