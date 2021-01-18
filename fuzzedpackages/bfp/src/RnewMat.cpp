#include "RnewMat.h"
#include "dataStructure.h"
#include <vector>
#include <numeric>
#include <algorithm>

using std::set;
using std::vector;
using std::accumulate;


ReturnMatrix getMultipleCols(const Matrix& M, const set<int>& s) // get different concatenated columns of matrix 
{
	
	
	// MATRIXSTORE(M, MStore)
	
	Matrix ret(M.Nrows(), s.size());

	set<int>::size_type cols = 1; // invariant: about to process column number cols
	for (set<int>::const_iterator i = s.begin(); i != s.end(); i++){
		ret.Column(cols++) = M.Column(*i);	
	}

	ret.Release(); return ret;
}


