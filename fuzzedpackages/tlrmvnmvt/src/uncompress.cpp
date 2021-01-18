#include <RcppEigen.h>
#include <vector>
#include "uncompress.h"

using namespace std;
using namespace Eigen;

Eigen::MatrixXd uncompress(const std::vector<Eigen::MatrixXd> &B, const
        std::vector<TLRNode> &UV, bool symm)
{
	int nb = B.size();
	int m = B[0].rows();
	int n = nb * m;
	MatrixXd M = Eigen::MatrixXd::Zero(n,n);
	// diag
	for(int i = 0; i < nb; i++)
		M.block(i*m,i*m,m,m) = B[i];
	// off-diag
	for(int j = 0; j < nb; j++)
		for (int i = j+1; i < nb; i++)
		{
			int idx = (2 * nb - 1 - j) * j / 2 + (i - 1 - j);
			int ncol = UV[idx].crtColNum;
			M.block(i*m,j*m,m,m) = UV[idx].U.block(0,0,m,ncol) * 
				UV[idx].V.block(0,0,m,ncol).transpose();
		}
	// symm
	if(symm)
		for(int i = 0; i < nb; i++)
			for(int j = i+1; j < nb; j++)
				M.block(i*m,j*m,m,m) = M.block(j*m,i*m,m,m).
					transpose();
	return M;
}

Eigen::MatrixXd uncompress(const std::vector<Eigen::MatrixXd> &B, const
        std::vector<Eigen::MatrixXd> &U, const std::vector<Eigen::MatrixXd> &V, 
	bool symm)
{
	int nb = B.size();
	int m = B[0].rows();
	int n = nb * m;
	MatrixXd M = Eigen::MatrixXd::Zero(n,n);
	// diag
	for(int i = 0; i < nb; i++)
		M.block(i*m,i*m,m,m) = B[i];
	// off-diag
	for(int j = 0; j < nb; j++)
		for (int i = j+1; i < nb; i++)
		{
			int idx = (2 * nb - 1 - j) * j / 2 + (i - 1 - j);
			M.block(i*m,j*m,m,m) = U[idx] * V[idx].transpose();
		}
	// symm
	if(symm)
		for(int i = 0; i < nb; i++)
			for(int j = i+1; j < nb; j++)
				M.block(i*m,j*m,m,m) = M.block(j*m,i*m,m,m).
					transpose();
	return M;
}

