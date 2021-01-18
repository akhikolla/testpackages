#ifdef _OPENMP
#include <omp.h>
#endif
#include "covariance.h"
#include "aca.h"
#include "tlr.h"

using namespace std;
using namespace Eigen;

void tlr_aca_covM(const Eigen::MatrixXd &geom, std::vector<Eigen::MatrixXd> &B,
        std::vector<TLRNode> &UV, std::function<double(double)> kernel, int m,
        const std::vector<int> &idxVec, double epsl, int allocSz)
{
        int n = geom.rows();
        int nb = n / m;
	int lastBlkDim = n % m;
        B.resize(nb + (lastBlkDim > 0));
        // diag blk
	#ifdef _OPENMP
        #pragma omp parallel for
	#endif
        for( int l = 0 ; l < nb ; l++ )
        {
                int i1 = l * m;
                B[l].resize( m , m );
                for( int i = 0 ; i < m ; i++ ) 
                        for( int j = 0 ; j <= i ; j++ )
                        {
                                double x = (geom.row(idxVec[i1+i]) - 
                                        geom.row(idxVec[i1+j])).norm();
                                B[l](i , j) = kernel( x );
                        }
        }
	if(lastBlkDim > 0)
	{
		B[nb] = MatrixXd::Zero(m, m);
		for( int i = 0 ; i < lastBlkDim ; i++ ) 
                        for( int j = 0 ; j <= i ; j++ )
                        {
				int i1 = nb * m;
                                double x = (geom.row(idxVec[i1+i]) - 
                                        geom.row(idxVec[i1+j])).norm();
                                B[nb](i , j) = kernel( x );
                        }
                B[nb].diagonal().segment(lastBlkDim, m-lastBlkDim).setOnes();
        }
        // off-diag part, store lower-tri only
	nb += (lastBlkDim > 0);
        int n_lrt = nb * (nb - 1) / 2;
	UV.resize(n_lrt);
	#ifdef _OPENMP
        #pragma omp parallel for 
	#endif
        for( int l = 0 ; l < n_lrt ; l++ )
        {
                int i = l % nb;
                int j = l / nb;
                if( i <= j )
                {
                        j = nb - j - 2;
                        i = nb - i - 1;
                }
                int blkIdx = (2 * nb - 1 - j) * j / 2 + ( i - 1 - j );
		UV[blkIdx].U.resize(m,allocSz);
		UV[blkIdx].V.resize(m,allocSz);
		UV[blkIdx].maxColNum = allocSz;
                function< double( int , int ) > kernel_aca =
                        [&i, &j, &kernel, &idxVec, &m, &n, &geom]
                        (int i2 , int j2)
                        {return i*m + i2 < n && j*m + j2 < n ? 
				kernel((geom.row(idxVec[i*m+i2]) - 
                                geom.row(idxVec[j*m+j2])).norm()) : 0.0;};
                aca(kernel_aca , m , epsl , UV[blkIdx]);
        }
}

Eigen::MatrixXd dense_covM(const Eigen::MatrixXd &geom, 
	std::function<double(double)> kernel)
{
	int n = geom.rows();
	MatrixXd covM(n,n);
	for(int i = 0; i < n; i++)
		for(int j = 0; j <= i; j++)
			covM(i,j) = kernel((geom.row(i) - geom.row(j)).norm());
	return covM;
}

void tlr_aca_covM(const Eigen::MatrixXd &covM, std::vector<Eigen::MatrixXd> &B,
        std::vector<TLRNode> &UV, int m, double epsl, int allocSz)
{
        int n = covM.rows();
        int nb = n / m;
	int lastBlkDim = n%m;
        B.resize(nb + (lastBlkDim > 0));
        // diag blk
	#ifdef _OPENMP
        #pragma omp parallel for
	#endif
        for( int l = 0 ; l < nb ; l++ )
		B[l] = covM.block(l*m,l*m,m,m);
	if(lastBlkDim > 0)
	{
		B[nb] = MatrixXd::Zero(m, m);
		B[nb].block(0,0,lastBlkDim,lastBlkDim) = covM.block(nb*m,nb*m,
			lastBlkDim,lastBlkDim);
		B[nb].diagonal().segment(lastBlkDim, m-lastBlkDim).setOnes();
	}
        // off-diag part, store lower-tri only
	nb += (lastBlkDim > 0);
        int n_lrt = nb * (nb - 1) / 2;
	UV.resize(n_lrt);
	const double *covMPtr = covM.data();
	#ifdef _OPENMP
        #pragma omp parallel for 
	#endif
        for( int l = 0 ; l < n_lrt ; l++ )
        {
                int i = l % nb;
                int j = l / nb;
                if( i <= j )
                {
                        j = nb - j - 2;
                        i = nb - i - 1;
                }
                int blkIdx = (2 * nb - 1 - j) * j / 2 + ( i - 1 - j );
		UV[blkIdx].U.resize(m,allocSz);
		UV[blkIdx].V.resize(m,allocSz);
		UV[blkIdx].maxColNum = allocSz;
                function<double(int , int)> kernel_aca =
                        [&i, &j, &m, &n, &covMPtr]
                        (int i2 , int j2)
                        {return i*m + i2 < n && j*m+j2 < n ?
				*(covMPtr + i*m + i2 + n*(j*m+j2)) : 0.0;};
                aca(kernel_aca , m , epsl , UV[blkIdx]);
        }
}
