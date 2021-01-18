#ifndef _QUIC_H_
#define _QUIC_H_

#include <algorithm>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <vector>
//#ifdef _OPENMP

//#endif
#ifdef __cplusplus
extern "C" {
#endif
#include "metis.h"
#include "omp.h"
#ifdef __cplusplus
}
#endif
#include <time.h>
using namespace std;


void NormalizeData(int p, int n, double *samples, double *samples_new, vector<int> &mapping);
double computeSij(const double *samples, long p, long n, long i, long j);

double innerproduct(vector<double> &x, vector<double> &y);
void vector_plus(vector<double> &x, vector<double> &y, vector<double> &z, double c);

class smat_t {
	public:
		long p;
		long nnz;
		vector<double> values;
		vector<long> row_ptr;
		vector<long> col_idx;
		int is_symmetric;
		vector<long> id_map;

		smat_t()
		{
			is_symmetric = 0;
		}


		smat_t(smat_t &X, smat_t &D, double alpha)
		{
			smat_t();
			is_symmetric = 0;
			p = X.p;
			nnz  = 0;
			values.resize(X.nnz+D.nnz);
			row_ptr.resize(p+1);
			col_idx.resize(X.nnz+D.nnz);
			for ( long i=0 ; i<p ; i++ )
			{
				row_ptr[i] = nnz;
				long idx1 = X.row_ptr[i];
				long idx2 = D.row_ptr[i];
				while ( (idx1 <X.row_ptr[i+1]) && (idx2<D.row_ptr[i+1]))
				{
					if ( X.col_idx[idx1] < D.col_idx[idx2] )
					{
						col_idx[nnz] = X.col_idx[idx1];
						values[nnz] = X.values[idx1];
						idx1++, nnz++;
					}
					else if ( X.col_idx[idx1] > D.col_idx[idx2])
					{
						col_idx[nnz] = D.col_idx[idx2];
						values[nnz] = alpha*D.values[idx2];
						idx2++, nnz++;
					}
					else
					{
						col_idx[nnz] = X.col_idx[idx1];
						values[nnz] = X.values[idx1]+alpha*D.values[idx2];
						idx1++, idx2++, nnz++;
					}
				}
				if ( idx1 < X.row_ptr[i+1] )
				{
					for ( ; idx1<X.row_ptr[i+1] ; idx1++, nnz++ )
					{
						col_idx[nnz] = X.col_idx[idx1];
						values[nnz] = X.values[idx1];
					}
				}
				else if ( idx2 < D.row_ptr[i+1])
				{
					for ( ; idx2<D.row_ptr[i+1] ; idx2++, nnz++ )
					{
						col_idx[nnz] = D.col_idx[idx2];
						values[nnz] = alpha*D.values[idx2];
					}
				}
			}
			row_ptr[p] = nnz;
		}

		smat_t(long p_) 
		{
			smat_t();
			p = p_;
			is_symmetric = 0;
		}

		smat_t(smat_t &X)
		{
			smat_t();
			p = X.p;
			nnz = X.nnz;
			values = X.values;
			row_ptr = X.row_ptr;
			col_idx = X.col_idx;
			is_symmetric = 0;
		}

		int identity(long p_)
		{
			p = p_;
			nnz = p;
			values.resize(p);
			row_ptr.resize(p+1);
			col_idx.resize(p);
			for ( long i=0 ; i<p ; i++)
			{
				values[i] =1;
				col_idx[i] = i;
				row_ptr[i] = i;
			}
			row_ptr[p] = p;
			return 1;
		}

		int reset() 
		{
			values.clear();
			col_idx.clear();
			row_ptr.clear();
			nnz = 0;
			return 1;
		}

		// Form the graph with node id in id_map. 
		// Cannot work when the id_map is not in increasing order. 
		void form_originalgraph()
		{
			int realp = id_map[p-1]+1;

			for ( long idx=0 ; idx<nnz ; idx++ )
				col_idx[idx] = id_map[col_idx[idx]];

			vector<long> tmp = row_ptr;
			row_ptr.resize(realp+1);
			for ( long i=0 ; i<=id_map[0] ; i++)
				row_ptr[i] = 0;
			long lastid = id_map[0];
			for ( long i=1 ; i<p ; i++ )
			{
				long nowid = id_map[i];
				for ( long ii=lastid+1 ; ii<=nowid ; ii++)
					row_ptr[ii] = tmp[i];
				lastid = nowid;
			}
			for ( long i=id_map[p-1]+1 ; i<=realp ; i++)
				row_ptr[i] = nnz;

			p = realp;
		}

		void dfs(long &ncluster, vector<long> &cluster_ind)
		{
			vector<long> stack (p);
			long top;
			ncluster = 0;
			cluster_ind.resize(p);
			for ( long i=0 ; i<p ; i++ )
				cluster_ind[i] = 0;
			for  ( long i=0 ; i<p ; i++ )
				if ( cluster_ind[i] == 0)
				{
					ncluster++;
					cluster_ind[i] = ncluster;
					stack[0] = i;
					top = 1;
					while (top >0 )
					{
						top--;
						long nowi = stack[top];
						for ( long idx = row_ptr[nowi] ; idx < row_ptr[nowi+1] ; idx++ )
						{
							long nowj = col_idx[idx];
							if ( cluster_ind[nowj] == 0 )
							{
								stack[top] = nowj;
								top++;
								cluster_ind[nowj] = ncluster;
							}
						}
					}
				}
		}

		void ComputeAx(vector<double> &x, vector<double> &Ax)
		{
			ComputeAx(x, Ax, p);
		}

		void ComputeAx_omp(vector<double> &x, vector<double> &Ax, long p_)
		{
			Ax.resize(p_);
			for ( long i=0 ; i<p_ ; i++ )
				Ax[i] = 0;
				for (long i=0 ; i<p_ ; i++ )
				{
					for ( long idx = row_ptr[i] ; idx < row_ptr[i+1] ; idx++)
					{
						long j = col_idx[idx];
						double v = values[idx];
						Ax[i] += v*x[j];
						if ( i!=j )
							Ax[j] += v*x[i];
					}
				}
		}


		void ComputeAx(vector<double> &x, vector<double> &Ax, long p_)
		{
			Ax.resize(p_);
			for ( long i=0 ; i<p_ ; i++ )
				Ax[i] = 0;
			if ( is_symmetric != 1 )
			{
				for (long i=0 ; i<p_ ; i++ )
					for ( long idx = row_ptr[i] ; idx < row_ptr[i+1] ; idx++)
					{
						long j = col_idx[idx];
						double v = values[idx];
						Ax[i] += v*x[j];
//						if ( is_symmetric != 1)
							if ( i!=j )
								Ax[j] += v*x[i];
					}
			}
			else
			{
				for (long i=0 ; i<p_ ; i++ )
				{
					double tmp=0;
					for ( long idx = row_ptr[i] ; idx < row_ptr[i+1] ; idx++)
					{
						long j = col_idx[idx];
						double v = values[idx];
						tmp += v*x[j];
					}
					Ax[i] = tmp;
				}
			}
		}

		int ComputeAinvb_omp(vector<double> &b, vector<double> &x, long p_, double tol = 1e-15)
		{
			vector<double>  r(p_), pp(p_), Ax_result(p_);
			int maxiter = 15;

			// Initial from x = 0
			for ( int i=0 ; i<p_ ; i++)
				x[i] = 0;
			r = b;
			pp = r;
			double r_norm = innerproduct(r,r);
			double initial = r_norm;
//			printf("initial rnorm: %lf\n", r_norm);
			if ( r_norm <1e-13)
				return 0;
			int iter;
			for ( iter = 0; iter<maxiter ; iter++)
			{
				ComputeAx_omp(pp, Ax_result, p_);
				double alpha = innerproduct(r, r)/innerproduct(Ax_result, pp);
				vector_plus(x, x, pp, alpha);
				vector_plus(r, r, Ax_result, (-1)*alpha);
				double r_norm_now = innerproduct(r,r);
//				printf("residual: %lf\n", r_norm_now);
				if ( r_norm_now < tol*initial)
					break;
				double beta = r_norm_now/r_norm;
				r_norm = r_norm_now;
				vector_plus(pp, r, pp, beta);
			}
			return (iter+1);
		}


		int ComputeAinvb(vector<double> &b, vector<double> &x, double tol = 1e-5)
		{
			return ComputeAinvb(b, x, p, tol);
		}

		int ComputeAinvb(vector<double> &b, vector<double> &x, long p_, double tol = 1e-15)
		{
			vector<double>  r(p_), pp(p_), Ax_result(p_);
			int maxiter = 15;

			// Initial from x = 0
			for ( int i=0 ; i<p_ ; i++)
				x[i] = 0;
			r = b;
			pp = r;
			double r_norm = innerproduct(r,r);
			double initial = r_norm;
//			printf("initial rnorm: %lf\n", r_norm);
			if ( r_norm <1e-13)
				return 0;
			int iter;
			for ( iter = 0; iter<maxiter ; iter++)
			{
				ComputeAx(pp, Ax_result, p_);
				double alpha = innerproduct(r, r)/innerproduct(Ax_result, pp);
				vector_plus(x, x, pp, alpha);
				vector_plus(r, r, Ax_result, (-1)*alpha);
				double r_norm_now = innerproduct(r,r);
//				printf("residual: %lf\n", r_norm_now);
				if ( r_norm_now < tol*initial)
					break;
				double beta = r_norm_now/r_norm;
				r_norm = r_norm_now;
				vector_plus(pp, r, pp, beta);
			}
			return (iter+1);
		}

		int ComputeLogdet_serial(double &result, double tol=1e-5)
		{
			double logdet = 0;
			// check if all diagonals are positive
			int pd = 1;
			long i;
			for ( i=0 ; i<p ; i++ ) 
			{
				if ( row_ptr[i+1]-1 < 0)
					break;
				if ( col_idx[row_ptr[i+1]-1] != i)
					break;
				if ( values[row_ptr[i+1]-1] < 0)
					break;
			}
			if ( i<p )
				return 0;
		
			logdet = log(values[0]);
			for ( long pend=1 ; pend<p ; pend++)
			{
				// compute
				double tmp = values[row_ptr[pend+1]-1];
				vector<double> b (pend,0);
				for ( long idx = row_ptr[pend] ; idx < row_ptr[pend+1]-1 ; idx++ )
					b[col_idx[idx]] = values[idx];

				if ( pend == 1 )
					tmp = tmp - b[0]*b[0]/values[0];
				else
				{
					vector<double> Ainvb (pend,0);
					ComputeAinvb_omp(b, Ainvb, pend, tol);
					tmp = tmp - innerproduct(b, Ainvb);
				}

				if ( tmp <= 0)
					return 0;
				logdet += log(tmp);
//				if ( pend %100 == 0)
//					printf("logdet iter %d: %lf\n", pend, logdet);
			}
			result = logdet;
			return 1;
		}

		int ComputeLogdet(double &result, double tol=1e-5)
		{
			double logdet = 0;
			// check if all diagonals are positive
			int pd = 1;
			long i;
			for ( i=0 ; i<p ; i++ ) 
			{
				if ( row_ptr[i+1]-1 < 0)
					break;
				if ( col_idx[row_ptr[i+1]-1] != i)
					break;
				if ( values[row_ptr[i+1]-1] < 0)
					break;
			}
			if ( i<p )
				return 0;
		
			logdet = 0;
			int errorflag = 0;
//#pragma omp parallel for shared(errorflag) reduction(+:logdet)
			for ( long pend=1 ; pend<p ; pend++)
			{
				if (errorflag==1) continue;
				// compute
				double tmp = values[row_ptr[pend+1]-1];
				vector<double> b (pend,0);
				for ( long idx = row_ptr[pend] ; idx < row_ptr[pend+1]-1 ; idx++ )
					b[col_idx[idx]] = values[idx];

				if ( pend == 1 )
					tmp = tmp - b[0]*b[0]/values[0];
				else
				{
					vector<double> Ainvb (pend,0);
					ComputeAinvb_omp(b, Ainvb, pend, tol);
					tmp = tmp - innerproduct(b, Ainvb);
				}

				if ( tmp <= tol)
					errorflag = 1;
//					return 0;
				logdet += log(tmp);
//				if ( pend %100 == 0)
//					printf("logdet iter %d: %lf\n", pend, logdet);
			}
			logdet += log(values[0]);
			result = logdet;
			if (errorflag == 1)
				return 0;
			else
				return 1;
		}
	
		double l1norm()
		{
			double result =0;
			for ( long i=0 ; i<p ; i++ )
				for ( long idx = row_ptr[i] ; idx<row_ptr[i+1] ; idx++ )
				{
					if ( col_idx[idx] == i )
						result += fabs(values[idx]);
					else
						result += fabs(values[idx])*2;
				}
			return result;
		}

		double ComputetrSX(const double *samples, long n)
		{
			double result = 0;
			const double *si, *sj;
			for ( long i=0 ; i<p ; i++ )
			{
				si = &(samples[n*i]);
				for ( long idx = row_ptr[i] ; idx<row_ptr[i+1] ; idx++ )
				{
					sj = &(samples[n*col_idx[idx]]);
					double tmp = 0;
					for ( long t=0 ; t<n ; t++ )
						tmp += si[t]*sj[t];
					if ( col_idx[idx] == i )
						result += values[idx]*tmp;
					else
						result += values[idx]*tmp*2;
				}
			}
			return result;
		}

		int copyfrom(smat_t &X)
		{
			p = X.p;
			nnz = X.nnz;
			values = X.values;
			row_ptr = X.row_ptr;
			col_idx = X.col_idx;
      return 0;
		}

		// make a symmetric matrix from lower triangular matrix
		int symmetricfrom(smat_t &X)
		{
			is_symmetric = 1;
			p=X.p;

			// Copy the idmap if there exists one
			if ( (unsigned)X.id_map.size() == X.p )
			{
				id_map.resize(X.p);
				for ( long i=0 ;i<X.p ; i++ )
					id_map[i] = X.id_map[i];
			}

			nnz = 0;
			row_ptr.resize(p+1,0);
			col_idx.resize(2*X.nnz);
			values.resize(2*X.nnz);
			vector<long> tmp(p);
			for ( long i=0 ; i<p ; i++ )
				tmp[i] = X.row_ptr[i];
			for (long i=0 ; i<p ; i++)
			{
				row_ptr[i] = nnz;
				for ( long idx = X.row_ptr[i] ; idx<X.row_ptr[i+1] ; idx++ )
				{
					col_idx[nnz] = X.col_idx[idx];
					values[nnz] = X.values[idx];
					nnz++;
				}
				for ( long j=i+1 ; j<p ; j++ )
				{
					if ( tmp[j] < X.row_ptr[j+1] )
					{
						if ( X.col_idx[tmp[j]] == i)
						{
							col_idx[nnz] = j;

							values[nnz] = X.values[tmp[j]];

							nnz++;
							tmp[j]++;
						}
					}
				}
			}
			row_ptr[p] = nnz;
      return 0;
		}

		void print(FILE *fp)
		{
			fprintf(fp, "p: %ld, nnz: %ld\n", p, nnz);
			for ( long i=0 ; i<p ; i++ )
			{
				for ( long idx = row_ptr[i] ; idx<row_ptr[i+1] ; idx++ )
				{
					if ( (unsigned)id_map.size() == p )
						fprintf(fp, "%ld %ld %f\n", id_map[i]+1, id_map[col_idx[idx]]+1, values[idx]);
					else
						fprintf(fp, "%ld %ld %f\n", i+1, col_idx[idx]+1, values[idx]);
				}
			}
		}


		void addblock(smat_t &A)
		{
			values.resize(nnz + A.nnz);
			for ( long i=0 ; i<A.nnz ; i++ )
				values[i+nnz] = A.values[i];
			col_idx.resize(nnz+A.nnz);
			for ( long i=0 ; i<A.nnz ; i++ )
				col_idx[i+nnz] = A.col_idx[i]+p;
			row_ptr.resize(p+A.p+1);
			for ( long i=0 ; i<=A.p ; i++ )
				row_ptr[i+p] = A.row_ptr[i]+nnz;
 			p += A.p;
			nnz += A.nnz;
		}

		double comp_within_ratio(vector<long> &block_ind)
		{
			double within_w= 0, total_w = 0;
			for ( long i=0 ; i<p ; i++ )
				for ( long idx = row_ptr[i] ; idx<row_ptr[i+1] ; idx++ )
				{
					if ( block_ind[i] == block_ind[col_idx[idx]])
						within_w += fabs(values[idx]);
					total_w += fabs(values[idx]);
				}
			return within_w/total_w;
		}
/*
		void random_graph(long p_, double rate_)
		{
			p = p_;
			nnz = 0;
			srand(0);
			col_idx.resize(int(p*p*rate_/2)+100);
			values.resize(int(p*p*rate_/2)+100);
			row_ptr.resize(p+1);
			for (long i=0 ; i<p ; i++ )
			{
				long subp = int((i+1)*rate_);
				row_ptr[i] = nnz;
				if ( subp > 0 )
				{
					vector<long> rowii (subpb);
					for ( long j=0 ; j<subp ; j++ )
						rowii[j] = rand()%p;
					sort(rowii.begin(), rowii.end());
	
					col_idx[nnz] = rowii[0];
					values[nnz] = rand()%2;
					if (values[nnz] == 0) values[nnz] = -1;
					values[nnz]*=0.1;
					nnz ++;

					for ( long j=1 ; j<subp ; j++ )
						if ( rowii[j] != rowii[j-1] )
						{
							col_idx[nnz] = rowii[j];
							values[nnz] = rand()%2;
							if ( values[nnz] == 0 ) values[nnz] = -1;
							values[nnz]*=0.1;
							nnz++;
						}
				}
				col_idx[nnz] = i;
				values[nnz] = 1;
				nnz++;
			}
			row_ptr[p] = nnz;
		}
*/
		void clustering(vector<long> &block_ind, long nblock)
		{
			int nvtxs = p ;
			int wgtflag = 0, numflag = 0,nparts = 1, options[11], edgecut=1  ;

			block_ind.resize(p);
			idxtype *xadj = (idxtype *)malloc(sizeof(idxtype)*(p+1));
			for ( long i=0 ; i<=p ; i++)
				xadj[i] = row_ptr[i];
			idxtype *adjncy = (idxtype *)malloc(sizeof(idxtype)*nnz);
			for ( long idx=0 ; idx<nnz ; idx++ )
				adjncy[idx] = col_idx[idx];
			idxtype *vwgt = (idxtype *)malloc(sizeof(idxtype)*p);
			for ( long i=0 ; i<p ; i++ )
				vwgt[i] = 1;
			idxtype *adjwgt = (idxtype *)malloc(sizeof(idxtype)*nnz);
			for ( long idx=0 ; idx<nnz ; idx++ )
				adjwgt[idx] = int(values[idx]*1000);
			nparts = nblock;
//idxtype *part = (idxtype *)malloc(sizeof(idxtype)*p);
			options[0] = 0;
      char passToIdxmalloc[] = "main: part";
			idxtype *part = idxmalloc(p, passToIdxmalloc);
			METIS_PartGraphKway(&nvtxs, xadj, adjncy, vwgt, adjwgt, &wgtflag, &numflag, &nparts, options, &edgecut, part); 
			for ( long i=0 ; i<p ; i++ )
				block_ind[i] = part[i];
			free(xadj);
			free(adjncy);
			free(vwgt);
			free(adjwgt);
			free(part);
		}
};

extern "C" {
void QUIC(int p, int n, double* samples, double lambda, double tol, int msg, int maxIter, int nblock, int numthreads, smat_t& X, vector<double> &objlist, vector<double> &timelist);
}

#endif
