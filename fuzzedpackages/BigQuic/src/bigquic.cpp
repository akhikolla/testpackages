#include <Rcpp.h>
// [[Rcpp::depends(Rcpp)]]
using namespace Rcpp;
// If your compiler is to old, just disable / remove the following line
// [[Rcpp::plugins(cpp11)]]
#include "bigquic.h"
#include <unistd.h>
//#include <random>
#define EPS (double(2.22E-16))

double innerproduct(vector<double> &x, vector<double> &y) {
   long p = x.size();
   double tmp = 0;
   for (long i = 0; i < p; i++)
      tmp += x[i] * y[i];
   return tmp;
}

void vector_plus(vector<double> &x, vector<double> &y, vector<double> &z, double c) {
   long p = x.size();
   for (long i = 0; i < p; i++)
      x[i] = y[i] + z[i] * c;
}

// A p by p sample covariance matrix cannot be stored in memory,
// so we have to compute it on the fly.

double computeSij(const double *samples, long p, long n, long i, long j) {
   const double *si = &(samples[n * i]), *sj = &(samples[n * j]);
   double result = 0;
   for (long t = 0; t < n; t++)
      result += si[t] * sj[t];
   return result;
}

unsigned long IsDiag(const smat_t &X) {
   for (long i = 0; i < X.p; i++)
      for (long idx = X.row_ptr[i]; idx < X.row_ptr[i + 1]; idx++)
         if (X.col_idx[idx] != i)
            return 0;
   return 1;
}

static inline void CoordinateDescentUpdate(
        unsigned long p, const double* const S, const double* const Lambda,
        const double* X, const double* W, double* U, double* D,
        unsigned long i, unsigned long j, double& normD, double& diffD) {
   unsigned long ip = i*p;
   unsigned long jp = j*p;
   unsigned long ij = ip + j;

   double a = W[ij] * W[ij];
   if (i != j)
      a += W[ip + i] * W[jp + j];
   double ainv = 1.0 / a; // multiplication is cheaper than division

   double b = S[ij] - W[ij];
   for (unsigned long k = 0; k < p; k++)
      b += W[ip + k] * U[k * p + j];

   double l = Lambda[ij] * ainv;
   double c = X[ij] + D[ij];
   double f = b*ainv;
   double mu;
   normD -= fabs(D[ij]);
   if (c > f) {
      mu = -f - l;
      if (c + mu < 0.0) {
         mu = -c;
         D[ij] = -X[ij];
      } else {
         D[ij] += mu;
      }
   } else {
      mu = -f + l;
      if (c + mu > 0.0) {
         mu = -c;
         D[ij] = -X[ij];
      } else {
         D[ij] += mu;
      }
   }
   diffD += fabs(mu);
   normD += fabs(D[ij]);
   if (mu != 0.0) {
      for (unsigned long k = 0; k < p; k++)
         U[ip + k] += mu * W[jp + k];
      if (i != j) {
         for (unsigned long k = 0; k < p; k++)
            U[jp + k] += mu * W[ip + k];
      }
   }
}

// Return the objective value.

static inline double DiagNewton(long p, long n, const double* samples,
        const double lambda, smat_t *X, smat_t *D, double &trgradgD, long nblock_more, vector<long> &block_more_ptr) {
   D->row_ptr.resize(p + 1);

   long total = 0; // for maintaining total amount of nonzero D
   double logdet = 0.0;
   double l1normX = 0.0;
   double trSX = 0.0;

   double timebegin; // = omp_get_wtime();

   // For storing columns of S
   vector<vector<double> > Sblock_more;


   // Computing updates
   for (long bi = 0; bi < nblock_more; bi++) {

      long subpi = block_more_ptr[bi + 1] - block_more_ptr[bi];
      Sblock_more.resize(subpi);

      // Compute the current portion of S
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for (long i = block_more_ptr[bi]; i < block_more_ptr[bi + 1]; i++) {
         vector<double> ei(p, 0);
         long nowi = i - block_more_ptr[bi];
         ei[i] = 1;
         Sblock_more[nowi].resize(p);
         for (long j = 0; j <= i; j++)
            Sblock_more[nowi][j] = computeSij(samples, p, n, i, j);
      }

      // Update
      for (long i = block_more_ptr[bi]; i < block_more_ptr[bi + 1]; i++) {
         D->row_ptr[i] = total;
         long nowi = i - block_more_ptr[bi];
         // i!=j
         for (int j = 0; j < i; j++) {
            double b = Sblock_more[nowi][j];
            double ainv = (X->values[i])*(X->values[j]);
            double l = lambda*ainv;
            double f = b*ainv;
            double mu;
            double Dij = 0;
            if (0 > f) {
               mu = -f - l;
               if (mu > 0)
                  Dij = mu;
            } else {
               mu = -f + l;
               if (mu < 0)
                  Dij = mu;
            }
            if (Dij != 0) {
               D->values.push_back(Dij);
               D->col_idx.push_back(j);
               total += 1;
               trgradgD += Dij * b * 2; // why?
            }
         }
         // j==i
         double c = X->values[i];
         logdet += log(c);
         l1normX += fabs(c) * lambda;
         double Sk = Sblock_more[nowi][i];
         trSX += c*Sk;
         double ainv = c*c;
         double b = Sk - 1.0 / c;
         double l = lambda*ainv;
         double f = b*ainv;
         double mu;
         double Dij = 0;
         if (c > f) {
            mu = -f - l;
            if (c + mu < 0.0)
               Dij = -c;
            else
               Dij = mu;
         } else {
            mu = -f + l;
            if (c + mu > 0.0)
               Dij = -c;
            else
               Dij = mu;
         }
         if (Dij != 0) {
            D->values.push_back(Dij);
            D->col_idx.push_back(i);
            total += 1;
            trgradgD += Dij*b;
         }
      }
   }

   D->row_ptr[p] = total;
   D->nnz = total;
   //	printf("diagnewton time %lf\n", omp_get_wtime()-timebegin);
   //fflush(stdout);

   double fX = -logdet + trSX + l1normX;
   return fX;
}

// Normalize the data so that X^TX is the correlation matrix (with all diagonals equal to 1)
// memory of samples_new should be allocated before calling this procedure

void NormalizeData(int p, int n, double *samples, double *samples_new, vector<int> &mapping) {
   // Eliminate the random variables with std<1e-10
   double threshold = 1e-10;

   // number of random variables after elimination
   int subp = 0;

   vector<double> mean(p);
   for (long i = 0; i < p; i++)
      mean[i] = 0;
   for (long i = 0; i < p; i++)
      for (long j = 0; j < n; j++)
         mean[i] += samples[i * n + j];
   for (long i = 0; i < p; i++)
      mean[i] = mean[i] / n;

   vector<double> std(p);
   for (long i = 0; i < p; i++)
      std[i] = 0;
   for (long i = 0; i < p; i++)
      for (long j = 0; j < n; j++) {
         long nownum = i * n + j;
         std[i] = std[i] + (samples[nownum] - mean[i])*(samples[nownum] - mean[i]);
      }
   for (int i = 0; i < p; i++) {
      std[i] = sqrt(std[i] / (n - 1));
      if (std[i] > threshold)
         subp++;
   }
   // Create a mapping to eliminate variables with std < 1e-10
   int ii = 0;
   mapping.resize(subp);
   for (long i = 0; i < p; i++)
      if (std[i] > threshold)
         mapping[ii++] = i;
   // sample_new is the sampels after eliminate variables with std<1e-10
   for (long i = 0; i < subp; i++) {
      long ori_i = mapping[i];
      for (long j = 0; j < n; j++)
         samples_new[i * n + j] = (samples[ori_i * n + j] - mean[ori_i]) / std[ori_i] / sqrt(n - 1);
   }
}

void QUIC(int p, int n, double* samples, double lambda, double tol, int msg, int maxIter, int nblock, int numthreads, smat_t& X, vector<double> &objlist, vector<double> &timelist) {

   //SEED SET IN R BEFORE CALLING... ALWAYS CHECK THIS FOR DEBUGGING!

   //srand(1);

   //Temp Disabled for debugging, SWITCH BACK TO MERSENNE TWISTER FOR FINAL
   //std::random_device rd;
   //std::mt19937 engine(rd());
   //SET SEED = 1 for testing same results
   //std::mt19937 engine(1);

   unsigned long maxNewtonIter = maxIter;
   double cdSweepTol = 0.0000000005;
   unsigned long max_lineiter = 20;
   double fX = 1e+15;
   double fX1 = 1e+15;
   double fXprev = 1e+15;
   double sigma = 0.001;

   int error_occur = 0;
   //	pair_t* activeSet;
   double l1normX = 0.0;
   double trSX = 0.0;
   double logdetX = 0.0;
   #ifdef _OPENMP
   double timeBegin = omp_get_wtime();
   #else
   double timeBegin = 0;
   #endif


   #ifdef _OPENMP
   int maxnumthreads = omp_get_max_threads();
   if (numthreads > maxnumthreads)
      numthreads = maxnumthreads;
   omp_set_num_threads(numthreads);
   #else
   int maxnumthreads = 1;
   #endif


   // Create block indicator why?
   // block_ptr is the cluster indicator for block coordinate descent
   vector<long> block_ptr(nblock + 1);
   for (long i = 0; i < nblock; i++)
      block_ptr[i] = long(p / nblock) * i;
   block_ptr[nblock] = p;

   // Block for temp W and S when computing gradient
   // This is different from block_ptr; the only goal for computing
   // gradient block by block is for parallelization.
   long block_more_size = 5000;
   long nblock_more = p / block_more_size;

   // nblock_more cannot be smaller than nblock (we can store at most
   // nblock*p elements in memory.
   if (nblock_more < nblock)
      nblock_more = nblock;
   vector<long> block_more_ptr(nblock_more + 1);
   for (long i = 0; i < nblock_more; i++)
      block_more_ptr[i] = long(p / nblock_more) * i;
   block_more_ptr[nblock_more] = p;

//   char hostname[40];
//   gethostname(hostname, 40);

   //if (msg >= 0) {
   if (msg >= 1) {
      Rcout << "Start BigQUIC solver" << endl;
      Rcout << "p: " << p << ", n: " << n << ", lambda:" << lambda << ", tol:" << tol << ", msg: " << msg << ", maxiter: " << maxIter << ", nblock:" << nblock << ", numthreads: " << numthreads << endl;
      //Rcout << "hostname:" << hostname << ", p: " << p << ", n: " << n << ", lambda:" << lambda << ", tol:" << tol << ", msg: " << msg << ", maxiter: " << maxIter << ", nblock:" << nblock << ", numthreads: " << numthreads << endl;
      //printf("Start BigQUIC solver\nhostname:%s, p: %d, n: %d, lambda:%lf, tol:%lf, msg: %d, maxiter: %d, nblock:%d, numthreads: %d\n", hostname, p, n, lambda, tol, msg, maxIter, nblock, numthreads);
      //fflush(stdout);
   }
   double quic_time = 0;
   double total_time_cd = 0;
   double total_time_Wcomp = 0;
   double total_time_Ucomp = 0;
   double total_time_FormActive = 0;
   double total_time_linesearch = 0;
   double total_time_gradient = 0;
   double total_scan_freeset = 0;

   unsigned long pathIdx = 0;
   unsigned long NewtonIter = 1;
   // Cache the diagonal of W, which will be reused a lot during CD
   vector<double> Wdiag(p);

   for (; NewtonIter <= maxNewtonIter; NewtonIter++) {
      Rcpp::checkUserInterrupt();
      double newton_loop_time_begin; // = omp_get_wtime();
      double trgradgD = 0.0; // maintain trgradG during coordinate descent
      double normD = 0.0;
      double diffD = 0.0;
      double subgrad = 1e+15;
      smat_t D(p);
      if (NewtonIter == 1 && IsDiag(X)) {
         if (msg >= 1) {
            Rcout << "Newton iteration 1." << endl;
            Rcout << "  X is a diagonal matrix." << endl;
            //printf("Newton iteration 1.\n");
            //printf("  X is a diagonal matrix.\n");
         }
         // Initial D
         D.reset();
         fX = DiagNewton(p, n, samples, lambda, &X, &D, trgradgD, nblock_more, block_more_ptr);
      }
      else {
         // Empty D
         D.row_ptr.resize(p + 1);
         D.nnz = 0;
         D.values.clear();
         D.col_idx.clear();

         long numActive = 0;
         subgrad = 0.0;

         vector<long> buffercolidx(p);
         // Store Sij, Xij, Wij on free set
         vector<double> Slist;
         vector<double> Xlist;
         vector<double> Wlist;


         // Tranfrom the matrix from lower triangular to symmetric (computing matrix vector multiplication is easier to parallel
         smat_t X_sym;
         X_sym.symmetricfrom(X);

         // Construct Free Set
         if (msg >= 1) {
            Rcout << "Constructing Free Set" << endl;
            //printf("Constructing Free Set\n");
            //fflush(stdout);
         }
         // Stopping tolerance for CG when computing gradient
         double g_tol = 1e-10;

         for (int bi = 0; bi < nblock_more; bi++) {
            // For parallel
            vector<vector<double> > Wiblock_more, Sblock_more;
            long subpi = block_more_ptr[bi + 1] - block_more_ptr[bi];
            Wiblock_more.resize(subpi);
            Sblock_more.resize(subpi);
            int iserr = 0;

            // Compute W_i for a range of i
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
            for (long i = block_more_ptr[bi]; i < block_more_ptr[bi + 1]; i++) {
               vector<double> ei(p, 0);
               long nowi = i - block_more_ptr[bi];
               ei[i] = 1;
               Wiblock_more[nowi].resize(p);
               X_sym.ComputeAinvb(ei, Wiblock_more[nowi], g_tol);
            }

            // Compute S_i for a range of i
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
            for (long i = block_more_ptr[bi]; i < block_more_ptr[bi + 1]; i++) {
               long nowi = i - block_more_ptr[bi];
               Sblock_more[nowi].resize(p);
               for (long j = 0; j <= i; j++)
                  Sblock_more[nowi][j] = computeSij(samples, p, n, i, j);
            }

            //				total_time_gradient += omp_get_wtime()-timebegin;

            //				timebegin = omp_get_wtime();
            // Identify free sets.
            for (long i = block_more_ptr[bi]; i < block_more_ptr[bi + 1]; i++) {
               long row_nnz = 0; // store number of free elements at this row/column
               long nowi = i - block_more_ptr[bi];

               D.row_ptr[i] = D.nnz;
               for (long j = 0, idx = X.row_ptr[i]; j <= i; j++) {
                  double Sij = Sblock_more[nowi][j];
                  double Wij = Wiblock_more[nowi][j];
                  double g = Sij - Wij;
                  if (i == j)
                     Wdiag[i] = Wij;

                  // Current Xij value (need to do this because X is
                  // not stored in a dense format)
                  double Xij = 0.0;
                  if (idx < X.row_ptr[i + 1])
                     if (X.col_idx[idx] == j) {
                        Xij = X.values[idx];
                        idx++;
                     }

                  if (Xij != 0.0 || (fabs(g) > lambda)) {
                     buffercolidx[row_nnz] = j;
                     row_nnz++;
                     numActive++;
                     Slist.push_back(Sij);
                     Xlist.push_back(Xij);
                     Wlist.push_back(Wij);
                     if (Xij > 0)
                        g += lambda;
                     else if (Xij < 0)
                        g -= lambda;
                     else
                        g = fabs(g) - lambda;
                     subgrad += fabs(g);
                  }
               }
               D.col_idx.resize(row_nnz + D.nnz);
               for (long j = 0; j < row_nnz; j++)
                  D.col_idx[D.nnz + j] = buffercolidx[j];
               D.nnz += row_nnz;
            }
            //				total_scan_freeset += (omp_get_wtime()-timebegin);
         }

         D.row_ptr[p] = D.nnz;
         D.values.resize(D.nnz, 0);
         // Finished constructing free set

         if (msg >= 1) {
            /*
                printf("Newton iteration %ld.\n", NewtonIter);
                printf("  Active set size = %ld.\n", numActive);
                printf("  sub-gradient = %e, l1-norm of X = %e.\n",
                      subgrad, l1normX);
             */
            Rcout << "Newton iteration " << NewtonIter << "." << endl;
            Rcout << "  Active set size = " << numActive << "." << endl;
            Rcout << "  sub-gradient = " << subgrad << ", l1-norm of X = " << l1normX << "." << endl;
         }

         // For random permutation
         vector<long> activeset_i(D.nnz);
         vector<long> activeset_j(D.nnz);
         vector<long> activeset_idx(D.nnz);


         vector<long> block_ind(p);
         vector<long> i_to_reali(p);
         vector<long> reali_to_i(p);
         for (long i = 0; i < p; i++)
            i_to_reali[i] = i;

         smat_t graph(D);
         for (long idx = 0; idx < graph.nnz; idx++) {
            double Sij = Slist[idx];
            double Xij = Xlist[idx];
            double Wij = Wlist[idx];
            double g = fabs(Sij - Wij);
            if (Xij == 0) {
               g -= lambda;
               if (g < 0)
                  g = 0;
            }
            graph.values[idx] = g;
         }

         smat_t symgraph;
         symgraph.symmetricfrom(graph);
         if (nblock > 1)
            symgraph.clustering(block_ind, nblock);
         else {
            for (long i = 0; i < p; i++)
               block_ind[i] = 0;
         }

         //			printf("within ratio: %lf\n", symgraph.comp_within_ratio(block_ind));

         // sorting according to block_ind (clustering result)
         vector<pair<long, long> > vp;
         for (long i = 0; i < p; i++)
            vp.push_back(make_pair(block_ind[i], i));
         sort(vp.begin(), vp.end());

         // Create the mapping between sorted index and original indexe
         for (long i = 0; i < p; i++)
            i_to_reali[i] = vp[i].second;
         for (long i = 0; i < p; i++)
            reali_to_i[vp[i].second] = i;

         vector<long> block_size(nblock, 0);
         for (long i = 0; i < p; i++)
            block_size[block_ind[i]]++;
         block_ptr[0] = 0;
         for (long i = 1; i <= nblock; i++)
            block_ptr[i] = block_ptr[i - 1] + block_size[i - 1];

         // coordinate descent
         for (unsigned long cdSweep = 1; cdSweep <= (1 + NewtonIter / 3);
                 cdSweep++) {
            diffD = 0.0;

            // Stopping condition for Hessian W\otimes W
            //				double nowtol = 1e-15;
            double nowtol = subgrad * subgrad * 0.01;
            if (nowtol > 1e-4)
               nowtol = 1e-4;
            //				double nowtol = 1e-5;

            // Generate random permutation of blocks
            vector<int> block_per(nblock);
            for (int i = 0; i < nblock; i++)
               block_per[i] = i;
            for (int i = 0; i < nblock; i++) {
               //std::uniform_real_distribution<double> dist(0.0, nblock - i);
               //Temp Disabled for debugging, SWITCH BACK TO MERSENNE TWISTER FOR FINAL
               //int j = i + dist(engine);
               //int j = i + rand() % (nblock - i);
               int j = i + R::runif(0,nblock - i - .00000000001);
               int tmp = block_per[i];
               block_per[i] = block_per[j];
               block_per[j] = tmp;
            }

            for (long bbi = 0; bbi < nblock; bbi++) {
               int bi = block_per[bbi];
               vector<vector<double> > Wiblock;
               long subpi = block_ptr[bi + 1] - block_ptr[bi];
               Wiblock.resize(subpi);
               int alliter = 0;

               // Hessian computation : w_i for i in bi
               //#pragma omp parallel for
               for (long i = block_ptr[bi]; i < block_ptr[bi + 1]; i++) {
                  vector<double> ei(p, 0);
                  long reali = i_to_reali[i];
                  ei[reali] = 1;
                  Wiblock[i - block_ptr[bi]].resize(p);

                  int nowiter = X_sym.ComputeAinvb(ei, Wiblock[i - block_ptr[bi]], nowtol);
                  alliter += nowiter;
               }
               //					total_time_Wcomp += omp_get_wtime()-timebegin;
               // Code for block permutation
               /*
                  vector<int> block_innerper (bi+1);
                  for ( int ii=0 ; ii<bi+1 ; ii++ )
                  block_innerper[ii] = ii;
                  for ( int ii=0 ; ii<bi+1-1 ; ii++ )
                  {
                  int jj = ii+rand()%(bi+1-ii);
                  int tmp = block_innerper[ii];
                  block_innerper[ii] = block_innerper[jj];
                  block_innerper[jj] = tmp;
                  }
                */

               //					for ( long bbj = 0; bbj<=bi ; bbj++)
               for (long bj = 0; bj <= bi; bj++) {
                  long ibegin = block_ptr[bi];
                  long jbegin = block_ptr[bj];

                  //						long bj = block_innerper[bbj];
                  if (msg >= 2) {
                     //printf("handling block i: %d - %d, j: %d-%d\n", block_ptr[bi], block_ptr[bi+1], block_ptr[bj], block_ptr[bj+1]);
                     //fflush(stdout);
                  }
                  // Compute W_j for j in bj
                  long subpj = block_ptr[bj + 1] - block_ptr[bj];
                  vector<vector<double> > Wjblock;
                  Wjblock.resize(subpj);

                  //						timebegin = omp_get_wtime();
                  // flag == 1: boundary nodes
                  // flag == 0: not boundary nodes, don't need to compute W_j
                  vector<int> flag(subpj, 0);
                  if (bi == bj) {
                     for (long i = 0; i < subpj; i++)
                        flag[i] = 1;
                  } else {
                     for (long i = 0; i < subpj; i++)
                        flag[i] = 0;
                  }
                  long subnnz = 0;

                  // Form active set in the current block
                  for (long reali = 0; reali < p; reali++) {
                     if (block_ind[reali] == bi) {
                        for (long idx = D.row_ptr[reali]; idx < D.row_ptr[reali + 1]; idx++) {
                           long realj = D.col_idx[idx];
                           if (block_ind[realj] == bj) {
                              activeset_i[subnnz] = reali;
                              activeset_j[subnnz] = realj;
                              activeset_idx[subnnz] = idx;
                              subnnz++;
                              flag[reali_to_i[realj] - block_ptr[bj]] = 1;
                           }
                        }
                     }

                     // Need this because D is upper triangular
                     if (bi != bj) {
                        if (block_ind[reali] == bj) {
                           for (long idx = D.row_ptr[reali]; idx < D.row_ptr[reali + 1]; idx++) {
                              long realj = D.col_idx[idx];
                              if (block_ind[realj] == bi) {
                                 activeset_i[subnnz] = realj;
                                 activeset_j[subnnz] = reali;
                                 activeset_idx[subnnz] = idx;
                                 subnnz++;
                                 flag[reali_to_i[reali] - block_ptr[bj]] = 1;
                              }
                           }
                        }
                     }
                  }


                  // Form local active set, a better way
                  /*						for ( long i=block_ptr[bi] ; i<block_ptr[bi+1] ; i++ )
                                    {
                                    long reali = i_to_reali[i];
                                    for ( long idx=D.row_ptr[reali] ; idx<D.row_ptr[reali+1] ; idx++)
                                    {
                                    long realj = D.col_idx[idx];
                                    if ( block_ind[realj] == bj)
                                    {
                                    activeset_i[subnnz] = reali;
                                    activeset_idx[subnnz] = idx;
                                    subnnz ++;
                                    flag[reali_to_i[realj]-block_ptr[bj]] = 1;
                                    }
                                    }
                                    }*/
                  //						total_time_FormActive += omp_get_wtime()-timebegin;

                  double tttbegin; //= omp_get_wtime();
                  // Form W for nodes in the j-th block
                  if (bi != bj) {
                     long num_boundary = 0;
                     //							timebegin = omp_get_wtime();

                     //#pragma omp parallel for
                     for (long j = block_ptr[bj]; j < block_ptr[bj + 1]; j++) {

                        if (flag[j - block_ptr[bj]] == 1) {
                           Wjblock[j - block_ptr[bj]].resize(p);
                           vector<double> ej(p, 0);
                           ej[i_to_reali[j]] = 1;
                           X_sym.ComputeAinvb(ej, Wjblock[j - block_ptr[bj]], nowtol);
                           //									X.ComputeAinvb(ej, Wjblock[j-block_ptr[bj]], nowtol); // Computed using assymmetric matrix
                           num_boundary++;
                        }
                     }

                     if (msg >= 2) {
                        //printf("number of column computation: %d,  subpj= %d\n", num_boundary, subpj);
                        //fflush(stdout);
                     }
                     //							total_time_Wcomp += omp_get_wtime()-timebegin;
                  }
                  //						printf("time for computing off-diagonal Wj: %lf\n", omp_get_wtime()-tttbegin);
                  //						timebegin = omp_get_wtime();
                  // Compute U=DW for a  block
                  vector<vector<double> > U(subpj);
                  //#pragma omp parallel for
                  for (long j = 0; j < subpj; j++)
                     if (flag[j] == 1) {
                        U[j].resize(p);
                        if (bi != bj)
                           D.ComputeAx(Wjblock[j], U[j]);
                        else
                           D.ComputeAx(Wiblock[j], U[j]);
                     }

                  //						total_time_Ucomp += omp_get_wtime()-timebegin;

                  // Compute Pij (see the coordinate update in the BigQUIC paper)
                  vector<double> activeset_pij(subnnz);
                  //#pragma omp parallel for
                  for (long ii = 0; ii < subnnz; ii++) {
                     long reali = activeset_i[ii];
                     long realj = activeset_j[ii];
                     long i = reali_to_i[reali];
                     long j = reali_to_i[realj];
                     double b = 0;
                     for (long k = 0; k < p; k++)
                        if (block_ind[k] != bi && block_ind[k] != bj)
                           b += Wiblock[i - ibegin][k] * U[j - jbegin][k];
                     activeset_pij[ii] = b;
                  }

                  //						timebegin = omp_get_wtime();

                  for (int inneriter = 0; inneriter < 1; inneriter++) {
                     // Random seed for debugging
                     //							srand(cdSweep*100+bi+bj);

                     // Random permutation for the active set (within a block)
                     for (long i = 0; i < subnnz; i++) {
                        //Temp Disabled for debugging, SWITCH BACK TO MERSENNE TWISTER FOR FINAL
                        //std::uniform_real_distribution<double> dist2(0.0, subnnz-i);
                        //long j = i + dist2(engine);
                        //long j = i + rand() % (subnnz - i);
                        long j = i + R::runif(0,subnnz - i);
                        long itmp = activeset_i[i];
                        long jtmp = activeset_j[i];
                        long idxtmp = activeset_idx[i];
                        double pij = activeset_pij[i];
                        activeset_i[i] = activeset_i[j];
                        activeset_idx[i] = activeset_idx[j];
                        activeset_j[i] = activeset_j[j];
                        activeset_i[j] = itmp;
                        activeset_idx[j] = idxtmp;
                        activeset_j[j] = jtmp;
                        activeset_pij[i] = activeset_pij[j];
                        activeset_pij[j] = pij;
                     }

                     // Coordinate descent updates
                     for (long ii = 0; ii < subnnz; ii++) {
                        long idx = activeset_idx[ii];
                        long reali = activeset_i[ii];
                        long realj = activeset_j[ii];
                        long i = reali_to_i[reali];
                        long j = reali_to_i[realj];
                        double Sij = Slist[idx];
                        double Xij = Xlist[idx];
                        double Wij = Wlist[idx];
                        double Dij = D.values[idx];

                        double a = Wij*Wij;
                        if (i != j)
                           a += Wdiag[reali] * Wdiag[realj];

                        double ainv = 1.0 / a; // multiplication is cheaper than division

                        double b = activeset_pij[ii];
                        double bb1 = 0, bb2 = 0;
                        //#pragma omp parallel for reduction(+:bb1)
                        for (long k = block_ptr[bi]; k < block_ptr[bi + 1]; k++) {
                           long realk = i_to_reali[k];
                           bb1 += Wiblock[i - ibegin][realk] * U[j - jbegin][realk];
                        }

                        b += bb1;
                        if (bi != bj) {
                           bb2 = 0;
                           //#pragma omp parallel for reduction(+:bb2)
                           for (long k = block_ptr[bj]; k < block_ptr[bj + 1]; k++) {
                              long realk = i_to_reali[k];
                              bb2 += Wiblock[i - ibegin][realk] * U[j - jbegin][realk];
                           }
                           b += bb2;
                        }

                        b += (Sij - Wij);

                        double l = lambda*ainv;
                        double c = Xij + D.values[idx];
                        double f = b*ainv;
                        double mu;
                        normD -= fabs(D.values[idx]);
                        if (c > f) {
                           mu = -f - l;
                           if (c + mu < 0.0) {
                              mu = -c;
                              D.values[idx] = -Xij;
                           } else {
                              D.values[idx] += mu;
                           }
                        } else {
                           mu = -f + l;
                           if (c + mu > 0.0) {
                              mu = -c;
                              D.values[idx] = -Xij;
                           } else {
                              D.values[idx] += mu;
                           }
                        }
                        diffD += fabs(mu);
                        normD += fabs(D.values[idx]);
                        if (i == j)
                           trgradgD += (Sij - Wij)*(D.values[idx] - Dij);
                        else
                           trgradgD += (Sij - Wij)*(D.values[idx] - Dij)*2;
                        if (mu != 0.0) {
                           //#pragma omp parallel for
                           for (long k = 0; k < subpj; k++)
                              if (flag[k] == 1) {
                                 if (bi != bj)
                                    U[k][reali] += mu * Wjblock[j - jbegin][i_to_reali[k + jbegin]];
                                 else
                                    U[k][reali] += mu * Wiblock[j - jbegin][i_to_reali[k + jbegin]];
                              }

                           if (i != j) {
                              //#pragma omp parallel for
                              for (long k = 0; k < subpj; k++)
                                 if (flag[k] == 1)
                                    U[k][realj] += mu * Wiblock[i - ibegin][i_to_reali[k + jbegin]];
                           }
                        }
                     }
                  }
                  //						total_time_cd += omp_get_wtime()-timebegin;
               }
            }
            if (msg >= 1) {
               //printf("  Coordinate descent sweep %ld. norm of D = %e, "
               //		"change in D = %e.\n", cdSweep, normD, diffD);
               //printf("diffD: %lf,  normD: %lf, cdSweepTol: %lf\n", diffD, normD, cdSweepTol);
            }
            if (diffD <= normD * cdSweepTol)
               break;

            if (diffD > 1e10) {
               error_occur = 1;
               break;
            }

         }
      }

      if (error_occur == 1) {
         //printf("ERROR Occurs!\n");
         double logdetX1tmp = 0;
         int flag_pd = X.ComputeLogdet(logdetX1tmp, 1e-10);
         if (!flag_pd) {
            //printf("ERROR occurs because X lack of positive definite");
         }
         //fflush(stdout);
         break;
      }

      // Line Search!

      double alpha = 1.0;
      double l1normXD = 0.0;
      double fX1prev = 1e+15;
      double timebegin; // = omp_get_wtime();
      for (unsigned long lineiter = 0; lineiter < max_lineiter;
              lineiter++) {
         // Generate X+alpha*D;  compute l1norm, trSX1
         smat_t X_alphaD(X, D, alpha);

         double l1normX1 = X_alphaD.l1norm() * lambda;
         double trSX1 = X_alphaD.ComputetrSX(samples, n);

         double logdetX1 = 0.0;
         int flag_pd = X_alphaD.ComputeLogdet(logdetX1, 1e-5);

         if (!flag_pd) {
            if (msg >= 1) {
               //printf("    Line search step size %e.  Lack of positive "
               //		"definiteness.\n", alpha);
            }
            alpha *= 0.5;
            continue;
         }

         fX1 = (trSX1 + l1normX1) - logdetX1;
         //			total_time_linesearch += omp_get_wtime()-timebegin;

         if (alpha == 1.0)
            l1normXD = l1normX1;
         if (fX1 <= fX + alpha * sigma * (trgradgD + l1normXD - l1normX) ||
                 normD == 0) {
            if (msg >= 1) {
               //printf("    Line search step size chosen: %e.\n", alpha);
            }
            fXprev = fX;
            fX = fX1;
            l1normX = l1normX1;
            logdetX = logdetX1;
            trSX = trSX1;
            X.copyfrom(X_alphaD);
            break;
         }
         if (msg >= 1) {
            //printf("    Line search step size %e.\n", alpha);
            //printf("      Objective value would not decrease sufficiently: "
            //		"%e.\n", fX1 - fX);
         }
         if (fX1prev < fX1) {
            fXprev = fX;
            l1normX = l1normX1;
            logdetX = logdetX1;
            trSX = trSX1;
            X.copyfrom(X_alphaD);
            break;
         }
         fX1prev = fX1;
         alpha *= 0.5;
      }
      //if (msg >=0){
      if (msg >= 1){
         //printf("Iter %ld: obj %lf time %lf\n", NewtonIter, fX1, omp_get_wtime()-timeBegin);
         #ifdef _OPENMP
         Rcout << "Iter " << NewtonIter << ": obj " << fX1 << " time " << omp_get_wtime() - timeBegin << endl;
         #else
               Rcout << "Iter " << NewtonIter << ": obj " << fX1 << " time " << "None, time was measured by OMP which is not available on your system" << endl;
         #endif
      }

      // compute W = inv(X):
      //		ptrdiff_t info;
      //		ptrdiff_t p0 = p;
      //		dpotri_((char*) "U", &p0, W, &p0, &info);
      //
      //		for (unsigned long i = 0; i < p; i++) {
      //			for (unsigned long j = 0; j <= i; j++) {
      //				double tmp = W[i*p+j];
      //				W[j*p+i] = tmp;
      //			}
      //		}
      //		for (unsigned long i = 0, k = 0; i < p; i++, k += p)
      //			for (unsigned long j = 0; j <= i; j++)
      //				X[k+j] += alpha*D[k+j];

      #ifdef _OPENMP
      quic_time = omp_get_wtime() - timeBegin;
      #else
          quic_time = 0;
      #endif
      timelist.push_back(quic_time);
      objlist.push_back(fX);

      // Check for convergence.
      if (subgrad * alpha >= l1normX * tol && (fabs((fX - fXprev) / fX) >= EPS))
         continue;
      //		if (mode =='P') {
      //			if (opt != NULL)
      //				opt[pathIdx] = fX;
      //			if (iter != NULL)
      //				iter[pathIdx] = NewtonIter;
      //			if (dGap != NULL) {
      //				double logdetW = projLogDet(p, S, W, U, Lambda);
      //				double gap = -logdetW - p - logdetX + trSX + l1normX;
      //				dGap[pathIdx] = gap;
      //			}
      //			for (unsigned long i = 0, k = 0; i < p; i++, k += p)
      //				for (unsigned long j = i+1; j < p; j++)
      //					X[k+j] = X[j*p+i];
      //			double elapsedTime = (clock() - timeBegin)/CLOCKS_PER_SEC;
      //			if (cputime != NULL)
      //				cputime[pathIdx] = elapsedTime;
      //			// Next lambda.
      //			pathIdx++;
      //			if (pathIdx == pathLen)
      //				break;
      //			if (msg > QUIC_MSG_NO)
      //				MSG("  New scaling value: %e\n", path[pathIdx]);
      //			unsigned long p2 = p*p;
      //			memcpy(X + p2, X, p2*sizeof(double));
      //			memcpy(W + p2, W, p2*sizeof(double));
      //			X += p2;
      //			W += p2;
      //			for (unsigned long i = 0; i < p*p; i++)
      //				Lambda[i] = Lambda0[i]*path[pathIdx];
      //			l1normX = l1normX/path[pathIdx-1]*path[pathIdx];
      //			continue;
      //		}

      break;
   }

   smat_t X_sym;
   X_sym.symmetricfrom(X);
   X.copyfrom(X_sym);

}

