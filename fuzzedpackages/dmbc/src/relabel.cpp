// relabel.cpp

#include "dmbc.h"
// [[Rcpp::depends(RcppProgress)]]
#include "progress.hpp"
#include "progbar.h"

// Relabeling algorithm as in Celeux et al. (2000)
void relabel_celeux(
	double* theta,
	double* z_chain,
	double* alpha_chain,
	double* eta_chain,
	double* sigma2_chain,
	double* lambda_chain,
	double* prob_chain,
	int* x_ind_chain,
	int init,
	int n,
	int p,
	int S,
	int M,
	int R,
	int G,
  int verbose){
	// estimate the G centers and variances using the first 'init' iterations
	double* theta_m = new double[R*G];
	double* theta_m_byrow = new double[R*G];
	double* theta_ms = new double[R*G];
	double* theta_v = new double[R*G];
	double* theta_v_byrow = new double[R*G];
	double* theta_m_old = new double[R*G];
	for(int r = 0; r < R; r++){
		for(int g = 0; g < G; g++){
			theta_m[r + R*g] = 0;
			theta_ms[r + R*g] = 0;
			for(int m = 0; m < init; m++){
				theta_m[r + R*g] += theta[m + M*r + M*R*g];
				theta_ms[r + R*g] += pow(theta[m + M*r + M*R*g], 2);
			}
			theta_m[r + R*g] /= init;
			theta_m_byrow[g + G*r] = theta_m[r + R*g];
			theta_ms[r + R*g] /= init;
			theta_v[r + R*g] = theta_ms[r + R*g] - pow(theta_m[r + R*g], 2);
			theta_v_byrow[g + G*r] = theta_v[r + R*g];
		}
	}

  ETAProgressBar pb;
  Progress prg(M, verbose, pb); // create the progress monitor

	// find all permutations of initial centers and variances estimates
	int J = factorial(G);
	int* all_perm = new int[G*J];
	permutations(all_perm, G, J, 1);
	double* theta_m_all = new double[R*G*J];
	double* theta_v_all = new double[R*G*J];
	for(int r = 0; r < R; r++){
		for(int g = 0; g < G; g++){
			for(int j = 0; j < J; j++){
				theta_m_all[r + R*g + R*G*j] = theta_m_byrow[(all_perm[g + G*j] - 1) + G*r];
				theta_v_all[r + R*g + R*G*j] = theta_v_byrow[(all_perm[g + G*j] - 1) + G*r];
			}
		}
	}

	int* j_star = new int[M];
	double* dist_theta = new double[R*G*J];
	double* csums = new double[J];
	bool perm_id = false;
	int* perm_tmp = new int[G];
	int* perm = new int[G];
	double* theta_new = new double[G];
	double* z_new = new double[n*p*G];
	double* alpha_new = new double[G];
	double* eta_new = new double[G];
	double* sigma2_new = new double[G];
	double* lambda_new = new double[G];
	double* prob_new = new double[S*G];
	double* x_ind_new = new double[S*G];
	for(int m = init; m < M; m++){
    if(!Progress::check_abort()){
      prg.increment(); // update progress

  		// calculate squared distances of each MCMC sample (from the (init + 1)-th one on) from the G group centers
  		for(int j = 0; j < J; j++){
  			for(int g = 0; g < G; g++){
  				for(int r = 0; r < R; r++){
  					dist_theta[r + R*g + R*G*j] = pow(theta[m + M*r + M*R*g] - theta_m_all[r + R*g + R*G*j], 2)/theta_v_all[r + R*g + R*G*j];
  				}
  			}
  		}
  		
  		for(int j = 0; j < J; j++){
  			csums[j] = 0;
  		}
  		colsums(csums, dist_theta, R*G, J);

  		// allocate the m-th MCMC sample to the closest group
  		which_min(&j_star[m], csums, J);
  		for(int g = 0; g < G; g++){
  			perm[g] = g + 1;
  			perm_tmp[g] = all_perm[g + (j_star[m] - 1)*G];
  		}
  		R_qsort_int_I(perm_tmp, perm, 1, G);
  		perm_id = true;
  		for(int g = 0; g < G; g++){
  			if(perm[g] != all_perm[g]){
  				perm_id = false;
  				break;
  			}
  		}
  		if(!perm_id){
  			for(int r = 0; r < R; r++){
  				for(int g = 0; g < G; g++){
  					theta_new[g] = theta[m + M*r + M*R*(perm[g] - 1)];
  				}
  				for(int g = 0; g < G; g++){
  					theta[m + M*r + M*R*g] = theta_new[g];
  				}
  			}
  		}

  		// update centers and variances
  		for(int r = 0; r < R; r++){
  			for(int g = 0; g < G; g++){
  				theta_m_old[r + R*g] = theta_m[r + R*g];
  				theta_m[r + R*g] = (m/(double)(m + 1))*theta_m[r + R*g] + (1/(double)(m + 1))*theta[m + M*r + M*R*g];
  				theta_v[r + R*g] = (m/(double)(m + 1))*theta_v[r + R*g] + (m/(double)(m + 1))*pow(theta_m_old[r + R*g] - theta_m[r + R*g], 2) + 
  					(1/(double)(m + 1))*pow(theta[m + M*r + M*R*g] - theta_m[r + R*g], 2);
  				theta_m_byrow[g + G*r] = theta_m[r + R*g];
  				theta_v_byrow[g + G*r] = theta_v[r + R*g];
  			}
  		}
  		for(int r = 0; r < R; r++){
  			for(int g = 0; g < G; g++){
  				for(int j = 0; j < J; j++){
  					theta_m_all[r + R*g + R*G*j] = theta_m_byrow[(all_perm[g + G*j] - 1) + G*r];
  					theta_v_all[r + R*g + R*G*j] = theta_v_byrow[(all_perm[g + G*j] - 1) + G*r];
  				}
  			}
  		}

  		// apply the same permutations to the "auxiliary" objects
  		if(!perm_id){
  			for(int g = 0; g < G; g++){
  				for(int i = 0; i < n; i++){
  					for(int j = 0; j < p; j++){
  						z_new[i + n*j + n*p*g] = z_chain[m + M*i + M*n*j + M*n*p*(perm[g] - 1)];
  					}
  				}
  				alpha_new[g] = alpha_chain[m + M*(perm[g] - 1)];
  				eta_new[g] = eta_chain[m + M*(perm[g] - 1)];
  				sigma2_new[g] = sigma2_chain[m + M*(perm[g] - 1)];
  				lambda_new[g] = lambda_chain[m + M*(perm[g] - 1)];
  				for(int s = 0; s < S; s++){
  					prob_new[s + S*g] = prob_chain[m + M*s + M*S*(perm[g] - 1)];
  					x_ind_new[s + S*g] = x_ind_chain[m + M*s + M*S*(perm[g] - 1)];
  				}
  			}
  			for(int g = 0; g < G; g++){
  				for(int i = 0; i < n; i++){
  					for(int j = 0; j < p; j++){
  						z_chain[m + M*i + M*n*j + M*n*p*g] = z_new[i + n*j + n*p*g];
  					}
  				}
  				alpha_chain[m + M*g] = alpha_new[g];
  				eta_chain[m + M*g] = eta_new[g];
  				sigma2_chain[m + M*g] = sigma2_new[g];
  				lambda_chain[m + M*g] = lambda_new[g];
  				for(int s = 0; s < S; s++){
  					prob_chain[m + M*s + M*S*g] = prob_new[s + S*g];
  					x_ind_chain[m + M*s + M*S*g] = x_ind_new[s + S*g];
  				}
  			}
  		}
      R_CheckUserInterrupt();
    }
	}

	delete[] theta_m;
	delete[] theta_m_byrow;
	delete[] theta_ms;
	delete[] theta_v;
	delete[] theta_v_byrow;
	delete[] theta_m_old;
	delete[] all_perm;
	delete[] theta_m_all;
	delete[] theta_v_all;
	delete[] j_star;
	delete[] dist_theta;
	delete[] csums;
	delete[] perm_tmp;
	delete[] perm;
	delete[] theta_new;
	delete[] z_new;
	delete[] alpha_new;
	delete[] eta_new;
	delete[] sigma2_new;
	delete[] lambda_new;
	delete[] prob_new;
	delete[] x_ind_new;
}

void pack_par(
	double* theta,
	const double* z,
	const double* alpha,
	const double* lambda,
	int n,
	int p,
	int M,
	int G){
	int r = n*(n - 1)/2;
	double* z_g = new double[n*p];
	double* delta = new double[r];
	double* alpha_plus_delta = new double[r];
	double* pi_g = new double[r];

	for(int g = 0; g < G; g++){
		for(int m = 0; m < M; m++){
			for(int j = 0; j < p; j++){
				for(int i = 0; i < n; i++){
					z_g[i + n*j] = z[m + M*i + M*n*j + M*n*p*g];
				}
			}
			dist(delta, z_g, n, p);
			for(int k = 0; k < r; k++){
				alpha_plus_delta[k] = alpha[m + M*g] + delta[k];
			}
			expit(pi_g, alpha_plus_delta, r);
			for(int k = 0; k < r; k++){
				theta[m + M*k + M*(r + 1)*g] = pi_g[k];
			}
			theta[m + M*r + M*(r + 1)*g] = lambda[m + M*g];
		}
	}

	delete[] z_g;
	delete[] delta;
	delete[] alpha_plus_delta;
	delete[] pi_g;
}
