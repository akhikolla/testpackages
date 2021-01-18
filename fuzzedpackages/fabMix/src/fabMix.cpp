#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
List compute_A_B_G_D_and_simulate_mu_Lambda( arma::mat SigmaINV, Rcpp::List suff_statistics, arma::mat OmegaINV, int K,  arma::vec priorConst1, arma::mat T_INV, arma::vec v_r) {
	int i, j, k; 
	int p = SigmaINV.n_cols;
	arma::mat sy = suff_statistics["sy"];
	int q = sy.n_cols;
	arma::vec cluster_size = suff_statistics["cluster_size"], lambdaMean(q), muMean(p);
	arma::mat A(p,p), B(p,1), G(q,q), D(p,q), mu(K,p), tmp(q,q), sx = suff_statistics["sx"], Y(1,q), M(1,p);
	arma::cube Lambda(p,q,K), Lambda_out(K,p,q), syy(q,q,K), sxy(p,q,K), syyOLD = suff_statistics["syy"], sxyOLD = suff_statistics["sxy"];
	double u;

	for(k = 0; k < K; k++){
		for(i = 0; i < q; i++){
			for(j = 0; j < q; j++){
				syy(i,j,k) = syyOLD(k,i,j);
			}
		}
		for(i = 0; i < p; i++){
			for(j = 0; j < q; j++){
				sxy(i,j,k) = sxyOLD(k,i,j);
			}
		}
	}

	Lambda = arma::zeros(p,q,K);
	for(k = 0; k < K; k++ ){
		A = arma::zeros(p,p);
		for(i = 0; i < p; i++){
			A(i, i) = 1.0/( cluster_size(k)*SigmaINV(i, i) + T_INV(i, i) );
			for(j = 0; j < v_r(i); j++){
				D(i, j) = sxy(i, j, k)*SigmaINV(i,i) - A(i,i)*SigmaINV(i,i)*SigmaINV(i,i)*sx(k,i)*sy(k,j);
			}
			u = A(i,i)*SigmaINV(i,i)*SigmaINV(i,i);
			tmp = OmegaINV( arma::span(0, v_r(i) - 1), arma::span(0, v_r(i) - 1) ) + syy.slice(k)(arma::span(0, v_r(i) - 1), arma::span(0, v_r(i) - 1))*SigmaINV(i,i)- u * sy( arma::span(k,k), arma::span(0,v_r(i)-1) ).t() * sy( arma::span(k,k), arma::span(0,v_r(i)-1) ); 
			G = inv(tmp( arma::span(0, v_r(i) - 1), arma::span(0, v_r(i) - 1) ));
			lambdaMean = G( arma::span(0, v_r(i) - 1), arma::span(0, v_r(i) - 1) ) * D( arma::span(i, i), arma::span(0, v_r(i) - 1) ).t(); 
			Y = arma::randn(1, v_r(i));
			Lambda.slice(k)( arma::span(i, i), arma::span(0, v_r(i) - 1) ) = repmat(lambdaMean,1,1).t() + Y * chol( G( arma::span(0, v_r(i) - 1), arma::span(0, v_r(i) - 1) ) );
			for(j = 0; j < q; j++){
				Lambda_out(k,i,j) = Lambda(i,j,k);
			}
		}
		B = SigmaINV * ( sx( arma::span(k, k), arma::span(0, p - 1) ).t() - Lambda.slice(k) * sy( arma::span(k, k), arma::span(0, q - 1) ).t() ) + priorConst1;
		muMean = A * B;
		M = arma::randn(1,p);
		mu(arma::span(k,k), arma::span(0,p-1) ) = repmat(muMean,1,1).t() + M*chol( A );
	}

	List result;
	result["Lambdas"] = Lambda_out;
	result["mu"] = mu;
	return(result);
}

// [[Rcpp::export]]
List compute_A_B_G_D_and_simulate_mu_Lambda_Sj( arma::cube SigmaINV, Rcpp::List suff_statistics, arma::mat OmegaINV, int K,  arma::vec priorConst1, arma::mat T_INV, arma::vec v_r) {
	int i, j, k; 
	int p = SigmaINV.n_cols;
	arma::mat sy = suff_statistics["sy"];
	int q = sy.n_cols;
	arma::vec cluster_size = suff_statistics["cluster_size"], lambdaMean(q), muMean(p);
	arma::mat A(p,p), B(p,1), G(q,q), D(p,q), mu(K,p), tmp(q,q), sx = suff_statistics["sx"], Y(1,q), M(1,p);
	arma::cube Lambda(p,q,K), Lambda_out(K,p,q), syy(q,q,K), sxy(p,q,K), syyOLD = suff_statistics["syy"], sxyOLD = suff_statistics["sxy"], SigmaINV_reshaped(p,p,K);
	double u;

	for(k = 0; k < K; k++){
		for(i = 0; i < q; i++){
			for(j = 0; j < q; j++){
				syy(i,j,k) = syyOLD(k,i,j);
			}
		}
		for(i = 0; i < p; i++){
			for(j = 0; j < q; j++){
				sxy(i,j,k) = sxyOLD(k,i,j);
			}
			for(j = 0; j < p; j++){
				SigmaINV_reshaped(i,j,k) = SigmaINV(k,i,j);
			}
		}
	}

	Lambda = arma::zeros(p,q,K);
	for(k = 0; k < K; k++ ){
		A = arma::zeros(p,p);
		for(i = 0; i < p; i++){
			A(i, i) = 1.0/( cluster_size(k)*SigmaINV_reshaped(i, i, k) + T_INV(i, i) );
			for(j = 0; j < v_r(i); j++){
				D(i, j) = sxy(i, j, k)*SigmaINV_reshaped(i,i,k) - A(i,i)*SigmaINV_reshaped(i,i,k)*SigmaINV_reshaped(i,i,k)*sx(k,i)*sy(k,j);
			}
			u = A(i,i)*SigmaINV_reshaped(i,i,k)*SigmaINV_reshaped(i,i,k);
			tmp = OmegaINV( arma::span(0, v_r(i) - 1), arma::span(0, v_r(i) - 1) ) + syy.slice(k)(arma::span(0, v_r(i) - 1), arma::span(0, v_r(i) - 1))*SigmaINV_reshaped(i,i,k)- u * sy( arma::span(k,k), arma::span(0,v_r(i)-1) ).t() * sy( arma::span(k,k), arma::span(0,v_r(i)-1) ); 
			G = inv(tmp( arma::span(0, v_r(i) - 1), arma::span(0, v_r(i) - 1) ));
			lambdaMean = G( arma::span(0, v_r(i) - 1), arma::span(0, v_r(i) - 1) ) * D( arma::span(i, i), arma::span(0, v_r(i) - 1) ).t(); 
			Y = arma::randn(1, v_r(i));
			Lambda.slice(k)( arma::span(i, i), arma::span(0, v_r(i) - 1) ) = repmat(lambdaMean,1,1).t() + Y * chol( G( arma::span(0, v_r(i) - 1), arma::span(0, v_r(i) - 1) ) );
			for(j = 0; j < q; j++){
				Lambda_out(k,i,j) = Lambda(i,j,k);
			}
		}
		B = SigmaINV_reshaped.slice(k) * ( sx( arma::span(k, k), arma::span(0, p - 1) ).t() - Lambda.slice(k) * sy( arma::span(k, k), arma::span(0, q - 1) ).t() ) + priorConst1;
		muMean = A * B;
		M = arma::randn(1,p);
		mu(arma::span(k,k), arma::span(0,p-1) ) = repmat(muMean,1,1).t() + M*chol( A );
	}

	List result;
	result["Lambdas"] = Lambda_out;
	result["mu"] = mu;
	return(result);
}


// [[Rcpp::export]]
List compute_A_B_G_D_and_simulate_mu_Lambda_CCU( arma::mat SigmaINV, Rcpp::List suff_statistics, arma::mat OmegaINV, int K,  arma::vec priorConst1, arma::mat T_INV, arma::vec v_r) {
	int i, j, k; 
	int p = SigmaINV.n_cols;
	arma::mat sy = suff_statistics["sy"];
	int q = sy.n_cols;
	arma::vec cluster_size = suff_statistics["cluster_size"], lambdaMean(q), muMean(p);
	arma::mat A(p,p), B(p,1), G(q,q), D(p,q), mu(K,p), tmp(q,q), sx = suff_statistics["sx"], Y(1,q), M(1,p);
	arma::cube Lambda(p,q,K), Lambda_out(K,p,q), syy(q,q,K), sxy(p,q,K), syyOLD = suff_statistics["syy"], sxyOLD = suff_statistics["sxy"];
	double u;

	for(k = 0; k < K; k++){
		for(i = 0; i < q; i++){
			for(j = 0; j < q; j++){
				syy(i,j,k) = syyOLD(k,i,j);
			}
		}
		for(i = 0; i < p; i++){
			for(j = 0; j < q; j++){
				sxy(i,j,k) = sxyOLD(k,i,j);
			}
		}
	}
//	Lambda
	Lambda = arma::zeros(p,q,K);
	D = arma::zeros(p,q);
	for(i = 0; i < p; i++){
		G = arma::zeros(q,q);
		for(k = 0; k < K; k++ ){
			G( arma::span(0, v_r(i) - 1), arma::span(0, v_r(i) - 1) ) = G( arma::span(0, v_r(i) - 1), arma::span(0, v_r(i) - 1) ) + syy.slice(k)(arma::span(0, v_r(i) - 1), arma::span(0, v_r(i) - 1))*SigmaINV(i,i);
			D.row(i) = D.row(i) + sxy.slice(k).row(i) * SigmaINV(i,i);
		}
		tmp = OmegaINV( arma::span(0, v_r(i) - 1), arma::span(0, v_r(i) - 1) ) + G( arma::span(0, v_r(i) - 1), arma::span(0, v_r(i) - 1) ); 
		G = inv(tmp( arma::span(0, v_r(i) - 1), arma::span(0, v_r(i) - 1) ));
		lambdaMean = G( arma::span(0, v_r(i) - 1), arma::span(0, v_r(i) - 1) ) * D( arma::span(i, i), arma::span(0, v_r(i) - 1) ).t(); 
		Y = arma::randn(1, v_r(i));
		k = 0;
		Lambda.slice(k)( arma::span(i, i), arma::span(0, v_r(i) - 1) ) = repmat(lambdaMean,1,1).t() + Y * chol( G( arma::span(0, v_r(i) - 1), arma::span(0, v_r(i) - 1) ) );
		for(j = 0; j < q; j++){
			Lambda_out(k,i,j) = Lambda(i,j,k);
		}
	}
//	mu
	for(k = 0; k < K; k++ ){
		A = arma::zeros(p,p);
		for(i = 0; i < p; i++){
			A(i, i) = 1.0/( cluster_size(k)*SigmaINV(i, i) + T_INV(i, i) );
		}
		B = SigmaINV * ( sx( arma::span(k, k), arma::span(0, p - 1) ).t() - Lambda.slice(0) * sy( arma::span(k, k), arma::span(0, q - 1) ).t() ) + priorConst1;
		muMean = A * B;
		M = arma::randn(1,p);
		mu(arma::span(k,k), arma::span(0,p-1) ) = repmat(muMean,1,1).t() + M*chol( A );
		if(k > 0){
			for(i = 0; i < p; i++){
				for(j = 0; j < q; j++){
					Lambda_out(k,i,j) = Lambda(i,j,0);
				}
			}
		}
	}

	List result;
	result["Lambdas"] = Lambda_out;
	result["mu"] = mu;
	return(result);
}



// [[Rcpp::export]]
List compute_A_B_G_D_and_simulate_mu_Lambda_CUU( arma::cube SigmaINV, Rcpp::List suff_statistics, arma::mat OmegaINV, int K,  arma::vec priorConst1, arma::mat T_INV, arma::vec v_r) {
	int i, j, k; 
	int p = SigmaINV.n_cols;
	arma::mat sy = suff_statistics["sy"];
	int q = sy.n_cols;
	arma::vec cluster_size = suff_statistics["cluster_size"], lambdaMean(q), muMean(p);
	arma::mat A(p,p), B(p,1), G(q,q), D(p,q), mu(K,p), tmp(q,q), sx = suff_statistics["sx"], Y(1,q), M(1,p);
	arma::cube Lambda(p,q,K), Lambda_out(K,p,q), syy(q,q,K), sxy(p,q,K), syyOLD = suff_statistics["syy"], sxyOLD = suff_statistics["sxy"], SigmaINV_reshaped(p,p,K);
	double u;

	for(k = 0; k < K; k++){
		for(i = 0; i < q; i++){
			for(j = 0; j < q; j++){
				syy(i,j,k) = syyOLD(k,i,j);
			}
		}
		for(i = 0; i < p; i++){
			for(j = 0; j < q; j++){
				sxy(i,j,k) = sxyOLD(k,i,j);
			}
			for(j = 0; j < p; j++){
				SigmaINV_reshaped(i,j,k) = SigmaINV(k,i,j);
			}
		}
	}
//	Lambda
	Lambda = arma::zeros(p,q,K);
	D = arma::zeros(p,q);
	for(i = 0; i < p; i++){
		G = arma::zeros(q,q);
		for(k = 0; k < K; k++ ){
			G( arma::span(0, v_r(i) - 1), arma::span(0, v_r(i) - 1) ) = G( arma::span(0, v_r(i) - 1), arma::span(0, v_r(i) - 1) ) + syy.slice(k)(arma::span(0, v_r(i) - 1), arma::span(0, v_r(i) - 1))*SigmaINV_reshaped(i,i,k);
			D.row(i) = D.row(i) + sxy.slice(k).row(i) * SigmaINV_reshaped(i,i,k);
		}
		tmp = OmegaINV( arma::span(0, v_r(i) - 1), arma::span(0, v_r(i) - 1) ) + G( arma::span(0, v_r(i) - 1), arma::span(0, v_r(i) - 1) ); 
		G = inv(tmp( arma::span(0, v_r(i) - 1), arma::span(0, v_r(i) - 1) ));
		lambdaMean = G( arma::span(0, v_r(i) - 1), arma::span(0, v_r(i) - 1) ) * D( arma::span(i, i), arma::span(0, v_r(i) - 1) ).t(); 
		Y = arma::randn(1, v_r(i));
		k = 0;
		Lambda.slice(k)( arma::span(i, i), arma::span(0, v_r(i) - 1) ) = repmat(lambdaMean,1,1).t() + Y * chol( G( arma::span(0, v_r(i) - 1), arma::span(0, v_r(i) - 1) ) );
		for(j = 0; j < q; j++){
			Lambda_out(k,i,j) = Lambda(i,j,k);
		}
	}
//	mu
	for(k = 0; k < K; k++ ){
		A = arma::zeros(p,p);
		for(i = 0; i < p; i++){
			A(i, i) = 1.0/( cluster_size(k)*SigmaINV_reshaped(i,i,k) + T_INV(i, i) );
		}
		B = SigmaINV_reshaped.slice(k) * ( sx( arma::span(k, k), arma::span(0, p - 1) ).t() - Lambda.slice(0) * sy( arma::span(k, k), arma::span(0, q - 1) ).t() ) + priorConst1;
		muMean = A * B;
		M = arma::randn(1,p);
		mu(arma::span(k,k), arma::span(0,p-1) ) = repmat(muMean,1,1).t() + M*chol( A );
		if(k > 0){
			for(i = 0; i < p; i++){
				for(j = 0; j < q; j++){
					Lambda_out(k,i,j) = Lambda(i,j,0);
				}
			}
		}
	}

	List result;
	result["Lambdas"] = Lambda_out;
	result["mu"] = mu;
	return(result);
}



