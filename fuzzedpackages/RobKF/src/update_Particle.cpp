#include <RcppEigen.h>
#include "Particle.h"
#include <iostream>

void update_Particle(struct Particle & Sampled_Particle, const Eigen::MatrixXd & A, const Eigen::MatrixXd & C, const Eigen::MatrixXd & Sigma_Inn, const Eigen::MatrixXd & Sigma_Add, std::list<Eigen::MatrixXd> Y)
{

	int ii;

	Eigen::MatrixXd mu_pred, mu_new, Sigma_pred, Sigma_new, Inn_Variance, Inn_Precision, Kalman_Gain, Innovation, V, W; 

	V = Eigen::MatrixXd::Identity(C.rows(),C.rows());
	W = Eigen::MatrixXd::Identity(C.cols(),C.cols());

	if (Sampled_Particle.anomaly_type == 1)
	{

		ii         = Sampled_Particle.anomaly_comp;
		W(ii,ii)  += 1/(Sampled_Particle.anomaly_strength);

	}

	if (Sampled_Particle.anomaly_type == 2)
	{

		ii         = Sampled_Particle.anomaly_comp;
		V(ii,ii)  += 1/(Sampled_Particle.anomaly_strength);

	}


	std::list<Eigen::MatrixXd >::iterator it_Y = Y.begin();

	advance(it_Y,Sampled_Particle.horizon-1);

	mu_pred    =  A * Sampled_Particle.ancestor->mu;
	Sigma_pred =  A * Sampled_Particle.ancestor->Sigma * A.transpose() + Sigma_Inn*W;

	Inn_Variance  = C*Sigma_pred*C.transpose() + Sigma_Add*V;
	Inn_Precision = Inn_Variance.inverse();
	Kalman_Gain   = Sigma_pred*C.transpose()*Inn_Precision;

	Innovation = *it_Y - C*mu_pred;

	mu_new    = mu_pred + Kalman_Gain*Innovation;
	Sigma_new = Sigma_pred - Kalman_Gain *C*Sigma_pred;

	for (ii = 0; ii < Sampled_Particle.horizon-1; ii ++)
	{
		
		it_Y--;
		mu_pred    =  A * mu_new;
		Sigma_pred =  A * Sigma_new * A.transpose() + Sigma_Inn;

		Inn_Variance  = C*Sigma_pred*C.transpose() + Sigma_Add;
		Inn_Precision = Inn_Variance.inverse();
		Kalman_Gain   = Sigma_pred*C.transpose()*Inn_Precision;

		Innovation = *it_Y - C*mu_pred;

		mu_new    = mu_pred + Kalman_Gain*Innovation;
		Sigma_new = Sigma_pred - Kalman_Gain *C*Sigma_pred;


	}

	Sampled_Particle.mu    = mu_new;
	Sampled_Particle.Sigma = Sigma_new;

};
