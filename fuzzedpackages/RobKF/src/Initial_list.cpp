#include "Particle.h"
#include <RcppEigen.h>

std::list < struct Particle > Initial_list(const Eigen::MatrixXd & mu_0, const Eigen::MatrixXd & Sigma_0, int Number_of_Particles)
{

	std::list < struct Particle > Output;

	struct Particle First_Particle;

	First_Particle.mu    = mu_0;
	First_Particle.Sigma = Sigma_0;

	First_Particle.log_prob = 0;

	First_Particle.horizon = -1;
	First_Particle.ancestor = NULL;

	First_Particle.anomaly_type = 0;
	First_Particle.anomaly_comp = -1;
	
	First_Particle.anomaly_strength = -1.0;   

	for (int ii = 0; ii < Number_of_Particles; ii++)
	{

		First_Particle.id = -1;
		First_Particle.position = ii;
		Output.push_back(First_Particle);

	}
	
	return Output;

};
