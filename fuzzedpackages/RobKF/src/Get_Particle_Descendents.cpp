#include "Particle.h"
#include <RcppEigen.h>
#include <iostream>

std::list < struct Particle > Get_Particle_Descendents(const struct Particle & Ancestor, const Eigen::MatrixXd  & relevant_Y, const std::vector < double > & sigma_tilde, const Eigen::MatrixXd & Sigma_Add,
 int Number, double s, const std::vector <double> & General_Weights_Add)
{

	std::list < struct Particle > Output;
	std::list < struct Particle > Additions;

	Eigen::MatrixXd Residual        = relevant_Y - Ancestor.obs_Pred;
	Eigen::MatrixXd Pre_Denominator = Ancestor.obs_Prec;
	Eigen::MatrixXd Pre_Numerator   = Pre_Denominator * Residual;
 
	const double log_likelihood = - 0.5 * ( (Residual.transpose() * Pre_Numerator ) .value() - log(Ancestor.obs_Prec.determinant()) );

	Additions = Get_Particle_Descendents_typical(Ancestor, log_likelihood);

	Output.splice(Output.end(),Additions);

	Additions = Get_Particle_Additive_Descendents(Ancestor, log_likelihood, Number, s, sigma_tilde, Sigma_Add, Pre_Numerator, Pre_Denominator, General_Weights_Add);

	Output.splice(Output.end(),Additions);

	return Output;


}
