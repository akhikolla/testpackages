#include "Particle.h"
#include <RcppEigen.h>
#include <R.h>

std::list < struct Particle > Get_Particle_Additive_Descendents(const struct Particle & Ancestor, const double & log_likelihood, const int & Number, const double & s, const std::vector<double> sigma_tilde, const Eigen::MatrixXd & Sigma_Add, 
const Eigen::MatrixXd & Pre_Numerator, const Eigen::MatrixXd & Pre_Denominator, const std::vector <double> General_Weights_Add)
{

	std::list < struct Particle > Output;
	std::list < struct Particle > Additions;

	double Numerator, Denominator, General_Weight;

	for (int ii = 0; ii < Sigma_Add.rows(); ii++)
	{

		Numerator   = Pre_Numerator(ii,0);
		Numerator   = Numerator*Numerator;
		Denominator = Pre_Denominator(ii,ii);

		// Move this higher up the function 
		General_Weight = General_Weights_Add[ii];

		Additions = Get_Particle_Descendents_V(Ancestor, ii, Number, log_likelihood, sigma_tilde[ii], Sigma_Add(ii,ii), Numerator, Denominator, s, General_Weight);

		Output.splice(Output.end(),Additions);

	}

	return Output;

};
