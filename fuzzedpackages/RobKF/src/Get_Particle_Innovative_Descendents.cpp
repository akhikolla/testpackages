#include "Particle.h"
#include <RcppEigen.h>
#include <iostream>

std::list < struct Particle > Get_Particle_Innovative_Descendents(const struct Particle & Ancestor, const Eigen::MatrixXd & relevant_Y, const Eigen::MatrixXd & C_Augmented, const std::vector < double > & sigma_hat, 
const Eigen::MatrixXd & Sigma_Inn, int Number, double s, const std::vector <double> & General_Weights, const std::vector <int> & Sample_From, const int & horizon)
{

	std::list < struct Particle > Output;
	std::list < struct Particle > Additions;

	Eigen::MatrixXd Residual        = relevant_Y - Ancestor.obs_Pred;
	Eigen::MatrixXd Pre_Denominator = C_Augmented.transpose() * Ancestor.obs_Prec * C_Augmented;
	Eigen::MatrixXd Pre_Numerator   = C_Augmented.transpose() * Ancestor.obs_Prec * Residual;

	double Numerator, Denominator, General_Weight, log_likelihood;

	log_likelihood = - 0.5 * ((Residual.transpose() * Ancestor.obs_Prec * Residual).value() - std::log(Ancestor.obs_Prec.determinant()) ) ;

	for (int ii = 0; ii < Sigma_Inn.rows(); ii++)
	{

		if (Sample_From[ii] == 1)
		{

			Numerator   = Pre_Numerator(ii,0);
			Numerator   = Numerator*Numerator;
			Denominator = Pre_Denominator(ii,ii);

			/*std::cout << std::endl;
			std::cout << "Innovation";
			std::cout << std::endl;*/
			

			// Move this higher up the function 
			General_Weight = General_Weights[ii];


			/*std::cout << "General_Weight";
			std::cout << std::endl;
			std::cout << General_Weight;
			std::cout << std::endl;
			std::cout << 0.5*log(sigma_hat[ii]);
			std::cout << std::endl;

			std::cout << log_likelihood;
			std::cout << std::endl;

			std::cout << "Numerator";
			std::cout << std::endl;
			std::cout << Numerator;
			std::cout << std::endl;

			std::cout << "Denominator";
			std::cout << std::endl;
			std::cout << Denominator;
			std::cout << std::endl;
			std::cout << C_Augmented;
			std::cout << std::endl;
			std::cout << Ancestor.obs_Prec;
			std::cout << std::endl;*/

			Additions = Get_Particle_Descendents_W(Ancestor, ii, Number, log_likelihood, sigma_hat[ii], Sigma_Inn(ii,ii), Numerator, Denominator, s, General_Weight, horizon);


			Output.splice(Output.end(),Additions);

		}


	}

	return Output;


}
