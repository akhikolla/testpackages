#include <RcppEigen.h>
#include <math.h>  

// [[Rcpp::export]]
Eigen::MatrixXd Sigma_Limit(Eigen::MatrixXd Sigma0, Eigen::MatrixXd C, Eigen::MatrixXd A, Eigen::MatrixXd Sigma_Inn, Eigen::MatrixXd Sigma_Add, double epsilon) 
{

	double squarediff;

	Eigen::MatrixXd Pred_Var, Pred_Obs_Var, K, current;

	do 
	{

		current = Sigma0;		

		Pred_Var     = A * current * A.transpose() + Sigma_Inn;
		Pred_Obs_Var = C * Pred_Var * C.transpose() + Sigma_Add;
		K            = Pred_Var * C.transpose() * Pred_Obs_Var.inverse();

		Sigma0 = Pred_Var - Pred_Var * C.transpose() *  Pred_Obs_Var.inverse() * C * Pred_Var; 

		squarediff = (Sigma0-current).squaredNorm();



	} while (sqrt(squarediff) > epsilon);


	                     
    	return Sigma0;

}


