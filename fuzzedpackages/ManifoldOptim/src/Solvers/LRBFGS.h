/*
This file defines the class of the limited-memory Riemannian BFGS method in [HGA2015]
	[HGA2015]: Wen Huang, K. A. Gallivan, and P.-A. Absil. A Broyden Class of Quasi-Newton Methods for Riemannian Optimization. 
	SIAM Journal on Optimization, 25(3):1660?685, 2015

Solvers --> SolversLS --> LRBFGS

---- WH
*/

#ifndef LRBFGS_H
#define LRBFGS_H

#include <cstring>
#include "SolversLS.h"
#include "def.h"

/*Define the namespace*/
namespace ROPTLIB{

	class LRBFGS : public SolversLS{
	public:
		/*The contructor of LRBFGS method. It calls the function Solvers::Initialization.
		INPUT : prob is the problem which defines the cost function, gradient and possible the action of Hessian
		and specifies the manifold of domain.
		initialx is the initial iterate.*/
		LRBFGS(const Problem *prob, const Variable *initialx);

		/*Destructor. Delete the arrays and vectors used in LRBFGS, i.e., series S and Y, and series RHO*/
		virtual ~LRBFGS();

		/*Check whether the parameters about LRBFGS are legal or not.*/
		virtual void CheckParams();

		/*Run the algorithm. New memory for S, Y and RHO. Then call SolversLS::Run*/
		virtual void Run();

		/*PARAMSMAP is defined in "def.h" and it is a map from string to double, i.e., std::map<std::string, double> .
		This function is used to set the parameters by the mapping*/
		virtual void SetParams(PARAMSMAP params);

		/*Call Solvers::SetProbX function and set up the temporary objects for LRBFGS algorithm.
		INPUT:	prob is the problem which defines the cost function, gradient and possible the action of Hessian
		and specifies the manifold of domain.
		initialx is the initial iterate.*/
		virtual void SetProbX(const Problem *prob, const Variable *initialx);

		/*Setting parameters (member variables) to be default values */
		virtual void SetDefaultParams();

		/* ===============public parameters below================= */

		/*specify whether the cost function is convex or not.
		If yes, then the initial Hessian approximation is a scalar times identity, where the scalar is to
		measure the magnitude of eigenvalues, otherwise, it is identity.
		Default: false*/
		bool isconvex;

		/*The same as \epsilon in [LF01, (3.2)]
			[LF01]: D.-H. Li and M. Fukushima. On the global convergence of the BFGS method for nonconvex unconstrained optimization problems.
			SIAM Journal on Optimization, 11(4):1054?064, 2001
			Default: 10^{-4}*/
		double nu;

		/*The same as \alpha in [LF01, (3.2)]
			[LF01]: D.-H. Li and M. Fukushima. On the global convergence of the BFGS method for nonconvex unconstrained optimization problems.
			SIAM Journal on Optimization, 11(4):1054?064, 2001
			Default: 1*/
		double mu;

		/*the number of pairs of s and y;  The same as \ell in [HGA2015, Algorithm 2]
			[HGA2015]: Wen Huang, K. A. Gallivan, and P.-A. Absil. A Broyden Class of Quasi-Newton Methods for Riemannian Optimization.
			SIAM Journal on Optimization, 25(3):1660?685, 2015.
			Default: 4*/
		integer LengthSY;
	protected:

		bool isupdated; /*Mark whether the inverse Hessian approximation is updated*/
		double betay, inpsy, inpss, inpyy;  /*betay: \|\xi\| / \|\mathcal{T}_{R_\xi} \xi\| in the locking condition;
											inpsy: g(s, y); inpss: g(s, s); inpyy: g(y, y);
											*/

		/*Compute the search direction. [HGA2015, Steps 3 to 14 in Algorithm 2]
			[HGA2015]: Wen Huang, K. A. Gallivan, and P.-A. Absil. A Broyden Class of Quasi-Newton Methods for Riemannian Optimization.
			SIAM Journal on Optimization, 25(3):1660?685, 2015.
			*/
		virtual void GetSearchDir();

		/*update the pairs of s and y. Add the latest one and remove the oldest one if necessary.
		transport them to the tangent space of x2*/
		virtual void UpdateData();

		/*Print information specific to LRBFGS*/
		virtual void PrintInfo();

		Vector *s, *y; /*the s and y of current step*/
		Vector **S, **Y; /*The stored pairs of s and y*/
		double *RHO; /*the sequence of 1 / g(s_k, y_k), where Hessian approximation k-th iteration is updated*/
		double rho, gamma; /*rho: 1/g(s, y) at current iteration, gamma: g(s, y) / g(y, y)*/
		integer Currentlength; /*The current length of array S, Y and RHO*/
		integer beginidx; /*The starting index of S, Y and RHO at current iteration*/

	private:
		/*new memory for the double array Vs, type Vector, with length l*/
		void NewVectors(Vector ** &Vs, integer l);

		/*delete memory for the double array Vs, type Vector, with length l*/
		void DeleteVectors(Vector ** &Vs, integer l);
	};
} /*end of ROPTLIB namespace*/
#endif // end of RBROYDENFAMILY_H
