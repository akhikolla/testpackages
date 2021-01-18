/*
This file defines the class of the Riemannian BFGS method in [RW0212]
	[RW0212]: W. Ring and B. Wirth. Optimization methods on Riemannian manifolds and their application to shape space. 
	SIAM Journal on Optimization, 22(2):596?27, January 2012.

Solvers --> SolversLS --> RWRBFGS

---- WH
*/
#ifndef RWRBFGS_H
#define RWRBFGS_H

#include "SolversLS.h"
#include "def.h"

/*Define the namespace*/
namespace ROPTLIB{

	class RWRBFGS : public SolversLS{
	public:
		/*The contructor of RBFGS method. It calls the function Solvers::Initialization.
		INPUT : prob is the problem which defines the cost function, gradient and possible the action of Hessian
		and specifies the manifold of domain.
		initialx is the initial iterate.
		initialH is the initial inverse Hessian approximation. If the input is nullptr, then the identity is used as the initial approximation*/
		RWRBFGS(const Problem *prob, const Variable *initialx, LinearOPE *initialH = nullptr);

		/*Destructor. Delete the vectors and Hessian approximation used in RBFGS, i.e., s and y, H and tildeH*/
		virtual ~RWRBFGS();

		/*Check whether the parameters about RWRBFGS are legal or not.*/
		virtual void CheckParams();

		/*PARAMSMAP is defined in "def.h" and it is a map from string to double, i.e., std::map<std::string, double> .
		This function is used to set the parameters by the mapping*/
		virtual void SetParams(PARAMSMAP params);

		/*Initialize the solvers by calling the "SetProbX" and "SetDefultParams" functions.
		INPUT:	prob is the problem which defines the cost function, gradient and possible the action of Hessian
		and specifies the manifold of domain.
		initialx is the initial iterate.
		initialH is the initial inverse Hessian approximation. If the input is nullptr, then the identity is used as the initial approximation*/
		virtual void Initialization(const Problem *prob, const Variable *initialx, LinearOPE *initialH = nullptr);

		/*Initialize the type of iterates x1, x2 and tangent vectors gf1, gf2, s, y, H and tildeH and obtian the problem and manifold information
		INPUT:	prob is the problem which defines the cost function, gradient and possible the action of Hessian
		and specifies the manifold of domain.
		initialx is the initial iterate.
		initialH is the initial inverse Hessian approximation. If the input is nullptr, then the identity is used as the initial approximation*/
		virtual void SetProbX(const Problem *prob, const Variable *initialx, LinearOPE *initialH = nullptr);

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
	protected:
		/*Compute the search direction. eta1 = H (-gf1) */
		virtual void GetSearchDir();

		/*Update the Hessian approximation if necessary*/
		virtual void UpdateData();

		/*Print information specific to RWRBFGS*/
		virtual void PrintInfo();

		bool isupdated; /*Mark whether the inverse Hessian approximation is updated*/
		double inpsy, inpss; /*inpsy: g(s, y); inpss: g(s, s);*/
		Vector *s, *y;/*the s and y of current step*/
		LinearOPE *H, *tildeH; /*The inverse Hessian approximations for current and next iterations respectively*/
	};
} /*end of ROPTLIB namespace*/

#endif // end of RWRBFGS_H
