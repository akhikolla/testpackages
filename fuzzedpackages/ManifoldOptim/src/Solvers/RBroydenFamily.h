/*
This file defines the class of the Riemannian Broyden method in [HGA2015]
		[HGA2015]: Wen Huang, K. A. Gallivan, and P.-A. Absil. A Broyden Class of Quasi-Newton Methods for Riemannian Optimization.
		SIAM Journal on Optimization, 25(3):1660?685, 2015

Solvers --> SolversLS --> RBroydenFamily

---- WH
*/

#ifndef RBROYDENFAMILY_H
#define RBROYDENFAMILY_H

#include "SolversLS.h"
#include "def.h"

/*Define the namespace*/
namespace ROPTLIB{

	class RBroydenFamily : public SolversLS{
	public:
		/*The contructor of RBroydenFamily method. It calls the function Solvers::Initialization.
		INPUT : prob is the problem which defines the cost function, gradient and possible the action of Hessian
		and specifies the manifold of domain.
		initialx is the initial iterate.
		initialH is the initial inverse Hessian approximation. If the input is nullptr, then the identity is used as the initial approximation*/
		RBroydenFamily(const Problem *prob, const Variable *initialx, LinearOPE *initialH = nullptr);

		/*Destructor. Delete the vectors and Hessian approximation used in RBroydenFamily, i.e., s and y, H and tildeH*/
		virtual ~RBroydenFamily();

		/*The phi coefficient in [HGA2015, (2.3)]
			[HGA2015]: Wen Huang, K. A. Gallivan, and P.-A. Absil. A Broyden Class of Quasi-Newton Methods for Riemannian Optimization.
			SIAM Journal on Optimization, 25(3):1660?685, 2015
			*/
		double Phi(Variable *x2, Vector *y, Vector *s, LinearOPE *tildeH, double inpsy, double yHy, Vector *u);

		/*Check whether the parameters about RBFGS are legal or not.*/
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
		bool isupdated; /*Mark whether the inverse Hessian approximation is updated*/
		double betay, phic, inpsy, inpss;	/*betay: \|\xi\| / \|\mathcal{T}_{R_\xi} \xi\| in the locking condition;
											phic: the coefficient (1-phic) BFGS + phi DFP
											inpsy: g(s, y); inpss: g(s, s);
											*/

		/*Compute the search direction. eta1 = H (-gf1) */
		virtual void GetSearchDir();

		/*Update the Hessian approximation if necessary*/
		virtual void UpdateData();

		/*Print information specific to RBroydenFamily*/
		virtual void PrintInfo();

		Vector *s, *y, *u;	/*the s, y and u of current step*/
		LinearOPE *H, *tildeH; /*The inverse Hessian approximations for current and next iterations respectively*/
	};
} /*end of ROPTLIB namespace*/
#endif // end of RBROYDENFAMILY_H
