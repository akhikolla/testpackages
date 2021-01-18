/*
This file defines the class of the limited-memory Riemannian trust-region symmetric rank one update method in [HAG2014]
[HAG2014]: W. Huang, P.-A. Absil, and K. A. Gallivan. A Riemannian symmetric rank-one trustregion method. 
		Mathematical Programming, 150(2):179?16, February 2015.

Solvers --> SolversTR --> LRTRSR1

---- WH
*/

#ifndef LRTRSR1_H
#define LRTRSR1_H

#include "ForDebug.h"
#include "SolversTR.h"
#include "def.h"

/*Define the namespace*/
namespace ROPTLIB{

	class LRTRSR1 : public SolversTR{
	public:
		/*The contructor of LRTRSR1 method. It calls the function Solvers::Initialization.
		INPUT : prob is the problem which defines the cost function, gradient and possible the action of Hessian
		and specifies the manifold of domain.
		initialx is the initial iterate.*/
		LRTRSR1(const Problem *prob, const Variable *initialx);

		/*Destructor. Delete the arrays and vectors used in LRTRSR1,
		i.e., vector series S and Y and YMGS, double series RHO, matrices PMGQ, SS, and YY, and permutation indices P (see below for details of the series)*/
		virtual ~LRTRSR1();

		/*Check whether the parameters about LRBFGS are legal or not.*/
		virtual void CheckParams();

		/*Run the algorithm. New memory for vector series S and Y and YMGS, double series RHO,
			matrices PMGQ, SS, and YY, and permutation indices P (see below for details of the series). Then call SolversTR::Run*/
		virtual void Run();

		/*PARAMSMAP is defined in "def.h" and it is a map from string to double, i.e., std::map<std::string, double> .
		This function is used to set the parameters by the mapping*/
		virtual void SetParams(PARAMSMAP params);

		/*Call Solvers::SetProbX function and set up the temporary objects for LRTRSR1 algorithm.
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

		/*the number of pairs of s and y;  The same as \ell in [HAG2014, Algorithm 2]
				[HAG2014]: W. Huang, P.-A. Absil, and K. A. Gallivan. A Riemannian symmetric rank-one trustregion method.
				Mathematical Programming, 150(2):179?16, February 2015.
				Default: 4*/
		integer LengthSY;

	protected:
		/*Compute result = H[Eta], where H is the Hessian or the Hessian approximation
			See details in [HAG2014, Section 4]
			[HAG2014]: W. Huang, P.-A. Absil, and K. A. Gallivan. A Riemannian symmetric rank-one trustregion method.
			Mathematical Programming, 150(2):179?16, February 2015.*/
		virtual void HessianEta(Vector *Eta, Vector *result);

		/*update the pairs of s and y. Add the latest one and remove the oldest one if necessary.*/
		virtual void UpdateData();

		/*This function is called when the candidate is accepted. It transports tangent vectors paris of s and y
		from the tangent space at x1 to the tangent space at x2.*/
		virtual void Acceptence();

		/*Print information specific to LRTRSR1*/
		virtual void PrintInfo();

		bool isupdated; /*Mark whether the inverse Hessian approximation is updated*/
		bool ischangedSandY; /*Mark whether S and Y is updated.*/
		double inpsy, inpss, inpyy;/*inpsy: g(s, y); inpss: g(s, s); inpyy: g(y, y);*/
		Vector *s, *y; /*the s and y of current step*/
		Vector **S, **Y, **YMGS; /*The stored pairs of s and y, and also Y - gamma S*/
		double *SS, *SY, *PMGQ; /*SS is S^\flat S which is the matrix Q defined in [HAG2014, (46)],
									SY is the matrix P defined in [HAG2014, (46)]
									PMGQ is the P - gamma Q matrix defined in [HAG2014, (46)] */
		integer *P;	/*the pemuation matrix when computing the LU decomposition for the matrix PMGQ*/

		double gamma;/*gamma: g(y, y) / g(s, y)*/
		integer Currentlength;	/*The current length of array S, Y and RHO*/
		integer beginidx;	/*The starting index of S, Y and RHO at current iteration*/
	private:
		/*new memory for the double array Vs, type Vector, with length l*/
		void NewVectors(Vector ** &Vs, integer l);

		/*delete memory for the double array Vs, type Vector, with length l*/
		void DeleteVectors(Vector ** &Vs, integer l);
	};
} /*end of ROPTLIB namespace*/
#endif // end of LRTRSR1_H
