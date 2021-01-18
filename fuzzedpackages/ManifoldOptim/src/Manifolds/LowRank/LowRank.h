/*
This file defines the class for the low-rank manifold R_r^{m times n}, which is represented by St(r, m) times R^{r times r} times St(r, n).
A tangent vector in T_x R_r^{m times n} is representated by 
etax = dot{U} D V^T + U dot{D} V^T + U D dot{V}^T, where dot{U}^T U = 0, dot{V}^T V = 0.
The used Riemannian metric is 
g(etax, xix) = trace(dot{U}_1^T dot{U}_2) + trace(dot{D}_1^T dot{D}_2) + trace(dot{V}_1^T dot{V}_2),
where etax = dot{U}_1 D V^T + U dot{D}_1 V^T + U D dot{V}_1^T and xix = dot{U}_2 D V^T + U dot{D}_2 V^T + U D dot{V}_2^T.

Manifold --> ProductManifold --> LowRank

---- WH
*/

#ifndef LOWRANK_H
#define LOWRANK_H

#include "ProductManifold.h"
#include "Stiefel.h"
#include "Euclidean.h"
#include "LowRankVariable.h"
#include "LowRankVector.h"

/*Define the namespace*/
namespace ROPTLIB{

	class LowRank : public ProductManifold{
	public:
		/*Construct the low rank manifold of m by n matrices with rank r.
		It is represented by St(r, m) times R^{r times r} times St(r, n), i.e.,
		X = U D V^T. U \in St(r, m), D \in R^{r times r} and V \in St(r, n).
		Note that D is not necessary a diagonal matrix.*/
		LowRank(integer m, integer n, integer r);

		/*Delete the manifold by deleting each component.*/
		~LowRank(void);

		/*Tangent vector etax = \dot(U) D V^T + U \dot{D} V^T + U D \dot{V}^T. Let \dot{U} = U Omega_U + U_perp K_U and \dot{V} = V Omega_V + V_perp K_V
		It follows that etax = U_perp K_U D V^T + U D K_V^T V_perp^T + U (Omega_U D + \dot{D} + D Omega_V^T) V^T
		The intrinsic representation would be given by vectorizing (A, B, C) := (K_U D, D K_V^T, Omega_U D + \dot{D} + D Omega_V^T) .*/
		virtual void ObtainIntr(Variable *x, Vector *etax, Vector *result) const;

		/*Compute the extrinsic approach by given (A, B, C) := (K_U D, D K_V^T, Omega_U D + \dot{D} + D Omega_V^T).
		\dot{U} = U_\perp A D^{-1}, \dot{D} = C, \dot{V} = V_\perp B^T D^{-T}.
		Note that ObtainExtr(ObtainIntr(\dot{U}, \dot{D}, \dot{V})) may not recover (\dot{U}, \dot{D}, \dot{V}).
		However, this is fine since \dot(U) D V^T + U \dot{D} V^T + U D \dot{V}^T do not change, which implies they represent a same tangent vector.
		*/
		virtual void ObtainExtr(Variable *x, Vector *intretax, Vector *result) const;

		/*Perform the retraction of each manifold*/
		virtual void Retraction(Variable *x, Vector *etax, Variable *result) const;

		/*Compute the tangent vector result satisfying
		g_y(\mathcal{T}_{R_etax}(xix), xiy) = g_x(xix, result) for all xix \in T_x M,
		where y = R_x(etax), xiy \in T_y M.
		The vector transport by differentiated retraction of the individual manifold is known.
		The closed form can be computed. (TODO: add details in notes?)*/
		virtual void coTangentVector(Variable *x, Vector *etax, Variable *y, Vector *xiy, Vector *result) const;

		/*Used in "coTangentVector". (TODO: add details in notes?)*/
		virtual void ExtrProjectionStiePerp(Variable *x, Vector *v, Vector *result) const;

		/*Perform the vector transport by differentiated retraction of each individual manifold.*/
		virtual void DiffRetraction(Variable *x, Vector *etax, Variable *y, Vector *xix, Vector *result, bool IsEtaXiSameDir) const;

		/*etax is in the ambient space. This function projects etax onto the tangent space of x, i.e., result = P_{T_x M} v*/
		virtual void ExtrProjection(Variable *x, Vector *v, Vector *result) const;

	protected:
		integer m; /*the number of row*/
		integer n; /*the number of column*/
		integer r; /*the rank of the matrix*/
	};
} /*end of ROPTLIB namespace*/
#endif // end of LOWRANK_H
