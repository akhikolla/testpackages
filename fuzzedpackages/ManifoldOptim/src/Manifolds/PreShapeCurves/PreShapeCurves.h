#ifndef PRESHAPECURVES_H
#define PRESHAPECURVES_H

#include "ForDebug.h"
#include "PSCVariable.h"
#include "PSCVector.h"
#include "Manifold.h"
#include "def.h"
#include "ElasticCurvesRO.h"
#include "PreShapePathStraighten.h"
#include "MyMatrix.h"

/*Define the namespace*/
namespace ROPTLIB{

	class PreShapeCurves : public Manifold{
	public:
		PreShapeCurves(integer r, integer c, integer n);
		virtual ~PreShapeCurves();
		virtual void EucHvToHv(Variable *x, Vector *etax, Vector *exix, Vector* result, const Problem *prob) const;
		virtual void CheckParams(void) const;

		virtual double Metric(Variable *x, Vector *etax, Vector *xix) const;
		virtual void Retraction(Variable *x, Vector *etax, Variable *result) const;
		virtual void EucGradToGrad(Variable *x, Vector *egf, Vector *gf, const Problem *prob) const;
		static double InnerProd_Q(const double *q1, const double *q2, integer innumP, integer indim);
		static void CovIntegral(const double *Dalpha, const double *alpha, integer innumC, integer innumP, integer indim, double *u);
		static void BackTrans(const double *u, const double *alpha, integer innumC, integer innumP, integer indim, double *utilde);
		static void GradVec(const double *utilde, const double *u, integer innumC, integer innumP, integer indim, double *w);

	private:
		integer numP; // number of points
		integer dim;  // dimension
		integer numC; // number of curves on a path
	};
} /*end of ROPTLIB namespace*/
#endif // end of EUCLIDEAN_H
