
#ifndef PRESHAPEPATHSTRAIGHTEN_H
#define PRESHAPEPATHSTRAIGHTEN_H

#include "PreShapeCurves.h"
#include "PSCVariable.h"
#include "PSCVector.h"
#include "Problem.h"
#include "SharedSpace.h"
#include "ElasticCurvesRO.h"
#include "def.h"
#include "MyMatrix.h"

/*Define the namespace*/
namespace ROPTLIB{

	//iniPath is the initial path connecting two end points
	//numP is the number of points on each curve
	//dim is the dimension of the space where the curve is in
	//numC is the number of curves on each path

	class PreShapePathStraighten : public Problem{
	public:
		PreShapePathStraighten(integer innumP, integer indim, integer innumC);
		virtual ~PreShapePathStraighten();
		virtual double f(Variable *x) const;    //cost function

		//virtual void RieGrad(Variable *x, Vector *gf) const;
		//virtual void RieHessianEta(Variable *x, Vector *etax, Vector *xix) const;

		virtual void EucGrad(Variable *x, Vector *egf) const;   //
		static void Item_1(const double *q, integer innumP, integer indim, double *q_c); //Projection onto closed curve space
		static void Item_2(const double *q, integer innumP, integer indim, double *w); //Projecting arbitrary points in L2 into tangent space at q
		static void Item_3(const double *w, const double *q1, const double *q2, integer innumP, integer indim, double *wbar); //Parallel Transport
		static double InnerProd_Q(const double *q1, const double *q2, integer innumP, integer indim);
		//virtual void EucHessianEta(Variable *x, Vector *etax, Vector *exix) const;

	private:
		integer numP;
		integer dim;
		integer numC;
	};
} /*end of ROPTLIB namespace*/
#endif // end of PRESHAPEPATHSTRAIGHTEN_H
