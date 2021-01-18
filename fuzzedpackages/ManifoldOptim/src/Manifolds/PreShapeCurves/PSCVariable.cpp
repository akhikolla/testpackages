
#include "PSCVariable.h"

/*Define the namespace*/
namespace ROPTLIB{

	PSCVariable::PSCVariable(integer r, integer c, integer n)
	{
		Element::Initialization(3, r, c, n);
	};

	PSCVariable *PSCVariable::ConstructEmpty(void) const
	{
		return new PSCVariable(size[0], size[1], size[2]);
	};

	void PSCVariable::RandInManifold(void)
	{
		Element::RandGaussian();
		OUTSTREAM << "PSCVariable::RandInManifold(). TODO" << std::endl;//--
	};

	void PSCVariable::Generate(double *initial, double *end)
	{
		double *q1, *q2;
		double theta, tau, coeff;
		q1 = initial;
		q2 = end;
		NewMemoryOnWrite();   //???????????
		integer numP = size[0];
		integer dim = size[1];
		integer numC = size[2];
		integer NXD = numP*dim;
		double *temp = new double[numP*dim];


		theta = std::acos(PreShapePathStraighten::InnerProd_Q(q1, q2, numP, dim));
		if (theta > 0.0001)
		{
			for (integer t = 0; t < numC; t++)
			{
				tau = static_cast<double> (t) / (numC - 1);
				for (integer i = 0; i < numP*dim; i++)
				{
					temp[i] = (std::sin((1.0 - tau)*theta)*q1[i] + std::sin(tau*theta)*q2[i]) / std::sin(theta);
				}
				PreShapePathStraighten::Item_1(temp, numP, dim, Space + t*numP*dim);
				coeff = 1.0 / std::sqrt(PreShapePathStraighten::InnerProd_Q(Space + t*numP*dim, Space + t*numP*dim, numP, dim));
				dscal_(&NXD, &coeff, Space + t*numP*dim, &GLOBAL::IONE);
			}
		}
		else
		{
			for (integer t = 0; t < numC; t++)
			{
				dcopy_(&NXD, q1, &GLOBAL::IONE, Space + t*numP*dim, &GLOBAL::IONE);
			}
		}

		delete[] temp;
	}
} /*end of ROPTLIB namespace*/
