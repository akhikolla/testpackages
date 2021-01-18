
#include "GrassVariable.h"

/*Define the namespace*/
namespace ROPTLIB{

	GrassVariable::GrassVariable(integer r, integer c, integer n)
	{
		Element::Initialization(3, r, c, n);
	};

	GrassVariable *GrassVariable::ConstructEmpty(void) const
	{
		return new GrassVariable(size[0], size[1], size[2]);
	};

	void GrassVariable::RandInManifold(void)
	{
		this->RandGaussian();
		double *xU = this->ObtainWriteEntireData();
		integer n = size[0], p = size[1];
		integer *jpvt = new integer[p];
		integer lwork = 2 * p + (1 + p) * 64, info; // 64 = INITIALBLOCKSIZE
		double *tau = new double[p + lwork];
		double *work = tau + p;
		for (integer i = 0; i < p; i++)
			jpvt[i] = 0;
		// QR decomposition for xU using householder reflectors, reflectors and R are stored in xU.
		// details: http://www.netlib.org/lapack/explore-html/db/de5/dgeqp3_8f.html
		dgeqp3_(&n, &p, xU, &n, jpvt, tau, work, &lwork, &info);
		if (info < 0)
			Rprintf("Error in qr decomposition!\n");
		// Compute the orthonormal matrix by using the reflectors defined in xU and tau
		// details: http://www.netlib.org/lapack/explore-html/d9/d1d/dorgqr_8f.html
		dorgqr_(&n, &p, &p, xU, &n, tau, work, &lwork, &info);
		if (info < 0)
			Rprintf("Error in forming Q matrix!\n");
		delete[] jpvt;
		delete[] tau;
	};
}; /*end of ROPTLIB namespace*/
