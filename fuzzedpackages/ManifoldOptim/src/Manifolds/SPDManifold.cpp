#include "SPDManifold.h"

/*Define the namespace*/
namespace ROPTLIB{

	SPDManifold::SPDManifold(integer inn)
	{
		n = inn;
		IsIntrApproach = true;
		HasHHR = false;
		UpdBetaAlone = false;
		HasLockCon = false;
		name.assign("SPDManifold");
		IntrinsicDim = n * (n + 1) / 2;
		ExtrinsicDim = n * n;
		EMPTYEXTR = new SPDVector(n, n);
		EMPTYINTR = new SPDVector(IntrinsicDim, 1);
	};

	SPDManifold::~SPDManifold(void)
	{
		delete EMPTYEXTR;
		delete EMPTYINTR;
	};

	void SPDManifold::CheckParams(void) const
	{
		Manifold::CheckParams();
		OUTSTREAM << name << " PARAMETERS:" << std::endl;
		OUTSTREAM << "row           :" << std::setw(15) << n << ",\t";
		OUTSTREAM << "col           :" << std::setw(15) << n << std::endl;
	};

	void SPDManifold::EucGradToGrad(Variable *x, Vector *egf, Vector *gf, const Problem *prob) const
	{
		const double *xptr = x->ObtainReadData();
		const double *egfptr = egf->ObtainReadData();
		double *gfptr = gf->ObtainWriteEntireData();
		integer dim = n;
		double *tmp = new double[dim * dim];
		dgemm_(GLOBAL::T, GLOBAL::N, &dim, &dim, &dim, &GLOBAL::DONE, const_cast<double *> (xptr), &dim, const_cast<double *> (egfptr), &dim,
			&GLOBAL::DZERO, tmp, &dim);
		dgemm_(GLOBAL::N, GLOBAL::N, &dim, &dim, &dim, &GLOBAL::DONE, const_cast<double *> (tmp), &dim, const_cast<double *> (xptr), &dim,
			&GLOBAL::DZERO, gfptr, &dim);
		delete[] tmp;
	};

	void SPDManifold::EucHvToHv(Variable *x, Vector *etax, Vector *exix, Vector* xix, const Problem *prob) const
	{
		OUTSTREAM << "warning:SPDManifold::EucHvToHv has not been done!" << std::endl;
		exix->CopyTo(xix);
	};

	void SPDManifold::CholeskyRepresentation(Variable *x) const
	{
		const double *xM = x->ObtainReadData();
		Variable *L = x->ConstructEmpty();
		SharedSpace *SharedL = new SharedSpace(L);
		double *LM = L->ObtainWriteEntireData();
		for (integer i = 0; i < n; i++)
		{
			for (integer j = i; j < n; j++)
			{
				LM[i + j * n] = 0;
				LM[j + i * n] = xM[j + i * n];
			}
		}

		integer info, N = n;
		dpotrf_(GLOBAL::L, &N, LM, &N, &info);
		x->AddToTempData("L", SharedL);
		if (info != 0)
		{
			OUTSTREAM << "Warning: SPDManifold::CholeskyRepresentation fails with info:" << info << "!" << std::endl;
		}
	};

	void SPDManifold::ExtrProjection(Variable *x, Vector *etax, Vector *result) const
	{
		const double *etaxTV = etax->ObtainReadData();
		double *resultTV = result->ObtainWriteEntireData();
		for (integer i = 0; i < n; i++)
		{
			resultTV[i + i * n] = etaxTV[i + i * n];
			for (integer j = i + 1; j < n; j++)
			{
				resultTV[i + j * n] = (etaxTV[i + j * n] + etaxTV[j + i * n]) * 0.5;
				resultTV[j + i * n] = resultTV[i + j * n];
			}
		}
	};

	void SPDManifold::ObtainIntr(Variable *x, Vector *etax, Vector *result) const
	{//L^{-1} etax L^{-T}, where x = L L^T
		if (!x->TempDataExist("L"))
		{
			CholeskyRepresentation(x);
		}
		const SharedSpace *SharedL = x->ObtainReadTempData("L");
		Variable *LElement = SharedL->GetSharedElement();
		const double *L = LElement->ObtainReadData();
		double *E = new double[n * n];
		integer length = n * n, N = n, info;
		dcopy_(&length, const_cast<double *> (etax->ObtainReadData()), &GLOBAL::IONE, E, &GLOBAL::IONE);
		/*Solve the linear system L X = E, i.e., X = L^{-1} E. The solution X is stored in E
		Details: http://www.netlib.org/lapack/explore-html/d6/d6f/dtrtrs_8f.html*/
		dtrtrs_(GLOBAL::L, GLOBAL::N, GLOBAL::N, &N, &N, const_cast<double *> (L), &N, E, &N, &info);
		if (info != 0)
		{
			OUTSTREAM << "warning: SPDManifold::ObtainIntr fails with info:" << info << "!" << std::endl;
		}

		/*E <-- E^T*/
		double tmp;
		for (integer i = 0; i < n; i++)
		{
			for (integer j = i + 1; j < n; j++)
			{
				tmp = E[i + j * n];
				E[i + j * n] = E[j + i * n];
				E[j + i * n] = tmp;
			}
		}
		/*Solve the linear system L X = E, i.e., X = L^{-1} E. The solution X is stored in E
		Details: http://www.netlib.org/lapack/explore-html/d6/d6f/dtrtrs_8f.html */
		dtrtrs_(GLOBAL::L, GLOBAL::N, GLOBAL::N, &N, &N, const_cast<double *> (L), &N, E, &N, &info);
		if (info != 0)
		{
			OUTSTREAM << "warning: SPDManifold::ObtainIntr fails with info:" << info << "!" << std::endl;
		}
		/*We don't have to do: E <-- E^T, since E is symmetric*/
		double *resultTV = result->ObtainWriteEntireData();
		integer idx = 0;
		double r2 = sqrt(2.0);
		for (integer i = 0; i < n; i++)
		{
			resultTV[idx] = E[i + i * n];
			idx++;
		}

		for (integer i = 0; i < n; i++)
		{
			for (integer j = i + 1; j < n; j++)
			{
				resultTV[idx] = E[j + i * n] * r2;
				idx++;
			}
		}

		delete[] E;
	};

	void SPDManifold::ObtainExtr(Variable *x, Vector *intretax, Vector *result) const
	{
		if (!x->TempDataExist("L"))
		{
			CholeskyRepresentation(x);
		}
		const SharedSpace *SharedL = x->ObtainReadTempData("L");
		Variable *LElement = SharedL->GetSharedElement();
		const double *L = LElement->ObtainReadData();
		const double *intretaxTV = intretax->ObtainReadData();
		double *resultTV = result->ObtainWriteEntireData();

		integer idx = 0;
		double r2 = sqrt(2.0);
		for (integer i = 0; i < n; i++)
		{
			resultTV[i + i * n] = intretaxTV[idx];
			idx++;
		}

		for (integer i = 0; i < n; i++)
		{
			for (integer j = i + 1; j < n; j++)
			{
				resultTV[j + i * n] = intretaxTV[idx] / r2;
				resultTV[i + j * n] = resultTV[j + i * n];
				idx++;
			}
		}

		double *E = new double[n * n];
		integer N = n;

		/*E <-- L resultTV */
		dgemm_(GLOBAL::N, GLOBAL::N, &N, &N, &N, &GLOBAL::DONE, const_cast<double *> (L), &N, resultTV, &N,
			&GLOBAL::DZERO, E, &N);
		/*resultTV <-- E L^T */
		dgemm_(GLOBAL::N, GLOBAL::T, &N, &N, &N, &GLOBAL::DONE, E, &N, const_cast<double *> (L), &N,
			&GLOBAL::DZERO, resultTV, &N);

		delete[] E;
	};

	void SPDManifold::Retraction(Variable *x, Vector *etax, Variable *result) const
	{ // assume intrinsic representation is used for the tangent vector etax.
		if (!x->TempDataExist("L"))
		{
			CholeskyRepresentation(x);
		}
		const SharedSpace *SharedL = x->ObtainReadTempData("L");
		Variable *LElement = SharedL->GetSharedElement();
		const double *L = LElement->ObtainReadData();

		/*Compute the extrinsic representation*/
		Vector *exetax = EMPTYEXTR->ConstructEmpty();
		ObtainExtr(x, etax, exetax);

		integer N = n, info;
		double *LiE = new double[N * N];
		const double *etaxTV = exetax->ObtainReadData();
		integer length = n * n;
		dcopy_(&length, const_cast<double *> (etaxTV), &GLOBAL::IONE, LiE, &GLOBAL::IONE);
		/*Solve the linear system L X = E, i.e., X = L^{-1} E. The solution X is stored in LiE
		Details: http://www.netlib.org/lapack/explore-html/d6/d6f/dtrtrs_8f.html */
		dtrtrs_(GLOBAL::L, GLOBAL::N, GLOBAL::N, &N, &N, const_cast<double *> (L), &N, LiE, &N, &info);
		if (info != 0)
		{
			OUTSTREAM << "warning: SPDManifold::Retraction fails with info:" << info << "!" << std::endl;
		}
		double *resultTV = result->ObtainWriteEntireData();
		/*Compute result = LiE^T LiE = E L^{-T} L^{-1} E*/
		dgemm_(GLOBAL::T, GLOBAL::N, &N, &N, &N, &GLOBAL::DONE, LiE, &N, LiE, &N, &GLOBAL::DZERO, resultTV, &N);
		delete[] LiE;
		double scalar = 0.5;
		/*result <-- 0.5 * result*/
		dscal_(&length, &scalar, resultTV, &GLOBAL::IONE);

		/*result <-- etax + result*/
		daxpy_(&length, &GLOBAL::DONE, const_cast<double *> (etaxTV), &GLOBAL::IONE, resultTV, &GLOBAL::IONE);

		const double *xM = x->ObtainReadData();
		/*result <-- x + result*/
		daxpy_(&length, &GLOBAL::DONE, const_cast<double *> (xM), &GLOBAL::IONE, resultTV, &GLOBAL::IONE);
		delete exetax;

		if (!result->TempDataExist("L"))
		{
			CholeskyRepresentation(result);
		}
	};

	void SPDManifold::DiffRetraction(Variable *x, Vector *etax, Variable *y, Vector *xix, Vector *result, bool IsEtaXiSameDir) const
	{ // xix + 0.5 xix x^{-1} etax + 0.5 etax x^{-1} xix
		if (!x->TempDataExist("L"))
		{
			CholeskyRepresentation(x);
		}
		const SharedSpace *SharedL = x->ObtainReadTempData("L");
		Variable *LElement = SharedL->GetSharedElement();
		const double *L = LElement->ObtainReadData();

		/*Compute the extrinsic representation*/
		Vector *exetax = EMPTYEXTR->ConstructEmpty();
		Vector *exxix = EMPTYEXTR->ConstructEmpty();
		ObtainExtr(x, etax, exetax);
		ObtainExtr(x, xix, exxix);
		double *LiE = new double[2 * n * n];
		double *LiX = LiE + n * n;
		const double *exetaxTV = exetax->ObtainReadData();
		const double *exxixTV = exxix->ObtainReadData();
		integer length = n * n, N = n, info;
		dcopy_(&length, const_cast<double*> (exetaxTV), &GLOBAL::IONE, LiE, &GLOBAL::IONE);
		dcopy_(&length, const_cast<double*> (exxixTV), &GLOBAL::IONE, LiX, &GLOBAL::IONE);
		delete exetax;

		/*Solve the linear system L X = Eta, i.e., X = L^{-1} Eta. The solution X is stored in LiE
		Details: http://www.netlib.org/lapack/explore-html/d6/d6f/dtrtrs_8f.html */
		dtrtrs_(GLOBAL::L, GLOBAL::N, GLOBAL::N, &N, &N, const_cast<double *> (L), &N, LiE, &N, &info);

		/*Solve the linear system L X = Xix, i.e., X = L^{-1} Xix. The solution X is stored in LiX
		Details: http://www.netlib.org/lapack/explore-html/d6/d6f/dtrtrs_8f.html */
		dtrtrs_(GLOBAL::L, GLOBAL::N, GLOBAL::N, &N, &N, const_cast<double *> (L), &N, LiX, &N, &info);

		Vector *exresult = EMPTYEXTR->ConstructEmpty();
		double *resultTV = exresult->ObtainWriteEntireData();
		/*resultTV <-- etax L^{-T} L^{-1} xix*/
		dgemm_(GLOBAL::T, GLOBAL::N, &N, &N, &N, &GLOBAL::DONE, LiE, &N, LiX, &N, &GLOBAL::DZERO, resultTV, &N);
		for (integer i = 0; i < n; i++)
		{
			for (integer j = i + 1; j < n; j++)
			{
				resultTV[j + i * n] = (resultTV[j + i * n] + resultTV[i + j * n]) * 0.5;
				resultTV[i + j * n] = resultTV[j + i * n];
			}
		}
		delete[] LiE;

		daxpy_(&length, &GLOBAL::DONE, const_cast<double*> (exxixTV), &GLOBAL::IONE, resultTV, &GLOBAL::IONE);
		delete exxix;
		ObtainIntr(y, exresult, result);
		delete exresult;

		if (IsEtaXiSameDir && (HasHHR || UpdBetaAlone))
		{
			const double *etaxTV = etax->ObtainReadData();
			const double *xixTV = xix->ObtainReadData();
			double EtatoXi = sqrt(Metric(x, etax, etax) / Metric(x, xix, xix));
			SharedSpace *beta = new SharedSpace(1, 3);
			double *betav = beta->ObtainWriteEntireData();
			betav[0] = sqrt(Metric(x, etax, etax) / Metric(x, result, result)) / EtatoXi;
			betav[1] = Metric(x, etax, etax);
			betav[2] = Metric(x, result, result) * EtatoXi * EtatoXi;
			etax->AddToTempData("beta", beta);

			if (HasHHR)
			{
				Vector *TReta = result->ConstructEmpty();
				result->CopyTo(TReta);
				ScaleTimesVector(x, betav[0] * EtatoXi, TReta, TReta);
				SharedSpace *SharedTReta = new SharedSpace(TReta);
				etax->AddToTempData("betaTReta", SharedTReta);
			}
		}
	};

	void SPDManifold::coTangentVector(Variable *x, Vector *etax, Variable *y, Vector *xiy, Vector *result) const
	{
		OUTSTREAM << "SPDManifold::coTangentVector has not been done!" << std::endl;
		xiy->CopyTo(result);
	};

	double SPDManifold::Beta(Variable *x, Vector *etax) const
	{
		if (!HasHHR && !UpdBetaAlone)
			return 1;

		if (!etax->TempDataExist("beta"))
		{ /*In case that beta is not computed, then compute it.*/
			Variable *y = x->ConstructEmpty();
			Vector *xiy = etax->ConstructEmpty();
			Retraction(x, etax, y);
			DiffRetraction(x, etax, y, etax, xiy, true);
			delete y;
			delete xiy;
		}

		/*If the beta has been computed in differentiated retraction, then obtain it.
		Beta should be almost always computed before.*/
		const SharedSpace *beta = etax->ObtainReadTempData("beta");
		const double *betav = beta->ObtainReadData();
		return betav[0];
	};
} /*end of ROPTLIB namespace*/
