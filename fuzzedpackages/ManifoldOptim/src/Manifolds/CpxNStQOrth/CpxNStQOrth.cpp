
#include "CpxNStQOrth.h"

/*Define the namespace*/
namespace ROPTLIB{

	CpxNStQOrth::CpxNStQOrth(integer r, integer c)
	{
		n = r;
		p = c;
		IsIntrApproach = true;
		HasHHR = false;
		UpdBetaAlone = false;
		name.assign("Complex Noncompact Stiefel manifold quotient unitary group");
		IntrinsicDim = 2 * r * c - c * c;
		ExtrinsicDim = 2 * r * c;
		EMPTYEXTR = new CSOVector(2 * r, c);
		EMPTYINTR = new CSOVector(IntrinsicDim);
	};

	CpxNStQOrth::~CpxNStQOrth(void)
	{
		delete EMPTYEXTR;
		delete EMPTYINTR;
	};

	double CpxNStQOrth::Metric(Variable *x, Vector *etax, Vector *xix) const
	{
		Vector *exetax = EMPTYEXTR->ConstructEmpty();
		Vector *exxix = EMPTYEXTR->ConstructEmpty();
		ObtainExtr(x, etax, exetax);
		ObtainExtr(x, xix, exxix);
		double result = Manifold::Metric(x, exetax, exxix);
		delete exetax;
		delete exxix;
		return result;
	};

	void CpxNStQOrth::CheckParams(void) const
	{
		Manifold::CheckParams();
		OUTSTREAM << name << " PARAMETERS:" << std::endl;
		if (p == 1)
			OUTSTREAM << "n           :" << std::setw(15) << n << std::endl;
		else
		{
			OUTSTREAM << "n           :" << std::setw(15) << n << ",\t";
			OUTSTREAM << "p           :" << std::setw(15) << p << std::endl;
		}
	};

	void CpxNStQOrth::EucGradToGrad(Variable *x, Vector *egf, Vector *gf, const Problem *prob) const
	{
		if (prob->GetUseHess())
		{
			Vector *segf = egf->ConstructEmpty();
			segf->NewMemoryOnWrite(); // I don't remember the reason. It seems to be required.
			egf->CopyTo(segf);
			SharedSpace *Sharedegf = new SharedSpace(segf);
			x->AddToTempData("EGrad", Sharedegf);
		}
		ExtrProjection(x, egf, gf);
	};

	void CpxNStQOrth::EucHvToHv(Variable *x, Vector *etax, Vector *exix, Vector* xix, const Problem *prob) const
	{
		const double *xM = x->ObtainReadData();
		const SharedSpace *Sharedegf = x->ObtainReadTempData("EGrad");
		Vector *segf = Sharedegf->GetSharedElement();
		const double *egf = segf->ObtainReadData();

		doublecomplex *XHX = new doublecomplex[2 * p * p];
		doublecomplex *XHM = XHX + p * p;

		// XHX <- XM^H XM
		Matrix MxM(xM, n, p), MXHX((double*)XHX, p, p), Megf(egf, n, p), MXHM((double*)XHM, p, p);
		Matrix::CGEMM(GLOBAL::ZONE, MxM, true, MxM, false, GLOBAL::ZZERO, MXHX);
		// XHM <- XM^H egf
		Matrix::CGEMM(GLOBAL::ZONE, MxM, true, Megf, false, GLOBAL::ZZERO, MXHM);

		// XHM <- XHM - XHM^H
		for (integer i = 0; i < p; i++)
		{
			for (integer j = i; j < p; j++)
			{
				XHM[i + j * p].r -= XHM[j + i * p].r;
				XHM[i + j * p].i += XHM[j + i * p].i;
				XHM[j + i * p].r = -XHM[i + j * p].r;
				XHM[j + i * p].i = XHM[i + j * p].i;
			}
		}
		Matrix::CSYL(MXHX, MXHX, MXHM);

		exix->CopyTo(xix);
		double *xixTV = xix->ObtainWritePartialData();
		const double *etaxTV = etax->ObtainReadData();
		Matrix MetaxTV(etaxTV, n, p), MxixTV(xixTV, n, p);
		Matrix::CGEMM(GLOBAL::ZNONE, MetaxTV, false, MXHM, false,
			GLOBAL::ZONE, MxixTV);

		delete[] XHX;

		ExtrProjection(x, xix, xix);
	};

	void CpxNStQOrth::Projection(Variable *x, Vector *v, Vector *result) const
	{
		if (IsIntrApproach)
			IntrProjection(x, v, result);
		else
			ExtrProjection(x, v, result);
	};

	void CpxNStQOrth::Retraction(Variable *x, Vector *etax, Variable *result) const
	{
		Vector *exetax = EMPTYEXTR->ConstructEmpty();
		ObtainExtr(x, etax, exetax);
		Manifold::Retraction(x, exetax, result);
		delete exetax;
	};

	void CpxNStQOrth::coTangentVector(Variable *x, Vector *etax, Variable *y, Vector *xiy, Vector *result) const
	{
		OUTSTREAM << "TODO" << std::endl;//----
		xiy->CopyTo(result);
	};

	void CpxNStQOrth::DiffRetraction(Variable *x, Vector *etax, Variable *y, Vector *xix, Vector *result, bool IsEtaXiSameDir) const
	{
		Vector *exxix = EMPTYEXTR->ConstructEmpty();
		Vector *exresult = EMPTYEXTR->ConstructEmpty();
		ObtainExtr(x, xix, exxix);
		ExtrProjection(y, exxix, exresult);
		ObtainIntr(y, exresult, result);
		delete exxix;
		delete exresult;
	};

	double CpxNStQOrth::Beta(Variable *x, Vector *etax) const
	{
		return 1;
	};


	void CpxNStQOrth::ComputeHHR(Variable *x) const
	{
		const double *xM = x->ObtainReadData();
		SharedSpace *HouseHolderResult = new SharedSpace(2, x->Getsize()[0], x->Getsize()[1]);
		double *ptrHHR = HouseHolderResult->ObtainWriteEntireData();
		SharedSpace *HHRTau = new SharedSpace(1, x->Getsize()[1] * 2);
		double *tau = HHRTau->ObtainWriteEntireData();

		// for complex, N is the size divided by 2
		integer N = x->Getsize()[0] / 2, P = x->Getsize()[1], inc = 1;
		// ptrHHR <- xM, details: http://www.netlib.org/lapack/explore-html/da/d6c/dcopy_8f.html
		integer Length = 2 * N * P;
		dcopy_(&Length, const_cast<double *> (xM), &inc, ptrHHR, &inc);
		integer *jpvt = new integer[P];
		integer info;
		integer lwork = -1;
		doublecomplex lworkopt;
		double *rwork = new double[P];
		// compute the size of space required in the dgeqp3
#ifndef MATLAB_MEX_FILE
		zgeqp3_(&N, &P, (doublecomplex *)ptrHHR, &N, jpvt, (doublecomplex *)tau, &lworkopt, &lwork, rwork, &info);
#else
		zgeqp3_(&N, &P, ptrHHR, &N, jpvt, tau, (double *) &lworkopt, &lwork, rwork, &info);
#endif
		lwork = static_cast<integer> (lworkopt.r);
		doublecomplex *work = new doublecomplex[lwork];
		for (integer i = 0; i < P; i++)
			jpvt[i] = i + 1;
		// QR decomposition for ptrHHR using Householder reflections. Householder reflectors and R are stored in ptrHHR.
		// details: http://www.netlib.org/lapack/explore-html/db/de5/dgeqp3_8f.html
#ifndef MATLAB_MEX_FILE
		zgeqp3_(&N, &P, (doublecomplex *)ptrHHR, &N, jpvt, (doublecomplex *)tau, work, &lwork, rwork, &info);
#else
		zgeqp3_(&N, &P, ptrHHR, &N, jpvt, tau, (double *) work, &lwork, rwork, &info);
#endif

		x->AddToTempData("HHR", HouseHolderResult);
		x->AddToTempData("HHRTau", HHRTau);
		if (info < 0)
			OUTSTREAM << "Error in qr decomposition!" << std::endl;
		for (integer i = 0; i < P; i++)
		{
			if (jpvt[i] != (i + 1))
				OUTSTREAM << "Error in qf retraction!" << std::endl;
		}
		delete[] jpvt;
		delete[] work;
		delete[] rwork;
	};

	void CpxNStQOrth::ObtainIntr(Variable *x, Vector *etax, Vector *result) const
	{
		if (!x->TempDataExist("HHR"))
		{
			ComputeHHR(x);
		}

		const double *xM = x->ObtainReadData();
		const double *etaxTV = etax->ObtainReadData();
		const SharedSpace *HHR = x->ObtainReadTempData("HHR");
		const SharedSpace *HHRTau = x->ObtainReadTempData("HHRTau");
		double *resultTV = result->ObtainWriteEntireData();
		const double *ptrHHR = HHR->ObtainReadData();
		const double *ptrHHRTau = HHRTau->ObtainReadData();

		// for complex, N is the size divided by 2
		integer N = x->Getsize()[0] / 2, P = x->Getsize()[1], inc = 1, Length = 2 * N * P;
		integer info;
		integer lwork = -1;
		doublecomplex lworkopt;
		doublecomplex *tempspace = new doublecomplex[n * p];
		// compute the size of space required in the dormqr
#ifndef MATLAB_MEX_FILE
		zunmqr_(GLOBAL::L, GLOBAL::C, &N, &P, &P, (doublecomplex *)(const_cast<double *> (ptrHHR)), &N,
			(doublecomplex *)(const_cast<double *> (ptrHHRTau)), tempspace, &N, &lworkopt, &lwork, &info);
#else
		zunmqr_(GLOBAL::L, GLOBAL::C, &N, &P, &P, (const_cast<double *> (ptrHHR)), &N,
			(const_cast<double *> (ptrHHRTau)), (double *) tempspace, &N, (double *) &lworkopt, &lwork, &info);
#endif
		lwork = static_cast<integer> (lworkopt.r);
		doublecomplex *work = new doublecomplex[lwork];
		// tempspace <- etaxTV, details: http://www.netlib.org/lapack/explore-html/da/d6c/dcopy_8f.html
		dcopy_(&Length, const_cast<double *> (etaxTV), &inc, (double *)tempspace, &inc);
		// tempspace <- Q^T * tempspace, where Q is the orthogonal matrix defined as the product elementary reflectors defined by ptrHHR and ptrHHRTau,
		// details: http://www.netlib.org/lapack/explore-html/da/d82/dormqr_8f.html
#ifndef MATLAB_MEX_FILE
		zunmqr_(GLOBAL::L, GLOBAL::C, &N, &P, &P, (doublecomplex *)(const_cast<double *> (ptrHHR)), &N,
			(doublecomplex *)(const_cast<double *> (ptrHHRTau)), tempspace, &N, work, &lwork, &info);
#else
		zunmqr_(GLOBAL::L, GLOBAL::C, &N, &P, &P, (const_cast<double *> (ptrHHR)), &N,
			(const_cast<double *> (ptrHHRTau)), (double *) tempspace, &N, (double *) work, &lwork, &info);
#endif
		doublecomplex sign = { 0, 0 };
		for (integer i = 0; i < p; i++)
		{
			sign.r = (((doublecomplex *)ptrHHR)[i + n * i].r >= 0) ? 1 : -1;
#ifndef MATLAB_MEX_FILE
			zscal_(&P, &sign, tempspace + i, &N);
#else
			zscal_(&P, (double *) &sign, (double *) (tempspace + i), &N);
#endif
		}
		Matrix MxM(xM, n, p), MetaxTV(etaxTV, n, p), Mtempspace((double *)tempspace, p, p, n);
		Matrix::CGEMM(GLOBAL::ZONE, MxM, true, MetaxTV, false,
			GLOBAL::ZZERO, Mtempspace);
		//OUTSTREAM << Matrix((double *)xM, 2 * n, p) << std::endl;//---
		//OUTSTREAM << Matrix((double *)etaxTV, 2 * n, p) << std::endl;//---
		//OUTSTREAM << Matrix((double *)tempspace, 2 * n, p) << std::endl;//---
		double r2 = sqrt(2.0);
		double factor = sqrt(Manifold::Metric(x, x, x));
		integer idx = 0;
		for (integer i = 0; i < p; i++)
		{
			resultTV[idx] = tempspace[i + i * n].r / factor;
			idx++;
		}
		for (integer i = 0; i < p; i++)
		{
			for (integer j = i + 1; j < p; j++)
			{
				resultTV[idx] = r2 * tempspace[j + i * n].r / factor;
				idx++;
				resultTV[idx] = r2 * tempspace[j + i * n].i / factor;
				idx++;
			}
		}

		for (integer i = 0; i < p; i++)
		{
			for (integer j = p; j < n; j++)
			{
				resultTV[idx] = tempspace[j + i * n].r;
				idx++;
				resultTV[idx] = tempspace[j + i * n].i;
				idx++;
			}
		}
		delete[] work;
		delete[] tempspace;
	};

	void CpxNStQOrth::ObtainExtr(Variable *x, Vector *intretax, Vector *result) const
	{
		if (!x->TempDataExist("HHR"))
		{
			ComputeHHR(x);
		}

		const double *xM = x->ObtainReadData();
		const SharedSpace *HHR = x->ObtainReadTempData("HHR");
		const SharedSpace *HHRTau = x->ObtainReadTempData("HHRTau");
		const double *ptrHHR = HHR->ObtainReadData();
		const double *ptrHHRTau = HHRTau->ObtainReadData();
		const double *intretaxTV = intretax->ObtainReadData();
		double *resultTV = result->ObtainWriteEntireData();

		integer N = x->Getsize()[0] / 2, P = x->Getsize()[1], inc = 1, Length = 2 * N * P;
		integer info;
		integer idx = 0;
		doublecomplex *S = new doublecomplex[p * p];
		double r2 = sqrt(2.0);
		double factor = sqrt(Manifold::Metric(x, x, x));
		for (integer i = 0; i < p; i++)
		{
			S[i + i * p].r = intretaxTV[idx] * factor;
			S[i + i * p].i = 0;
			((doublecomplex *)resultTV)[i + i * n].r = 0;
			((doublecomplex *)resultTV)[i + i * n].i = 0;
			idx++;
		}
		for (integer i = 0; i < p; i++)
		{
			for (integer j = i + 1; j < p; j++)
			{
				S[j + i * p].r = intretaxTV[idx] / r2 * factor;
				S[i + j * p].r = intretaxTV[idx] / r2 * factor;
				((doublecomplex *)resultTV)[j + i * n].r = 0;
				((doublecomplex *)resultTV)[i + j * n].r = 0;
				idx++;
				S[j + i * p].i = intretaxTV[idx] / r2 * factor;
				S[i + j * p].i = -intretaxTV[idx] / r2 * factor;
				((doublecomplex *)resultTV)[j + i * n].i = 0;
				((doublecomplex *)resultTV)[i + j * n].i = 0;
				idx++;
			}
		}
		for (integer i = 0; i < p; i++)
		{
			for (integer j = p; j < n; j++)
			{
				((doublecomplex *)resultTV)[j + i * n].r = intretaxTV[idx];
				idx++;
				((doublecomplex *)resultTV)[j + i * n].i = intretaxTV[idx];
				idx++;
			}
		}

		doublecomplex sign = { 0, 0 };
		for (integer i = 0; i < p; i++)
		{
			sign.r = (((doublecomplex *)ptrHHR)[i + n * i].r >= 0) ? 1 : -1;
			// result(i, :) <- sign * result(i, :), details: http://www.netlib.org/lapack/explore-html/d4/dd0/dscal_8f.html
#ifndef MATLAB_MEX_FILE
			zscal_(&P, &sign, ((doublecomplex *)resultTV) + i, &N);
#else
			zscal_(&P, (double *) &sign, resultTV + 2 * i, &N);
#endif
		}
		integer lwork = -1;
		doublecomplex lworkopt;
		// compute the size of space required in the dormqr
#ifndef MATLAB_MEX_FILE
		zunmqr_(GLOBAL::L, GLOBAL::N, &N, &P, &P, (doublecomplex *)(const_cast<double *> (ptrHHR)), &N,
			(doublecomplex *)(const_cast<double *> (ptrHHRTau)), (doublecomplex *)resultTV, &N, &lworkopt, &lwork, &info);
#else
		zunmqr_(GLOBAL::L, GLOBAL::N, &N, &P, &P, (const_cast<double *> (ptrHHR)), &N,
			(const_cast<double *> (ptrHHRTau)), resultTV, &N, (double *) &lworkopt, &lwork, &info);
#endif
		lwork = static_cast<integer> (lworkopt.r);

		doublecomplex *work = new doublecomplex[lwork];
		// resultTV <- Q * resultTV, where Q is the orthogonal matrix defined as the product elementary reflectors defined by ptrHHR and ptrHHRTau,
		// details: http://www.netlib.org/lapack/explore-html/da/d82/dormqr_8f.html
#ifndef MATLAB_MEX_FILE
		zunmqr_(GLOBAL::L, GLOBAL::N, &N, &P, &P, (doublecomplex *)(const_cast<double *> (ptrHHR)), &N,
			(doublecomplex *)(const_cast<double *> (ptrHHRTau)), (doublecomplex *)resultTV, &N, work, &lwork, &info);
#else
		zunmqr_(GLOBAL::L, GLOBAL::N, &N, &P, &P, (const_cast<double *> (ptrHHR)), &N,
			(const_cast<double *> (ptrHHRTau)), resultTV, &N, (double *) work, &lwork, &info);
#endif
		delete[] work;

#ifndef MATLAB_MEX_FILE
		// solve for (R^H)^{-1} S
		ztrtrs_(GLOBAL::U, GLOBAL::C, GLOBAL::N, &P, &P, (doublecomplex *)ptrHHR, &N, S, &P, &info);
		// solve for R^{-1} S, this S is now equal to (X^H X)^{-1} S
		ztrtrs_(GLOBAL::U, GLOBAL::N, GLOBAL::N, &P, &P, (doublecomplex *)ptrHHR, &N, S, &P, &info);
#else
		// solve for (R^H)^{-1} S
		ztrtrs_(GLOBAL::U, GLOBAL::C, GLOBAL::N, &P, &P, const_cast<double *> (ptrHHR), &N, (double *) S, &P, &info);
		// solve for R^{-1} S, this S is now equal to (X^H X)^{-1} S
		ztrtrs_(GLOBAL::U, GLOBAL::N, GLOBAL::N, &P, &P, const_cast<double *> (ptrHHR), &N, (double *)S, &P, &info);
#endif
		Matrix MxM(xM, n, p), MS((double *)S, p, p), Mresult(resultTV, n, p);
		Matrix::CGEMM(GLOBAL::ZONE, MxM, false, MS, false, GLOBAL::ZONE, Mresult);
		delete[] S;
	};

	void CpxNStQOrth::IntrProjection(Variable *x, Vector *v, Vector *result) const
	{
		v->CopyTo(result);
	};

	void CpxNStQOrth::ExtrProjection(Variable *x, Vector *v, Vector *result) const
	{
		const double *xM = x->ObtainReadData();
		const double *V = v->ObtainReadData();

		doublecomplex *XHX = new doublecomplex[2 * p * p];
		doublecomplex *XHV = XHX + p * p;

		Matrix MxM(xM, n, p), MXHX((double*)XHX, p, p), MV(V, n, p), MXHV((double*)XHV, p, p);
		// XHX <- XM^H XM
		Matrix::CGEMM(GLOBAL::ZONE, MxM, true, MxM, false, GLOBAL::ZZERO, MXHX);
		// XHV <- XM^H V
		Matrix::CGEMM(GLOBAL::ZONE, MxM, true, MV, false, GLOBAL::ZZERO, MXHV);

		// XHV <- XHV - XHV^H
		for (integer i = 0; i < p; i++)
		{
			for (integer j = i; j < p; j++)
			{
				XHV[i + j * p].r -= XHV[j + i * p].r;
				XHV[i + j * p].i += XHV[j + i * p].i;
				XHV[j + i * p].r = -XHV[i + j * p].r;
				XHV[j + i * p].i = XHV[i + j * p].i;
			}
		}
		Matrix::CSYL(MXHX, MXHX, MXHV);
		v->CopyTo(result);
		double *resultTV = result->ObtainWritePartialData();
		Matrix MresultTV(resultTV, n, p);
		Matrix::CGEMM(GLOBAL::ZNONE, MxM, false, MXHV, false, GLOBAL::ZONE, MresultTV);

		delete[] XHX;
	};
} /*end of ROPTLIB namespace*/
