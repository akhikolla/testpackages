
#include "Manifold.h"

/*Define the namespace*/
namespace ROPTLIB{

	Manifold::~Manifold(void)
	{
	};

	double Manifold::Metric(Variable *x, Vector *etax, Vector *xix) const
	{
		const double *v1 = etax->ObtainReadData();
		const double *v2 = xix->ObtainReadData();

		integer inc = 1, N = etax->Getlength();
		// output v1^T v2, details: http://www.netlib.org/lapack/explore-html/d5/df6/ddot_8f.html
		return ddot_(&N, const_cast<double *> (v1), &inc, const_cast<double *> (v2), &inc);
	};

	void Manifold::LinearOPEEta(Variable *x, LinearOPE *Hx, Vector *etax, Vector *result) const
	{
		if (etax == result)
		{
			OUTSTREAM << "The arguments of etax and result should not be the same!" << std::endl;
		}
		integer ell = Hx->Getsize()[0];
		const double *v = etax->ObtainReadData();
		const double *M = Hx->ObtainReadData();
		double *resultTV = result->ObtainWriteEntireData();

		char *transn = const_cast<char *> ("n");
		double one = 1, zero = 0;
		integer inc = 1, N = ell;
		// resultTV <- M * v; details: http://www.netlib.org/lapack/explore-html/dc/da8/dgemv_8f.html
		dgemv_(transn, &N, &N, &one, const_cast<double *> (M), &N, const_cast<double *> (v), &inc, &zero, resultTV, &inc);
	};

	void Manifold::ScaleTimesVector(Variable *x, double scalar, Vector *etax, Vector *result) const
	{
		etax->CopyTo(result);
		double *resultTV = result->ObtainWritePartialData();
		integer N = etax->Getlength(), inc = 1;
		// resultTV <- scalar * resultTV, details: http://www.netlib.org/lapack/explore-html/d4/dd0/dscal_8f.html
		dscal_(&N, &scalar, resultTV, &inc);
	};

	void Manifold::VectorAddVector(Variable *x, Vector *etax, Vector *xix, Vector *result) const
	{
		VectorLinearCombination(x, 1.0, etax, 1.0, xix, result);
	};

	void Manifold::VectorMinusVector(Variable *x, Vector *etax, Vector *xix, Vector *result) const
	{
		VectorLinearCombination(x, 1.0, etax, -1.0, xix, result);
	};

	void Manifold::scalarVectorAddVector(Variable *x, double scalar, Vector *etax, Vector *xix, Vector *result) const
	{
		VectorLinearCombination(x, scalar, etax, 1.0, xix, result);
	};

	void Manifold::scalarVectorMinusVector(Variable *x, double scalar, Vector *etax, Vector *xix, Vector *result) const
	{
		VectorLinearCombination(x, scalar, etax, -1.0, xix, result);
	};

	void Manifold::VectorLinearCombination(Variable *x, double scalar1, Vector *etax, double scalar2, Vector *xix, Vector *result) const
	{
		const double *etaxTV = etax->ObtainReadData();
		const double *xixTV = xix->ObtainReadData();
		double *resultTV = result->ObtainWriteEntireData();

		integer N1 = etax->Getlength(), N2 = xix->Getlength(), N3 = result->Getlength(), inc = 1;
		integer N = (N1 > N2) ? N2 : N1;
		N = (N > N3) ? N3 : N;
		if (resultTV == etaxTV)
		{
			// resultTV <- scalar1 * resultTV, details: http://www.netlib.org/lapack/explore-html/d4/dd0/dscal_8f.html
			dscal_(&N, &scalar1, resultTV, &inc);
			// resultTV <- scalar2 * xixTV + resultTV, details: http://www.netlib.org/lapack/explore-html/d9/dcd/daxpy_8f.html
			daxpy_(&N, &scalar2, const_cast<double *> (xixTV), &inc, resultTV, &inc);
		}
		else
			if (resultTV == xixTV)
			{
				// resultTV <- scalar2 * resultTV, details: http://www.netlib.org/lapack/explore-html/d4/dd0/dscal_8f.html
				dscal_(&N, &scalar2, resultTV, &inc);
				// resultTV <- scalar1 * etaxTV + resultTV, details: http://www.netlib.org/lapack/explore-html/d9/dcd/daxpy_8f.html
				daxpy_(&N, &scalar1, const_cast<double *> (etaxTV), &inc, resultTV, &inc);
			}
			else
			{
				// resultTV <- etaxTV, details: http://www.netlib.org/lapack/explore-html/da/d6c/dcopy_8f.html
				dcopy_(&N, const_cast<double *> (etaxTV), &inc, resultTV, &inc);
				// resultTV <- scalar1 * resultTV, details: http://www.netlib.org/lapack/explore-html/d4/dd0/dscal_8f.html
				dscal_(&N, &scalar1, resultTV, &inc);
				// resultTV <- scalar2 * xixTV + resultTV, details: http://www.netlib.org/lapack/explore-html/d9/dcd/daxpy_8f.html
				daxpy_(&N, &scalar2, const_cast<double *> (xixTV), &inc, resultTV, &inc);
			}
	};

	void Manifold::Projection(Variable *x, Vector *v, Vector *result) const
	{
		v->CopyTo(result);
	};

	void Manifold::RandomTangentVectors(Variable *x, integer N, Vector **result_arr) const // Be careful
	{
	};

	void Manifold::Retraction(Variable *x, Vector *etax, Variable *result) const
	{
		const double *v = etax->ObtainReadData();
		const double *xM = x->ObtainReadData();
		double *resultM = result->ObtainWriteEntireData();

		integer inc = 1, N = x->Getlength();
		double one = 1;
		// resultM <- xM, details: http://www.netlib.org/lapack/explore-html/da/d6c/dcopy_8f.html
		if (resultM != xM)
			dcopy_(&N, const_cast<double *> (xM), &inc, resultM, &inc);
		// resultM <- v + resultTV, details: http://www.netlib.org/lapack/explore-html/d9/dcd/daxpy_8f.html
		daxpy_(&N, &one, const_cast<double *> (v), &inc, resultM, &inc);
	};

	void Manifold::coTangentVector(Variable *x, Vector *etax, Variable *y, Vector *xiy, Vector *result) const
	{
		xiy->CopyTo(result);
	};

	void Manifold::DiffRetraction(Variable *x, Vector *etax, Variable *y, Vector *xix, Vector *result, bool IsEtaXiSameDir) const
	{
		xix->CopyTo(result);
	};

	double Manifold::Beta(Variable *x, Vector *etax) const
	{
		return 1;
	};

	void Manifold::VectorTransport(Variable *x, Vector *etax, Variable *y, Vector *xix, Vector *result) const
	{
		if (HasHHR)
			LCVectorTransport(x, etax, y, xix, result);
		else
			xix->CopyTo(result);
	};

	void Manifold::InverseVectorTransport(Variable *x, Vector *etax, Variable *y, Vector *xiy, Vector *result) const
	{
		if (HasHHR)
			LCInverseVectorTransport(x, etax, y, xiy, result);
		else
			xiy->CopyTo(result);
	};

	void Manifold::HInvTran(Variable *x, Vector *etax, Variable *y, LinearOPE *Hx, integer start, integer end, LinearOPE *result) const
	{
		if (HasHHR)
			LCHInvTran(x, etax, y, Hx, start, end, result);
		else
			Hx->CopyTo(result);
	};

	void Manifold::TranH(Variable *x, Vector *etax, Variable *y, LinearOPE *Hx, integer start, integer end, LinearOPE *result) const
	{
		if (HasHHR)
			LCTranH(x, etax, y, Hx, start, end, result);
		else
			Hx->CopyTo(result);
	};

	void Manifold::TranHInvTran(Variable *x, Vector *etax, Variable *y, LinearOPE *Hx, LinearOPE *result) const
	{
		if (HasHHR)
			LCTranHInvTran(x, etax, y, Hx, result);
		else
			Hx->CopyTo(result);
	};

	void Manifold::HaddScaledRank1OPE(Variable *x, LinearOPE *Hx, double scalar, Vector *etax, Vector *xix, LinearOPE *result) const
	{
		const double *veta = etax->ObtainReadData();
		const double *vxi = xix->ObtainReadData();
		Hx->CopyTo(result);
		double *resultL = result->ObtainWritePartialData();
		integer ell = Hx->Getsize()[0], N = ell, inc = 1;
		// resultL <- scalar * veta * vxi^T + resultL, details: http://www.netlib.org/lapack/explore-html/dc/da8/dger_8f.html
		dger_(&N, &N, &scalar, const_cast<double *> (veta), &inc, const_cast<double *> (vxi), &inc, resultL, &N);
	};

	void Manifold::ObtainEtaxFlat(Variable *x, Vector *etax, Vector *etaxflat) const
	{
		etax->CopyTo(etaxflat);
	};

	void Manifold::ObtainIntr(Variable *x, Vector *etax, Vector *result) const
	{
		etax->CopyTo(result);
	};

	void Manifold::ObtainExtr(Variable *x, Vector *intretax, Vector *result) const
	{
		intretax->CopyTo(result);
	};

	void Manifold::Obtainnu1nu2forLC(Variable *x, Vector *etax, Variable *y) const
	{
		Vector *eps1 = etax->ConstructEmpty();
		Vector *nu1 = etax->ConstructEmpty();
		Vector *nu2 = etax->ConstructEmpty();
		if (!etax->TempDataExist("beta") || !etax->TempDataExist("betaTReta"))
		{
			DiffRetraction(x, etax, y, etax, eps1, true);
		}
		HasHHR = false; VectorTransport(x, etax, y, etax, eps1); HasHHR = true;
		const SharedSpace *TReta = etax->ObtainReadTempData("betaTReta");
		Vector *TRetaVector = TReta->GetSharedElement();
		SharedSpace *Sharedtau1tau2 = new SharedSpace(1, 2);
		SharedSpace *Sharednu1 = new SharedSpace(nu1);
		SharedSpace *Sharednu2 = new SharedSpace(nu2);
		double *tau1tau2 = Sharedtau1tau2->ObtainWriteEntireData();
		ScaleTimesVector(x, 2.0, eps1, nu1);
		VectorLinearCombination(x, -1.0, eps1, -1.0, TRetaVector, nu2);
		tau1tau2[0] = 2.0 / Metric(x, nu1, nu1);
		tau1tau2[1] = 2.0 / Metric(x, nu2, nu2);

		etax->AddToTempData("tau1tau2", Sharedtau1tau2);
		etax->AddToTempData("nu1", Sharednu1);
		etax->AddToTempData("nu2", Sharednu2);
		delete eps1;
	};

	void Manifold::LCVectorTransport(Variable *x, Vector *etax, Variable *y, Vector *xix, Vector *result) const
	{
		if (!etax->TempDataExist("nu1nu2"))
		{
			Obtainnu1nu2forLC(x, etax, y);
		}
		HasHHR = false; VectorTransport(x, etax, y, xix, result); HasHHR = true;
		const SharedSpace *Sharedtau1tau2 = etax->ObtainReadTempData("tau1tau2");
		const double *tau1tau2 = Sharedtau1tau2->ObtainReadData();
		const SharedSpace *Sharednu1 = etax->ObtainReadTempData("nu1");
		Vector *nu1 = Sharednu1->GetSharedElement();
		const SharedSpace *Sharednu2 = etax->ObtainReadTempData("nu2");
		Vector *nu2 = Sharednu2->GetSharedElement();
		double temp = 0;
		temp = -Metric(x, result, nu1);
		scalarVectorAddVector(x, temp * tau1tau2[0], nu1, result, result);
		temp = -Metric(x, result, nu2);
		scalarVectorAddVector(x, temp * tau1tau2[1], nu2, result, result);
	};

	void Manifold::LCInverseVectorTransport(Variable *x, Vector *etax, Variable *y, Vector *xiy, Vector *result) const
	{
		if (!etax->TempDataExist("nu1nu2"))
		{
			Obtainnu1nu2forLC(x, etax, y);
		}
		const SharedSpace *Sharedtau1tau2 = etax->ObtainReadTempData("tau1tau2");
		const double *tau1tau2 = Sharedtau1tau2->ObtainReadData();
		const SharedSpace *Sharednu1 = etax->ObtainReadTempData("nu1");
		Vector *nu1 = Sharednu1->GetSharedElement();
		const SharedSpace *Sharednu2 = etax->ObtainReadTempData("nu2");
		Vector *nu2 = Sharednu2->GetSharedElement();
		double temp = 0;
		temp = -Metric(x, xiy, nu2);
		scalarVectorAddVector(x, temp * tau1tau2[1], nu2, xiy, result);
		temp = -Metric(x, result, nu1);
		scalarVectorAddVector(x, temp * tau1tau2[0], nu1, result, result);
		HasHHR = false; InverseVectorTransport(x, etax, y, result, result); HasHHR = true;
	};

	void Manifold::LCHInvTran(Variable *x, Vector *etax, Variable *y, LinearOPE *Hx, integer start, integer end, LinearOPE *result) const
	{
		if (!etax->TempDataExist("nu1nu2"))
		{
			Obtainnu1nu2forLC(x, etax, y);
		}
		const SharedSpace *Sharedtau1tau2 = etax->ObtainReadTempData("tau1tau2");
		const double *tau1tau2 = Sharedtau1tau2->ObtainReadData();
		const SharedSpace *Sharednu1 = etax->ObtainReadTempData("nu1");
		Vector *nu1 = Sharednu1->GetSharedElement();
		const SharedSpace *Sharednu2 = etax->ObtainReadTempData("nu2");
		Vector *nu2 = Sharednu2->GetSharedElement();

		double temp = 0;
		const double *nu1TV = nu1->ObtainReadData();
		const double *nu2TV = nu2->ObtainReadData();

		//double *nu1TV = new double[nu1->Getlength() * 2];
		//double *nu2TV = nu1TV + nu1->Getlength();
		//nu1->CopytoArray(nu1TV);
		//nu2->CopytoArray(nu2TV);
		HasHHR = false; HInvTran(x, etax, y, Hx, start, end, result); HasHHR = true;
		double *resultTV = result->ObtainWritePartialData();
		char *sider = const_cast<char *> ("r");
		integer ell = Hx->Getsize()[0], length = etax->Getlength();
		double *work = new double[ell];

		// resultTV(:, start : start + length - 1) <- resultTV(:, start : start + length - 1) * (I - tau1tau2(0) * nu1TV * nu1TV^T),
		// details: http://www.netlib.org/lapack/explore-html/db/d10/dlarfx_8f.html
		dlarfx_(sider, &ell, &length, const_cast<double *> (nu1TV), const_cast<double *> (tau1tau2), resultTV + start * ell, &ell, work);
		// resultTV(:, start : start + length - 1) <- resultTV(:, start : start + length - 1) * (I - tau1tau2(1) * nu2TV * nu2TV^T),
		// details: http://www.netlib.org/lapack/explore-html/db/d10/dlarfx_8f.html
		dlarfx_(sider, &ell, &length, const_cast<double *> (nu2TV), const_cast<double *> (tau1tau2 + 1), resultTV + start * ell, &ell, work);
		delete[] work;
		//delete[] nu1TV;
	};

	void Manifold::LCTranH(Variable *x, Vector *etax, Variable *y, LinearOPE *Hx, integer start, integer end, LinearOPE *result) const
	{
		if (!etax->TempDataExist("nu1nu2"))
		{
			Obtainnu1nu2forLC(x, etax, y);
		}
		const SharedSpace *Sharedtau1tau2 = etax->ObtainReadTempData("tau1tau2");
		const double *tau1tau2 = Sharedtau1tau2->ObtainReadData();
		const SharedSpace *Sharednu1 = etax->ObtainReadTempData("nu1");
		Vector *nu1 = Sharednu1->GetSharedElement();
		const SharedSpace *Sharednu2 = etax->ObtainReadTempData("nu2");
		Vector *nu2 = Sharednu2->GetSharedElement();
		double temp = 0;
		const double *nu1TV = nu1->ObtainReadData();
		const double *nu2TV = nu2->ObtainReadData();
		//double *nu1TV = new double[nu1->Getlength() * 2];
		//double *nu2TV = nu1TV + nu1->Getlength();
		//nu1->CopytoArray(nu1TV);
		//nu2->CopytoArray(nu2TV);
		HasHHR = false; TranH(x, etax, y, Hx, start, end, result); HasHHR = true;
		double *resultTV = result->ObtainWritePartialData();

		char *sidel = const_cast<char *> ("l");
		integer ell = Hx->Getsize()[0], length = etax->Getlength();
		double *work = new double[ell];
		// resultTV(start : start + length - 1, :) <- (I - tau1tau2(0) * nu1TV * nu1TV^T) * resultTV(start : start + length - 1, :),
		// details: http://www.netlib.org/lapack/explore-html/db/d10/dlarfx_8f.html
		dlarfx_(sidel, &length, &ell, const_cast<double *> (nu1TV), const_cast<double *> (tau1tau2), resultTV + start, &ell, work);
		// resultTV(start : start + length - 1, :) <- (I - tau1tau2(1) * nu2TV * nu2TV^T) * resultTV(start : start + length - 1, :),
		// details: http://www.netlib.org/lapack/explore-html/db/d10/dlarfx_8f.html
		dlarfx_(sidel, &length, &ell, const_cast<double *> (nu2TV), const_cast<double *> (tau1tau2 + 1), resultTV + start, &ell, work);
		delete[] work;
		//delete[] nu1TV;
	};

	void Manifold::LCTranHInvTran(Variable *x, Vector *etax, Variable *y, LinearOPE *Hx, LinearOPE *result) const
	{
		if (!etax->TempDataExist("nu1nu2"))
		{
			Obtainnu1nu2forLC(x, etax, y);
		}
		const SharedSpace *Sharedtau1tau2 = etax->ObtainReadTempData("tau1tau2");
		const double *tau1tau2 = Sharedtau1tau2->ObtainReadData();
		const SharedSpace *Sharednu1 = etax->ObtainReadTempData("nu1");
		Vector *nu1 = Sharednu1->GetSharedElement();
		const SharedSpace *Sharednu2 = etax->ObtainReadTempData("nu2");
		Vector *nu2 = Sharednu2->GetSharedElement();
		double temp = 0;
		const double *nu1TV = nu1->ObtainReadData();
		const double *nu2TV = nu2->ObtainReadData();
		HasHHR = false; TranHInvTran(x, etax, y, Hx, result); HasHHR = true;
		double *resultTV = result->ObtainWritePartialData();
		char *sidel = const_cast<char *> ("l"), *sider = const_cast<char *> ("r");
		integer ell = Hx->Getsize()[0], length = etax->Getlength();
		double *work = new double[ell];
		// resultTV <- resultTV * (I - tau1tau2(0) * nu1TV * nu1TV^T),
		// details: http://www.netlib.org/lapack/explore-html/db/d10/dlarfx_8f.html
		dlarfx_(sider, &ell, &length, const_cast<double *> (nu1TV), const_cast<double *> (tau1tau2), resultTV, &ell, work);
		// resultTV <- resultTV * (I - tau1tau2(1) * nu2TV * nu2TV^T),
		// details: http://www.netlib.org/lapack/explore-html/db/d10/dlarfx_8f.html
		dlarfx_(sider, &ell, &length, const_cast<double *> (nu2TV), const_cast<double *> (tau1tau2 + 1), resultTV, &ell, work);
		// resultTV <- (I - tau1tau2(0) * nu1TV * nu1TV^T) * resultTV,
		// details: http://www.netlib.org/lapack/explore-html/db/d10/dlarfx_8f.html
		dlarfx_(sidel, &length, &ell, const_cast<double *> (nu1TV), const_cast<double *> (tau1tau2), resultTV, &ell, work);
		// resultTV <- (I - tau1tau2(1) * nu2TV * nu2TV^T) * resultTV,
		// details: http://www.netlib.org/lapack/explore-html/db/d10/dlarfx_8f.html
		dlarfx_(sidel, &length, &ell, const_cast<double *> (nu2TV), const_cast<double *> (tau1tau2 + 1), resultTV, &ell, work);
		delete[] work;
		//delete[] nu1TV;
	};

	void Manifold::IntrProjection(Variable *x, Vector *v, Vector *result) const
	{
		v->CopyTo(result);
	};

	void Manifold::ExtrProjection(Variable *x, Vector *v, Vector *result) const
	{
		v->CopyTo(result);
	};

	void Manifold::CheckParams(void) const
	{
		OUTSTREAM << "GENERAL PARAMETERS:" << std::endl;
		OUTSTREAM << "name          :" << std::setw(15) << name << ",\t";
		OUTSTREAM << "IsIntrApproach:" << std::setw(15) << IsIntrApproach << std::endl;
		OUTSTREAM << "IntrinsicDim  :" << std::setw(15) << IntrinsicDim << ",\t";
		OUTSTREAM << "ExtrinsicDim  :" << std::setw(15) << ExtrinsicDim << std::endl;
		OUTSTREAM << "HasHHR        :" << std::setw(15) << HasHHR << ",\t";
		OUTSTREAM << "UpdBetaAlone  :" << std::setw(15) << UpdBetaAlone << std::endl;
		OUTSTREAM << "HasLockCon    :" << std::setw(15) << HasLockCon << std::endl;
	};

	void Manifold::CheckIntrExtr(Variable *x) const
	{
		OUTSTREAM << "==============Check Intrinsic/Extrinsic transform=========" << std::endl;
		Vector *exetax = EMPTYEXTR->ConstructEmpty();
		Vector *inetax = EMPTYINTR->ConstructEmpty();

		x->Print("x");
		exetax->RandGaussian();
		ExtrProjection(x, exetax, exetax);
		exetax->Print("exetax1");
		ObtainIntr(x, exetax, inetax);
		OUTSTREAM << "extr inp:" << Manifold::Metric(x, exetax, exetax) << ", intr inp:" << Manifold::Metric(x, inetax, inetax) << std::endl;
		inetax->Print("inetax1");
		ObtainExtr(x, inetax, exetax);
		exetax->Print("exetax2");
		ObtainIntr(x, exetax, inetax);
		inetax->Print("inetax2");
		OUTSTREAM << "exeta1 and inetax1 should approximately equal exetax2 and inetax2 respectively!" << std::endl;

		delete exetax;
		delete inetax;
	};

	void Manifold::CheckRetraction(Variable *x) const
	{
		OUTSTREAM << "==============Check Retraction=========" << std::endl;
		Vector *etax, *FDetax;
		etax = EMPTYEXTR->ConstructEmpty();
		FDetax = EMPTYEXTR->ConstructEmpty();
		etax->RandGaussian();
		ExtrProjection(x, etax, etax);
		//x->Print("x:");
		etax->Print("etax:");
		double eps = 1e-5;
		Variable *y = x->ConstructEmpty();
		ScaleTimesVector(x, eps, etax, etax);
		if (IsIntrApproach)
		{
			Vector *inetax = EMPTYINTR->ConstructEmpty();
			ObtainIntr(x, etax, inetax);
			Retraction(x, inetax, y);
			delete inetax;
		}
		else
		{
			Retraction(x, etax, y);
		}
		VectorMinusVector(x, y, x, FDetax);
		ScaleTimesVector(x, 1.0 / eps, FDetax, FDetax);
		FDetax->Print("FDetax:");

		OUTSTREAM << "etax should approximately equal FDetax = (R(eps etax)-R(etax))/eps!" << std::endl;
		delete etax;
		delete FDetax;
		delete y;
	};

	void Manifold::CheckDiffRetraction(Variable *x, bool IsEtaXiSameDir) const
	{
		OUTSTREAM << "==============Check Differentiated Retraction=========" << std::endl;
		Vector *etax, *xix, *zetax;
		etax = EMPTYEXTR->ConstructEmpty();
		xix = EMPTYEXTR->ConstructEmpty();
		zetax = EMPTYEXTR->ConstructEmpty();
		etax->RandGaussian();
		ExtrProjection(x, etax, etax);
		if (IsEtaXiSameDir)
		{
			etax->CopyTo(xix);
		}
		else
		{
			xix->RandGaussian();
			ExtrProjection(x, xix, xix);
		}
		Variable *y = x->ConstructEmpty();
		if (IsIntrApproach)
		{
			Vector *inetax = EMPTYINTR->ConstructEmpty();
			Vector *inxix = EMPTYINTR->ConstructEmpty();
			Vector *inzetax = EMPTYINTR->ConstructEmpty();
			ObtainIntr(x, etax, inetax);
			ObtainIntr(x, xix, inxix);
			Retraction(x, inetax, y);
			DiffRetraction(x, inetax, y, inxix, inzetax, IsEtaXiSameDir);
			ObtainExtr(y, inzetax, zetax);
			delete inetax;
			delete inxix;
			delete inzetax;
		}
		else
		{
			Retraction(x, etax, y);
			DiffRetraction(x, etax, y, xix, zetax, IsEtaXiSameDir);
		}
		y->Print("y:");
		zetax->Print("zetax:");
		Variable *yeps = x->ConstructEmpty();
		double eps = 1e-5;
		scalarVectorAddVector(x, eps, xix, etax, etax);
		if (IsIntrApproach)
		{
			Vector *inetax = EMPTYINTR->ConstructEmpty();
			ObtainIntr(x, etax, inetax);
			Retraction(x, inetax, yeps);
			delete inetax;
		}
		else
		{
			Retraction(x, etax, yeps);
		}
		VectorMinusVector(x, yeps, y, zetax);
		ScaleTimesVector(x, 1.0 / eps, zetax, zetax);
		ExtrProjection(y, zetax, zetax);
		zetax->Print("FDzetax:");
		OUTSTREAM << "zetax = T_{R_etax} xix should approximately equal FDzetax = (R(etax+eps xix) - R(etax))/eps!" << std::endl;

		delete etax;
		delete xix;
		delete zetax;
		delete yeps;
		delete y;
	};

	void Manifold::CheckLockingCondition(Variable *x) const
	{
		OUTSTREAM << "==============Check Locking Condition=========" << std::endl;
		Vector *etax, *xix, *zetax;
		etax = EMPTYEXTR->ConstructEmpty();
		xix = EMPTYEXTR->ConstructEmpty();
		zetax = EMPTYEXTR->ConstructEmpty();
		etax->RandGaussian();
		ExtrProjection(x, etax, etax);
		ScaleTimesVector(x, genrand_real2() + 0.5, etax, xix);
		Variable *y = x->ConstructEmpty();
		if (IsIntrApproach)
		{
			Vector *inetax = EMPTYINTR->ConstructEmpty();
			Vector *inxix = EMPTYINTR->ConstructEmpty();
			Vector *inzetax = EMPTYINTR->ConstructEmpty();
			ObtainIntr(x, etax, inetax);
			ObtainIntr(x, xix, inxix);
			Retraction(x, inetax, y);
			DiffRetraction(x, inetax, y, inxix, inzetax, true);
			if (inetax->TempDataExist("beta"))
			{
				const SharedSpace *beta = inetax->ObtainReadTempData("beta");
				const double *betav = beta->ObtainReadData();
				OUTSTREAM << "beta = |etax| / |T_{etax} etax|:" << betav[0] << std::endl;
			}
			else
			{
				OUTSTREAM << "beta:" << 1 << std::endl;
			}
			OUTSTREAM << "|xix| / |T_{etax} xix|:" << sqrt(Metric(x, inxix, inxix) / Metric(x, inzetax, inzetax)) << std::endl;
			ScaleTimesVector(x, sqrt(Metric(x, inxix, inxix) / Metric(x, inzetax, inzetax)),
				inzetax, inzetax);
			ObtainExtr(y, inzetax, zetax);
			zetax->Print("Beta DiffRetraction zetax:");
			VectorTransport(x, inetax, y, inxix, inzetax);
			ObtainExtr(y, inzetax, zetax);
			zetax->Print("Vector Transport zetax:");
			delete inetax;
			delete inxix;
			delete inzetax;
		}
		else
		{
			Retraction(x, etax, y);
			DiffRetraction(x, etax, y, xix, zetax, true);
			if (etax->TempDataExist("beta"))
			{
				const SharedSpace *beta = etax->ObtainReadTempData("beta");
				const double *betav = beta->ObtainReadData();
				OUTSTREAM << "beta = |etax| / |T_{etax} etax|:" << betav[0] << std::endl;
			}
			else
			{
				OUTSTREAM << "beta:" << 1 << std::endl;
			}
			OUTSTREAM << "|xix| / |T_{etax} xix|:" << sqrt(Metric(x, xix, xix) / Metric(x, zetax, zetax)) << std::endl;
			ScaleTimesVector(x, sqrt(Metric(x, xix, xix) / Metric(x, zetax, zetax)),
				zetax, zetax);
			zetax->Print("Beta DiffRetraction zetax:");
			VectorTransport(x, etax, y, xix, zetax);
			zetax->Print("Vector Transport zetax:");
		}
		OUTSTREAM << "Beta DiffRetraction zetax should approximately equal Vector Transport zetax!" << std::endl;

		delete etax;
		delete xix;
		delete zetax;
		delete y;
	};

	void Manifold::CheckcoTangentVector(Variable *x) const
	{
		OUTSTREAM << "==============Check CoTangentVector=========" << std::endl;
		Vector *etax, *xix, *zetay, *xiy, *zetax;
		etax = EMPTYEXTR->ConstructEmpty();
		xix = EMPTYEXTR->ConstructEmpty();
		zetay = EMPTYEXTR->ConstructEmpty();
		zetax = EMPTYEXTR->ConstructEmpty();
		xiy = EMPTYEXTR->ConstructEmpty();
		etax->RandGaussian();
		ExtrProjection(x, etax, etax);

		xix->RandGaussian();
		ExtrProjection(x, xix, xix);

		Variable *y = x->ConstructEmpty();
		if (IsIntrApproach)
		{
			Vector *inetax = EMPTYINTR->ConstructEmpty();
			Vector *inxix = EMPTYINTR->ConstructEmpty();
			Vector *inzetay = EMPTYINTR->ConstructEmpty();
			Vector *inxiy = EMPTYINTR->ConstructEmpty();
			Vector *inzetax = EMPTYINTR->ConstructEmpty();
			ObtainIntr(x, etax, inetax);
			ObtainIntr(x, xix, inxix);
			Retraction(x, inetax, y);
			DiffRetraction(x, inetax, y, inxix, inzetay, false);
			ObtainExtr(y, inzetay, zetay);

			xiy->RandGaussian();
			ExtrProjection(y, xiy, xiy);
			ObtainIntr(y, xiy, inxiy);
			OUTSTREAM << "<xiy, T_{R_{eta}} xix>:" << Metric(y, inxiy, inzetay) << std::endl;

			coTangentVector(x, inetax, y, inxiy, inzetax);
			ObtainExtr(x, inzetax, zetax);
			OUTSTREAM << "C(x, etax, xiy) [xix]:" << Metric(x, inzetax, inxix) << std::endl;
			delete inetax;
			delete inxix;
			delete inzetay;
			delete inxiy;
			delete inzetax;
		}
		else
		{
			Retraction(x, etax, y);
			DiffRetraction(x, etax, y, xix, zetay, false);
			xiy->RandGaussian();
			ExtrProjection(y, xiy, xiy);
			ScaleTimesVector(y, sqrt(Metric(y, xiy, xiy)), xiy, xiy);
			OUTSTREAM << "<xiy, T_{R_{eta}} xix>:" << Metric(y, xiy, zetay) << std::endl;
			coTangentVector(x, etax, y, xiy, zetax);
			OUTSTREAM << "C(x, etax, xiy) [xix]:" << Metric(x, zetax, xix) << std::endl;
		}
		OUTSTREAM << "<xiy, T_{R_{eta}} xix> should approximately equal C(x, etax, xiy) [xix]!" << std::endl;


		delete etax;
		delete xix;
		delete zetay;
		delete zetax;
		delete xiy;
		delete y;
	};

	void Manifold::CheckIsometryofVectorTransport(Variable *x) const
	{
		OUTSTREAM << "==============Check Isometry of the Vector Transport=========" << std::endl;
		Vector *etax, *xix, *zetay;
		etax = EMPTYEXTR->ConstructEmpty();
		xix = EMPTYEXTR->ConstructEmpty();
		zetay = EMPTYEXTR->ConstructEmpty();
		etax->RandGaussian();
		ExtrProjection(x, etax, etax);

		xix->RandGaussian();
		ExtrProjection(x, xix, xix);

		Variable *y = x->ConstructEmpty();
		if (IsIntrApproach)
		{
			Vector *inetax = EMPTYINTR->ConstructEmpty();
			Vector *inxix = EMPTYINTR->ConstructEmpty();
			Vector *inzetay = EMPTYINTR->ConstructEmpty();
			ObtainIntr(x, etax, inetax);
			ObtainIntr(x, xix, inxix);
			Retraction(x, inetax, y);
			VectorTransport(x, inetax, y, inxix, inzetay);
			OUTSTREAM << "Before vector transport:" << Metric(x, inxix, inxix)
				<< ", After vector transport:" << Metric(y, inzetay, inzetay) << std::endl;
			delete inetax;
			delete inxix;
			delete inzetay;
		}
		else
		{
			Retraction(x, etax, y);
			VectorTransport(x, etax, y, xix, zetay);
			OUTSTREAM << "Before vector transport:" << Metric(x, xix, xix)
				<< ", After vector transport:" << Metric(y, zetay, zetay) << std::endl;
		}
		OUTSTREAM << "|xix| (Before vector transport) should approximately equal |T_{R_etax} xix| (After vector transport)" << std::endl;

		delete etax;
		delete xix;
		delete zetay;
		delete y;
	};

	void Manifold::CheckIsometryofInvVectorTransport(Variable *x) const
	{
		OUTSTREAM << "==============Check Isometry of the Inverse Vector Transport=========" << std::endl;
		Vector *etax, *xix, *zetay;
		etax = EMPTYEXTR->ConstructEmpty();
		xix = EMPTYEXTR->ConstructEmpty();
		zetay = EMPTYEXTR->ConstructEmpty();

		etax->RandGaussian();
		ExtrProjection(x, etax, etax);

		Variable *y = x->ConstructEmpty();
		if (IsIntrApproach)
		{
			Vector *inetax = EMPTYINTR->ConstructEmpty();
			Vector *inxix = EMPTYINTR->ConstructEmpty();
			Vector *inzetay = EMPTYINTR->ConstructEmpty();
			ObtainIntr(x, etax, inetax);
			Retraction(x, inetax, y);
			zetay->RandGaussian();
			ExtrProjection(y, zetay, zetay);
			ScaleTimesVector(y, sqrt(Metric(y, zetay, zetay)), zetay, zetay);
			ObtainIntr(y, zetay, inzetay);

			InverseVectorTransport(x, inetax, y, inzetay, inxix);
			OUTSTREAM << "Before inverse vector transport:" << Metric(y, inzetay, inzetay)
				<< ", After inverse vector transport:" << Metric(x, inxix, inxix) << std::endl;
			delete inetax;
			delete inxix;
			delete inzetay;
		}
		else
		{
			Retraction(x, etax, y);
			zetay->RandGaussian();
			ExtrProjection(x, zetay, zetay);
			InverseVectorTransport(x, etax, y, zetay, xix);
			OUTSTREAM << "Before inverse vector transport:" << Metric(y, zetay, zetay)
				<< ", After inverse vector transport:" << Metric(x, xix, xix) << std::endl;
		}
		OUTSTREAM << "|zetay| (Before inverse vector transport) should approximately equal |T_{R_etax}^{-1} zetay| (After inverse vector transport)" << std::endl;

		delete etax;
		delete xix;
		delete zetay;
		delete y;
	};

	void Manifold::CheckVecTranComposeInverseVecTran(Variable *x) const
	{
		OUTSTREAM << "==============Check Vector Transport Compose Inverse Vector Transport=========" << std::endl;
		Vector *etax, *xix, *zetay;
		etax = EMPTYEXTR->ConstructEmpty();
		xix = EMPTYEXTR->ConstructEmpty();
		zetay = EMPTYEXTR->ConstructEmpty();

		etax->RandGaussian();
		ExtrProjection(x, etax, etax);
		xix->RandGaussian();
		ExtrProjection(x, xix, xix);

		Variable *y = x->ConstructEmpty();
		if (IsIntrApproach)
		{
			Vector *inetax = EMPTYINTR->ConstructEmpty();
			Vector *inxix = EMPTYINTR->ConstructEmpty();
			Vector *inzetay = EMPTYINTR->ConstructEmpty();
			ObtainIntr(x, etax, inetax);
			Retraction(x, inetax, y);
			ObtainIntr(x, xix, inxix);
			xix->Print("xix:");
			VectorTransport(x, inetax, y, inxix, inzetay);
			InverseVectorTransport(x, inetax, y, inzetay, inxix);
			ObtainExtr(x, inxix, xix);
			xix->Print("T^{-1} ciric T xix:");
			OUTSTREAM << "xix and T^{-1} ciric T xix should be similar!" << std::endl;
			delete inetax;
			delete inxix;
			delete inzetay;
		}
		else
		{
			Retraction(x, etax, y);
			xix->Print("xix:");
			VectorTransport(x, etax, y, xix, zetay);
			InverseVectorTransport(x, etax, y, zetay, xix);
			xix->Print("T^{-1} ciric T xix:");
			OUTSTREAM << "xix and T^{-1} ciric T xix should be similar!" << std::endl;
		}
		delete etax;
		delete xix;
		delete zetay;
		delete y;
	};

	void Manifold::CheckTranHInvTran(Variable *x) const
	{
		OUTSTREAM << "==============Check Transport of a Hessian approximation=========" << std::endl;
		Vector *etax;
		Variable *y;
		LinearOPE *Hx, *result;

		etax = EMPTYEXTR->ConstructEmpty();
		etax->RandGaussian();
		ExtrProjection(x, etax, etax);

		y = x->ConstructEmpty();
		if (IsIntrApproach)
		{
			Vector *inetax = EMPTYINTR->ConstructEmpty();
			ObtainIntr(x, etax, inetax);
			Retraction(x, inetax, y);
			Hx = new LinearOPE(EMPTYINTR->Getlength());
			Hx->ScaledIdOPE();
			Hx->Print("Hx before:");
			result = new LinearOPE(EMPTYINTR->Getlength());
			TranHInvTran(x, inetax, y, Hx, result);

			result->Print("Hx after:");
			delete inetax;
		}
		else
		{
			Hx = new LinearOPE(EMPTYEXTR->Getlength());
			Hx->ScaledIdOPE();
			Hx->Print("Hx before:");
			result = new LinearOPE(EMPTYEXTR->Getlength());
			Retraction(x, etax, y);
			Vector *zetay1 = EMPTYEXTR->ConstructEmpty();
			Vector *zetay2 = EMPTYEXTR->ConstructEmpty();
			zetay1->RandGaussian();
			ExtrProjection(y, zetay1, zetay1);
			TranHInvTran(x, etax, y, Hx, result);
			result->Print("Hx after:");
			zetay1->Print("zetay:");
			LinearOPEEta(y, result, zetay1, zetay2);
			zetay2->Print("Hx zetay:");
			delete zetay1;
			delete zetay2;
		}

		delete etax;
		delete y;
		delete Hx;
		delete result;
	};

	void Manifold::CheckHaddScaledRank1OPE(Variable *x) const
	{
		OUTSTREAM << "==============Check Rank one Update to a Hessian Approximation=========" << std::endl;
		LinearOPE *Hx, *result;
		double scalar = 1.0;
		Vector *etax, *xix;
		etax = EMPTYEXTR->ConstructEmpty();
		etax->RandGaussian();
		ExtrProjection(x, etax, etax);

		xix = EMPTYEXTR->ConstructEmpty();
		xix->RandGaussian();
		ExtrProjection(x, xix, xix);

		if (IsIntrApproach)
		{
			Vector *inetax = EMPTYINTR->ConstructEmpty();
			Vector *inxix = EMPTYINTR->ConstructEmpty();
			ObtainIntr(x, etax, inetax);
			ObtainIntr(x, xix, inxix);
			Hx = new LinearOPE(EMPTYINTR->Getlength());
			Hx->ScaledIdOPE();
			Hx->Print("Hx before:");
			result = new LinearOPE(EMPTYINTR->Getlength());
			HaddScaledRank1OPE(x, Hx, scalar, inetax, inxix, result);
			inetax->Print("etax:");
			inxix->Print("xix:");
			result->Print("Hx after:");
			delete inetax;
			delete inxix;
		}
		else
		{
			Hx = new LinearOPE(EMPTYEXTR->Getlength());
			Hx->ScaledIdOPE();
			Hx->Print("Hx before:");
			result = new LinearOPE(EMPTYEXTR->Getlength());
			HaddScaledRank1OPE(x, Hx, scalar, etax, xix, result);
			etax->Print("etax:");
			xix->Print("xix:");
			result->Print("Hx after:");
		}
		delete Hx;
		delete result;
		delete etax;
		delete xix;
	};

	void Manifold::SetParams(PARAMSMAP params)
	{
		PARAMSMAP::iterator iter;
		for (iter = params.begin(); iter != params.end(); iter++)
		{
			if (iter->first == static_cast<std::string> ("HasHHR"))
			{
				SetHasHHR(((static_cast<integer> (iter->second)) != 0));
			}
		}
	};
} /*end of ROPTLIB namespace*/
