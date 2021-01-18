
#include "LRBFGS.h"

/*Define the namespace*/
namespace ROPTLIB{

	LRBFGS::LRBFGS(const Problem *prob, const Variable *initialx, const Variable *insoln)
	{
		Initialization(prob, initialx, insoln);
	};

	void LRBFGS::SetProbX(const Problem *prob, const Variable *initialx, const Variable *insoln)
	{
		SolversLS::SetProbX(prob, initialx, insoln);
		const Vector *EMPTYETA;
		if (prob->GetDomain()->GetIsIntrinsic())
			EMPTYETA = prob->GetDomain()->GetEMPTYINTR();
		else
			EMPTYETA = prob->GetDomain()->GetEMPTYEXTR();
		s = EMPTYETA->ConstructEmpty();
		y = EMPTYETA->ConstructEmpty();
		Py = EMPTYETA->ConstructEmpty();
		prob->SetUseGrad(true);
		prob->SetUseHess(false);
	};

	void LRBFGS::SetDefaultParams()
	{
		SolversLS::SetDefaultParams();
		isconvex = false;
		BBratio = 1;
		Num_pre_BB = 0;
		nu = 1e-4;
		mu = 1;
		LengthSY = 4;
		S = nullptr;
		Y = nullptr;
		Currentlength = 0;
		beginidx = 0;
		RHO = nullptr;
		gamma = 1;
		InitSteptype = QUADINTMOD;
		SolverName.assign("LRBFGS");
	};

	LRBFGS::~LRBFGS(void)
	{
		delete s;
		delete y;
		delete Py;
		DeleteVectors(S, LengthSY);
		DeleteVectors(Y, LengthSY);
		if (RHO != nullptr)
			delete[] RHO;
	};

	void LRBFGS::Run(void)
	{
		DeleteVectors(S, LengthSY);
		NewVectors(S, LengthSY);
		DeleteVectors(Y, LengthSY);
		NewVectors(Y, LengthSY);
		if (RHO != nullptr)
			delete[] RHO;
		RHO = new double[LengthSY];
		SolversLS::Run();
	};

	void LRBFGS::CheckParams(void)
	{
		SolversLS::CheckParams();
		char YES[] = "YES";
		char NO[] = "NO";
		char *status;

		Rprintf("LRBFGS METHOD PARAMETERS:\n");
		status = (nu >= 0 && nu < 1) ? YES : NO;
		Rprintf("nu            :%15g[%s],\t", nu, status);
		status = (mu >= 0) ? YES : NO;
		Rprintf("mu            :%15g[%s]\n", mu, status);
		status = YES;
		Rprintf("isconvex      :%15d[%s],\t", isconvex, status);
		status = (LengthSY >= 0) ? YES : NO;
		Rprintf("LengthSY      :%15d[%s]\n", LengthSY, status);
	};

	void LRBFGS::GetSearchDir(void)
	{
		HvLRBFGS(gf1, eta1);
		Mani->ScaleTimesVector(x1, -1.0, eta1, eta1);
	};

	void LRBFGS::UpdateData(void)
	{
		UpdateDataLRBFGS();
	};

	void LRBFGS::PrintInfo(void)
	{
		Rprintf("\n\tbetay:%.3e,rho:%.3e,gamma:%.3e,inpss:%.3e,inpsy:%.3e,inpyy:%.3e,IsUpdateHessian:%d,", betay, rho, gamma, inpss, inpsy, inpyy, isupdated);
		Rprintf("\n");
	};
}; /*end of ROPTLIB namespace*/
