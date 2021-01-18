
#include "RWRBFGS.h"

/*Define the namespace*/
namespace ROPTLIB{

	RWRBFGS::RWRBFGS(const Problem *prob, const Variable *initialx, LinearOPE *initialH)
	{
		Initialization(prob, initialx, initialH);
	};

	void RWRBFGS::Initialization(const Problem *prob, const Variable *initialx, LinearOPE *initialH)
	{
		SetProbX(prob, initialx, initialH);
		SetDefaultParams();
	};

	void RWRBFGS::SetProbX(const Problem *prob, const Variable *initialx, LinearOPE *initialH)
	{
		SolversLS::SetProbX(prob, initialx);
		const Vector *EMPTYETA;
		if (prob->GetDomain()->GetIsIntrinsic())
			EMPTYETA = prob->GetDomain()->GetEMPTYINTR();
		else
			EMPTYETA = prob->GetDomain()->GetEMPTYEXTR();
		bool initHisnull = (initialH == nullptr);
		if (initHisnull)
		{
			if (prob->GetDomain()->GetIsIntrinsic())
			{
				initialH = new LinearOPE(prob->GetDomain()->GetEMPTYINTR()->Getlength());
			}
			else
			{
				initialH = new LinearOPE(prob->GetDomain()->GetEMPTYEXTR()->Getlength());
			}
			initialH->ScaledIdOPE();
		}

		H = initialH->ConstructEmpty();
		tildeH = initialH->ConstructEmpty();
		initialH->CopyTo(H);
		s = EMPTYETA->ConstructEmpty();
		y = EMPTYETA->ConstructEmpty();

		prob->SetUseGrad(true);
		prob->SetUseHess(false);
		if (initHisnull)
			delete initialH;
	};

	void RWRBFGS::SetDefaultParams()
	{
		SolversLS::SetDefaultParams();
		isconvex = false;
		nu = 1e-4;
		mu = 1;
		InitSteptype = QUADINTMOD;
		SolverName.assign("RWRBFGS");
	};

	RWRBFGS::~RWRBFGS(void)
	{
		delete s;
		delete y;
		delete H;
		delete tildeH;
	};

	void RWRBFGS::CheckParams(void)
	{
		SolversLS::CheckParams();
		char YES[] = "YES";
		char NO[] = "NO";
		char *status;

		OUTSTREAM << "RWRBFGS METHOD PARAMETERS:" << std::endl;
		status = (nu >= 0 && nu < 1) ? YES : NO;
		OUTSTREAM << "nu            :" << std::setw(15) << nu << "[" << status << "],\t";
		status = (mu >= 0) ? YES : NO;
		OUTSTREAM << "mu            :" << std::setw(15) << mu << "[" << status << "]" << std::endl;
		status = YES;
		OUTSTREAM << "isconvex      :" << std::setw(15) << isconvex << "[" << status << "]" << std::endl;
	};

	void RWRBFGS::GetSearchDir(void)
	{
		Mani->LinearOPEEta(x1, H, gf1, eta1); nH++;
		Mani->ScaleTimesVector(x1, -1.0, eta1, eta1);
	};

	void RWRBFGS::UpdateData(void)
	{
		eta2->CopyTo(s);
		Mani->coTangentVector(x1, eta2, x2, gf2, y); nV++;
		Mani->VectorMinusVector(x1, y, gf1, y);
		inpsy = Mani->Metric(x1, s, y);
		if (isconvex && iter == 1 && inpsy > 0)
			H->ScaledIdOPE(inpsy / Mani->Metric(x1, y, y));
		inpss = Mani->Metric(x1, s, s);
		if (inpsy / inpss >= nu * pow(ngf, mu) && inpss > std::numeric_limits<double>::epsilon()
			&& inpsy > std::numeric_limits<double>::epsilon())
		{
			Mani->LinearOPEEta(x2, H, y, zeta); // zeta = H y
			Mani->HaddScaledRank1OPE(x2, H, -1.0 / inpsy, s, zeta, H);
			Mani->LinearOPEEta(x2, H, y, zeta); // zeta = H y

			Mani->HaddScaledRank1OPE(x2, H, -1.0 / inpsy, zeta, s, H);
			Mani->HaddScaledRank1OPE(x2, H, 1.0 / inpsy, s, s, H);
			Mani->TranHInvTran(x1, eta2, x2, H, tildeH);
			tildeH->CopyTo(H);
			isupdated = true;
		}
		else
		{
			isupdated = false;
			Mani->TranHInvTran(x1, eta2, x2, H, tildeH);
			tildeH->CopyTo(H);
		}
	};

	void RWRBFGS::PrintInfo(void)
	{
		Rprintf("\n\tinpss:%.3e,inpsy:%.3e,IsUpdateHessian:%d,", inpss, inpsy, isupdated);
		Rprintf("\n");
	};

	void RWRBFGS::SetParams(PARAMSMAP params)
	{
		SolversLS::SetParams(params);
		PARAMSMAP::iterator iter;
		for (iter = params.begin(); iter != params.end(); iter++)
		{
			if (iter->first == static_cast<std::string> ("isconvex"))
			{
				isconvex = ((static_cast<integer> (iter->second)) != 0);
			}
			else
				if (iter->first == static_cast<std::string> ("nu"))
				{
					nu = iter->second;
				}
				else
					if (iter->first == static_cast<std::string> ("mu"))
					{
						mu = iter->second;
					}
		}
	};
} /*end of ROPTLIB namespace*/
