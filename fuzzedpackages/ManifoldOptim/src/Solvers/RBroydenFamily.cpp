
#include "RBroydenFamily.h"

/*Define the namespace*/
namespace ROPTLIB{

	RBroydenFamily::RBroydenFamily(const Problem *prob, const Variable *initialx, LinearOPE *initialH)
	{
		Initialization(prob, initialx, initialH);
	};

	void RBroydenFamily::Initialization(const Problem *prob, const Variable *initialx, LinearOPE *initialH)
	{
		SetProbX(prob, initialx, initialH);
		SetDefaultParams();
	};

	void RBroydenFamily::SetProbX(const Problem *prob, const Variable *initialx, LinearOPE *initialH)
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
		u = EMPTYETA->ConstructEmpty();
		if (initHisnull)
			delete initialH;
		prob->SetUseGrad(true);
		prob->SetUseHess(false);
	};

	void RBroydenFamily::SetDefaultParams()
	{
		SolversLS::SetDefaultParams();
		isconvex = false;
		nu = 1e-4;
		mu = 1;
		InitSteptype = QUADINTMOD;
		SolverName.assign("RBroydenFamily");
	};

	RBroydenFamily::~RBroydenFamily(void)
	{
		delete s;
		delete y;
		delete u;
		delete H;
		delete tildeH;
	};

	double RBroydenFamily::Phi(Variable *x2, Vector *y, Vector *s, LinearOPE *tildeH, double inpsy, double yHy, Vector *u)
	{
		return 1;
	};

	void RBroydenFamily::CheckParams(void)
	{
		SolversLS::CheckParams();
		char YES[] = "YES";
		char NO[] = "NO";
		char *status;

		OUTSTREAM << "RBROYDENFAMILY METHOD PARAMETERS:" << std::endl;
		status = (nu >= 0 && nu < 1) ? YES : NO;
		OUTSTREAM << "nu            :" << std::setw(15) << nu << "[" << status << "],\t";
		status = (mu >= 0) ? YES : NO;
		OUTSTREAM << "mu            :" << std::setw(15) << mu << "[" << status << "]" << std::endl;
		status = YES;
		OUTSTREAM << "isconvex      :" << std::setw(15) << isconvex << "[" << status << "]" << std::endl;
	};

	void RBroydenFamily::GetSearchDir(void)
	{
		Mani->LinearOPEEta(x1, H, gf1, eta1);
		Mani->ScaleTimesVector(x1, -1.0, eta1, eta1);
	};

	void RBroydenFamily::UpdateData(void)
	{
		double yHy;
		Mani->VectorTransport(x1, eta2, x2, eta2, s); nV++;
		Mani->VectorTransport(x1, eta2, x2, gf1, zeta); nVp++;
		betay = Mani->Beta(x1, eta2);
		//Mani->VectorLinearCombination(x2, 1.0 / betay, gf2, -1.0, zeta, y);
		Mani->scalarVectorMinusVector(x2, 1.0 / betay, gf2, zeta, y);
		inpsy = Mani->Metric(x2, s, y);
		if (isconvex && iter == 1 && inpsy > 0)
			H->ScaledIdOPE(inpsy / Mani->Metric(x2, y, y));
		Mani->TranHInvTran(x1, eta2, x2, H, tildeH);
		inpss = Mani->Metric(x2, s, s);
		if (inpsy / inpss >= nu * pow(ngf, mu) && (ngf / ngf0 < 1e-3 ||
			(inpss > std::numeric_limits<double>::epsilon() && inpsy > std::numeric_limits<double>::epsilon())))
		{
			Mani->LinearOPEEta(x2, tildeH, y, zeta); // zeta = tildeH y
			yHy = Mani->Metric(x2, y, zeta);
			Mani->VectorLinearCombination(x2, 1.0 / inpsy, s, -1.0 / yHy, zeta, u);
			phic = Phi(x2, y, s, tildeH, inpsy, yHy, u);
			Mani->HaddScaledRank1OPE(x2, tildeH, -1.0 / yHy, zeta, zeta, H);
			Mani->HaddScaledRank1OPE(x2, H, 1.0 / inpsy, s, s, H);
			Mani->HaddScaledRank1OPE(x2, H, phic * yHy, u, u, H);
			isupdated = true;
		}
		else
		{
			isupdated = false;
			tildeH->CopyTo(H);
		}
	};

	void RBroydenFamily::PrintInfo(void)
	{
		Rprintf("\n\tbetay:%.3e,Phic:%.3e,inpss:%.3e,inpsy:%.3e,IsUpdateHessian:%d,", betay, phic, inpss, inpsy, isupdated);
		Rprintf("\n");
	};

	void RBroydenFamily::SetParams(PARAMSMAP params)
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
