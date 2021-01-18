
#include "RTRSR1.h"

/*Define the namespace*/
namespace ROPTLIB{

	RTRSR1::RTRSR1(const Problem *prob, const Variable *initialx, LinearOPE *initialB)
	{
		Initialization(prob, initialx, initialB);
	};

	void RTRSR1::Initialization(const Problem *prob, const Variable *initialx, LinearOPE *initialB)
	{
		SetProbX(prob, initialx, initialB);
		SetDefaultParams();
	};

	void RTRSR1::SetProbX(const Problem *prob, const Variable *initialx, LinearOPE *initialB)
	{
		SolversTR::SetProbX(prob, initialx);
		const Vector *EMPTYETA;
		if (prob->GetDomain()->GetIsIntrinsic())
			EMPTYETA = prob->GetDomain()->GetEMPTYINTR();
		else
			EMPTYETA = prob->GetDomain()->GetEMPTYEXTR();
		bool initBisnull = (initialB == nullptr);
		if (initBisnull)
		{
			if (prob->GetDomain()->GetIsIntrinsic())
			{
				initialB = new LinearOPE(prob->GetDomain()->GetEMPTYINTR()->Getlength());
			}
			else
			{
				initialB = new LinearOPE(prob->GetDomain()->GetEMPTYEXTR()->Getlength());
			}
			initialB->ScaledIdOPE();
		}
		B = initialB->ConstructEmpty();
		tildeB = initialB->ConstructEmpty();
		initialB->CopyTo(B);
		s = EMPTYETA->ConstructEmpty();
		y = EMPTYETA->ConstructEmpty();
		prob->SetUseGrad(true);
		prob->SetUseHess(false);

		if (initBisnull)
			delete initialB;
	};

	void RTRSR1::SetDefaultParams()
	{
		SolversTR::SetDefaultParams();
		theta = 0.1;
		kappa = 0.1;
		SolverName.assign("RTRSR1");
		isconvex = false;
	};

	RTRSR1::~RTRSR1(void)
	{
		delete s;
		delete y;
		delete B;
		delete tildeB;
	};

	void RTRSR1::HessianEta(Vector *Eta, Vector *result)
	{
		Mani->LinearOPEEta(x1, B, Eta, result);
	};

	void RTRSR1::UpdateData(void)
	{
		double denorminator, norm2ymBs;
		double mintolsq = std::numeric_limits<double>::epsilon();
		Prob->Grad(x2, gf2); ng++;
		eta2->CopyTo(s);
		Mani->InverseVectorTransport(x1, eta2, x2, gf2, eta1); nV++;
		Mani->VectorMinusVector(x1, eta1, gf1, y);
		Mani->VectorMinusVector(x1, y, zeta, zeta);
		denorminator = Mani->Metric(x1, s, zeta);
		if (isconvex && iter == 1 && Mani->Metric(x1, s, y) > 0)
			B->ScaledIdOPE(Mani->Metric(x1, y, y) / Mani->Metric(x1, s, y));
		inpss = Mani->Metric(x1, s, s);
		norm2ymBs = Mani->Metric(x1, zeta, zeta);
		if (denorminator * denorminator >= mintolsq * inpss * norm2ymBs && (norm2ymBs >= mintolsq || ngf / ngf0 < 1e-3))
		{
			Mani->HaddScaledRank1OPE(x1, B, 1.0 / denorminator, zeta, zeta, B);
			isupdated = true;
		}
		else
		{
			isupdated = false;
		}
	};

	void RTRSR1::Acceptence(void)
	{
		Mani->TranHInvTran(x1, eta2, x2, B, tildeB);
		tildeB->CopyTo(B);
	};

	void RTRSR1::CheckParams(void)
	{
		SolversTR::CheckParams();
		char YES[] = "YES";
		char NO[] = "NO";
		char *status;

		OUTSTREAM << "RTRSR1 METHOD PARAMETERS:" << std::endl;
		status = YES;
		OUTSTREAM << "isconvex      :" << std::setw(15) << isconvex << "[" << status << "]" << std::endl;
	};

	void RTRSR1::PrintInfo(void)
	{
		Rprintf("\n\tinpss:%.3e,IsUpdateHessian:%d,", inpss, isupdated);
		Rprintf("\n");
	};

	void RTRSR1::SetParams(PARAMSMAP params)
	{
		SolversTR::SetParams(params);
		PARAMSMAP::iterator iter;
		for (iter = params.begin(); iter != params.end(); iter++)
		{
			if (iter->first == static_cast<std::string> ("isconvex"))
			{
				isconvex = ((static_cast<integer> (iter->second)) != 0);
			}
		}
	};
} /*end of ROPTLIB namespace*/
