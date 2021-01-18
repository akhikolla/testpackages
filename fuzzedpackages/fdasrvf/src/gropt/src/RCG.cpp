
#include "RCG.h"

/*Define the namespace*/
namespace ROPTLIB{

	RCG::RCG(const Problem *prob, const Variable *initialx, const Variable *insoln)
	{
		Initialization(prob, initialx, insoln);
	};

	void RCG::SetProbX(const Problem *prob, const Variable *initialx, const Variable *insoln)
	{
		SolversLS::SetProbX(prob, initialx, insoln);
		prob->SetUseGrad(true);
		prob->SetUseHess(false);
	};

	void RCG::SetDefaultParams()
	{
		SolversLS::SetDefaultParams();
		sigma = 0;
		RCGmethod = HESTENES_STIEFEL;
		ManDim = std::numeric_limits<integer>::max();
		InitSteptype = BBSTEP;

		SolverName.assign("RCG");

		RCGmethodSetnames = new std::string[RCGMETHODSLENGTH];
		RCGmethodSetnames[FLETCHER_REEVES].assign("FLETCHER_REEVES");
		RCGmethodSetnames[POLAK_RIBIERE_MOD].assign("POLAK_RIBIERE_MOD");
		RCGmethodSetnames[HESTENES_STIEFEL].assign("HESTENES_STIEFEL");
		RCGmethodSetnames[FR_PR].assign("FR_PR");
		RCGmethodSetnames[DAI_YUAN].assign("DAI_YUAN");
		RCGmethodSetnames[HAGER_ZHANG].assign("HAGER_ZHANG");
	};

	RCG::~RCG(void)
	{
		delete[] RCGmethodSetnames;
	};

	void RCG::CheckParams(void)
	{
		SolversLS::CheckParams();
		char YES[] = "YES";
		char NO[] = "NO";
		char *status;

		Rprintf("RCG METHOD PARAMETERS:\n");
		status = (ManDim >= 0 && ManDim <= std::numeric_limits<integer>::max()) ? YES : NO;
		Rprintf("ManDim        :%15d[%s],\t", ManDim, status);
		status = (RCGmethod >= 0 && RCGmethod <= RCGMETHODSLENGTH) ? YES : NO;
		Rprintf("RCGmethod     :%15s[%s]\n", RCGmethodSetnames[RCGmethod].c_str(), status);
	};

	void RCG::PrintInfo(void)
	{
		if (iter % ManDim == 0 || Mani->Metric(x1, eta1, gf1) >= -std::numeric_limits<double>::epsilon()) // restart and safeguard
			Rprintf("\n\tsigma:%.3e,Reset search direction to the negative gradient,", sigma);
		else
			Rprintf("\n\tsigma:%.3e,", sigma);
		Rprintf("\n");
	};

	void RCG::GetSearchDir(void)
	{
		PreConditioner(x1, gf1, Pgf1);
		if (iter % ManDim == 0 || Mani->Metric(x1, eta1, gf1) / ngf / ngf >= -std::sqrt(std::numeric_limits<double>::epsilon())) // restart and safeguard
		{
			Mani->ScaleTimesVector(x1, -1.0, Pgf1, eta1);
		}
	};

	void RCG::UpdateData(void)
	{
		if (iter % ManDim != 0)
		{
			PreConditioner(x2, gf2, Pgf2);
			if (RCGmethod == FLETCHER_REEVES)
			{
				sigma = Mani->Metric(x2, gf2, Pgf2) / Mani->Metric(x1, gf1, Pgf1);
				Mani->VectorTransport(x1, eta2, x2, eta1, zeta); nV++;
			}
			else
				if (RCGmethod == POLAK_RIBIERE_MOD)
				{
					Mani->VectorTransport(x1, eta2, x2, gf1, zeta); nV++;
					Mani->VectorMinusVector(x2, gf2, zeta, zeta);
					sigma = Mani->Metric(x2, zeta, Pgf2) / Mani->Metric(x1, gf1, Pgf1);
					if (LineSearch_LS == STRONGWOLFE && sigma <= 0)
						sigma = 0;
					else
					{
						Mani->VectorTransport(x1, eta2, x2, eta1, zeta); nVp++;
					}
				}
				else
					if (RCGmethod == HESTENES_STIEFEL)
					{
						double numerator, denominator;
						Mani->VectorTransport(x1, eta2, x2, eta1, zeta); nV++;
						Mani->VectorTransport(x1, eta2, x2, gf1, eta1); nVp++;
						Mani->VectorMinusVector(x2, gf2, eta1, eta1);
						numerator = Mani->Metric(x2, eta1, Pgf2);
						denominator = Mani->Metric(x2, zeta, eta1);
						sigma = numerator / denominator;
					}
					else
						if (RCGmethod == FR_PR)
						{
							double sigmaFR, sigmaPR;
							Mani->VectorTransport(x1, eta2, x2, gf1, zeta); nV++;
							Mani->VectorMinusVector(x2, gf2, zeta, zeta);
							sigmaPR = Mani->Metric(x2, zeta, Pgf2) / Mani->Metric(x1, gf1, Pgf1);
							sigmaFR = Mani->Metric(x2, gf2, Pgf2) / Mani->Metric(x1, gf1, Pgf1);
							if (sigmaPR < -sigmaFR)
								sigma = -sigmaFR;
							else
								if (sigmaPR > sigmaFR)
									sigma = sigmaFR;
								else
									sigma = sigmaPR;
							Mani->VectorTransport(x1, eta2, x2, eta1, zeta); nVp++;
						}
						else
							if (RCGmethod == DAI_YUAN)
							{
								Mani->VectorTransport(x1, eta2, x2, eta1, zeta); nV++;

								Mani->VectorTransport(x1, eta2, x2, gf1, eta1); nVp++;
								Mani->VectorMinusVector(x2, gf2, eta1, eta1);
								sigma = Mani->Metric(x2, gf2, Pgf2) / Mani->Metric(x2, zeta, eta1);
							}
							else
								if (RCGmethod == HAGER_ZHANG)
								{
									double temp1, temp2;
									Mani->VectorTransport(x1, eta2, x2, eta1, zeta); nV++;

									Mani->VectorTransport(x1, eta2, x2, gf1, eta1); nVp++;
									Mani->VectorMinusVector(x2, gf2, eta1, eta1);
									temp1 = Mani->Metric(x2, eta1, zeta);
									temp2 = -2.0 * Mani->Metric(x2, eta1, eta1) / temp1;
									Mani->scalarVectorAddVector(x2, temp2, zeta, eta1, eta1);
									sigma = Mani->Metric(x2, eta1, Pgf2) / temp1;
								}
			Mani->scalarVectorMinusVector(x2, sigma, zeta, gf2, eta1);
		}
	};

	void RCG::SetParams(PARAMSMAP params)
	{
		SolversLS::SetParams(params);
		PARAMSMAP::iterator iter;
		for (iter = params.begin(); iter != params.end(); iter++)
		{
			if (iter->first == static_cast<std::string> ("ManDim"))
			{
				ManDim = static_cast<integer> (iter->second);
			}
			else
				if (iter->first == static_cast<std::string> ("RCGmethod"))
				{
					RCGmethod = static_cast<RCGmethods> (static_cast<integer> (iter->second));
				}
		}
	};
}; /*end of ROPTLIB namespace*/
