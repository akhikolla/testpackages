
#include "LRBFGS.h"

/*Define the namespace*/
namespace ROPTLIB{

	LRBFGS::LRBFGS(const Problem *prob, const Variable *initialx)
	{
		Initialization(prob, initialx);
	};

	void LRBFGS::SetProbX(const Problem *prob, const Variable *initialx)
	{
		SolversLS::SetProbX(prob, initialx);
		const Vector *EMPTYETA;
		if (prob->GetDomain()->GetIsIntrinsic())
			EMPTYETA = prob->GetDomain()->GetEMPTYINTR();
		else
			EMPTYETA = prob->GetDomain()->GetEMPTYEXTR();
		s = EMPTYETA->ConstructEmpty();
		y = EMPTYETA->ConstructEmpty();
		prob->SetUseGrad(true);
		prob->SetUseHess(false);
	};

	void LRBFGS::SetDefaultParams()
	{
		SolversLS::SetDefaultParams();
		isconvex = false;
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

	void LRBFGS::NewVectors(Vector ** &Vs, integer l)
	{
		Vs = new Vector *[l];
		for (integer i = 0; i < l; i++)
			Vs[i] = eta1->ConstructEmpty();
	};

	void LRBFGS::DeleteVectors(Vector ** &Vs, integer l)
	{
		if (Vs != nullptr)
		{
			for (integer i = 0; i < l; i++)
				delete Vs[i];
			delete[] Vs;
		}
	};

	void LRBFGS::CheckParams(void)
	{
		SolversLS::CheckParams();
		char YES[] = "YES";
		char NO[] = "NO";
		char *status;

		OUTSTREAM << "LRBFGS METHOD PARAMETERS:" << std::endl;
		status = (nu >= 0 && nu < 1) ? YES : NO;
		OUTSTREAM << "nu            :" << std::setw(15) << nu << "[" << status << "],\t";
		status = (mu >= 0) ? YES : NO;
		OUTSTREAM << "mu            :" << std::setw(15) << mu << "[" << status << "]" << std::endl;
		status = YES;
		OUTSTREAM << "isconvex      :" << std::setw(15) << isconvex << "[" << status << "],\t";
		status = (LengthSY >= 0) ? YES : NO;
		OUTSTREAM << "LengthSY      :" << std::setw(15) << LengthSY << "[" << status << "]" << std::endl;
	};

	void LRBFGS::GetSearchDir(void)
	{
		double *xi = new double[Currentlength];
		double omega;
		integer idx;
		gf1->CopyTo(eta1);
		for (integer i = Currentlength - 1; i >= 0; i--)
		{
			idx = (beginidx + i) % LengthSY;
			xi[idx] = RHO[idx] * Mani->Metric(x1, S[idx], eta1);
			Mani->scalarVectorAddVector(x1, -xi[idx], Y[idx], eta1, eta1);
		}
		Mani->ScaleTimesVector(x1, gamma, eta1, eta1);
		for (integer i = 0; i < Currentlength; i++)
		{
			idx = (beginidx + i) % LengthSY;
			omega = RHO[idx] * Mani->Metric(x1, Y[idx], eta1);
			Mani->scalarVectorAddVector(x1, xi[idx] - omega, S[idx], eta1, eta1);
		}
		Mani->ScaleTimesVector(x1, -1.0, eta1, eta1);

		delete[] xi;
	};

	void LRBFGS::UpdateData(void)
	{
		Mani->VectorTransport(x1, eta2, x2, eta2, s); nV++;
		Mani->VectorTransport(x1, eta2, x2, gf1, zeta); nVp++;
		betay = Mani->Beta(x1, eta2);
		Mani->scalarVectorMinusVector(x2, 1.0 / betay, gf2, zeta, y);
		inpsy = Mani->Metric(x2, s, y);
		inpss = Mani->Metric(x2, s, s);
		inpyy = Mani->Metric(x2, y, y);
		rho = 1.0 / inpsy;
		if (inpsy / inpss >= nu * pow(ngf, mu) && (ngf / ngf0 < 1e-3 ||
			(inpss > std::numeric_limits<double>::epsilon() && inpsy > std::numeric_limits<double>::epsilon())))
		{
			gamma = inpsy / inpyy; /*Suggested in NW2006*/
			//gamma = inpss / inpsy; /*BB stepsize*/
			if (Currentlength < LengthSY)
			{
				y->CopyTo(Y[Currentlength]);
				s->CopyTo(S[Currentlength]);
				RHO[Currentlength] = rho;
				for (integer i = 0; i < Currentlength; i++)
				{
					Mani->VectorTransport(x1, eta2, x2, Y[i], Y[i]); nVp++;
					Mani->VectorTransport(x1, eta2, x2, S[i], S[i]); nVp++;
				}
				Currentlength++;
			}
			else
				if (LengthSY > 0)
				{
					integer idx;
					y->CopyTo(Y[beginidx]);
					s->CopyTo(S[beginidx]);
					RHO[beginidx] = rho;
					++beginidx; beginidx = (beginidx) % LengthSY;
					for (integer i = beginidx; i < beginidx + LengthSY - 1; i++)
					{
						idx = i % LengthSY;
						Mani->VectorTransport(x1, eta2, x2, Y[idx], Y[idx]); nVp++;
						Mani->VectorTransport(x1, eta2, x2, S[idx], S[idx]); nVp++;
					}
				}
			isupdated = true;
		}
		else
		{
			for (integer i = 0; i < Currentlength; i++)
			{
				Mani->VectorTransport(x1, eta2, x2, Y[i], Y[i]); nVp++;
				Mani->VectorTransport(x1, eta2, x2, S[i], S[i]); nVp++;
			}
			isupdated = false;
		}
	};

	void LRBFGS::PrintInfo(void)
	{
		Rprintf("\n\tbetay:%.3e,rho:%.3e,gamma:%.3e,inpss:%.3e,inpsy:%.3e,IsUpdateHessian:%d,", betay, rho, gamma, inpss, inpsy, isupdated);
		Rprintf("\n");
	};

	void LRBFGS::SetParams(PARAMSMAP params)
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
					else
						if (iter->first == static_cast<std::string> ("LengthSY"))
						{
							LengthSY = static_cast<integer> (iter->second);
						}
		}
	};
} /*end of ROPTLIB namespace*/
