
#include "LRTRSR1.h"

/*Define the namespace*/
namespace ROPTLIB{

	LRTRSR1::LRTRSR1(const Problem *prob, const Variable *initialx)
	{
		Initialization(prob, initialx);
	};

	void LRTRSR1::SetProbX(const Problem *prob, const Variable *initialx)
	{
		SolversTR::SetProbX(prob, initialx);
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

	void LRTRSR1::SetDefaultParams()
	{
		SolversTR::SetDefaultParams();
		theta = 0.1;
		kappa = 0.1;
		isconvex = false;
		LengthSY = 4;
		S = nullptr;
		Y = nullptr;
		YMGS = nullptr;
		inpss = 0;
		inpsy = 0;
		inpyy = 0;
		Currentlength = 0;
		beginidx = 0;
		SS = nullptr;
		SY = nullptr;
		PMGQ = nullptr;
		P = nullptr;
		gamma = 1;
		ischangedSandY = true;
		SolverName.assign("LRTRSR1");
	};

	LRTRSR1::~LRTRSR1(void)
	{
		delete s;
		delete y;
		DeleteVectors(S, LengthSY);
		DeleteVectors(Y, LengthSY);
		DeleteVectors(YMGS, LengthSY);
		if (SS != nullptr)
			delete[] SS;
		if (SY != nullptr)
			delete[] SY;
		if (PMGQ != nullptr)
			delete[] PMGQ;
		if (P != nullptr)
			delete[] P;
	};

	void LRTRSR1::Run(void)
	{
		DeleteVectors(S, LengthSY);
		NewVectors(S, LengthSY);
		DeleteVectors(Y, LengthSY);
		NewVectors(Y, LengthSY);
		DeleteVectors(YMGS, LengthSY);
		NewVectors(YMGS, LengthSY);
		if (SS != nullptr)
			delete[] SS;
		SS = new double[LengthSY * LengthSY];
		if (SY != nullptr)
			delete[] SY;
		SY = new double[LengthSY * LengthSY];
		if (PMGQ != nullptr)
			delete[] PMGQ;
		PMGQ = new double[LengthSY * LengthSY];
		if (P != nullptr)
			delete[] P;
		P = new integer[LengthSY];
		SolversTR::Run();
	};

	void LRTRSR1::NewVectors(Vector ** &Vs, integer l)
	{
		Vs = new Vector *[l];
		for (integer i = 0; i < l; i++)
			Vs[i] = eta1->ConstructEmpty();
	};

	void LRTRSR1::DeleteVectors(Vector ** &Vs, integer l)
	{
		if (Vs != nullptr)
		{
			for (integer i = 0; i < l; i++)
				delete Vs[i];
			delete[] Vs;
		}
	};

	void LRTRSR1::CheckParams(void)
	{
		SolversTR::CheckParams();
		char YES[] = "YES";
		char NO[] = "NO";
		char *status;

		OUTSTREAM << "LRTRSR1 METHOD PARAMETERS:" << std::endl;
		status = YES;
		OUTSTREAM << "isconvex      :" << std::setw(15) << isconvex << "[" << status << "],\t";
		status = (LengthSY >= 0) ? YES : NO;
		OUTSTREAM << "LengthSY      :" << std::setw(15) << LengthSY << "[" << status << "]" << std::endl;
	};

	void LRTRSR1::HessianEta(Vector *Eta, Vector *result)
	{
		/* This function makes use of SS, SY and gamma to evaluate the action of Hessian approximation [HAG2014, (64)].
			[HAG2014]: W. Huang, P.-A. Absil, and K. A. Gallivan. A Riemannian symmetric rank-one trustregion method.
			Mathematical Programming, 150(2):179?16, February 2015.

			SS is the Q in (46), SY is the P in (46), PMGQ is the P - gamma Q in (46).
			*/
		integer idx;
		double *v = new double[Currentlength];

		/*if S and Y has been updated in function: UpdateData(void), then PMGQ is recomputed and
			the LU decomposition of PMGQ is also recomputed. The LU decomposition is used to evaluate
			the action of PMGQ^{-1}. */
		if (ischangedSandY)
		{
			for (integer i = 0; i < Currentlength; i++)
			{
				idx = (i + beginidx) % LengthSY;
				Mani->scalarVectorAddVector(x1, -gamma, S[idx], Y[idx], YMGS[i]);
			}
			for (integer i = 0; i < Currentlength; i++)
			{
				for (integer j = 0; j < Currentlength; j++)
				{
					PMGQ[i + j * Currentlength] = SY[i + j * LengthSY] - gamma * SS[i + j * LengthSY];
				}
			}
			if (Currentlength > 0)
			{
				// compute LU
				integer info, CurLen = Currentlength;
				// LU decomposion for PMGQ, PMGQ = P * L * U, L and U are stored in PMGQ, the permutation matrix is in P
				// details: http://www.netlib.org/lapack/explore-html/d3/d6a/dgetrf_8f.html
				dgetrf_(&CurLen, &CurLen, PMGQ, &CurLen, P, &info);
				ischangedSandY = false;
			}
		}

		for (integer i = 0; i < Currentlength; i++)
			v[i] = Mani->Metric(x1, YMGS[i], Eta);

		if (Currentlength > 0)
		{
			char *trans = const_cast<char *> ("n");
			integer info, one = 1, CurLen = Currentlength;
			// solve linear system: PMGQ * X = v using the LU decomposition results from dgetrf, then solution is stored in v.
			// details: http://www.netlib.org/lapack/explore-html/d6/d49/dgetrs_8f.html
			dgetrs_(trans, &CurLen, &one, PMGQ, &CurLen, P, v, &CurLen, &info);
		}

		Mani->ScaleTimesVector(x1, gamma, Eta, result);
		for (integer i = 0; i < Currentlength; i++)
		{
			Mani->scalarVectorAddVector(x1, v[i], YMGS[i], result, result);
		}

		delete[] v;
	};

	void LRTRSR1::UpdateData(void)
	{
		double denorminator, norm2ymBs;
		double mintolsq = std::numeric_limits<double>::epsilon();
		double mintol = sqrt(mintolsq);
		Prob->Grad(x2, gf2); ng++;
		eta2->CopyTo(s);
		Mani->InverseVectorTransport(x1, eta2, x2, gf2, eta1); nV++;
		Mani->VectorMinusVector(x1, eta1, gf1, y);
		Mani->VectorMinusVector(x1, y, zeta, zeta);
		denorminator = Mani->Metric(x1, s, zeta);
		inpss = Mani->Metric(x1, s, s);
		norm2ymBs = Mani->Metric(x1, zeta, zeta);
		if (iter == 0) // This is for the robustness when the cost function is quadratic 
		{			   // and its Hessian is identity everywhere.
			inpsy = Mani->Metric(x1, s, y);
			inpyy = Mani->Metric(x1, y, y);
			gamma = inpyy / inpsy;
		}
		if (denorminator * denorminator >= mintolsq * inpss * norm2ymBs && (norm2ymBs >= mintolsq || ngf / ngf0 < 1e-3)
			&& (iter != 0 || fabs(gamma - inpsy / inpss) > mintol)) // This is for the robustness when the cost 
			// function is quadratic and its Hessian is identity everywhere.
		{
			inpsy = Mani->Metric(x1, s, y);
			inpyy = Mani->Metric(x1, y, y);
			gamma = inpyy / inpsy;

			/*if s and y are accepted, then S and Y need to be updated. It follows that the matrices SY and SS need to be update.*/
			if (Currentlength < LengthSY)
			{
				s->CopyTo(S[Currentlength]);
				y->CopyTo(Y[Currentlength]);
				SS[Currentlength + Currentlength * LengthSY] = Mani->Metric(x1, S[Currentlength], S[Currentlength]);
				SY[Currentlength + Currentlength * LengthSY] = Mani->Metric(x1, S[Currentlength], Y[Currentlength]);
				for (integer i = 0; i < Currentlength; i++)
				{
					SS[Currentlength + i * LengthSY] = Mani->Metric(x1, S[Currentlength], S[i]);
					SS[i + Currentlength * LengthSY] = SS[Currentlength + i * LengthSY];
					SY[Currentlength + i * LengthSY] = Mani->Metric(x1, S[Currentlength], Y[i]);
					SY[i + Currentlength * LengthSY] = SY[Currentlength + i * LengthSY];
				}
				Currentlength++;
			}
			else
			{
				s->CopyTo(S[beginidx]);
				y->CopyTo(Y[beginidx]);
				for (integer i = 0; i < LengthSY - 1; i++)
				{
					for (integer j = 0; j < LengthSY - 1; j++)
					{
						SS[i + j * LengthSY] = SS[i + 1 + (j + 1) * LengthSY];
						SY[i + j * LengthSY] = SY[i + 1 + (j + 1) * LengthSY];
					}
				}
				SS[LengthSY * LengthSY - 1] = Mani->Metric(x1, S[beginidx], S[beginidx]);
				SY[LengthSY * LengthSY - 1] = Mani->Metric(x1, S[beginidx], Y[beginidx]);
				integer idx = 0;
				for (integer i = 0; i < LengthSY - 1; i++)
				{
					idx = (i + beginidx + 1) % LengthSY;
					SS[i + (LengthSY - 1) * LengthSY] = Mani->Metric(x1, S[idx], S[beginidx]);
					SS[LengthSY - 1 + i * LengthSY] = SS[i + (LengthSY - 1) * LengthSY];
					SY[i + (LengthSY - 1) * LengthSY] = Mani->Metric(x1, Y[idx], S[beginidx]);
					SY[LengthSY - 1 + i * LengthSY] = SY[i + (LengthSY - 1) * LengthSY];
				}
				++beginidx; beginidx = (beginidx) % LengthSY;
			}
			isupdated = true;
			ischangedSandY = true;
		}
		else
		{
			isupdated = false;
		}
	};

	void LRTRSR1::Acceptence(void)
	{
		for (integer i = 0; i < Currentlength; i++)
		{
			Mani->VectorTransport(x1, eta2, x2, S[i], S[i]);
			Mani->VectorTransport(x1, eta2, x2, Y[i], Y[i]);
		}
		ischangedSandY = true;
	};

	void LRTRSR1::PrintInfo(void)
	{
		Rprintf("\n\tgamma:%.3e,inpss:%.3e,inpsy:%.3e,inpyy:%.3e,IsUpdateHessian:%d,", gamma, inpss, inpsy, inpyy, isupdated);
		Rprintf("\n");
	};

	void LRTRSR1::SetParams(PARAMSMAP params)
	{
		SolversTR::SetParams(params);
		PARAMSMAP::iterator iter;
		for (iter = params.begin(); iter != params.end(); iter++)
		{
			if (iter->first == static_cast<std::string> ("isconvex"))
			{
				isconvex = ((static_cast<integer> (iter->second)) != 0);
			}
			else
				if (iter->first == static_cast<std::string> ("LengthSY"))
				{
					LengthSY = static_cast<integer> (iter->second);
				}
		}
	};
} /*end of ROPTLIB namespace*/
