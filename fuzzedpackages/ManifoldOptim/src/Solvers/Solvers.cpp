
#include "Solvers.h"

/*Define the namespace*/
namespace ROPTLIB{

	void Solvers::OutPutResults(Variable *inx1, double &inf1, double &inngf0, double &inngf, integer &initer,
		integer &innf, integer &inng, integer &innR, integer &innV, integer &innVp, double &inComTime,
		double *intimeSeries, double *infunSeries, double *ingradSeries, integer &inlengthSeries)
	{
		inx1 = x1; inf1 = f1; inngf0 = ngf0; inngf = ngf; initer = iter;
		innf = nf; inng = ng; innR = nR; innV = nV; innVp = nVp; inComTime = ComTime;
		intimeSeries = timeSeries; infunSeries = funSeries; ingradSeries = gradSeries; inlengthSeries = lengthSeries;
	};

	void Solvers::PrintGenInfo(void)
	{
		if (nV == 0)
			Rprintf("i:%d,f:%.3e,df/f:%.3e,|gf|:%.3e,time:%.2e,nf:%d,ng:%d,nR:%d,", iter, f2,
			((f1 - f2) / f2), ngf, static_cast<double>(getTickCount() - starttime) / CLK_PS, nf, ng, nR);
		else
			Rprintf("i:%d,f:%.3e,df/f:%.3e,|gf|:%.3e,time:%.2e,nf:%d,ng:%d,nR:%d,nV(nVp):%d(%d),", iter, f2,
			((f1 - f2) / f2), ngf, static_cast<double>(getTickCount() - starttime) / CLK_PS, nf, ng, nR, nV, nVp);
	};

	void Solvers::PrintInfo(void)
	{
		Rprintf("\n");
	};

	bool Solvers::IsStopped(void)
	{
		if (static_cast<double>(getTickCount() - starttime) / CLK_PS > TimeBound)
			return true;

		if (StopPtr != nullptr)
		{
			if (Prob->GetDomain()->GetIsIntrinsic())
			{
				const double * gf2space = gf2->GetSpace();
				if (gf2space != nullptr)
				{
					Vector *exgf2 = Prob->GetDomain()->GetEMPTYEXTR()->ConstructEmpty();
					Prob->GetDomain()->ObtainExtr(x2, gf2, exgf2);
					bool flag = StopPtr(x2, exgf2, f2, ngf, ngf0);
					delete exgf2;
					return flag;
				}
				return false;
			}
			else
			{
				return StopPtr(x2, gf2, f2, ngf, ngf0);
			}
		}

		if (Stop_Criterion == FUN_REL)
		{
			//if (fabs((f1 - f2) / f1) < Tolerance)
			//	numfref++;
			//else
			//	numfref = 0;
			//if (numfref > 4)
			//	return true;
			//return false;
			return fabs((f1 - f2) / f1) < Tolerance;
		}
		else
			if (Stop_Criterion == GRAD_F)
				return ngf < Tolerance;
			else
				if (Stop_Criterion == GRAD_F_0)
					return (ngf / ngf0) < Tolerance;

		OUTSTREAM << "Error: Stopping Criterion is not specefic!" << std::endl;
		return true;
	};

	void Solvers::CheckParams(void)
	{
		//status = (Stop_Criterion >= 0 && Stop_Criterion <= (FUN_REL | GRAD_F | GRAD_F_0 | TIMESEC)) ? YES : NO;
		//OUTSTREAM << "Stop_Criterion:" << std::setw(15);
		//for (integer i = 0; i < STOPCRITLENGTH; i++)
		//{
		//	if (static_cast<integer> (pow(2, i)) & Stop_Criterion)
		//	{
		//		OUTSTREAM << "(" << STOPCRITnames[i] << ")";
		//	}
		//}
		//OUTSTREAM << "[" << status << "]" << std::endl;

		std::string STOPCRITnames[STOPCRITLENGTH] = { "FUN_REL", "GRAD_F", "GRAD_F_0" };
		std::string DEBUGnames[DEBUGLENGTH] = { "NOOUTPUT", "FINALRESULT", "ITERRESULT", "DETAILED" };
		char YES[] = "YES";
		char NO[] = "NO";
		char *status;
		OUTSTREAM << "GENERAL PARAMETERS:" << std::endl;
		status = (Stop_Criterion >= 0 && Stop_Criterion < STOPCRITLENGTH) ? YES : NO;
		OUTSTREAM << "Stop_Criterion:" << std::setw(15) << STOPCRITnames[Stop_Criterion] << "[" << status << "],\t";
		status = (Tolerance > 0) ? YES : NO;
		OUTSTREAM << "Tolerance     :" << std::setw(15) << Tolerance << "[" << status << "]" << std::endl;
		status = (Max_Iteration > 0 && Max_Iteration >= Min_Iteration) ? YES : NO;
		OUTSTREAM << "Max_Iteration :" << std::setw(15) << Max_Iteration << "[" << status << "]" << ",\t";
		status = (Min_Iteration >= 0 && Min_Iteration <= Max_Iteration) ? YES : NO;
		OUTSTREAM << "Min_Iteration :" << std::setw(15) << Min_Iteration << "[" << status << "]" << std::endl;
		status = (OutputGap > 0) ? YES : NO;
		OUTSTREAM << "OutputGap     :" << std::setw(15) << OutputGap << "[" << status << "]" << ",\t";
		status = (Debug >= 0 && Debug < DEBUGLENGTH) ? YES : NO;
		OUTSTREAM << "DEBUG         :" << std::setw(15) << DEBUGnames[Debug] << "[" << status << "]" << std::endl;
	};

	void Solvers::Run(void)
	{
		if (Debug >= ITERRESULT)
		{
			if (timeSeries != nullptr)
				delete[] timeSeries;
			timeSeries = new double[1 + Max_Iteration];
			if (funSeries == nullptr)
				delete[] funSeries;
			funSeries = new double[1 + Max_Iteration];
			if (gradSeries == nullptr)
				delete[] gradSeries;
			gradSeries = new double[1 + Max_Iteration];
		}
		if (Debug >= FINALRESULT)
			Rprintf("=========================%s=========================\n", SolverName.c_str());
	};

	void Solvers::Initialization(const Problem *prob, const Variable *initialx)
	{
		SetProbX(prob, initialx);
		SetDefaultParams();
	};

	void Solvers::SetProbX(const Problem *prob, const Variable *initialx)
	{
		const Vector *EMPTYETA;
		if (prob->GetDomain()->GetIsIntrinsic())
			EMPTYETA = prob->GetDomain()->GetEMPTYINTR();
		else
			EMPTYETA = prob->GetDomain()->GetEMPTYEXTR();

		Mani = prob->GetDomain();
		Prob = prob;
		x1 = initialx->ConstructEmpty();

		initialx->CopyTo(x1);
		x2 = initialx->ConstructEmpty();
		gf1 = EMPTYETA->ConstructEmpty();
		gf2 = EMPTYETA->ConstructEmpty();
	};

	void Solvers::SetDefaultParams()
	{
		nf = 0; ng = 0; nV = 0; nVp = 0; nR = 0; nH = 0; lengthSeries = 0;
		timeSeries = nullptr; funSeries = nullptr; gradSeries = nullptr;
		StopPtr = nullptr;

		Stop_Criterion = GRAD_F_0;
		TimeBound = 60 * 60 * 24 * 365;//one year;
		Tolerance = 1e-6;
		Max_Iteration = 500;
		Min_Iteration = 0;
		OutputGap = 1;
		Debug = ITERRESULT;

		//numfref = 0;
	};

	Solvers::~Solvers(void)
	{
		delete x1;
		delete x2;
		delete gf1;
		delete gf2;
		if (Debug >= ITERRESULT)
		{
			if (timeSeries != nullptr)
				delete[] timeSeries;
			if (funSeries != nullptr)
				delete[] funSeries;
			if (gradSeries != nullptr)
				delete[] gradSeries;
		}
	};

	void Solvers::SetParams(PARAMSMAP params)
	{
		PARAMSMAP::iterator iter;
		for (iter = params.begin(); iter != params.end(); iter++)
		{
			if (iter->first == static_cast<std::string> ("Stop_Criterion"))
			{
				Stop_Criterion = static_cast<StopCrit> (static_cast<integer> (iter->second));
			}
			else
			if (iter->first == static_cast<std::string> ("Tolerance"))
			{
				Tolerance = iter->second;
			}
			else
			if (iter->first == static_cast<std::string> ("TimeBound"))
			{
				TimeBound = iter->second;
			}
			else
			if (iter->first == static_cast<std::string> ("Max_Iteration"))
			{
				Max_Iteration = static_cast<integer> (iter->second);
			}
			else
			if (iter->first == static_cast<std::string> ("Min_Iteration"))
			{
				Min_Iteration = static_cast<integer> (iter->second);
			}
			else
			if (iter->first == static_cast<std::string> ("OutputGap"))
			{
				OutputGap = static_cast<integer> (iter->second);
			}
			else
			if (iter->first == static_cast<std::string> ("DEBUG"))
			{
				Debug = static_cast<DEBUGINFO> (static_cast<integer> (iter->second));
			}
		}
	};
} /*end of ROPTLIB namespace*/
