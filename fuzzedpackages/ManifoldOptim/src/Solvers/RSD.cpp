
#include "RSD.h"

/*Define the namespace*/
namespace ROPTLIB{

	RSD::RSD(const Problem *prob, const Variable *initialx)
	{
		Initialization(prob, initialx);
	};

	void RSD::SetProbX(const Problem *prob, const Variable *initialx)
	{
		SolversLS::SetProbX(prob, initialx);
		prob->SetUseGrad(true);
		prob->SetUseHess(false);
	};

	void RSD::SetDefaultParams()
	{
		SolversLS::SetDefaultParams();
		InitSteptype = BBSTEP;
		SolverName.assign("RSD");
	};

	void RSD::GetSearchDir(void)
	{
		Mani->ScaleTimesVector(x1, -1, gf1, eta1);
	};
} /*end of ROPTLIB namespace*/
