
#include "RTRNewton.h"

/*Define the namespace*/
namespace ROPTLIB{

	RTRNewton::RTRNewton(const Problem *prob, const Variable *initialx)
	{
		Initialization(prob, initialx);
	};

	void RTRNewton::SetProbX(const Problem *prob, const Variable *initialx)
	{
		SolversTR::SetProbX(prob, initialx);
		prob->SetUseGrad(true);
		prob->SetUseHess(true);
	};

	void RTRNewton::SetDefaultParams()
	{
		SolversTR::SetDefaultParams();
		SolverName.assign("RTRNewton");
	};

	void RTRNewton::HessianEta(Vector *Eta, Vector *result)
	{
		Prob->HessianEta(x1, Eta, result);
	};
} /*end of ROPTLIB namespace*/
