
#ifndef MRANKADAPTIVE_H
#define MRANKADAPTIVE_H

#include "SolversLS.h"
#include "def.h"

/*Define the namespace*/
namespace ROPTLIB{

	class MRankAdaptive : public SolversLS{
	public:
		MRankAdaptive(const Problem *prob, const Variable *initialx, Solvers *inSolver);
		virtual void SetProbX(const Problem *prob, const Variable *initialx);
		virtual void SetDefaultParams();
		virtual void Run();

		// parameters
		double Reps1, Reps2, Reps3;
		double Rca, Rcr;
		double Rtau2;
		Solvers *innerSolver;
	protected:
		virtual void GetSearchDir();
		virtual void InitialStepSize();
	};

	namespace MRA{
		double RDelta;
		bool InnerStop(Variable *x, Vector *gf, double f, double ngf, double ngf0);
	};
} /*end of ROPTLIB namespace*/
#endif // end of RSD_H
