#ifndef _REPLACE_H

#define _REPLACE_H 1

#include "individual.h"
#include "populationreal.h"
#include "signal.h"

namespace realea {

class IReplace : public IReset {
    public:
	IReplace(void);
	virtual ~IReplace(void) {}
	virtual void reset(void);
	virtual unsigned getCandidate(PopulationRealPtr pop, tIndividualRealPtr newind)=0;
	virtual bool mustBeReplace(tIndividualRealPtr old, tIndividualRealPtr newind);

    private:
	unsigned m_total, m_success;
};

/**
 * @class ReplaceWorst
 *
 * Replace the worst individual in the population
 */
class ReplaceWorst : public IReplace {
    public:
	virtual unsigned getCandidate(PopulationRealPtr pop, tIndividualRealPtr newind);
}; 

/**
 * @class ReplaceWorst
 *
 * Replace the more similar individual in the population
 */
class ReplaceDC: public IReplace {
    public:
	virtual unsigned getCandidate(PopulationRealPtr pop, tIndividualRealPtr newind);
}; 


}

#endif
