#ifndef _SSGA_H

#define _SSGA_H 1

#include "iea.h"
#include "populationreal.h"
#include "cross.h"
#include "mutation.h"
#include "select.h"
#include "replace.h"

namespace realea {

/**
 * @class SSGA 
 *
 * @brief Implements the Steady-by-State Genetic Algorithm (SSGA).
 *
 */
class SSGA : public ICrossEAlgorithm {
public:
    /**
     * Constructor.
     *
     * @param random random number generator
     */
    SSGA(Random *random);
    ~SSGA(void);
    unsigned init(void);
    unsigned realApply(tChromosomeReal &sol, tFitness &fitness);

    unsigned getDefaultPopsize(void) {
	return 60;
     }


    /**
     * Set the parents selection method
     *
     * @param sel criterion.
     */
    void setSelect(ISelect *sel);
    /**
     * Set the Replacement Strategy
     * @param replace Replacement Strategy
     */
    void setReplacement(IReplace *replace);
    /**
     * Set the mutation criterion
     *
     * @param mutation mutation method
     */
    void setMutation(IMutation *mutation);
    /**
     * Set the mutation rate
     * @param prob_mutation probability of mutation of new individuals (between 0 and 1)
     */
    void setMutationRate(double prob_mutation);

    void setProblem(Problem *problem);

private:
    /** 
     * Crossover operator 
     *
     * @param mom position of the 'mother'
     * @param dad position of the 'father' individual
     *
     * @param crom Chromosome, output  
     *
     */
    void cross(unsigned mom, unsigned dad, tChromosomeReal &crom);

    unsigned m_initEval;
    int m_pmut;

    ISelect *m_select;
    IReplace *m_replace;
    Mutation *m_mutation;
    IMutation *m_imutation;
};

}

#endif
