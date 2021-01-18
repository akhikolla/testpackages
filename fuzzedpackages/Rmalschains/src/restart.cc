#include "restart.h"

using namespace realea;

void RestartBest::apply(PopulationReal *pop_alg, Problem *problem, IEvalInd *initeval) {
    pop_alg->reset(problem->getDomain(),pop_alg->getBest());
    pop_alg->eval(initeval);
}

RestartReduce::RestartReduce(double scale) : RestartBest() {
    m_scale = scale;
}

void RestartReduce::apply(PopulationReal *pop_alg, Problem *problem, IEvalInd *initeval) {
    DomainRealPtr domain = problem->getDomain();
    tIndividualReal *best = pop_alg->getInd(pop_alg->getBest());
    domain->setDomainCenter(best->sol(), m_scale);

    RestartBest::apply(pop_alg, problem, initeval);
}

