/**
 * Copyright 2008, Daniel Molina Cabrera <danimolina@gmail.com>
 * 
 * This file is part of software Realea
 * 
 * Realea is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * Realea is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "de.h"
#include "random.h"
#include <cassert>

#define CopyVector(a,b) memcpy((a),(b),nDim*sizeof(double))
#define Element(a,b,c)  (a->getInd(b).gen(c))

using namespace realea;

DE::DE(Random *random) : ClassEAlgorithm(random), m_CR(-1), m_F(-1) {
    m_crossover = 'e';
}

void DE::setF(double F) {
   assert(F >= 0 && F <= 1); 
   m_F = F;
}

void DE::setCR(double CR) {
    assert(CR >= 0 && CR <= 1);
    m_CR = CR;
}

void DE::setCrossoverBin(void) {
    m_crossover = 'b';
}

void DE::setCrossoverExp(void) {
    m_crossover = 'e';
}


DE::~DE(void) {
}

unsigned DE::init(void) {
    // Comprueba los elementos
    if (m_F == -1) {
	throw new ConfigException("DE::F");
    }

    if (m_CR == -1) {
	throw new ConfigException("DE::CR");
    }

    // Primero inicio la población
    m_pop->reset(m_problem->getDomain());
    // Inicio los distintos elementos: Running, Cross ...
    reset();
    // Evaluo la población
    m_pop->eval(m_init_eval);
    return m_running->numEval();
}

unsigned DE::realApply(tChromosomeReal &sol, tFitness &fitness) {
    tChromosomeReal crom(m_problem->getDimension());
    tFitness best_fit;
    unsigned better;

    unsigned popsize = m_pop->size();
    m_running->reset();

   // While is not finish
   while (!m_running->isFinish()) {
        if (m_stat)
	    m_stat->newGeneration();

	// Selecciono los individuos a cruzar
	for (unsigned ind = 0; ind < popsize && !m_running->isFinish(); ++ind) {
	   // Aplico el cruce sobre el individuo
	    cross(m_pop, ind, crom);
	    // Genero el nuevo individuo
	    tIndividualRealPtr newind = m_pop->getInstance(crom);
	    // Lo evalúo
	    m_new_eval->eval(newind);

	    if (newind->isBetter(m_pop->getInd(ind))) 
		m_pop->replace(ind, newind);
	    else 
		delete newind;

	} // De modificar cada individuo 

	better = m_pop->getBest();
	best_fit = m_pop->getInd(better)->perf();

	if (m_stat)
	   m_stat->endGeneration(best_fit);

   } // De running

   // Obtengo el mejor
   unsigned pos = m_pop->getBest();
   tIndividualRealPtr best= m_pop->getInd(pos);

   tChromosomeReal bestsol= best->sol();
   copy(bestsol.begin(), bestsol.end(), sol.begin());
   fitness = best->perf();
   return m_running->numEval();
}

void DE::cross(PopulationReal *pop, unsigned pos, tChromosomeReal &crom) {
    if (m_crossover == 'b') {
	crossBin(pop, pos, crom);
    }
    else if (m_crossover == 'e') {
	crossExp(pop, pos, crom);
    }

}

void DE::crossExp(PopulationReal *pop, unsigned pos, tChromosomeReal &crom) {
    int r1, r2, r3;
    tIndividualReal *I1, *I2, *I3;
    int n, nDim;
    int popsize = pop->size();
    int *sample = new int[popsize];

    // Muestreo r1 y r2
    initSample(sample, popsize);
    // Evito que salga el mismo individuo
    sample[pos] = popsize-1;
    --popsize;
    r1 = m_random->getSample(sample, &popsize);
    I1 = pop->getInd(r1);
    r2 = m_random->getSample(sample, &popsize);
    I2 = pop->getInd(r2);
    r3 = m_random->getSample(sample, &popsize);
    I3 = pop->getInd(r3);
    delete[] sample;

    // Obtengo n
    nDim = pop->ndim();
    n = m_random->randint(0, nDim-1);

    tChromosomeReal origin = pop->getInd(pos)->sol();
    copy(origin.begin(), origin.end(), crom.begin());

    for (int i = 0; m_random->rand() < m_CR && (i < nDim); i++) {
		crom[n] = I1->gen(n)
							+ m_F * (I2->gen(n)
							- I3->gen(n));
		n = (n + 1) % nDim;
    }

    // Compruebo que no se salga
    m_problem->getDomain()->clip(crom);
}



void DE::crossBin(PopulationReal *pop, unsigned pos, tChromosomeReal &crom) {
    int r1, r2, r3;
    tIndividualReal *I1, *I2, *I3;
    int n, nDim;
    int popsize = pop->size();
    int *sample = new int[popsize];

    // Muestreo r1 y r2
    initSample(sample, popsize);
    // Evito que salga el mismo individuo
    sample[pos] = popsize-1;
    --popsize;
    r1 = m_random->getSample(sample, &popsize);
    I1 = pop->getInd(r1);
    r2 = m_random->getSample(sample, &popsize);
    I2 = pop->getInd(r2);
    r3 = m_random->getSample(sample, &popsize);
    I3 = pop->getInd(r3);
    delete[] sample;

    // Obtengo n
    nDim = pop->ndim();
    n = m_random->randint(0, nDim-1);

    tChromosomeReal origin = pop->getInd(pos)->sol();
    copy(origin.begin(), origin.end(), crom.begin());

    for (int i = 0; i < nDim; i++) {
	if (m_random->rand() < m_CR) {
		crom[n] = I1->gen(n)
							+ m_F * (I2->gen(n)
							- I3->gen(n));
	}
	n = (n + 1) % nDim;
    }

    // Compruebo que no se salga
    m_problem->getDomain()->clip(crom);
}

void DE::setCross(ICrossBinaryPtr cross) {
    throw new ConfigException("DE::cross can not be changed");
}
