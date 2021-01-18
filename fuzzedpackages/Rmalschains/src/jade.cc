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

#include "jade.h"
#include "populationreal.h"
#include "random.h"
#include "debug.h"
#include <cassert>
#include <cmath>
#include <deque>
#include <vector>

#define CopyVector(a,b) memcpy((a),(b),nDim*sizeof(double))
#define Element(a,b,c)  (a->getInd(b).gen(c))

static bool DEBUG=false;

using namespace realea;

JADE::JADE(Random *random) : ClassEAlgorithm(random), m_CR(-1), m_F(-1) {
    m_crossover = 'e';
    m_meanF = 0.5;
    m_popReductions = 0;
}

void JADE::setDebug(void) {
    DEBUG = true;
}


JADE::~JADE(void) {
}

unsigned JADE::init(void) {


    // Primero inicio la población
    m_pop->reset(m_problem->getDomain());
    // Inicio los distintos elementos: Running, Cross ...
    currentPopReduction = 1;
    reset();
    m_G = 0;
    // Evaluo la población
    m_pop->eval(m_init_eval);
    m_meanCR = 0.5;
    m_meanF = 0.5;
    p =  0.05;
    c = 0.1;
    return m_running->numEval();
}

unsigned JADE::realApply(tChromosomeReal &sol, tFitness &fitness) {
    tChromosomeReal crom(m_problem->getDimension());
    tFitness best_fit;
    unsigned better, currentEval;
    double std = 0.1;
    int successful=0;
    double sF,sF2, sCr;



    unsigned popsize = m_pop->size();
    unsigned MaxFES = m_running->maxEval();
    int countdown = currentPopReduction*MaxFES/(m_popReductions+1);
    m_running->reset();

   // While is not finish
   while (!m_running->isFinish()) {
        if (m_stat)
	    m_stat->newGeneration();

        if(DEBUG) print_info("m_G = %d\n" ,m_G);

        sF = sF2 = sCr= 0;
        successful = 0;
            if(DEBUG){
                print_info("mean_F = %f\n" ,m_meanF);
                print_info("mean_CR = %f\n" ,m_meanCR);
            }


	// Selecciono los individuos a cruzar
	for (unsigned ind = 0; ind < popsize && !m_running->isFinish(); ++ind) {
                        
            do{
                m_F = m_random->normal(std)+m_meanF;
            }while (m_F <= 0 || m_F > 1);

            do{
                m_CR = m_random->normal(std)+m_meanCR;
            }while (m_CR <= 0 || m_CR > 1);



	   // Aplico el cruce sobre el individuo
	    cross(m_pop, ind, crom);
	    // Genero el nuevo individuo
	    tIndividualRealPtr newind = m_pop->getInstance(crom);
	    // Lo evalúo
	    m_new_eval->eval(newind);
            
            if (newind->isBetter(m_pop->getInd(ind)))
            {
                m_archive.push_back(m_pop->getInd(ind));
                m_pop->replaceWithoutDeleting(ind, newind);
                sF += m_F;
                sF2 += pow(m_F,2);
                sCr += m_CR;
                successful+=1;                
            }
            else
            {
                delete newind;
            }

	} // De modificar cada individuo

        if (successful>0)
        {
            m_meanF = (1-c)*m_meanF + c*sF2/sF;
            m_meanCR = (1-c)*m_meanCR + c*sCr/successful;
        }

        while (m_archive.size()>popsize)
        {
	    m_archive.erase(m_archive.begin()+m_random->randint(0,m_archive.size()-1));
        }
        
	better = m_pop->getBest();
	best_fit = m_pop->getInd(better)->perf();

	if (m_stat)
	   m_stat->endGeneration(best_fit);

        currentEval = m_running->numEval();

        //print_info("diffs : %d \n", countdown-currentEval);
        //print_info("evals : %d \n", currentEval);
        //print_info("maxfses : %d \n", MaxFES);

        if (m_popReductions > 0 && (countdown-currentEval)<=0 && popsize > 10 && currentEval < MaxFES-1) {
            //print_info("evals : %d \n", currentEval);
           //print_info("pop reduction : %d\n", currentPopReduction );
            //print_info("tot reduction : %d\n", m_popReductions );
            //print_info("pop size : %d\n", popsize );
           currentPopReduction++;
           countdown = currentPopReduction*MaxFES/(m_popReductions+1);
	   m_pop->reduceHalf();
	   popsize = m_pop->size();
	}

        m_G++;
   } // De running
    
   // Obtengo el mejor
   unsigned pos = m_pop->getBest();
   tIndividualRealPtr best= m_pop->getInd(pos);

   tChromosomeReal bestsol= best->sol();
   copy(bestsol.begin(), bestsol.end(), sol.begin());
   fitness = best->perf();
   return m_running->numEval();
}

void JADE::cross(PopulationReal *pop, unsigned pos, tChromosomeReal &crom) {
    int r1, r2, r3;
    tIndividualReal *frompBests, *fromPop, *fromPopAndArchive, *current;
    int nDim;
    int popsize = pop->size();
    int archivesize = m_archive.size();
    unsigned bestsToConsider = floor(popsize*p + 0.5); //round(popsize*p);
    
    vector<unsigned> pbests = pop->getBests(bestsToConsider);

    int index = m_random->randint(0,bestsToConsider-1);
    r1 = pbests.at(index);
    frompBests = pop->getInd(r1);


    do{
        r2 = m_random->randint(0,popsize-1);
    }while (r2 == r1);

    fromPop = pop->getInd(r2);

    do{
        r3 = m_random->randint(0,popsize+archivesize-1);
    }while (r3 == r1 || r3 == r2);
    
    if(r3 >= popsize)
    {
        fromPopAndArchive = m_archive.at(r3-popsize);
    }
    else
        fromPopAndArchive = pop->getInd(r3);

    
    nDim = pop->ndim();
    
    tChromosomeReal origin = pop->getInd(pos)->sol();
    copy(origin.begin(), origin.end(), crom.begin());

    current = pop->getInd(pos);

    int jrand = m_random->randint(0,nDim-1);


    for (int i = 0; i < nDim; i++) {
	if (m_random->rand() < m_CR || i == jrand) {
            crom[i] = current->gen(i)
                            + m_F * (frompBests->gen(i) - current->gen(i))
                            + m_F * (fromPop->gen(i) - fromPopAndArchive->gen(i));
	}
    }

    // Compruebo que no se salga
    m_problem->getDomain()->clip(crom);
    
}





void JADE::setCross(ICrossBinaryPtr cross) {
    throw new ConfigException("JADE::cross can not be changed");
}



void JADE::setPopReductions(unsigned num) {
    m_popReductions = num;
}

