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

#include "sade.h"
#include "populationreal.h"
#include "random.h"
#include "debug.h"
#include <cassert>

#define CopyVector(a,b) memcpy((a),(b),nDim*sizeof(double))
#define Element(a,b,c)  (a->getInd(b).gen(c))

static bool DEBUG=false;

using namespace realea;

SADE::SADE(Random *random) : ClassEAlgorithm(random), m_CR(-1), m_F(-1) {
    m_crossover = 'e';
    m_meanF = 0.3;
    m_popReductions = 0;
}

void SADE::setDebug(void) {
    DEBUG = true;
}

SADE::~SADE(void) {
}

static double minValues(int value1, int value2, int value3) {
    if (value1 <= value2) {
        return (value1 <= value3 ? value1 : value3);
    }
    else
        return (value2 <= value3 ? value2 : value3);
}

unsigned SADE::init(void) {
    for(int k=0;k<SADE_NB_OF_STRATEGIES;k++)
    {
        for(int g=0; g<SADE_LP; g++)
        {
            failure_memory[k][g] = 1;
            success_memory[k][g] = 1;
            CR_memory[k][g] = 0.5;
            
        }
    }

    for(int k=0;k<SADE_NB_OF_STRATEGIES;k++)
    {
        CRmk[k]=0.5;
        strategy_prob[k]=(double)1/SADE_NB_OF_STRATEGIES;
    }

    // Primero inicio la población
    m_pop->reset(m_problem->getDomain());
    // Inicio los distintos elementos: Running, Cross ...
    reset();
    m_G = 0;
    // Evaluo la población
    m_pop->eval(m_init_eval);
    return m_running->numEval();
}

unsigned SADE::realApply(tChromosomeReal &sol, tFitness &fitness) {
    tChromosomeReal crom(m_problem->getDimension());
    tFitness best_fit;
    unsigned better;
    double std_F;
    int successful=0, failed=0;

    unsigned popsize = m_pop->size();
    unsigned MaxFES = m_running->maxEval();
    m_running->reset();
    

   // While is not finish
   while (!m_running->isFinish()) {
        if (m_stat)
	    m_stat->newGeneration();

        print_info("m_G = %d\n" ,m_G);


        //set the median of the CR parameter according to history
        //if the learning period is finished
        if(m_G > SADE_LP)
        {
            setStrategyProb();
            for(int k=0;k<SADE_NB_OF_STRATEGIES;k++)
            {
                failure_memory[k][m_G%SADE_LP]=0;
                success_memory[k][m_G%SADE_LP]=0;
                CRmk[k]=0;
                for(int g=0;g<SADE_LP;g++)
                    CRmk[k]+=CR_memory[k][g];
                CRmk[k]=CRmk[k]/SADE_LP;
            }
        }

        //set for each strategy its CR

        

        for(int k=0;k<SADE_NB_OF_STRATEGIES;k++)
        {
            CRk[k]= m_random->normal(0.1)+CRmk[k];
            while(CRk[k]<0 || CRk[k]>1 )
            {
                CRk[k]= m_random->normal(0.1)+CRmk[k];
            }
        }

        
        if (DEBUG)
        {
            printCRmk();
            printCRk();
            printStrategyProb();
            printSuccessMemory();
            printFailureMemory();
        }


	// Selecciono los individuos a cruzar
	for (unsigned ind = 0; ind < popsize && !m_running->isFinish(); ++ind) {
            
            m_K = m_random->rand();
            std_F = minValues(0.3, 1-m_meanF, m_meanF);
            m_F = m_random->normal(std_F)+m_meanF;

            int strategy = selectStrategy();

            m_CR = CRk[strategy];
            
	   // Aplico el cruce sobre el individuo
	    cross(m_pop, ind, crom,strategy);
	    // Genero el nuevo individuo
	    tIndividualRealPtr newind = m_pop->getInstance(crom);
	    // Lo evalúo
	    m_new_eval->eval(newind);


	    if (newind->isBetter(m_pop->getInd(ind)))
            {
                m_pop->replace(ind, newind);
                success_memory[strategy][m_G%SADE_LP]+=1;
                CR_memory[strategy][m_G%SADE_LP] = CRk[strategy];
                successful+=1;
                 
            }
	    else
            {
                failure_memory[strategy][m_G%SADE_LP]+=1;
		delete newind;
                failed+=1;
            }

	} // De modificar cada individuo 

        print_info(" success = %f",(double)successful/((double)(successful+failed)));
	better = m_pop->getBest();
	best_fit = m_pop->getInd(better)->perf();

	if (m_stat)
	   m_stat->endGeneration(best_fit);

	if (m_popReductions > 0 && m_G % (MaxFES / (m_popReductions+1)) ==
(MaxFES/(m_popReductions+1) - 1) && popsize > 10 && m_G < MaxFES-1) {
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

void SADE::cross(PopulationReal *pop, unsigned pos, tChromosomeReal &crom, int strategy) {
    switch (strategy)
    {
            case 0 : crossRand1Bin(pop,pos,crom);break;
            case 1 : crossRand2Bin(pop,pos,crom);break;
            case 2 : crossRandToBest2Bin(pop,pos,crom);break;
            case 3 : crossRand1Bin(pop,pos,crom);break;
        default : crossRand1Bin(pop,pos,crom);break;
    }
}


void SADE::crossRand1Bin(PopulationReal *pop, unsigned pos, tChromosomeReal &crom) {
    int r1, r2, r3;
    tIndividualReal *I1, *I2, *I3;
    int nDim;
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

    tChromosomeReal origin = pop->getInd(pos)->sol();
    copy(origin.begin(), origin.end(), crom.begin());

    for (int i = 0; i < nDim; i++) {
	if (m_random->rand() < m_CR) {
            crom[i] = I1->gen(i) + m_F * (I2->gen(i) - I3->gen(i));
	}
    }

    // Compruebo que no se salga
    m_problem->getDomain()->clip(crom);

}

void SADE::crossRand2Bin(PopulationReal *pop, unsigned pos, tChromosomeReal &crom) {
    int r1, r2, r3, r4, r5;
    tIndividualReal *I1, *I2, *I3, *I4, *I5;
    int nDim;
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
    r4 = m_random->getSample(sample, &popsize);
    I4 = pop->getInd(r4);
    r5 = m_random->getSample(sample, &popsize);
    I5 = pop->getInd(r5);
    delete[] sample;

    // Obtengo n
    nDim = pop->ndim();

    tChromosomeReal origin = pop->getInd(pos)->sol();
    copy(origin.begin(), origin.end(), crom.begin());

    
    for (int i = 0; i < nDim; i++) {
	if (m_random->rand() < m_CR) {
            crom[i] = I1->gen(i)
                            + m_F * (I2->gen(i) - I3->gen(i))
                            + m_F * (I4->gen(i) - I5->gen(i));
                
                
	}
    }

    // Compruebo que no se salga
    m_problem->getDomain()->clip(crom);
}

void SADE::crossRandToBest2Bin(PopulationReal *pop, unsigned pos, tChromosomeReal &crom) {
    int r1, r2, r3, r4;
    tIndividualReal *I1, *I2, *I3, *I4, *best, *current;
    int nDim;
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
    r4 = m_random->getSample(sample, &popsize);
    delete[] sample;
    // TODO: Revisar
    I4 = pop->getInd(r4);
    current = pop->getInd(pos);
    best = pop->getInd(pop->getBest());

    // Obtengo n
    nDim = pop->ndim();

    tChromosomeReal origin = pop->getInd(pos)->sol();
    copy(origin.begin(), origin.end(), crom.begin());

    for (int i = 0; i < nDim; i++) {
	if (m_random->rand() < m_CR) {
           crom[i] = current->gen(i)
                            + m_F * (best->gen(i) - current->gen(i))
                            + m_F * (I1->gen(i) - I2->gen(i))
                            + m_F * (I3->gen(i) - I4->gen(i));
                
	}
    }
    // Compruebo que no se salga
    m_problem->getDomain()->clip(crom);
}

void SADE::crossCurrentToRand1(PopulationReal *pop, unsigned pos, tChromosomeReal &crom) {
    int r1, r2, r3;
    tIndividualReal *I1, *I2, *I3, *current;
    int nDim;
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
    current = pop->getInd(pos);
    delete[] sample;

    // Obtengo n
    nDim = pop->ndim();

    tChromosomeReal origin = pop->getInd(pos)->sol();
    copy(origin.begin(), origin.end(), crom.begin());


    for (int i = 0; i < nDim; i++) {

        crom[i] = current->gen(i)
                + m_F * (I2->gen(i) - I3->gen(i))
                + m_K * (I1->gen(i) - current->gen(i));
        
    }
    // Compruebo que no se salga
    m_problem->getDomain()->clip(crom);
}



void SADE::setCross(ICrossBinaryPtr cross) {
    throw new ConfigException("SADE::cross can not be changed");
}


void SADE::setStrategyProb(){

    double sum_of_prob=0;
    int sum_of_failures;
    int sum_of_success;

    for(int k=0;k<SADE_NB_OF_STRATEGIES;k++)
    {
        sum_of_failures = sum_of_success = 0;
        for(int g=0; g<SADE_LP; g++)
        {
            sum_of_failures += failure_memory[k][g];
            sum_of_success += success_memory[k][g];
        }
        strategy_prob[k]=((double)sum_of_success/(double)(sum_of_failures+sum_of_success))+SADE_EPS;
        sum_of_prob += strategy_prob[k];
    }
    
    for(int k=0;k<SADE_NB_OF_STRATEGIES;k++)
    {
        strategy_prob[k]=strategy_prob[k]/sum_of_prob;
    }
}



void SADE::setAverageF(double meanF) {
     assert(meanF >= 0 && meanF <= 1);
     m_meanF = meanF;
}
int SADE::selectStrategy()
{
    double r = m_random->rand();
    double aux =0;
    bool end = false;
    int selected=0;

    for(int i=0; i< SADE_NB_OF_STRATEGIES && !end; i++){
        aux += strategy_prob[i];
        if(r <= aux){
             selected = i;
             end  =true;
        }
    }
    return selected;
}



void SADE::printStrategyProb()
{
    print_info("strategy prob = ");
    for(int k=0;k<SADE_NB_OF_STRATEGIES;k++)
    {
         print_info("%f ",strategy_prob[k]);
    }
     print_info("\n");
    
}



void SADE::printCRmk()
{
     print_info("CRmks = ");
    for(int k=0;k<SADE_NB_OF_STRATEGIES;k++)
    {
         print_info("%f ",CRmk[k]);
    }
     print_info("\n");
}

void SADE::printCRk()
{
     print_info("CRs = ");
    for(int k=0;k<SADE_NB_OF_STRATEGIES;k++)
    {
         print_info("%f ",CRk[k]);
    }
     print_info("\n");
}

void SADE::printSuccessMemory()
{
    print_info("success memory\n");
    for(int i=0;i<SADE_LP;i++)
    {
        for (int k=0;k<SADE_NB_OF_STRATEGIES;k++)
        {
            print_info("%d ",success_memory[k][i]);
        }
        print_info("\n");
    }
}

void SADE::printFailureMemory()
{
    print_error("failure memory\n");
    for(int i=0;i<SADE_LP;i++)
    {
        for (int k=0;k<SADE_NB_OF_STRATEGIES;k++)
        {
            print_error("%d ",failure_memory[k][i]);
        }
        print_error("\n");
    }
}

void SADE::printCRMemory()
{
    for(int i=0;i<SADE_LP;i++)
    {
        for (int k=0;k<SADE_NB_OF_STRATEGIES;k++)
        {
            print_info("%f ",CR_memory[k][i]);
        }
        print_info("\n");
    }
}

void SADE::setPopReductions(unsigned num) {
    m_popReductions = num;
}

