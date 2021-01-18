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
#include "populationpso.h"
#include "problem.h"
#include <cassert>
#include <cmath>

using namespace realea;

static ConfigPSO *m_config=NULL;

ConfigPSO::ConfigPSO(DomainRealPtr domain, double inergymin, double inergymax) : m_vmax(domain->getDimension() ) {
     unsigned size = domain->getDimension();
     double min, max;

     if (inergymin < 0.2 || inergymin > 1.2) {
	throw new ConfigException("ConfigPSO::ConfigPSO inergymin");
     }
     if (inergymax < inergymin || inergymax > 1.2) {
	throw new ConfigException("ConfigPSO::ConfigPSO inergymax");
     }


     for (unsigned i = 0; i < size; i++) {
        domain->getValues(i, &min,&max);
	m_vmax[i] = (max-min)/2.0;
     }

     m_c1 = 2.8;
     m_c2 = 1.3;
     assert(m_c1 + m_c2 >= 4);
     double p = m_c1+m_c2;
     double k = 1.0;
     m_x = 2*k/fabs(2-p-sqrt(p*p-4*p));
     m_wmin = inergymin;
     m_wmax = inergymax;
}

double ConfigPSO::w(double ratio) {
   assert(ratio >= 0 && ratio <= 1);
   return m_wmax+ratio*(m_wmin-m_wmax);
}


void realea::setGlobalConfigPSO(ConfigPSO *config) {
    assert(m_config == NULL);
    m_config = config;
}

void realea::delGlobalConfigPSO(void) {
    if (m_config != NULL) {
	delete m_config;
	m_config = NULL;
    }
}

void tIndividualPSO::initPSOParams(Random *random) {
    double vmaxi;
    assert(m_config);
    unsigned size = m_sol.size();
    m_current = m_sol;

    for (unsigned i = 0; i < size; i++) {
        vmaxi = m_config->vmax(i);
	m_velocity[i] = -vmaxi+2*random->rand()*vmaxi;
    }
}

tIndividualPSO::tIndividualPSO(const tChromosomeReal &initialPos, Random *random) : tIndividualReal(initialPos), m_velocity(initialPos.size() ) {
   initPSOParams(random);
}

tIndividualPSO::tIndividualPSO(const tChromosomeReal &initialPos, double fitness, Random *random) : tIndividualReal(initialPos, fitness) {
   initPSOParams(random);
}

void tIndividualPSO::move(tChromosomeReal &best, double ratio) {
    unsigned size = m_sol.size();
    double vmaxi,w;
    unsigned i;

    w = m_config->w(ratio);

    for (i = 0; i < size; ++i) {
	m_velocity[i] = m_config->x()*
				    (m_velocity[i]*w
				     + m_config->c1()*(m_sol[i]-m_current[i]) 
				     + m_config->c2()*(best[i]-m_current[i]));

	vmaxi = m_config->vmax(i);

	if (m_velocity[i] > vmaxi)
	    m_velocity[i] = vmaxi;

        m_current[i] = m_current[i]+m_velocity[i];
    }
}


  
tIndividualReal* PopulationPSO::getInstance(tChromosomeReal &crom) {
    return new tIndividualPSO(crom, m_random);
}

tIndividualReal* PopulationPSO::getInstance(tChromosomeReal &crom, double fitness) {
    return new tIndividualPSO(crom, fitness, m_random);
}

void tIndividualPSO::setCurrentFitness(double newfit) {
    assert(newfit < perf());
    change(m_current, newfit);
}

PopulationPSO::PopulationPSO(Random *random,unsigned int max, unsigned int pob) :
    PopulationReal(random, max, pob) {
}

void PopulationPSO::move(IEvalInd *evalInd, Running *running) {
    unsigned pos_best;
    tIndividualPSO *ind;
    tChromosomeReal best;
    double oldfit,newfit;
    unsigned i, popsize;

    pos_best = getBest();
    best = getInd(pos_best)->sol();

    popsize = size();

    for (i = 0; i < popsize && !running->isFinish(); ++i) {
	ind = (tIndividualPSO*) getInd(i);
	oldfit = ind->perf();
	ind->move(best,running->ratio());
	newfit = evalInd->eval(ind->current());

	if (running->isBetter(newfit, oldfit) ) {
	    change(i, ind->current(), newfit);
	    notifyObservers(i);
	}

    }

}
