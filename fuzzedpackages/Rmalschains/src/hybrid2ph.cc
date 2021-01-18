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
#include "hybrid2ph.h"
#include <cassert>
#include <cmath>

using namespace realea;
using namespace realea::internal;

void Hybrid2Ph::setEffortRatio(double ratio) {

	if (ratio == 1)
	   throw new string("Hybrid2Phd::effortRatio is not valide");

	m_effort = ratio;
}

unsigned Hybrid2Ph::init(void) {
   if (m_effort >= 0) {
      setIntensity((unsigned) ceil(m_running->maxEval()*m_effort));
   }

   initLs();

   return m_alg->init();
}

void Hybrid2Ph::setMaxEval(unsigned maxeval) {
   if (m_effort >= 0) {
      setIntensity((unsigned) ceil(maxeval*m_effort));
   }

   // Reduce the intensity
   assert(maxeval >= m_intensity);
   m_alg->setMaxEval(maxeval-m_intensity);
}

void Hybrid2Ph::setRunning(Running *running) {
   ProxyEA::setRunning(running);
   // Reduce the intensity
   m_alg->setRunning(m_running->getSubRunning(m_running->maxEval()-m_intensity));
}

unsigned Hybrid2Ph::realApply(tChromosomeReal &sol, tFitness &fitness) {
	unsigned num;

	num = m_alg->realApply(sol, fitness);

	// Apply the LS to the best one with the rest of intensity
	ILSParameters *params = m_ls->getInitOptions(sol);

        num += 
	   m_ls->apply(params, sol, fitness, m_intensity);

	delete params;

	return num;
}
