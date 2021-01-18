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

#include "hybrid.h"
#include <cassert>

using namespace realea;
using namespace realea::internal;

Hybrid::Hybrid(IEAlgorithm *alg, ILocalSearch *ls) : 
	ProxyEA(alg), m_ls(ls), m_random(NULL) {
	m_intensity = 0;
}


/**
 * init the LS method
 */
void Hybrid::initLs(void) {
    if (!m_random) {
       m_random = m_alg->getRandom();
    }

    m_ls->setPopulation(m_alg->getPop());
    m_ls->setProblem(m_problem);
    m_ls->setRunning(m_running);
    m_ls->setRandom(m_random);
    m_ls->setEval(m_eval);
}


void ProxyEA::setNewEval(IEval*eval) {
    m_alg->setNewEval(eval);
    m_eval = (IEvalInd *) new EvalRunning(eval, m_running);
}

