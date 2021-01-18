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

#include "populationjderand.h"
#include <cmath>

using namespace realea;

PopulationRealJDERand::PopulationRealJDERand(Random *random, unsigned int max, unsigned int pob) :
    PopulationReal(random,max, pob) {
}

void PopulationRealJDERand::restart(DomainRealPtr domain) {
    reset(domain, getBest());
}

tIndividualReal *PopulationRealJDERand::getInstance(tChromosomeReal &crom) {
    return new tIndividualRealJDERand(crom);
}

tIndividualReal *PopulationRealJDERand::getInstance(tChromosomeReal &crom, double fitness) {
    return new tIndividualRealJDERand(crom, fitness);
}

tIndividualRealJDERand::tIndividualRealJDERand(const tChromosomeReal &com) :
    tIndividualReal(com), 
    m_F(), m_CR() {
}

tIndividualRealJDERand::tIndividualRealJDERand(const tChromosomeReal &com, double fitness) :
    tIndividualReal(com, fitness), 
    m_F(), m_CR() {
}

double tIndividualRealJDERand::getF(string strategy) {
    map<string,double>::iterator it = m_F.find(strategy);

    if (it == m_F.end()) {
        return 0.5;
    }
    else
        return it->second;

}

double tIndividualRealJDERand::getCR(string strategy) {
    map<string,double>::iterator it = m_CR.find(strategy);

    if (it == m_CR.end()) {
        return 0.9;
    }
    else
        return it->second;

}

void tIndividualRealJDERand::setF(string strategy, double F) {
    m_F[strategy] = F;
}

void tIndividualRealJDERand::setCR(string strategy, double CR) {
    m_CR[strategy] = CR;
}

tIndividualRealJDERand::~tIndividualRealJDERand(void) {
}
