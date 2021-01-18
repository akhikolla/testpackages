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
#ifndef _MAIN_UTIL_H 

#define _MAIN_UTIL_H 1

#include "hybrid.h"
#include "iea.h"
#include "ilocalsearch.h"
#include "random.h"
#include <string>
#include <cstdlib>

using std::string;
using namespace realea;

void setVerbose(void);

double string_to_double( const std::string& s );

/**
 * Get the IEA using the string name 
 *
 * @param alg algorithm name. valide values: ssga|chc|de
 *
 * @param random random generator
 *
 * @return the IEA initialised
 */
IEA *get_EA(string alg, Random &random);

/**
 * Name the EA names using the sep separator
 *
 * @param sep separator (usually ' ', ',' or '|')
 */
string get_EANames(string sep);

void set_InitVerbose(void);

/**
 * Get the LS using the name
 *
 * @param arg_ls algorithm name. Valide values: sw|cmaes
 * @param random random generator
 *
 * @return the LS method
 */
ILocalSearch* get_LS(string arg_ls, DomainRealPtr domain=NULL, Random *random=NULL);
/**
 * Name the LS names using the sep separator
 *
 * @param sep separator (usually ' ', ',' or '|')
 */
string get_LSNames(string sep);

/**
 * Set the LS effort.
 *
 * @param hybrid hybrid algorithm
 *
 * @param effort. value "" does not do anything
 */
void set_Effort(Hybrid *hybrid, string effort); 

/**
 * Set the Maximum number
 *
 * @param ea EA 
 *
 * @param maxeval. value "" does not do anything
 */
void set_MaxEval(IEA *ea, string maxeval); 

#endif
