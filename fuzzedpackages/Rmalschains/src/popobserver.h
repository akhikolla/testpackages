#ifndef _POPULATION_OBSERVER_H 
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

#define _POPULATION_OBSERVER_H 1

#include "real.h"

namespace realea {

/**
 * @class PopulationObserver
 *
 * @brief This class allows to store an information for each individual
 *
 */
class PopulationObserver {
public:
	/**
	 * It is called when the population is restarted completely
	 */
	virtual void reset(void)=0;
	/**
	 * It is called when a specific individual is changed
	 *
	 * @param id identification of the individual changed
	 */
	virtual void notifyChange(unsigned id)=0;
	/**
 	 * It is called when the id are changed.
 	 * 
 	 * @param olid old identification id
 	 * @param newid new identification id
 	 * @see removeWorses
 	 */
	virtual void changeId(unsigned oldid, unsigned newid)=0;
	virtual ~PopulationObserver(void) {}
};

}

#endif
