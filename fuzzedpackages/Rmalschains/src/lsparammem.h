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
#ifndef _LS_PARAMETERS_MEMORY_H 

#define _LS_PARAMETERS_MEMORY_H 1

#include "ilocalsearch.h"
#include "popobserver.h"

namespace realea {

namespace internal {

/**
 * @class LSParamsMemory
 *
 * @brief This class store the LS parameters to each individual of a population
 */
class LSParametersMemory : public PopulationObserver {
public:
	/**
	 * It remove the memory when it is called
	 */
	void reset(void);

	/**
	 * It removes the LS parameter
	 *
	 * @param id identification individual
	 */
	void remove(unsigned id);
	/**
	 * It remove the memory of the individual when it is changed
	 * @param id identification individual
	 */
	void notifyChange(unsigned id);
	/**
	 * It update the ids of the memory elements, swapping them
	 *
	 * @param oldir old identification
	 * @param newid new identification
	 */
	void changeId(unsigned oldid, unsigned newid);
public:
	/**
	 * Insert the ILSParameters into the memory.
	 *
	 * If there exists the id, it is replaced, othercase it is inserted.
	 *
	 * @param id identification
	 * @param params parameter values to store
	 */
	void store(unsigned id, ILSParameters *params);
	/**
	 * Recover the memory
	 */
	ILSParameters *recover(unsigned id);
	
	/**
	 * Constructor
	 *
	 * @param tam initial population size
	 */
	LSParametersMemory(unsigned tam);

	/**
	 * Destructor. 
	 *
	 * Frees the structure and LS Parameters
	 */
	~LSParametersMemory(void);

private:
	vector<ILSParameters*> m_params;
};

	typedef vector<ILSParameters*> LSMemory;

}

}


#endif
