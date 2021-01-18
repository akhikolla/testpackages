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

#ifndef _SIGNALABLE

#define _SIGNALABLE 1

#include <list>
#include <cstddef>
using namespace std;

class NotifyObserver {
public:
	virtual void notify(void) {
	}

	virtual ~NotifyObserver(void) {}
};

class NotifyEvalObserver {
public:
	virtual void notify(double fitness)=0; 

	virtual ~NotifyEvalObserver(void) {}
};


class IReset {
public:
   virtual ~IReset(void) {}
   virtual void reset(void){}
   virtual void clear(void){}
};

class Resetable : public IReset {
private:
	list<IReset*> *m_observers;
	
public:
	Resetable(void): m_observers(NULL) {}
	virtual ~Resetable(void);
	virtual void realReset(void) {}
	virtual void realClear(void){}
	/**
	 * Permite añadir un elemento a la lista
	 */
	void appendSignal(IReset *obj);

	/**
	 * Permite eliminar un elemento de la lista
	 */
	void clearSignal(void);

	void reset(void);
	void clear(void);
};

#endif
