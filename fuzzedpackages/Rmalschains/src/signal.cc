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

#include "signal.h"

Resetable::~Resetable(void) {
    if (m_observers)
	delete m_observers;
}

void Resetable::appendSignal(IReset *obj) {
    if (m_observers == NULL) {
	m_observers = new list<IReset*>;
    }
    m_observers->push_back(obj);
}

void Resetable::clearSignal(void) {
   if (m_observers)
	m_observers->clear();
}

void Resetable::reset(void) {
	list<IReset*>::iterator item;

	if (m_observers == NULL)
	  return;

	for (item = m_observers->begin(); item != m_observers->end();
			item++) {
		(*item)->reset();
	}	

	realReset();
}

void Resetable::clear(void) {
	list<IReset *>::iterator item;

	if (m_observers == NULL)
	  return;

	for (item = m_observers->begin(); item != m_observers->end();
			item++) {
		(*item)->clear();
	}	

	realClear();
}
