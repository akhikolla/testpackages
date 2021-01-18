#include "lsparammem.h"
#include <cassert>

using namespace realea;
using namespace realea::internal;

LSParametersMemory::LSParametersMemory(unsigned tam) : m_params(tam) {
}

LSParametersMemory::~LSParametersMemory(void) {
	 reset();
}

void LSParametersMemory::store(unsigned id, ILSParameters *params) {
	ILSParameters *previous;

	if (id > m_params.size()) {
		throw ConfigException("LSParametersMemory::Size");
	}	

	previous = m_params[id];

	if (params != previous) {

	   if (previous != NULL) {
	      delete previous;
	   }

	   m_params[id] = params;
	}
}

ILSParameters *LSParametersMemory::recover(unsigned id) {
	 LSMemory::iterator pos;

	 if (id > m_params.size()) {
		  throw ConfigException("ILSParameters::recover"); 
	 }
	 
	 return m_params[id]; 
}

void LSParametersMemory::remove(unsigned id) {
	if (m_params[id] != NULL) {
	    delete m_params[id];
	    m_params[id] = NULL;
	}
}

void LSParametersMemory::reset(void) {
	 LSMemory::iterator elem;

	 for (elem = m_params.begin(); elem != m_params.end(); ++elem) {
		  if (*elem != NULL) {
		     delete (*elem);
		     *elem = NULL;
		  }
	 }

}


void LSParametersMemory::notifyChange(unsigned id) {
	 ILSParameters *param;
	 assert(id <= m_params.size());

	 param = m_params[id];

	 if (param != NULL) {
	    delete param;
	    m_params[id] = NULL;
	 }

}

void LSParametersMemory::changeId(unsigned oldir, unsigned newid) {
	 if (oldir != newid) {
		  swap(m_params[oldir], m_params[newid]);
	 }
}
