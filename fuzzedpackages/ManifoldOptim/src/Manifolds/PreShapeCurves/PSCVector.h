#ifndef PSCVECTOR_H
#define PSCVECTOR_H

#include "Element.h"
#include <new>
#include <iostream>
#include "def.h"

/*Define the namespace*/
namespace ROPTLIB{

	class PSCVector : public Element{
	public:
		PSCVector(integer r, integer c, integer n);
		virtual PSCVector *ConstructEmpty(void) const;
	};
} /*end of ROPTLIB namespace*/

#endif // end of EUCVECTOR_H
