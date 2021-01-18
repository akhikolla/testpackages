#ifndef PSCVARIABLE_H
#define PSCVARIABLE_H

#include "Element.h"
#include <new>
#include <iostream>
#include "def.h"
#include "PreShapePathStraighten.h"

/*Define the namespace*/
namespace ROPTLIB{

	class PSCVariable : public Element{
	public:
		// r is the number of points on a curve
		// l is the dimension of a space the curves are in
		// n is the number of curves on a path
		PSCVariable(integer r, integer l, integer n);
		void Generate(double *initial, double *end);
		virtual PSCVariable *ConstructEmpty(void) const;
		virtual void RandInManifold();
	};
} /*end of ROPTLIB namespace*/
#endif // end of EUCVARIABLE_H
