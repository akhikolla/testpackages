/*
This file defines the class of a point on the low-rank manifold R_r^{m times n}, which is represented by 
St(r, m) times R^{r times r} times St(r, n).

SmartSpace --> ProductElement --> LowRankVariable

---- WH
*/

#ifndef LOWRANKVARIABLE_H
#define LOWRANKVARIABLE_H

#include "ProductElement.h"
#include "StieVariable.h"
#include "EucVariable.h"

/*Define the namespace*/
namespace ROPTLIB{

	class LowRankVariable : public ProductElement{
	public:
		/*Construct an empty variable on the low-rank manifold R_r^{m times n} with only size information.*/
		LowRankVariable(integer m, integer n, integer r);

		/*Destruct by deleting all variables*/
		virtual ~LowRankVariable(void);

		/*Create an object of LowRankVariable with same size as this LowRankVariable.*/
		virtual LowRankVariable *ConstructEmpty(void) const;
	};
} /*end of ROPTLIB namespace*/
#endif // end of OBLIQUEVARIABLE_H
