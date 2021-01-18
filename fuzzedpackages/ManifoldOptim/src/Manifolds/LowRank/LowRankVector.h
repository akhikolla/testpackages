/*
This file defines the class of a point on the tangent space of low-rank manifold R_r^{m times n}

SmartSpace --> ProductElement --> LowRankVector

---- WH
*/

#ifndef LOWRANKVECTOR_H
#define LOWRANKVECTOR_H

#include "ProductElement.h"
#include "StieVector.h"
#include "EucVector.h"

/*Define the namespace*/
namespace ROPTLIB{

	class LowRankVector : public ProductElement{
	public:
		/*Construct an empty vector on the sphere of L^2([0, 1], R) with only size information.
		The representation of the tangent vector is
		R^{Ur times Uc} times R^{Drc times Drc} times R^{Vr times Vc} */
		LowRankVector(integer Ur, integer Uc, integer Drc, integer Vr, integer Vc);

		/*Destruct by deleting variables*/
		virtual ~LowRankVector(void);

		/*Create an object of LowRankVector with same size as this LowRankVector.*/
		virtual LowRankVector *ConstructEmpty(void) const;
	};
} /*end of ROPTLIB namespace*/
#endif // end of OBLIQUEVECTOR_H
