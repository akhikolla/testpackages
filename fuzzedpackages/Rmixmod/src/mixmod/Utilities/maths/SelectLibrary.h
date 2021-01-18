/***************************************************************************
                             SRC/mixmod/Utilities/maths/SelectLibrary.h  description
    copyright            : (C) MIXMOD Team - 2001-2016
    email                : contact@mixmod.org
 ***************************************************************************/

/***************************************************************************
    This file is part of MIXMOD
    
    MIXMOD is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    MIXMOD is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with MIXMOD.  If not, see <http://www.gnu.org/licenses/>.

    All informations available on : http://www.mixmod.org                                                                                               
***************************************************************************/
// Route to interface from any mathematical library to our code.
// Only the following include would change if we try another library.


#ifndef XEMmathLib
#define XEMmathLib 1 // default is Eigen
#endif

#if XEMmathLib == 0
#include "mixmod/Utilities/maths/NEWMAT.h" // Old library, still available

#elif XEMmathLib == 1
#include "mixmod/Utilities/maths/Eigen.h" // Default library, only required one

#elif XEMmathLib == 2
#include "mixmod/Utilities/maths/GSL.h" // Optional.

#elif XEMmathLib == 3
#include "mixmod/Utilities/maths/ITpp.h" // Optional.

#elif XEMmathLib == 4
#include "mixmod/Utilities/maths/Armadillo.h" // Optional.

#endif


/*
enum Storage_type {

	row = 0,         // stored by row
	column = 1,      // stored by colmun
};*/




// Current setting: NEWMAT [seg fault with Eigen, TOFIX]
//#include "mixmod/Utilities/maths/NEWMAT.h"
//#include "mixmod/Utilities/maths/Eigen.h"
