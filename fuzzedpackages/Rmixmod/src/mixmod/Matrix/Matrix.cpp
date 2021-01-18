/***************************************************************************
                             SRC/mixmod/Matrix/Matrix.cpp  description
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
#include "mixmod/Matrix/Matrix.h"
#include "mixmod/Utilities/maths/SelectLibrary.h"

namespace XEM {

//------------
// Constructor
//------------
Matrix::Matrix() {
	_s_pbDimension = 0;
}

Matrix::Matrix(int64_t pbDimension) {
	_s_pbDimension = pbDimension;
}

Matrix::Matrix(Matrix * A) {
	_s_pbDimension = A->getPbDimension();
}

//----------
//Destructor
//----------
Matrix::~Matrix() {
}

void Matrix::edit(std::ostream& flux, std::string before) {

	int64_t i;
	int64_t j;
	double** store = storeToArray();
	for (i = 0; i < _s_pbDimension; i++) {
		flux << '\t' << '\t' << '\t' << '\t';
		for (j = 0; j < _s_pbDimension; j++)
			putDoubleInStream(flux, store[i][j], " ");
		flux << '\n';
	}

	for (i = 0; i < _s_pbDimension; i++) {
		delete [] store[i];
		store[i] = NULL;
	}
	delete [] store;
	store = NULL;

	/*  flux << before << flush;
		  for(j=0;j<i;j++){
		flux << "0.000000 " << flush;
	  }
	  flux << _store << " "<< flush;
	  for(j=i+1;j<_s_pbDimension;j++){
		flux << "0.000000 " << flush;
	  }
	  flux << endl;
	}*/
}

}
