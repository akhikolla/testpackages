/***************************************************************************
                             SRC/mixmod/Kernel/Model/ModelType.cpp  description
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
#include "mixmod/Kernel/Model/ModelType.h"

namespace XEM {

//------------
// Constructor
//------------
ModelType::ModelType() {
	_nameModel = defaultGaussianModelName;
	_nbSubDimensionFree = 0;
	_tabSubDimensionFree = NULL;
	_subDimensionEqual = 0;
}

ModelType::ModelType(ModelName name, int64_t nbSubDimensionFree) {
	_nbSubDimensionFree = nbSubDimensionFree;
	_tabSubDimensionFree = NULL;
	_subDimensionEqual = 0;

	_nameModel = name;
}

ModelType::ModelType(const ModelType & iModelType) {
	_nameModel = iModelType._nameModel;
	_subDimensionEqual = iModelType._subDimensionEqual;
	_nbSubDimensionFree = iModelType._nbSubDimensionFree;
	if ((_nbSubDimensionFree != 0) && (iModelType._tabSubDimensionFree)) {
		_tabSubDimensionFree = new int64_t[_nbSubDimensionFree];
		recopyTab(iModelType._tabSubDimensionFree, _tabSubDimensionFree, _nbSubDimensionFree);
	}
	else {
		_tabSubDimensionFree = NULL;
	}
}

ModelType::~ModelType() {
	if (_tabSubDimensionFree) {
		//delete [] _tabSubDimensionFree; //TODO [bauder]: there is something wrong with this variable, 
		                                  //maybe reaffected to something already deleted.
		                                  //Should be fixed cleanly by tracking its affectations.
		                                  //Test after debug by running example da3.
		_tabSubDimensionFree = NULL;
	}
}

/// Comparison operator
bool ModelType::operator ==(const ModelType & modelType) const {

	if (_nameModel != modelType.getModelName()) return false;
	if (_nbSubDimensionFree != modelType._nbSubDimensionFree) return false;
	if (_subDimensionEqual != modelType._subDimensionEqual) return false;
	if (_tabSubDimensionFree) {
		for (int64_t i = 0; i < _nbSubDimensionFree; i++) {
			if (_tabSubDimensionFree[i] != modelType.getTabSubDimensionFreeI(i)) return false;
		}
	}
	return true;
}

void ModelType::input(std::ifstream & fi, int64_t nbCluster) {

	_nbSubDimensionFree = nbCluster;
	std::string keyWord = "";
	int64_t dim;
	std::string a = "";

	fi >> a;
	if (a.compare("Gaussian_p_L_I") == 0) {
		_nameModel = Gaussian_p_L_I;
	}
	else if (a.compare("Gaussian_p_Lk_I") == 0) {
		_nameModel = Gaussian_p_Lk_I;
	}
	else if (a.compare("Gaussian_p_L_B") == 0) {
		_nameModel = Gaussian_p_L_B;
	}
	else if (a.compare("Gaussian_p_Lk_B") == 0) {
		_nameModel = Gaussian_p_Lk_B;
	}
	else if (a.compare("Gaussian_p_L_Bk") == 0) {
		_nameModel = Gaussian_p_L_Bk;
	}
	else if (a.compare("Gaussian_p_Lk_Bk") == 0) {
		_nameModel = Gaussian_p_Lk_Bk;
	}
	else if (a.compare("Gaussian_p_L_C") == 0) {
		_nameModel = Gaussian_p_L_C;
	}
	else if (a.compare("Gaussian_p_Lk_C") == 0) {
		_nameModel = Gaussian_p_Lk_C;
	}
	else if (a.compare("Gaussian_p_L_D_Ak_D") == 0) {
		_nameModel = Gaussian_p_L_D_Ak_D;
	}
	else if (a.compare("Gaussian_p_Lk_D_Ak_D") == 0) {
		_nameModel = Gaussian_p_Lk_D_Ak_D;
	}
	else if (a.compare("Gaussian_p_L_Dk_A_Dk") == 0) {
		_nameModel = Gaussian_p_L_Dk_A_Dk;
	}
	else if (a.compare("Gaussian_p_Lk_Dk_A_Dk") == 0) {
		_nameModel = Gaussian_p_Lk_Dk_A_Dk;
	}
	else if (a.compare("Gaussian_p_L_Ck") == 0) {
		_nameModel = Gaussian_p_L_Ck;
	}
	else if (a.compare("Gaussian_p_Lk_Ck") == 0) {
		_nameModel = Gaussian_p_Lk_Ck;
	}
	else if (a.compare("Gaussian_pk_L_I") == 0) {
		_nameModel = Gaussian_pk_L_I;
	}
	else if (a.compare("Gaussian_pk_Lk_I") == 0) {
		_nameModel = Gaussian_pk_Lk_I;
	}
	else if (a.compare("Gaussian_pk_L_B") == 0) {
		_nameModel = Gaussian_pk_L_B;
	}
	else if (a.compare("Gaussian_pk_Lk_B") == 0) {
		_nameModel = Gaussian_pk_Lk_B;
	}
	else if (a.compare("Gaussian_pk_L_Bk") == 0) {
		_nameModel = Gaussian_pk_L_Bk;
	}
	else if (a.compare("Gaussian_pk_Lk_Bk") == 0) {
		_nameModel = Gaussian_pk_Lk_Bk;
	}
	else if (a.compare("Gaussian_pk_L_C") == 0) {
		_nameModel = Gaussian_pk_L_C;
	}
	else if (a.compare("Gaussian_pk_Lk_C") == 0) {
		_nameModel = Gaussian_pk_Lk_C;
	}
	else if (a.compare("Gaussian_pk_L_D_Ak_D") == 0) {
		_nameModel = Gaussian_pk_L_D_Ak_D;
	}
	else if (a.compare("Gaussian_pk_Lk_D_Ak_D") == 0) {
		_nameModel = Gaussian_pk_Lk_D_Ak_D;
	}
	else if (a.compare("Gaussian_pk_Lk_Dk_A_Dk") == 0) {
		_nameModel = Gaussian_pk_Lk_Dk_A_Dk;
	}
	else if (a.compare("Gaussian_pk_L_Dk_A_Dk") == 0) {
		_nameModel = Gaussian_pk_L_Dk_A_Dk;
	}
	else if (a.compare("Gaussian_pk_L_Ck") == 0) {
		_nameModel = Gaussian_pk_L_Ck;
	}
	else if (a.compare("Gaussian_pk_Lk_Ck") == 0) {
		_nameModel = Gaussian_pk_Lk_Ck;
	}

	// Binary models
	else if (a.compare("Binary_p_E") == 0) {
		_nameModel = Binary_p_E;
	}
	else if (a.compare("Binary_p_Ek") == 0) {
		_nameModel = Binary_p_Ek;
	}
	else if (a.compare("Binary_p_Ej") == 0) {
		_nameModel = Binary_p_Ej;
	}
	else if (a.compare("Binary_p_Ekj") == 0) {
		_nameModel = Binary_p_Ekj;
	}
	else if (a.compare("Binary_p_Ekjh") == 0) {
		_nameModel = Binary_p_Ekjh;
	}
	else if (a.compare("Binary_pk_E") == 0) {
		_nameModel = Binary_pk_E;
	}
	else if (a.compare("Binary_pk_Ek") == 0) {
		_nameModel = Binary_pk_Ek;
	}
	else if (a.compare("Binary_pk_Ej") == 0) {
		_nameModel = Binary_pk_Ej;
	}
	else if (a.compare("Binary_pk_Ekj") == 0) {
		_nameModel = Binary_pk_Ekj;
	}
	else if (a.compare("Binary_pk_Ekjh") == 0) {
		_nameModel = Binary_pk_Ekjh;
	}

	//HDDA models
	else if (a.compare("Gaussian_HD_pk_AkjBkQkD") == 0) {
		_nameModel = Gaussian_HD_pk_AkjBkQkD;
		fi >> keyWord;
		ConvertBigtoLowString(keyWord);
		if (keyWord.compare("subdimensionequal") == 0) {
			fi >> dim;
			_subDimensionEqual = dim;
		}
		else {
			THROW(InputException, wrongSubDimension);
		}
	}
	else if (a.compare("Gaussian_HD_pk_AkjBQkD") == 0) {
		_nameModel = Gaussian_HD_pk_AkjBQkD;
		fi >> keyWord;
		ConvertBigtoLowString(keyWord);
		if (keyWord.compare("subdimensionequal") == 0) {
			fi >> dim;
			_subDimensionEqual = dim;
		}
		else {
			THROW(InputException, wrongSubDimension);
		}
	}
	else if (a.compare("Gaussian_HD_pk_AjBkQkD") == 0) {
		_nameModel = Gaussian_HD_pk_AjBkQkD;
		fi >> keyWord;
		ConvertBigtoLowString(keyWord);
		if (keyWord.compare("subdimensionequal") == 0) {
			fi >> dim;
			_subDimensionEqual = dim;
		}
		else {
			THROW(InputException, wrongSubDimension);
		}
	}
	else if (a.compare("Gaussian_HD_pk_AjBQkD") == 0) {
		_nameModel = Gaussian_HD_pk_AjBQkD;
		fi >> keyWord;
		ConvertBigtoLowString(keyWord);
		if (keyWord.compare("subdimensionequal") == 0) {
			fi >> dim;
			_subDimensionEqual = dim;
		}
		else {
			THROW(InputException, wrongSubDimension);
		}
	}
	else if (a.compare("Gaussian_HD_pk_AkBkQkD") == 0) {
		_nameModel = Gaussian_HD_pk_AkBkQkD;
		fi >> keyWord;
		ConvertBigtoLowString(keyWord);
		if (keyWord.compare("subdimensionequal") == 0) {
			fi >> dim;
			_subDimensionEqual = dim;
		}
		else {
			THROW(InputException, wrongSubDimension);
		}
	}
	else if (a.compare("Gaussian_HD_pk_AkBQkD") == 0) {
		_nameModel = Gaussian_HD_pk_AkBQkD;
		fi >> keyWord;
		ConvertBigtoLowString(keyWord);
		if (keyWord.compare("subdimensionequal") == 0) {
			fi >> dim;
			_subDimensionEqual = dim;
		}
		else {
			THROW(InputException, wrongSubDimension);
		}
	}
	else if (a.compare("Gaussian_HD_pk_AkjBkQkDk") == 0) {
		_nameModel = Gaussian_HD_pk_AkjBkQkDk;
		fi >> keyWord;
		ConvertBigtoLowString(keyWord);
		if (keyWord.compare("subdimensionfree") == 0) {
			_tabSubDimensionFree = new int64_t[_nbSubDimensionFree];
			for (int64_t k = 0; k < _nbSubDimensionFree; k++) {
				fi >> dim;
				_tabSubDimensionFree[k] = dim;
			}
		}
		else {
			THROW(InputException, wrongSubDimension);
		}
	}
	else if (a.compare("Gaussian_HD_pk_AkBkQkDk") == 0) {
		_nameModel = Gaussian_HD_pk_AkBkQkDk;
		fi >> keyWord;
		ConvertBigtoLowString(keyWord);
		if (keyWord.compare("subdimensionfree") == 0) {
			_tabSubDimensionFree = new int64_t[_nbSubDimensionFree];
			for (int64_t k = 0; k < _nbSubDimensionFree; k++) {
				fi >> dim;
				_tabSubDimensionFree[k] = dim;
			}
		}
		else {
			THROW(InputException, wrongSubDimension);
		}
	}
	else if (a.compare("Gaussian_HD_p_AkjBkQkD") == 0) {
		_nameModel = Gaussian_HD_p_AkjBkQkD;
		fi >> keyWord;
		ConvertBigtoLowString(keyWord);
		if (keyWord.compare("subdimensionequal") == 0) {
			fi >> dim;
			_subDimensionEqual = dim;
		}
		else {
			THROW(InputException, wrongSubDimension);
		}
	}
	else if (a.compare("Gaussian_HD_p_AkjBQkD") == 0) {
		_nameModel = Gaussian_HD_p_AkjBQkD;
		fi >> keyWord;
		ConvertBigtoLowString(keyWord);
		if (keyWord.compare("subdimensionequal") == 0) {
			fi >> dim;
			_subDimensionEqual = dim;
		}
		else {
			THROW(InputException, wrongSubDimension);
		}
	}
	else if (a.compare("Gaussian_HD_p_AjBkQkD") == 0) {
		_nameModel = Gaussian_HD_p_AjBkQkD;
		fi >> keyWord;
		ConvertBigtoLowString(keyWord);
		if (keyWord.compare("subdimensionequal") == 0) {
			fi >> dim;
			_subDimensionEqual = dim;
		}
		else {
			THROW(InputException, wrongSubDimension);
		}
	}
	else if (a.compare("Gaussian_HD_p_AjBQkD") == 0) {
		_nameModel = Gaussian_HD_p_AjBQkD;
		fi >> keyWord;
		ConvertBigtoLowString(keyWord);
		if (keyWord.compare("subdimensionequal") == 0) {
			fi >> dim;
			_subDimensionEqual = dim;
		}
		else {
			THROW(InputException, wrongSubDimension);
		}
	}
	else if (a.compare("Gaussian_HD_p_AkBkQkD") == 0) {
		_nameModel = Gaussian_HD_p_AkBkQkD;
		fi >> keyWord;
		ConvertBigtoLowString(keyWord);
		if (keyWord.compare("subdimensionequal") == 0) {
			fi >> dim;
			_subDimensionEqual = dim;
		}
		else {
			THROW(InputException, wrongSubDimension);
		}
	}
	else if (a.compare("Gaussian_HD_p_AkBQkD") == 0) {
		_nameModel = Gaussian_HD_p_AkBQkD;
		fi >> keyWord;
		ConvertBigtoLowString(keyWord);
		if (keyWord.compare("subdimensionequal") == 0) {
			fi >> dim;
			_subDimensionEqual = dim;
		}
		else {
			THROW(InputException, wrongSubDimension);
		}
	}
	else if (a.compare("Gaussian_HD_p_AkjBkQkDk") == 0) {
		_nameModel = Gaussian_HD_p_AkjBkQkDk;
		fi >> keyWord;
		ConvertBigtoLowString(keyWord);
		if (keyWord.compare("subdimensionfree") == 0) {
			_tabSubDimensionFree = new int64_t[_nbSubDimensionFree];
			int64_t k;
			for (k = 0; k < _nbSubDimensionFree; k++) {
				fi >> dim;
				_tabSubDimensionFree[k] = dim;
			}
		}
		else {
			THROW(InputException, wrongSubDimension);
		}
	}
	else if (a.compare("Gaussian_HD_p_AkBkQkDk") == 0) {
		_nameModel = Gaussian_HD_p_AkBkQkDk;
		fi >> keyWord;
		ConvertBigtoLowString(keyWord);
		if (keyWord.compare("subdimensionfree") == 0) {
			_tabSubDimensionFree = new int64_t [_nbSubDimensionFree];
			for (int64_t k = 0; k < _nbSubDimensionFree; k++) {
				fi >> dim;
				_tabSubDimensionFree[k] = dim;
			}
		}
		else {
			THROW(InputException, wrongSubDimension);
		}
	}

	else {
		THROW(InputException, wrongModelType);
	}
	
	//delete[] a;
	//delete[] keyWord;
}

// <<
std::ostream & operator<<(std::ostream & fo, ModelType & modelType) {
	std::string name = ModelNameToString(modelType._nameModel);
	fo << name << endl;
	if (modelType._subDimensionEqual != 0) {
		fo << "subDimensionEqual = " << modelType._subDimensionEqual << endl;
		;
	}
	if (modelType._nbSubDimensionFree != 0 && modelType._tabSubDimensionFree) {
		fo << "subDimensionFree : " << endl;
		for (int64_t i = 0; i < modelType._nbSubDimensionFree; i++) {
			fo << modelType._tabSubDimensionFree[i] << endl;
		}
		fo << endl;
	}
	return fo;
}

//-------------------------
// print ModelType Shortcut
//-------------------------
void ModelType::printShortcut(std::ostream & flux)const {

	switch (_nameModel) {
	case (Gaussian_p_L_B):
	case (Gaussian_p_Lk_B):
	case (Gaussian_p_L_Bk):
	case (Gaussian_p_Lk_Bk):
	case (Gaussian_pk_L_B):
	case (Gaussian_pk_Lk_B):
	case (Gaussian_pk_L_Bk):
	case (Gaussian_pk_Lk_Bk):
		flux << "D" << flush;
		break;

	case (Gaussian_p_L_C):
	case (Gaussian_p_Lk_C):
	case (Gaussian_p_L_D_Ak_D):
	case (Gaussian_p_Lk_D_Ak_D):
	case (Gaussian_p_L_Dk_A_Dk):
	case (Gaussian_p_Lk_Dk_A_Dk):
	case (Gaussian_p_L_Ck):
	case (Gaussian_p_Lk_Ck):
	case (Gaussian_pk_L_C):
	case (Gaussian_pk_Lk_C):
	case (Gaussian_pk_L_D_Ak_D):
	case (Gaussian_pk_Lk_D_Ak_D):
	case (Gaussian_pk_L_Dk_A_Dk):
	case (Gaussian_pk_Lk_Dk_A_Dk):
	case (Gaussian_pk_L_Ck):
	case (Gaussian_pk_Lk_Ck):
		flux << "G" << flush;
		break;

	case (Gaussian_p_L_I):
	case (Gaussian_p_Lk_I):
	case (Gaussian_pk_L_I):
	case (Gaussian_pk_Lk_I):
		flux << "S" << flush;
		break;

	// Binary models
	case (Binary_p_E):
	case (Binary_p_Ek):
	case (Binary_p_Ej):
	case (Binary_p_Ekj):
	case (Binary_p_Ekjh):
	case (Binary_pk_E):
	case (Binary_pk_Ek):
	case (Binary_pk_Ej):
	case (Binary_pk_Ekj):
	case (Binary_pk_Ekjh):
		flux << "B" << flush;
		break;

	case (Gaussian_HD_pk_AkjBkQkD):
	case (Gaussian_HD_pk_AkjBkQkDk):
	case (Gaussian_HD_pk_AkjBQkD):
	case (Gaussian_HD_pk_AjBkQkD):
	case (Gaussian_HD_pk_AjBQkD):
	case (Gaussian_HD_pk_AkBkQkD):
	case (Gaussian_HD_pk_AkBkQkDk):
	case (Gaussian_HD_pk_AkBQkD):
	case (Gaussian_HD_p_AkjBkQkD):
	case (Gaussian_HD_p_AkjBkQkDk):
	case (Gaussian_HD_p_AkjBQkD):
	case (Gaussian_HD_p_AjBkQkD):
	case (Gaussian_HD_p_AjBQkD):
	case (Gaussian_HD_p_AkBkQkD):
	case (Gaussian_HD_p_AkBkQkDk):
	case (Gaussian_HD_p_AkBQkD):
		flux << "H" << flush;
		break;
	default:
		THROW(OtherException, internalMixmodError);
	}
}

//---------------------
// print out Model Type
//---------------------
void ModelType::print(std::ostream & flux) const {

	switch (_nameModel) {

	// Gaussian models
	case (Gaussian_p_L_B):
		flux << "p_L_B         ";
		break;
	case (Gaussian_p_Lk_B):
		flux << "p_Lk_B        ";
		break;
	case (Gaussian_p_L_Bk):
		flux << "p_L_Bk        ";
		break;
	case (Gaussian_p_Lk_Bk):
		flux << "p_Lk_Bk       ";
		break;
	case (Gaussian_pk_L_B):
		flux << "pk_L_B        ";
		break;
	case (Gaussian_pk_Lk_B):
		flux << "pk_Lk_B       ";
		break;
	case (Gaussian_pk_L_Bk):
		flux << "pk_L_Bk       ";
		break;
	case (Gaussian_pk_Lk_Bk):
		flux << "pk_Lk_Bk      ";
		break;
	case (Gaussian_p_L_C):
		flux << "p_L_C         ";
		break;
	case (Gaussian_p_Lk_C):
		flux << "p_Lk_C        ";
		break;
	case (Gaussian_p_L_D_Ak_D):
		flux << "p_L_D_Ak_D    ";
		break;
	case (Gaussian_p_Lk_D_Ak_D):
		flux << "p_Lk_D_Ak_D   ";
		break;
	case (Gaussian_p_L_Dk_A_Dk):
		flux << "p_L_Dk_A_Dk   ";
		break;
	case (Gaussian_p_Lk_Dk_A_Dk):
		flux << "p_Lk_Dk_A_Dk  ";
		break;
	case (Gaussian_p_L_Ck):
		flux << "p_L_Ck        ";
		break;
	case (Gaussian_p_Lk_Ck):
		flux << "p_Lk_Ck       ";
		break;
	case (Gaussian_pk_L_C):
		flux << "pk_L_C        ";
		break;
	case (Gaussian_pk_Lk_C):
		flux << "pk_Lk_C       ";
		break;
	case (Gaussian_pk_L_D_Ak_D):
		flux << "pk_L_D_Ak_D   ";
		break;
	case (Gaussian_pk_Lk_D_Ak_D):
		flux << "pk_Lk_D_Ak_D  ";
		break;
	case (Gaussian_pk_L_Dk_A_Dk):
		flux << "pk_L_Dk_A_Dk  ";
		break;
	case (Gaussian_pk_Lk_Dk_A_Dk):
		flux << "pk_Lk_Dk_A_Dk ";
		break;
	case (Gaussian_pk_L_Ck):
		flux << "pk_L_Ck       ";
		break;
	case (Gaussian_pk_Lk_Ck):
		flux << "pk_Lk_Ck      ";
		break;
	case (Gaussian_p_L_I):
		flux << "p_L_I         ";
		break;
	case (Gaussian_p_Lk_I):
		flux << "p_Lk_I        ";
		break;
	case (Gaussian_pk_L_I):
		flux << "pk_L_I        ";
		break;
	case (Gaussian_pk_Lk_I):
		flux << "pk_Lk_I       ";
		break;

	// Binary models
	case (Binary_p_E):
		flux << "Binary_p_E    ";
		break;
	case (Binary_p_Ek):
		flux << "Binary_p_Ek   ";
		break;
	case (Binary_p_Ej):
		flux << "Binary_p_Ej   ";
		break;
	case (Binary_p_Ekj):
		flux << "Binary_p_Ekj  ";
		break;
	case (Binary_p_Ekjh):
		flux << "Binary_p_Ekjh  ";
		break;
	case (Binary_pk_E):
		flux << "Binary_pk_E   ";
		break;
	case (Binary_pk_Ek):
		flux << "Binary_pk_Ek  ";
		break;
	case (Binary_pk_Ej):
		flux << "Binary_pk_Ej  ";
		break;
	case (Binary_pk_Ekj):
		flux << "Binary_pk_Ekj ";
		break;
	case (Binary_pk_Ekjh):
		flux << "Binary_pk_Ekjh ";
		break;

	//HDDA models
	case (Gaussian_HD_pk_AkjBkQkD):
		flux << "HD_pk_AkjBkQkD   ";
		break;
	case (Gaussian_HD_pk_AkjBkQkDk):
		flux << "HD_pk_AkjBkQkDk   ";
		break;
	case (Gaussian_HD_pk_AkjBQkD):
		flux << "HD_pk_AkjBQkD    ";
		break;
	case (Gaussian_HD_pk_AjBkQkD):
		flux << "HD_pk_AjBkQkD    ";
		break;
	case (Gaussian_HD_pk_AjBQkD):
		flux << "HD_pk_AjBQkD     ";
		break;
	case (Gaussian_HD_pk_AkBkQkD):
		flux << "HD_pk_AkBkQkD    ";
		break;
	case (Gaussian_HD_pk_AkBkQkDk):
		flux << "HD_pk_AkBkQkDk    ";
		break;
	case (Gaussian_HD_pk_AkBQkD):
		flux << "HD_pk_AkBQkD     ";
		break;
	case (Gaussian_HD_p_AkjBkQkD):
		flux << "HD_p_AkjBkQkD    ";
		break;
	case (Gaussian_HD_p_AkjBkQkDk):
		flux << "HD_p_AkjBkQkDk     ";
		break;
	case (Gaussian_HD_p_AkjBQkD):
		flux << "HD_p_AkjBQkD     ";
		break;
	case (Gaussian_HD_p_AjBkQkD):
		flux << "HD_p_AjBkQkD     ";
		break;
	case (Gaussian_HD_p_AjBQkD):
		flux << "HD_p_AjBQkD      ";
		break;
	case (Gaussian_HD_p_AkBkQkD):
		flux << "HD_p_AkBkQkD     ";
		break;
	case (Gaussian_HD_p_AkBkQkDk):
		flux << "HD_p_AkBkQkDk     ";
		break;
	case (Gaussian_HD_p_AkBQkD):
		flux << "HD_p_AkBQkD      ";
		break;
	case (UNKNOWN_MODEL_NAME):
		flux << "UNKNOWN_MODEL_NAME";
		break;
	default:
		THROW(OtherException, internalMixmodError);
	}
	flux << flush;
}

//---------------
/// editModelType
//----------------
void ModelType::edit(std::ostream & oFile) {

	oFile << "\t\t\tModel Type : ";

	switch (_nameModel) {

	// Gaussian models
	case (Gaussian_p_L_B):
		oFile << "Gaussian Diagonal Model : ";
		oFile << "p_L_B";
		break;
	case (Gaussian_p_Lk_B):
		oFile << "p_Lk_B";
		break;
		oFile << "Gaussian Diagonal Model : ";
	case (Gaussian_p_L_Bk):
		oFile << "Gaussian Diagonal Model : ";
		oFile << "p_L_Bk";
		break;
	case (Gaussian_p_Lk_Bk):
		oFile << "Gaussian Diagonal Model : ";
		oFile << "p_Lk_Bk";
		break;
	case (Gaussian_pk_L_B):
		oFile << "Gaussian Diagonal Model : ";
		oFile << " pk_L_B";
		break;
	case (Gaussian_pk_Lk_B):
		oFile << "Gaussian Diagonal Model : ";
		oFile << " pk_Lk_B";
		break;
	case (Gaussian_pk_L_Bk):
		oFile << "Gaussian Diagonal Model : ";
		oFile << "pk_L_Bk";
		break;
	case (Gaussian_pk_Lk_Bk):
		oFile << "Gaussian Diagonal Model : ";
		oFile << "pk_Lk_Bk";
		break;
	case (Gaussian_p_L_I):
		oFile << "Gaussian Spherical Model : ";
		oFile << "p_L_I";
		break;
	case (Gaussian_p_Lk_I):
		oFile << "Gaussian Spherical Model : ";
		oFile << "p_Lk_I";
		break;
	case (Gaussian_pk_L_I):
		oFile << "Gaussian Spherical Model : ";
		oFile << "pk_L_I";
		break;
	case (Gaussian_pk_Lk_I):
		oFile << "Gaussian Spherical Model : ";
		oFile << "pk_Lk_I";
		break;
	case (Gaussian_p_L_C):
		oFile << "Gaussian Ellipsoidal Model : ";
		oFile << "p_L_C";
		break;
	case (Gaussian_p_Lk_C):
		oFile << "Gaussian Ellipsoidal Model : ";
		oFile << "p_Lk_C";
		break;
	case (Gaussian_p_L_D_Ak_D):
		oFile << "Gaussian Ellipsoidal Model : ";
		oFile << "p_L_D_Ak_D";
		break;
	case (Gaussian_p_Lk_D_Ak_D):
		oFile << "Gaussian Ellipsoidal Model : ";
		oFile << "p_Lk_D_Ak_D";
		break;
	case (Gaussian_p_L_Dk_A_Dk):
		oFile << "Gaussian Ellipsoidal Model : ";
		oFile << "p_L_Dk_A_Dk";
		break;
	case (Gaussian_p_Lk_Dk_A_Dk):
		oFile << "Gaussian Ellipsoidal Model : ";
		oFile << "p_Lk_Dk_A_Dk";
		break;
	case (Gaussian_p_L_Ck):
		oFile << "Gaussian Ellipsoidal Model : ";
		oFile << "p_L_Ck";
		break;
	case (Gaussian_p_Lk_Ck):
		oFile << "Gaussian Ellipsoidal Model : ";
		oFile << "p_Lk_Ck";
		break;
	case (Gaussian_pk_L_C):
		oFile << "Gaussian Ellipsoidal Model : ";
		oFile << "pk_L_C";
		break;
	case (Gaussian_pk_Lk_C):
		oFile << "Gaussian Ellipsoidal Model : ";
		oFile << "pk_Lk_C";
		break;
	case (Gaussian_pk_L_D_Ak_D):
		oFile << "Gaussian Ellipsoidal Model : ";
		oFile << "pk_L_D_Ak_D";
		break;
	case (Gaussian_pk_Lk_D_Ak_D):
		oFile << "Gaussian Ellipsoidal Model : ";
		oFile << "pk_Lk_D_Ak_D";
		break;
	case (Gaussian_pk_L_Dk_A_Dk):
		oFile << "Gaussian Ellipsoidal Model : ";
		oFile << "pk_L_Dk_A_Dk";
		break;
	case (Gaussian_pk_Lk_Dk_A_Dk):
		oFile << "Gaussian Ellipsoidal Model : ";
		oFile << "pk_Lk_Dk_A_Dk";
		break;
	case (Gaussian_pk_L_Ck):
		oFile << "Gaussian Ellipsoidal Model : ";
		oFile << "pk_L_Ck";
		break;
	case (Gaussian_pk_Lk_Ck):
		oFile << "Gaussian Ellipsoidal Model : ";
		oFile << "pk_Lk_Ck";
		break;

		// Binary models //
	case (Binary_p_E):
		oFile << "Binary Model : ";
		oFile << "Binary_p_E";
		break;
	case (Binary_p_Ek):
		oFile << "Binary Model : ";
		oFile << "Binary_p_Ek";
		break;
	case (Binary_p_Ej):
		oFile << "Binary Model : ";
		oFile << "Binary_p_Ej";
		break;
	case (Binary_p_Ekj):
		oFile << "Binary Model : ";
		oFile << "Binary_p_Ekj";
		break;
	case (Binary_p_Ekjh):
		oFile << "Binary Model : ";
		oFile << "Binary_p_Ekjh";
		break;
	case (Binary_pk_E):
		oFile << "Binary Model : ";
		oFile << "Binary_pk_E";
		break;
	case (Binary_pk_Ek):
		oFile << "Binary Model : ";
		oFile << "Binary_pk_Ek";
		break;
	case (Binary_pk_Ej):
		oFile << "Binary Model : ";
		oFile << "Binary_pk_Ej";
		break;
	case (Binary_pk_Ekj):
		oFile << "Binary Model : ";
		oFile << "Binary_pk_Ekj";
		break;
	case (Binary_pk_Ekjh):
		oFile << "Binary Model : ";
		oFile << "Binary_pk_Ekjh";
		break;

	case (Gaussian_HD_pk_AkjBkQkD):
		oFile << "HD Model : ";
		oFile << "Gaussian_HD_pk_AkjBkQkD";
		break;
	case (Gaussian_HD_pk_AkjBkQkDk):
		oFile << "HD Model : ";
		oFile << "Gaussian_HD_pk_AkjBkQkDk";
		break;
	case (Gaussian_HD_pk_AkjBQkD):
		oFile << "HD Model : ";
		oFile << "Gaussian_HD_pk_AkjBQkD";
		break;
	case (Gaussian_HD_pk_AjBkQkD):
		oFile << "HD Model : ";
		oFile << "Gaussian_HD_pk_AjBkQkD";
		break;
	case (Gaussian_HD_pk_AjBQkD):
		oFile << "HD Model : ";
		oFile << "Gaussian_HD_pk_AjBQkD";
		break;
	case (Gaussian_HD_pk_AkBkQkD):
		oFile << "HD Model : ";
		oFile << "Gaussian_HD_pk_AkBkQkD";
		break;
	case (Gaussian_HD_pk_AkBkQkDk):
		oFile << "HD Model : ";
		oFile << "Gaussian_HD_pk_AkBkQkDk";
		break;
	case (Gaussian_HD_pk_AkBQkD):
		oFile << "HD Model : ";
		oFile << "Gaussian_HD_pk_AkBQkD";
		break;
	case (Gaussian_HD_p_AkjBkQkD):
		oFile << "HD Model : ";
		oFile << "Gaussian_HD_p_AkjBkQkD";
		break;
	case (Gaussian_HD_p_AkjBkQkDk):
		oFile << "HD Model : ";
		oFile << "Gaussian_HD_p_AkjBkQkDk";
		break;
	case (Gaussian_HD_p_AkjBQkD):
		oFile << "HD Model : ";
		oFile << "Gaussian_HD_p_AkjBQkD";
		break;
	case (Gaussian_HD_p_AjBkQkD):
		oFile << "HD Model : ";
		oFile << "Gaussian_HD_p_AjBkQkD";
		break;
	case (Gaussian_HD_p_AjBQkD):
		oFile << "HD Model : ";
		oFile << "Gaussian_HD_p_AjBQkD";
		break;
	case (Gaussian_HD_p_AkBkQkD):
		oFile << "HD Model : ";
		oFile << "Gaussian_HD_p_AkBkQkD";
		break;
	case (Gaussian_HD_p_AkBkQkDk):
		oFile << "HD Model : ";
		oFile << "Gaussian_HD_p_AkBkQkDk";
		break;
	case (Gaussian_HD_p_AkBQkD):
		oFile << "HD Model : ";
		oFile << "Gaussian_HD_p_AkBQkD";
		break;

	default:
		oFile << "Model Type Error";
	}
	oFile << endl;
	oFile << "\t\t\t----------" << endl << endl;
}

void ModelType::setTabSubDimensionFree(int64_t iTabSubDimensionFree, int64_t position) {
	if (!isHD(_nameModel) || !isFreeSubDimension(_nameModel)) THROW(InputException,wrongModelInSetSubDimensionFree);	if (position>=0 && position<_nbSubDimensionFree){
    if (_tabSubDimensionFree == NULL) {
      _tabSubDimensionFree = new int64_t[_nbSubDimensionFree];
    }
    _tabSubDimensionFree[position] = iTabSubDimensionFree;
	}
	else{
	  THROW(InputException, wrongModelPositionInSetSubDimensionFree);
	}
}

void ModelType::setSubDimensionEqual(int64_t iSubDimensionEqual) {
  if (!isHD(_nameModel) || isFreeSubDimension(_nameModel)) THROW (InputException,wrongModelInSetSubDimensionEqual);
	_subDimensionEqual = iSubDimensionEqual;
}




}
