/***************************************************************************
                             SRC/mixmod/Utilities/Util.cpp  description
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
#include "mixmod/Utilities/Util.h"
#include "mixmod/Utilities/Random.h"
#include <ctype.h>
#include <set>

namespace XEM {

int VERBOSE = 0;
int MASSICCC = 0;
IoMode IOMODE = IoMode::NUMERIC;

void putDoubleInStream(std::ostream& output, double value, std::string appendChars)
{
	if (IOMODE == IoMode::BINARY) {
		uint64_t tmp;
		memcpy(&tmp, &value, sizeof(value));
		output << hex << tmp << appendChars << endl;
		//output << hexfloat <<_proba[i][k] << "\t"; //Writing is ok but reading isn't...
	}
	else
		output << value << appendChars << endl;
}

double getDoubleFromStream(std::istream& input)
{
	double value;
	if (IOMODE == IoMode::BINARY) {
		stringstream stream;
		uint64_t tmp;
		input >> hex >> tmp;
		//input >> hexfloat >> value; //Reading is not OK...
		memcpy(&value, &tmp, sizeof(tmp));
	}
	else
		input >> value;
	return value;
}

double powAndCheckIfNotNull(double a, double b, const Exception & errorType) {
	double res;
	res = pow(a, b);
	if (res == 0.0) throw errorType;
	return res;
}

//-----------------------
// return the nearest int64_t
//-----------------------
int64_t Round(double d) {
	int64_t res = (int64_t) (d + 0.5);
	return res;
}

//------------------------------------
// Convert big char of str in low char
//------------------------------------
void ConvertBigtoLowString(std::string & str) {
	for (unsigned int i = 0; i < str.length(); i++) {
		if (isupper(str[i]))
			str[i] = tolower(str[i]);
	}
}

//XEMCriterionNameToString
std::string ModelNameToString(const ModelName & modelName) {
	std::string res;

	switch (modelName) {
	case (UNKNOWN_MODEL_NAME):
		res = "UNKNOWN_MODEL_TYPE";
		break;
	// Gaussian models
	case (Gaussian_p_L_B):
		res = "Gaussian_p_L_B";
		break;
	case (Gaussian_p_Lk_B):
		res = "Gaussian_p_Lk_B";
		break;
	case (Gaussian_p_L_Bk):
		res = "Gaussian_p_L_Bk";
		break;
	case (Gaussian_p_Lk_Bk):
		res = "Gaussian_p_Lk_Bk";
		break;
	case (Gaussian_pk_L_B):
		res = "Gaussian_pk_L_B";
		break;
	case (Gaussian_pk_Lk_B):
		res = "Gaussian_pk_Lk_B";
		break;
	case (Gaussian_pk_L_Bk):
		res = "Gaussian_pk_L_Bk";
		break;
	case (Gaussian_pk_Lk_Bk):
		res = "Gaussian_pk_Lk_Bk";
		break;
	case (Gaussian_p_L_I):
		res = "Gaussian_p_L_I";
		break;
	case (Gaussian_p_Lk_I):
		res = "Gaussian_p_Lk_I";
		break;
	case (Gaussian_pk_L_I):
		res = "Gaussian_pk_L_I";
		break;
	case (Gaussian_pk_Lk_I):
		res = "Gaussian_pk_Lk_I";
		break;
	case (Gaussian_p_L_C):
		res = "Gaussian_p_L_C";
		break;
	case (Gaussian_p_Lk_C):
		res = "Gaussian_p_Lk_C";
		break;
	case (Gaussian_p_L_D_Ak_D):
		res = "Gaussian_p_L_D_Ak_D";
		break;
	case (Gaussian_p_Lk_D_Ak_D):
		res = "Gaussian_p_Lk_D_Ak_D";
		break;
	case (Gaussian_p_L_Dk_A_Dk):
		res = "Gaussian_p_L_Dk_A_Dk";
		break;
	case (Gaussian_p_Lk_Dk_A_Dk):
		res = "Gaussian_p_Lk_Dk_A_Dk";
		break;
	case (Gaussian_p_L_Ck):
		res = "Gaussian_p_L_Ck";
		break;
	case (Gaussian_p_Lk_Ck):
		res = "Gaussian_p_Lk_Ck";
		break;
	case (Gaussian_pk_L_C):
		res = "Gaussian_pk_L_C";
		break;
	case (Gaussian_pk_Lk_C):
		res = "Gaussian_pk_Lk_C";
		break;
	case (Gaussian_pk_L_D_Ak_D):
		res = "Gaussian_pk_L_D_Ak_D";
		break;
	case (Gaussian_pk_Lk_D_Ak_D):
		res = "Gaussian_pk_Lk_D_Ak_D";
		break;
	case (Gaussian_pk_L_Dk_A_Dk):
		res = "Gaussian_pk_L_Dk_A_Dk";
		break;
	case (Gaussian_pk_Lk_Dk_A_Dk):
		res = "Gaussian_pk_Lk_Dk_A_Dk";
		break;
	case (Gaussian_pk_L_Ck):
		res = "Gaussian_pk_L_Ck";
		break;
	case (Gaussian_pk_Lk_Ck):
		res = "Gaussian_pk_Lk_Ck";
		break;

	// Binary models
	case (Binary_p_E):
		res = "Binary_p_E";
		break;
	case (Binary_p_Ek):
		res = "Binary_p_Ek";
		break;
	case (Binary_p_Ej):
		res = "Binary_p_Ej";
		break;
	case (Binary_p_Ekj):
		res = "Binary_p_Ekj";
		break;
	case (Binary_p_Ekjh):
		res = "Binary_p_Ekjh";
		break;
	case (Binary_pk_E):
		res = "Binary_pk_E";
		break;
	case (Binary_pk_Ek):
		res = "Binary_pk_Ek";
		break;
	case (Binary_pk_Ej):
		res = "Binary_pk_Ej";
		break;
	case (Binary_pk_Ekj):
		res = "Binary_pk_Ekj";
		break;
	case (Binary_pk_Ekjh):
		res = "Binary_pk_Ekjh";
		break;

	case (Gaussian_HD_pk_AkjBkQkD):
		res = "Gaussian_HD_pk_AkjBkQkD";
		break;
	case (Gaussian_HD_pk_AkjBkQkDk):
		res = "Gaussian_HD_pk_AkjBkQkDk";
		break;
	case (Gaussian_HD_pk_AkjBQkD):
		res = "Gaussian_HD_pk_AkjBQkD";
		break;
	case (Gaussian_HD_pk_AjBkQkD):
		res = "Gaussian_HD_pk_AjBkQkD";
		break;
	case (Gaussian_HD_pk_AjBQkD):
		res = "Gaussian_HD_pk_AjBQkD";
		break;
	case (Gaussian_HD_pk_AkBkQkD):
		res = "Gaussian_HD_pk_AkBkQkD";
		break;
	case (Gaussian_HD_pk_AkBkQkDk):
		res = "Gaussian_HD_pk_AkBkQkDk";
		break;
	case (Gaussian_HD_pk_AkBQkD):
		res = "Gaussian_HD_pk_AkBQkD";
		break;
	case (Gaussian_HD_p_AkjBkQkD):
		res = "Gaussian_HD_p_AkjBkQkD";
		break;
	case (Gaussian_HD_p_AkjBkQkDk):
		res = "Gaussian_HD_p_AkjBkQkDk";
		break;
	case (Gaussian_HD_p_AkjBQkD):
		res = "Gaussian_HD_p_AkjBQkD";
		break;
	case (Gaussian_HD_p_AjBkQkD):
		res = "Gaussian_HD_p_AjBkQkD";
		break;
	case (Gaussian_HD_p_AjBQkD):
		res = "Gaussian_HD_p_AjBQkD";
		break;
	case (Gaussian_HD_p_AkBkQkD):
		res = "Gaussian_HD_p_AkBkQkD";
		break;
	case (Gaussian_HD_p_AkBkQkDk):
		res = "Gaussian_HD_p_AkBkQkDk";
		break;
	case (Gaussian_HD_p_AkBQkD):
		res = "Gaussian_HD_p_AkBQkD";
		break;
	case (Heterogeneous_pk_E_L_B):
		res = "Heterogeneous_pk_E_L_B";
		break;
	case (Heterogeneous_pk_E_Lk_B):
		res = "Heterogeneous_pk_E_Lk_B";
		break;
	case (Heterogeneous_pk_E_L_Bk):
		res = "Heterogeneous_pk_E_L_Bk";
		break;
	case (Heterogeneous_pk_E_Lk_Bk):
		res = "Heterogeneous_pk_E_Lk_Bk";
		break;
	case (Heterogeneous_pk_Ek_L_B):
		res = "Heterogeneous_pk_Ek_L_B";
		break;
	case (Heterogeneous_pk_Ek_Lk_B):
		res = "Heterogeneous_pk_Ek_Lk_B";
		break;
	case (Heterogeneous_pk_Ek_L_Bk):
		res = "Heterogeneous_pk_Ek_L_Bk";
		break;
	case (Heterogeneous_pk_Ek_Lk_Bk):
		res = "Heterogeneous_pk_Ek_Lk_Bk";
		break;
	case (Heterogeneous_pk_Ej_L_B):
		res = "Heterogeneous_pk_Ej_L_B";
		break;
	case (Heterogeneous_pk_Ej_Lk_B):
		res = "Heterogeneous_pk_Ej_Lk_B";
		break;
	case (Heterogeneous_pk_Ej_L_Bk):
		res = "Heterogeneous_pk_Ej_L_Bk";
		break;
	case (Heterogeneous_pk_Ej_Lk_Bk):
		res = "Heterogeneous_pk_Ej_Lk_Bk";
		break;
	case (Heterogeneous_pk_Ekj_L_B):
		res = "Heterogeneous_pk_Ekj_L_B";
		break;
	case (Heterogeneous_pk_Ekj_Lk_B):
		res = "Heterogeneous_pk_Ekj_Lk_B";
		break;
	case (Heterogeneous_pk_Ekj_L_Bk):
		res = "Heterogeneous_pk_Ekj_L_Bk";
		break;
	case (Heterogeneous_pk_Ekj_Lk_Bk):
		res = "Heterogeneous_pk_Ekj_Lk_Bk";
		break;
	case (Heterogeneous_pk_Ekjh_L_B):
		res = "Heterogeneous_pk_Ekjh_L_B";
		break;
	case (Heterogeneous_pk_Ekjh_Lk_B):
		res = "Heterogeneous_pk_Ekjh_Lk_B";
		break;
	case (Heterogeneous_pk_Ekjh_L_Bk):
		res = "Heterogeneous_pk_Ekjh_L_Bk";
		break;
	case (Heterogeneous_pk_Ekjh_Lk_Bk):
		res = "Heterogeneous_pk_Ekjh_Lk_Bk";
		break;
	case (Heterogeneous_p_E_L_B):
		res = "Heterogeneous_p_E_L_B";
		break;
	case (Heterogeneous_p_E_Lk_B):
		res = "Heterogeneous_p_E_Lk_B";
		break;
	case (Heterogeneous_p_E_L_Bk):
		res = "Heterogeneous_p_E_L_Bk";
		break;
	case (Heterogeneous_p_E_Lk_Bk):
		res = "Heterogeneous_p_E_Lk_Bk";
		break;
	case (Heterogeneous_p_Ek_L_B):
		res = "Heterogeneous_p_Ek_L_B";
		break;
	case (Heterogeneous_p_Ek_Lk_B):
		res = "Heterogeneous_p_Ek_Lk_B";
		break;
	case (Heterogeneous_p_Ek_L_Bk):
		res = "Heterogeneous_p_Ek_L_Bk";
		break;
	case (Heterogeneous_p_Ek_Lk_Bk):
		res = "Heterogeneous_p_Ek_Lk_Bk";
		break;
	case (Heterogeneous_p_Ej_L_B):
		res = "Heterogeneous_p_Ej_L_B";
		break;
	case (Heterogeneous_p_Ej_Lk_B):
		res = "Heterogeneous_p_Ej_Lk_B";
		break;
	case (Heterogeneous_p_Ej_L_Bk):
		res = "Heterogeneous_p_Ej_L_Bk";
		break;
	case (Heterogeneous_p_Ej_Lk_Bk):
		res = "Heterogeneous_p_Ej_Lk_Bk";
		break;
	case (Heterogeneous_p_Ekj_L_B):
		res = "Heterogeneous_p_Ekj_L_B";
		break;
	case (Heterogeneous_p_Ekj_Lk_B):
		res = "Heterogeneous_p_Ekj_Lk_B";
		break;
	case (Heterogeneous_p_Ekj_L_Bk):
		res = "Heterogeneous_p_Ekj_L_Bk";
		break;
	case (Heterogeneous_p_Ekj_Lk_Bk):
		res = "Heterogeneous_p_Ekj_Lk_Bk";
		break;
	case (Heterogeneous_p_Ekjh_L_B):
		res = "Heterogeneous_p_Ekjh_L_B";
		break;
	case (Heterogeneous_p_Ekjh_Lk_B):
		res = "Heterogeneous_p_Ekjh_Lk_B";
		break;
	case (Heterogeneous_p_Ekjh_L_Bk):
		res = "Heterogeneous_p_Ekjh_L_Bk";
		break;
	case (Heterogeneous_p_Ekjh_Lk_Bk):
		res = "Heterogeneous_p_Ekjh_Lk_Bk";
		break;
	default:
		THROW(InputException, wrongModelType);
	}

	return res;
}

// StringToModelName
ModelName StringToModelName(const std::string & strModelName) {
	ModelName res = UNKNOWN_MODEL_NAME;

	if (strModelName.compare("UNKNOWN_MODEL_NAME") == 0)
		res = UNKNOWN_MODEL_NAME;

	// Gaussian models
	if (strModelName.compare("Gaussian_p_L_B") == 0)
		res = Gaussian_p_L_B;
	if (strModelName.compare("Gaussian_p_Lk_B") == 0)
		res = Gaussian_p_Lk_B;
	if (strModelName.compare("Gaussian_p_L_Bk") == 0)
		res = Gaussian_p_L_Bk;
	if (strModelName.compare("Gaussian_p_Lk_Bk") == 0)
		res = Gaussian_p_Lk_Bk;
	if (strModelName.compare("Gaussian_pk_L_B") == 0)
		res = Gaussian_pk_L_B;
	if (strModelName.compare("Gaussian_pk_Lk_B") == 0)
		res = Gaussian_pk_Lk_B;
	if (strModelName.compare("Gaussian_pk_L_Bk") == 0)
		res = Gaussian_pk_L_Bk;
	if (strModelName.compare("Gaussian_pk_Lk_Bk") == 0)
		res = Gaussian_pk_Lk_Bk;
	if (strModelName.compare("Gaussian_p_L_I") == 0)
		res = Gaussian_p_L_I;
	if (strModelName.compare("Gaussian_p_Lk_I") == 0)
		res = Gaussian_p_Lk_I;
	if (strModelName.compare("Gaussian_pk_L_I") == 0)
		res = Gaussian_pk_L_I;
	if (strModelName.compare("Gaussian_pk_Lk_I") == 0)
		res = Gaussian_pk_Lk_I;
	if (strModelName.compare("Gaussian_p_L_C") == 0)
		res = Gaussian_p_L_C;
	if (strModelName.compare("Gaussian_p_Lk_C") == 0)
		res = Gaussian_p_Lk_C;
	if (strModelName.compare("Gaussian_p_L_D_Ak_D") == 0)
		res = Gaussian_p_L_D_Ak_D;
	if (strModelName.compare("Gaussian_p_Lk_D_Ak_D") == 0)
		res = Gaussian_p_Lk_D_Ak_D;
	if (strModelName.compare("Gaussian_p_L_Dk_A_Dk") == 0)
		res = Gaussian_p_L_Dk_A_Dk;
	if (strModelName.compare("Gaussian_p_Lk_Dk_A_Dk") == 0)
		res = Gaussian_p_Lk_Dk_A_Dk;
	if (strModelName.compare("Gaussian_p_L_Ck") == 0)
		res = Gaussian_p_L_Ck;
	if (strModelName.compare("Gaussian_p_Lk_Ck") == 0)
		res = Gaussian_p_Lk_Ck;
	if (strModelName.compare("Gaussian_pk_L_C") == 0)
		res = Gaussian_pk_L_C;
	if (strModelName.compare("Gaussian_pk_Lk_C") == 0)
		res = Gaussian_pk_Lk_C;
	if (strModelName.compare("Gaussian_pk_L_D_Ak_D") == 0)
		res = Gaussian_pk_L_D_Ak_D;
	if (strModelName.compare("Gaussian_pk_Lk_D_Ak_D") == 0)
		res = Gaussian_pk_Lk_D_Ak_D;
	if (strModelName.compare("Gaussian_pk_L_Dk_A_Dk") == 0)
		res = Gaussian_pk_L_Dk_A_Dk;
	if (strModelName.compare("Gaussian_pk_Lk_Dk_A_Dk") == 0)
		res = Gaussian_pk_Lk_Dk_A_Dk;
	if (strModelName.compare("Gaussian_pk_L_Ck") == 0)
		res = Gaussian_pk_L_Ck;
	if (strModelName.compare("Gaussian_pk_Lk_Ck") == 0)
		res = Gaussian_pk_Lk_Ck;

	// Binary models
	if (strModelName.compare("Binary_p_E") == 0)
		res = Binary_p_E;
	if (strModelName.compare("Binary_p_Ek") == 0)
		res = Binary_p_Ek;
	if (strModelName.compare("Binary_p_Ej") == 0)
		res = Binary_p_Ej;
	if (strModelName.compare("Binary_p_Ekj") == 0)
		res = Binary_p_Ekj;
	if (strModelName.compare("Binary_p_Ekjh") == 0)
		res = Binary_p_Ekjh;
	if (strModelName.compare("Binary_pk_E") == 0)
		res = Binary_pk_E;
	if (strModelName.compare("Binary_pk_Ek") == 0)
		res = Binary_pk_Ek;
	if (strModelName.compare("Binary_pk_Ej") == 0)
		res = Binary_pk_Ej;
	if (strModelName.compare("Binary_pk_Ekj") == 0)
		res = Binary_pk_Ekj;
	if (strModelName.compare("Binary_pk_Ekjh") == 0)
		res = Binary_pk_Ekjh;

	//HDModel
	if (strModelName.compare("Gaussian_HD_pk_AkjBkQkD") == 0)
		res = Gaussian_HD_pk_AkjBkQkD;
	if (strModelName.compare("Gaussian_HD_pk_AkjBkQkDk") == 0)
		res = Gaussian_HD_pk_AkjBkQkDk;
	if (strModelName.compare("Gaussian_HD_pk_AkjBQkD") == 0)
		res = Gaussian_HD_pk_AkjBQkD;
	if (strModelName.compare("Gaussian_HD_pk_AjBkQkD") == 0)
		res = Gaussian_HD_pk_AjBkQkD;
	if (strModelName.compare("Gaussian_HD_pk_AjBQkD") == 0)
		res = Gaussian_HD_pk_AjBQkD;
	if (strModelName.compare("Gaussian_HD_pk_AkBkQkD") == 0)
		res = Gaussian_HD_pk_AkBkQkD;
	if (strModelName.compare("Gaussian_HD_pk_AkBkQkDk") == 0)
		res = Gaussian_HD_pk_AkBkQkDk;
	if (strModelName.compare("Gaussian_HD_pk_AkBQkD") == 0)
		res = Gaussian_HD_pk_AkBQkD;
	if (strModelName.compare("Gaussian_HD_p_AkjBkQkD") == 0)
		res = Gaussian_HD_p_AkjBkQkD;
	if (strModelName.compare("Gaussian_HD_p_AkjBkQkDk") == 0)
		res = Gaussian_HD_p_AkjBkQkDk;
	if (strModelName.compare("Gaussian_HD_p_AkjBQkD") == 0)
		res = Gaussian_HD_p_AkjBQkD;
	if (strModelName.compare("Gaussian_HD_p_AjBkQkD") == 0)
		res = Gaussian_HD_p_AjBkQkD;
	if (strModelName.compare("Gaussian_HD_p_AjBQkD") == 0)
		res = Gaussian_HD_p_AjBQkD;
	if (strModelName.compare("Gaussian_HD_p_AkBkQkD") == 0)
		res = Gaussian_HD_p_AkBkQkD;
	if (strModelName.compare("Gaussian_HD_p_AkBkQkDk") == 0)
		res = Gaussian_HD_p_AkBkQkDk;
	if (strModelName.compare("Gaussian_HD_p_AkBQkD") == 0)
		res = Gaussian_HD_p_AkBQkD;
	if (strModelName.compare("Heterogeneous_pk_E_L_B") == 0)
		res = Heterogeneous_pk_E_L_B;

	if (strModelName.compare("Heterogeneous_pk_E_Lk_B") == 0)
		res = Heterogeneous_pk_E_Lk_B;

	if (strModelName.compare("Heterogeneous_pk_E_L_Bk") == 0)
		res = Heterogeneous_pk_E_L_Bk;

	if (strModelName.compare("Heterogeneous_pk_E_Lk_Bk") == 0)
		res = Heterogeneous_pk_E_Lk_Bk;

	if (strModelName.compare("Heterogeneous_pk_Ek_L_B") == 0)
		res = Heterogeneous_pk_Ek_L_B;

	if (strModelName.compare("Heterogeneous_pk_Ek_Lk_B") == 0)
		res = Heterogeneous_pk_Ek_Lk_B;

	if (strModelName.compare("Heterogeneous_pk_Ek_L_Bk") == 0)
		res = Heterogeneous_pk_Ek_L_Bk;

	if (strModelName.compare("Heterogeneous_pk_Ek_Lk_Bk") == 0)
		res = Heterogeneous_pk_Ek_Lk_Bk;

	if (strModelName.compare("Heterogeneous_pk_Ej_L_B") == 0)
		res = Heterogeneous_pk_Ej_L_B;

	if (strModelName.compare("Heterogeneous_pk_Ej_Lk_B") == 0)
		res = Heterogeneous_pk_Ej_Lk_B;

	if (strModelName.compare("Heterogeneous_pk_Ej_L_Bk") == 0)
		res = Heterogeneous_pk_Ej_L_Bk;

	if (strModelName.compare("Heterogeneous_pk_Ej_Lk_Bk") == 0)
		res = Heterogeneous_pk_Ej_Lk_Bk;

	if (strModelName.compare("Heterogeneous_pk_Ekj_L_B") == 0)
		res = Heterogeneous_pk_Ekj_L_B;

	if (strModelName.compare("Heterogeneous_pk_Ekj_Lk_B") == 0)
		res = Heterogeneous_pk_Ekj_Lk_B;

	if (strModelName.compare("Heterogeneous_pk_Ekj_L_Bk") == 0)
		res = Heterogeneous_pk_Ekj_L_Bk;

	if (strModelName.compare("Heterogeneous_pk_Ekj_Lk_Bk") == 0)
		res = Heterogeneous_pk_Ekj_Lk_Bk;

	if (strModelName.compare("Heterogeneous_pk_Ekjh_L_B") == 0)
		res = Heterogeneous_pk_Ekjh_L_B;

	if (strModelName.compare("Heterogeneous_pk_Ekjh_Lk_B") == 0)
		res = Heterogeneous_pk_Ekjh_Lk_B;

	if (strModelName.compare("Heterogeneous_pk_Ekjh_L_Bk") == 0)
		res = Heterogeneous_pk_Ekjh_L_Bk;

	if (strModelName.compare("Heterogeneous_pk_Ekjh_Lk_Bk") == 0)
		res = Heterogeneous_pk_Ekjh_Lk_Bk;

	if (strModelName.compare("Heterogeneous_p_E_L_B") == 0)
		res = Heterogeneous_p_E_L_B;

	if (strModelName.compare("Heterogeneous_p_E_Lk_B") == 0)
		res = Heterogeneous_p_E_Lk_B;

	if (strModelName.compare("Heterogeneous_p_E_L_Bk") == 0)
		res = Heterogeneous_p_E_L_Bk;

	if (strModelName.compare("Heterogeneous_p_E_Lk_Bk") == 0)
		res = Heterogeneous_p_E_Lk_Bk;

	if (strModelName.compare("Heterogeneous_p_Ek_L_B") == 0)
		res = Heterogeneous_p_Ek_L_B;

	if (strModelName.compare("Heterogeneous_p_Ek_Lk_B") == 0)
		res = Heterogeneous_p_Ek_Lk_B;

	if (strModelName.compare("Heterogeneous_p_Ek_L_Bk") == 0)
		res = Heterogeneous_p_Ek_L_Bk;

	if (strModelName.compare("Heterogeneous_p_Ek_Lk_Bk") == 0)
		res = Heterogeneous_p_Ek_Lk_Bk;

	if (strModelName.compare("Heterogeneous_p_Ej_L_B") == 0)
		res = Heterogeneous_p_Ej_L_B;

	if (strModelName.compare("Heterogeneous_p_Ej_Lk_B") == 0)
		res = Heterogeneous_p_Ej_Lk_B;

	if (strModelName.compare("Heterogeneous_p_Ej_L_Bk") == 0)
		res = Heterogeneous_p_Ej_L_Bk;

	if (strModelName.compare("Heterogeneous_p_Ej_Lk_Bk") == 0)
		res = Heterogeneous_p_Ej_Lk_Bk;

	if (strModelName.compare("Heterogeneous_p_Ekj_L_B") == 0)
		res = Heterogeneous_p_Ekj_L_B;

	if (strModelName.compare("Heterogeneous_p_Ekj_Lk_B") == 0)
		res = Heterogeneous_p_Ekj_Lk_B;

	if (strModelName.compare("Heterogeneous_p_Ekj_L_Bk") == 0)
		res = Heterogeneous_p_Ekj_L_Bk;

	if (strModelName.compare("Heterogeneous_p_Ekj_Lk_Bk") == 0)
		res = Heterogeneous_p_Ekj_Lk_Bk;

	if (strModelName.compare("Heterogeneous_p_Ekjh_L_B") == 0)
		res = Heterogeneous_p_Ekjh_L_B;

	if (strModelName.compare("Heterogeneous_p_Ekjh_Lk_B") == 0)
		res = Heterogeneous_p_Ekjh_Lk_B;

	if (strModelName.compare("Heterogeneous_p_Ekjh_L_Bk") == 0)
		res = Heterogeneous_p_Ekjh_L_Bk;

	if (strModelName.compare("Heterogeneous_p_Ekjh_Lk_Bk") == 0)
		res = Heterogeneous_p_Ekjh_Lk_Bk;

	//Heterogeneous models
	if (res == UNKNOWN_MODEL_NAME)
		THROW(InputException, wrongModelType);

	return res;
}

ModelName getHeterogeneousModelName(const ModelName binaryName, const ModelName gaussianName) {
	string g_model = ModelNameToString(gaussianName);
	string b_model = ModelNameToString(binaryName);
	if (hasFreeProportion(gaussianName) && hasFreeProportion(binaryName)) {
		const string h_model = "Heterogeneous_" + b_model.substr(7) + g_model.substr(11);
		return StringToModelName(h_model);
	}
	else   if (!hasFreeProportion(gaussianName)&&!hasFreeProportion(binaryName)) {
		const string h_model = "Heterogeneous_" + b_model.substr(7) + g_model.substr(10);
		return StringToModelName(h_model);
	}

	THROW(InputException, badInputType);
}

ModelName getGaussianModelNamefromHeterogeneous(const ModelName HeterogeneousName) {
	ModelName res;
	switch (HeterogeneousName) {
	case (Heterogeneous_pk_E_L_B):
		res = Gaussian_pk_L_B;
		break;
	case (Heterogeneous_pk_E_Lk_B):
		res = Gaussian_pk_Lk_B;
		break;
	case (Heterogeneous_pk_E_L_Bk):
		res = Gaussian_pk_L_Bk;
		break;
	case (Heterogeneous_pk_E_Lk_Bk):
		res = Gaussian_pk_Lk_Bk;
		break;
	case (Heterogeneous_pk_Ek_L_B):
		res = Gaussian_pk_L_B;
		break;
	case (Heterogeneous_pk_Ek_Lk_B):
		res = Gaussian_pk_Lk_B;
		break;
	case (Heterogeneous_pk_Ek_L_Bk):
		res = Gaussian_pk_L_Bk;
		break;
	case (Heterogeneous_pk_Ek_Lk_Bk):
		res = Gaussian_pk_Lk_Bk;
		break;
	case (Heterogeneous_pk_Ej_L_B):
		res = Gaussian_pk_L_B;
		break;
	case (Heterogeneous_pk_Ej_Lk_B):
		res = Gaussian_pk_Lk_B;
		break;
	case (Heterogeneous_pk_Ej_L_Bk):
		res = Gaussian_pk_L_Bk;
		break;
	case (Heterogeneous_pk_Ej_Lk_Bk):
		res = Gaussian_pk_Lk_Bk;
		break;
	case (Heterogeneous_pk_Ekj_L_B):
		res = Gaussian_pk_L_B;
		break;
	case (Heterogeneous_pk_Ekj_Lk_B):
		res = Gaussian_pk_Lk_B;
		break;
	case (Heterogeneous_pk_Ekj_L_Bk):
		res = Gaussian_pk_L_Bk;
		break;
	case (Heterogeneous_pk_Ekj_Lk_Bk):
		res = Gaussian_pk_Lk_Bk;
		break;
	case (Heterogeneous_pk_Ekjh_L_B):
		res = Gaussian_pk_L_B;
		break;
	case (Heterogeneous_pk_Ekjh_Lk_B):
		res = Gaussian_pk_Lk_B;
		break;
	case (Heterogeneous_pk_Ekjh_L_Bk):
		res = Gaussian_pk_L_Bk;
		break;
	case (Heterogeneous_pk_Ekjh_Lk_Bk):
		res = Gaussian_pk_Lk_Bk;
		break;
	case (Heterogeneous_p_E_L_B):
		res = Gaussian_p_L_B;
		break;
	case (Heterogeneous_p_E_Lk_B):
		res = Gaussian_p_Lk_B;
		break;
	case (Heterogeneous_p_E_L_Bk):
		res = Gaussian_p_L_Bk;
		break;
	case (Heterogeneous_p_E_Lk_Bk):
		res = Gaussian_p_Lk_Bk;
		break;
	case (Heterogeneous_p_Ek_L_B):
		res = Gaussian_p_L_B;
		break;
	case (Heterogeneous_p_Ek_Lk_B):
		res = Gaussian_p_Lk_B;
		break;
	case (Heterogeneous_p_Ek_L_Bk):
		res = Gaussian_p_L_Bk;
		break;
	case (Heterogeneous_p_Ek_Lk_Bk):
		res = Gaussian_p_Lk_Bk;
		break;
	case (Heterogeneous_p_Ej_L_B):
		res = Gaussian_p_L_B;
		break;
	case (Heterogeneous_p_Ej_Lk_B):
		res = Gaussian_p_Lk_B;
		break;
	case (Heterogeneous_p_Ej_L_Bk):
		res = Gaussian_p_L_Bk;
		break;
	case (Heterogeneous_p_Ej_Lk_Bk):
		res = Gaussian_p_Lk_Bk;
		break;
	case (Heterogeneous_p_Ekj_L_B):
		res = Gaussian_p_L_B;
		break;
	case (Heterogeneous_p_Ekj_Lk_B):
		res = Gaussian_p_Lk_B;
		break;
	case (Heterogeneous_p_Ekj_L_Bk):
		res = Gaussian_p_L_Bk;
		break;
	case (Heterogeneous_p_Ekj_Lk_Bk):
		res = Gaussian_p_Lk_Bk;
		break;
	case (Heterogeneous_p_Ekjh_L_B):
		res = Gaussian_p_L_B;
		break;
	case (Heterogeneous_p_Ekjh_Lk_B):
		res = Gaussian_p_Lk_B;
		break;
	case (Heterogeneous_p_Ekjh_L_Bk):
		res = Gaussian_p_L_Bk;
		break;
	case (Heterogeneous_p_Ekjh_Lk_Bk):
		res = Gaussian_p_Lk_Bk;
		break;
	default:
		THROW(InputException, badInputType);
	}

	return res;
}

ModelName getBinaryModelNamefromHeterogeneous(const ModelName HeterogeneousName) {
	ModelName res;
	switch (HeterogeneousName) {
	case (Heterogeneous_pk_E_L_B):
		res = Binary_pk_E;
		break;
	case (Heterogeneous_pk_E_Lk_B):
		res = Binary_pk_E;
		break;
	case (Heterogeneous_pk_E_L_Bk):
		res = Binary_pk_E;
		break;
	case (Heterogeneous_pk_E_Lk_Bk):
		res = Binary_pk_E;
		break;
	case (Heterogeneous_pk_Ek_L_B):
		res = Binary_pk_Ek;
		break;
	case (Heterogeneous_pk_Ek_Lk_B):
		res = Binary_pk_Ek;
		break;
	case (Heterogeneous_pk_Ek_L_Bk):
		res = Binary_pk_Ek;
		break;
	case (Heterogeneous_pk_Ek_Lk_Bk):
		res = Binary_pk_Ek;
		break;
	case (Heterogeneous_pk_Ej_L_B):
		res = Binary_pk_Ej;
		break;
	case (Heterogeneous_pk_Ej_Lk_B):
		res = Binary_pk_Ej;
		break;
	case (Heterogeneous_pk_Ej_L_Bk):
		res = Binary_pk_Ej;
		break;
	case (Heterogeneous_pk_Ej_Lk_Bk):
		res = Binary_pk_Ej;
		break;
	case (Heterogeneous_pk_Ekj_L_B):
		res = Binary_pk_Ekj;
		break;
	case (Heterogeneous_pk_Ekj_Lk_B):
		res = Binary_pk_Ekj;
		break;
	case (Heterogeneous_pk_Ekj_L_Bk):
		res = Binary_pk_Ekj;
		break;
	case (Heterogeneous_pk_Ekj_Lk_Bk):
		res = Binary_pk_Ekj;
		break;
	case (Heterogeneous_pk_Ekjh_L_B):
		res = Binary_pk_Ekjh;
		break;
	case (Heterogeneous_pk_Ekjh_Lk_B):
		res = Binary_pk_Ekjh;
		break;
	case (Heterogeneous_pk_Ekjh_L_Bk):
		res = Binary_pk_Ekjh;
		break;
	case (Heterogeneous_pk_Ekjh_Lk_Bk):
		res = Binary_pk_Ekjh;
		break;
	case (Heterogeneous_p_E_L_B):
		res = Binary_p_E;
		break;
	case (Heterogeneous_p_E_Lk_B):
		res = Binary_p_E;
		break;
	case (Heterogeneous_p_E_L_Bk):
		res = Binary_p_E;
		break;
	case (Heterogeneous_p_E_Lk_Bk):
		res = Binary_p_E;
		break;
	case (Heterogeneous_p_Ek_L_B):
		res = Binary_p_Ek;
		break;
	case (Heterogeneous_p_Ek_Lk_B):
		res = Binary_p_Ek;
		break;
	case (Heterogeneous_p_Ek_L_Bk):
		res = Binary_p_Ek;
		break;
	case (Heterogeneous_p_Ek_Lk_Bk):
		res = Binary_p_Ek;
		break;
	case (Heterogeneous_p_Ej_L_B):
		res = Binary_p_Ej;
		break;
	case (Heterogeneous_p_Ej_Lk_B):
		res = Binary_p_Ej;
		break;
	case (Heterogeneous_p_Ej_L_Bk):
		res = Binary_p_Ej;
		break;
	case (Heterogeneous_p_Ej_Lk_Bk):
		res = Binary_p_Ej;
		break;
	case (Heterogeneous_p_Ekj_L_B):
		res = Binary_p_Ekj;
		break;
	case (Heterogeneous_p_Ekj_Lk_B):
		res = Binary_p_Ekj;
		break;
	case (Heterogeneous_p_Ekj_L_Bk):
		res = Binary_p_Ekj;
		break;
	case (Heterogeneous_p_Ekj_Lk_Bk):
		res = Binary_p_Ekj;
		break;
	case (Heterogeneous_p_Ekjh_L_B):
		res = Binary_p_Ekjh;
		break;
	case (Heterogeneous_p_Ekjh_Lk_B):
		res = Binary_p_Ekjh;
		break;
	case (Heterogeneous_p_Ekjh_L_Bk):
		res = Binary_p_Ekjh;
		break;
	case (Heterogeneous_p_Ekjh_Lk_Bk):
		res = Binary_p_Ekjh;
		break;
	default:
		THROW(InputException, badInputType);
	}

	return res;
}
// edit modelName (debug)
void edit(const ModelName & modelName) {
	cout << ModelNameToString(modelName);
}

bool isKeyword(std::string& name) {
	// List of keywords defined as a static variable, since its only usage is here.
	static std::set<std::string> keywords;
	keywords.insert("NbLines");
	keywords.insert("PbDimension");
	keywords.insert("NbNbCluster");
	keywords.insert("ListNbCluster");
	keywords.insert("NbModality");
	keywords.insert("NbCriterion");
	keywords.insert("ListCriterion");
	keywords.insert("NbModel");
	keywords.insert("ListModel");
	keywords.insert("subDimensionEqual");
	keywords.insert("subDimensionFree");
	keywords.insert("NbStrategy");
	keywords.insert("InitType");
	keywords.insert("InitFile");
	keywords.insert("NbAlgorithm");
	keywords.insert("Algorithm");
	keywords.insert("PartitionFile");
	keywords.insert("DataFile");
	keywords.insert("WeightFile");
	keywords.insert("NbCVBlocks");
	keywords.insert("CVinitBlocks");
	keywords.insert("NbDCVBlocks");
	keywords.insert("DCVinitBlocks");
	keywords.insert("SizeKeyword");
	return (keywords.find(name) != keywords.end());
}

//CriterionNameToString
std::string CriterionNameToString(const CriterionName & criterionName) {
	std::string res;
	switch (criterionName) {
	case UNKNOWN_CRITERION_NAME:
		res = "UNKNOWN_CRITERION_NAME" ;
		break ;
	case BIC:
		res = "BIC";
		break;
	case ICL:
		res = "ICL";
		break;
	case NEC:
		res = "NEC";
		break;
	case CV:
		res = "CV";
		break;
	case DCV:
		res = "DCV";
		break;
	}
	return res;
}

// StringtoCriterionName
CriterionName StringtoCriterionName(const std::string & str) {
	CriterionName res = UNKNOWN_CRITERION_NAME;
	if (str.compare("UNKNOWN_CRITERION_NAME") == 0)
		res = UNKNOWN_CRITERION_NAME;
	if (str.compare("BIC") == 0)
		res = BIC;
	if (str.compare("ICL") == 0)
		res = ICL;
	if (str.compare("NEC") == 0)
		res = NEC;
	if (str.compare("CV") == 0)
		res = CV;
	if (str.compare("DCV") == 0)
		res = DCV;

	if (res == UNKNOWN_CRITERION_NAME)
		THROW(InputException, badCriterion);

	return res;
}

// edit CriterionName (debug)
void edit(const CriterionName & criterionName) {
	cout << CriterionNameToString(criterionName);
}

//AlgoNameToString
std::string AlgoNameToString(const AlgoName & typeAlgo) {
	std::string res;
	switch (typeAlgo) {
	case UNKNOWN_ALGO_NAME:
		res = "UNKNOWN_ALGO_NAME";
		break;
	case EM:
		res = "EM";
		break;
	case CEM:
		res = "CEM";
		break;
	case SEM:
		res = "SEM";
		break;
	case MAP:
		res = "MAP";
		break;
	case M:
		res = "M";
		break;
	}
	return res;
}

//StringToAlgoName
AlgoName StringToAlgoName(const std::string & str) {
	AlgoName res = UNKNOWN_ALGO_NAME;
	if (str.compare("UNKNOWN_ALGO_NAME") == 0)
		res = UNKNOWN_ALGO_NAME;
	if (str.compare("EM") == 0)
		res = EM;
	if (str.compare("CEM") == 0)
		res = CEM;
	if (str.compare("SEM") == 0)
		res = SEM;
	if (str.compare("MAP") == 0)
		res = MAP;
	if (str.compare("M") == 0)
		res = M;
	if (res == UNKNOWN_ALGO_NAME)
		THROW(InputException, badAlgo);

	return res;
}

//AlgoNameToString
std::string AlgoStopNameToString(const AlgoStopName & algoStopName) {
	std::string res;
	switch (algoStopName) {
	case NO_STOP_NAME:
		res = "NO_STOP_NAME";
		break;
	case NBITERATION:
		res = "NBITERATION";
		break;
	case EPSILON:
		res = "EPSILON";
		break;
	case NBITERATION_EPSILON:
		res = "NBITERATION_EPSILON";
		break;
	}
	return res;
}

//StringToAlgoStopName
AlgoStopName StringToAlgoStopName(const std::string & str) {
	AlgoStopName res = NO_STOP_NAME;
	if (str.compare("NO_STOP_NAME") == 0)
		res = NO_STOP_NAME;
	if (str.compare("NBITERATION") == 0)
		res = NBITERATION;
	if (str.compare("EPSILON") == 0)
		res = EPSILON;
	if (str.compare("NBITERATION_EPSILON") == 0)
		res = NBITERATION_EPSILON;
	if (res == NO_STOP_NAME)
		THROW(InputException, badAlgoStop);
	return res;
}

// edit typeAlgo (debug)
void edit(const AlgoName & typeAlgo) {
	cout << AlgoNameToString(typeAlgo);
}

//FormatNumericToString
std::string FormatNumericFileToString(const FormatNumeric::FormatNumericFile & formatNumericFile) {
	std::string res;
	switch (formatNumericFile) {
	case FormatNumeric::txt:
		res = "txt";
		break;
	case FormatNumeric::hdf5:
		res = "hdf5";
		break;
	case FormatNumeric::XML:
		res = "XML";
		break;
	}
	return res;
}

//StringToFormatFile
FormatNumeric::FormatNumericFile StringToFormatNumericFile(const std::string & strFormatNumericFile) {
	FormatNumeric::FormatNumericFile res ;
	if (strFormatNumericFile.compare("txt") == 0) {
		res = FormatNumeric::txt;
	}
	else if (strFormatNumericFile.compare("hdf5") == 0) {
		res = FormatNumeric::hdf5;
	}
	else if (strFormatNumericFile.compare("XML") == 0) {
		res = FormatNumeric::XML;
	}
	else {
		THROW(OtherException, badFormat);
	}
	return res;
}

//TypePartitionToString
std::string TypePartitionToString(const TypePartition::TypePartition & typePartition) {
	std::string res;
	switch (typePartition) {
	case TypePartition::UNKNOWN_PARTITION:
		res = "UNKNOWN_PARTITION";
		break;
	case TypePartition::label:
		res = "label";
		break;
	case TypePartition::partition:
		res = "partition";
		break;
	}
	return res;
}

//StringToTypePartition
TypePartition::TypePartition StringToTypePartition(const std::string & strTypePartition) {
	TypePartition::TypePartition res = TypePartition::UNKNOWN_PARTITION;
	if (strTypePartition.compare("UNKNOWN_PARTITION") == 0)
		res = TypePartition::UNKNOWN_PARTITION;
	if (strTypePartition.compare("label") == 0)
		res = TypePartition::label;
	if (strTypePartition.compare("partition") == 0)
		res = TypePartition::partition;

	return res;
}

// StrategyInitNameToString
std::string StrategyInitNameToString(const StrategyInitName & strategyInitName) {
	std::string res;
	switch (strategyInitName) {
	case RANDOM:
		res = "RANDOM";
		break;
	case CEM_INIT:
		res = "CEM_INIT";
		break;
	case SEM_MAX:
		res = "SEM_MAX";
		break;
	case SMALL_EM:
		res = "SMALL_EM";
		break;
	case USER:
		res = "USER";
		break;
	case USER_PARTITION:
		res = "USER_PARTITION";
		break;
	}
	return res;
}
// StrategyInitNameToString for 3rd party env (R, Py, etc.)
std::string StrategyInitNameToStringApp(const StrategyInitName & strategyInitName) {
	std::string res;
	switch (strategyInitName) {
	case RANDOM:
		res = "random";
		break;
	case CEM_INIT:
		res = "CEM";
		break;
	case SEM_MAX:
		res = "SEMMax";
		break;
	case SMALL_EM:
		res = "smallEM";
		break;
	case USER:
		res = "parameter";
		break;
	case USER_PARTITION:
		res = "label";
		break;
	}
	return res;
}

//StringToStrategyInitName
StrategyInitName StringToStrategyInitName(const std::string & str) {

	StrategyInitName res;
	if (str.compare("RANDOM") == 0)
		res = RANDOM;
	if (str.compare("CEM_INIT") == 0)
		res = CEM_INIT;
	if (str.compare("SEM_MAX") == 0)
		res = SEM_MAX;
	if (str.compare("SMALL_EM") == 0)
		res = SMALL_EM;
	if (str.compare("PARAMETER") == 0)
		res = USER;
	if (str.compare("PARTITION") == 0)
		res = USER_PARTITION;

	return res;
}

// edit strategyInitName (debug)
void edit(const StrategyInitName & strategyInitName) {
	cout << StrategyInitNameToString(strategyInitName);
}

//edit algoStopName (debug)
void edit(const AlgoStopName & algoStopName) {
	cout << AlgoStopNameToString(algoStopName);
}

///printAlgoType
void printTypeAlgo(std::ostream & flux, const AlgoName & typeAlgo) {
	if (typeAlgo == EM)
		flux << "EM" << endl;
	else if (typeAlgo == CEM)
		flux << "CEM" << endl;
	else if (typeAlgo == SEM)
		flux << "SEM" << endl;
	else if (typeAlgo == MAP)
		flux << "MAP" << endl;
	else if (typeAlgo == M)
		flux << "M" << endl;
}

//-------------------------------
// is modelName a spherical Model
//-------------------------------
bool isSpherical(ModelName modelName) {
	bool res = false;
	if ((modelName == Gaussian_p_L_I)
			|| (modelName == Gaussian_p_Lk_I)
			|| (modelName == Gaussian_pk_L_I)
			|| (modelName == Gaussian_pk_Lk_I)) {
		res = true;
	}
	return res;
}

//------------------------------
// is modelName a diagonal Model
//------------------------------
bool isDiagonal(ModelName modelName) {
	bool res = false;
	if (   (modelName == Gaussian_p_L_B)
			|| (modelName == Gaussian_p_Lk_B)
			|| (modelName == Gaussian_p_L_Bk)
			|| (modelName == Gaussian_p_Lk_Bk)
			|| (modelName == Gaussian_pk_L_B)
			|| (modelName == Gaussian_pk_Lk_B)
			|| (modelName == Gaussian_pk_L_Bk)
			|| (modelName == Gaussian_pk_Lk_Bk)
			) {
		res = true;
	}
	return res;
}

//-----------------------------
// is modelName a general Model
//-----------------------------
bool isGeneral(ModelName modelName) {
	bool res = false;
	if ((modelName == Gaussian_p_L_C) ||
			(modelName == Gaussian_p_Lk_C) ||
			(modelName == Gaussian_p_L_D_Ak_D) ||
			(modelName == Gaussian_p_Lk_D_Ak_D) ||
			(modelName == Gaussian_p_L_Dk_A_Dk) ||
			(modelName == Gaussian_p_Lk_Dk_A_Dk) ||
			(modelName == Gaussian_p_L_Ck) ||
			(modelName == Gaussian_p_Lk_Ck) ||
			(modelName == Gaussian_pk_L_C) ||
			(modelName == Gaussian_pk_Lk_C) ||
			(modelName == Gaussian_pk_L_D_Ak_D) ||
			(modelName == Gaussian_pk_Lk_D_Ak_D) ||
			(modelName == Gaussian_pk_L_Dk_A_Dk) ||
			(modelName == Gaussian_pk_Lk_Dk_A_Dk) ||
			(modelName == Gaussian_pk_L_Ck) ||
			(modelName == Gaussian_pk_Lk_Ck)) {
		res =   true;
	}
	return res;
}

//-----------------------------------------
// is modelName a EDDA (Classical Gaussian)
//-----------------------------------------
bool isEDDA(ModelName modelName) {
	bool res = false;
	if (isSpherical(modelName) || isDiagonal(modelName) || isGeneral(modelName)) {
		res = true;
	}
	return res;
}

//---------------
// HD (HD or HDk)
//---------------
bool isHD(ModelName modelName) {
	bool res = false;
	if (  (modelName == Gaussian_HD_p_AkjBkQkDk)
			|| (modelName == Gaussian_HD_p_AkBkQkDk)
			|| (modelName == Gaussian_HD_p_AkjBkQkD)
			|| (modelName == Gaussian_HD_p_AjBkQkD)
			|| (modelName == Gaussian_HD_p_AkjBQkD)
			|| (modelName == Gaussian_HD_p_AjBQkD)
			|| (modelName == Gaussian_HD_p_AkBkQkD)
			|| (modelName == Gaussian_HD_p_AkBQkD)
			|| (modelName == Gaussian_HD_pk_AkjBkQkDk)
			|| (modelName == Gaussian_HD_pk_AkBkQkDk)
			|| (modelName == Gaussian_HD_pk_AkjBkQkD)
			|| (modelName == Gaussian_HD_pk_AjBkQkD)
			|| (modelName == Gaussian_HD_pk_AkjBQkD)
			|| (modelName == Gaussian_HD_pk_AjBQkD)
			|| (modelName == Gaussian_HD_pk_AkBkQkD)
			|| (modelName == Gaussian_HD_pk_AkBQkD)) {
		res = true;
	}
	return res;
}

bool isFreeSubDimension(ModelName modelName) {
	bool res = false;
	if (   (modelName == Gaussian_HD_p_AkjBkQkDk)
			|| (modelName == Gaussian_HD_p_AkBkQkDk)
			|| (modelName == Gaussian_HD_p_AkjBkQkD)
			|| (modelName == Gaussian_HD_p_AjBkQkD)
			|| (modelName == Gaussian_HD_p_AkBkQkD)
			|| (modelName == Gaussian_HD_pk_AkjBkQkDk)
			|| (modelName == Gaussian_HD_pk_AkBkQkDk)
			|| (modelName == Gaussian_HD_pk_AkjBkQkD)
			|| (modelName == Gaussian_HD_pk_AjBkQkD)
			|| (modelName == Gaussian_HD_pk_AkBkQkD)
		) {
		res = true;
	}
	return res;
}

bool isBinary(ModelName modelName) {
	bool res = false;
	if ( (modelName == Binary_p_E) ||
			(modelName == Binary_pk_E) ||
			(modelName == Binary_p_Ej) ||
			(modelName == Binary_pk_Ej) ||
			(modelName == Binary_p_Ek) ||
			(modelName == Binary_pk_Ek) ||
			(modelName == Binary_p_Ekj) ||
			(modelName == Binary_pk_Ekj) ||
			(modelName == Binary_p_Ekjh) ||
			(modelName == Binary_pk_Ekjh) ) {
		res = true;
	}
	return res;
}

bool isHeterogeneous(ModelName modelname) {
	if (
			modelname == Heterogeneous_p_E_L_B ||
			modelname == Heterogeneous_p_E_L_Bk ||
			modelname == Heterogeneous_p_E_Lk_B ||
			modelname == Heterogeneous_p_E_Lk_Bk ||
			modelname == Heterogeneous_p_Ej_L_B ||
			modelname == Heterogeneous_p_Ej_L_Bk ||
			modelname == Heterogeneous_p_Ej_Lk_B ||
			modelname == Heterogeneous_p_Ej_Lk_Bk ||
			modelname == Heterogeneous_p_Ek_L_B ||
			modelname == Heterogeneous_p_Ek_L_Bk ||
			modelname == Heterogeneous_p_Ek_Lk_B ||
			modelname == Heterogeneous_p_Ek_Lk_Bk ||
			modelname == Heterogeneous_p_Ekj_L_B ||
			modelname == Heterogeneous_p_Ekj_L_Bk ||
			modelname == Heterogeneous_p_Ekj_Lk_B ||
			modelname == Heterogeneous_p_Ekj_Lk_Bk ||
			modelname == Heterogeneous_p_Ekjh_L_B ||
			modelname == Heterogeneous_p_Ekjh_L_Bk ||
			modelname == Heterogeneous_p_Ekjh_Lk_B ||
			modelname == Heterogeneous_p_Ekjh_Lk_Bk ||
			modelname == Heterogeneous_pk_E_L_B ||
			modelname == Heterogeneous_pk_E_L_Bk ||
			modelname == Heterogeneous_pk_E_Lk_B ||
			modelname == Heterogeneous_pk_E_Lk_Bk ||
			modelname == Heterogeneous_pk_Ej_L_B ||
			modelname == Heterogeneous_pk_Ej_L_Bk ||
			modelname == Heterogeneous_pk_Ej_Lk_B ||
			modelname == Heterogeneous_pk_Ej_Lk_Bk ||
			modelname == Heterogeneous_pk_Ek_L_B ||
			modelname == Heterogeneous_pk_Ek_L_Bk ||
			modelname == Heterogeneous_pk_Ek_Lk_B ||
			modelname == Heterogeneous_pk_Ek_Lk_Bk ||
			modelname == Heterogeneous_pk_Ekj_L_B ||
			modelname == Heterogeneous_pk_Ekj_L_Bk ||
			modelname == Heterogeneous_pk_Ekj_Lk_B ||
			modelname == Heterogeneous_pk_Ekj_Lk_Bk ||
			modelname == Heterogeneous_pk_Ekjh_L_B ||
			modelname == Heterogeneous_pk_Ekjh_L_Bk ||
			modelname == Heterogeneous_pk_Ekjh_Lk_B ||
			modelname == Heterogeneous_pk_Ekjh_Lk_Bk
			)
		return true;
	else
		return false;
}

ModelGenre getModelGenre(ModelName name) {
	if (isBinary(name)) return QualitativeModel;

	if (isHeterogeneous(name)) return HeterogeneousModel;

	return QuantitativeModel;
}

//---------------------------------
// is modelType has free proportion
//---------------------------------

bool hasFreeProportion(ModelName modelName) {
	bool res = true;
	if (		 (modelName == Gaussian_p_L_I)
			|| (modelName == Gaussian_p_Lk_I)
			|| (modelName == Gaussian_p_L_B)
			|| (modelName == Gaussian_p_Lk_B)
			|| (modelName == Gaussian_p_L_Bk)
			|| (modelName == Gaussian_p_Lk_Bk)
			|| (modelName == Gaussian_p_L_C)
			|| (modelName == Gaussian_p_Lk_C)
			|| (modelName == Gaussian_p_L_D_Ak_D)
			|| (modelName == Gaussian_p_Lk_D_Ak_D)
			|| (modelName == Gaussian_p_L_Dk_A_Dk)
			|| (modelName == Gaussian_p_Lk_Dk_A_Dk)
			|| (modelName == Gaussian_p_L_Ck)
			|| (modelName == Gaussian_p_Lk_Ck)
			|| (modelName == Binary_p_E)
			|| (modelName == Binary_p_Ej)
			|| (modelName == Binary_p_Ek)
			|| (modelName == Binary_p_Ekj)
			|| (modelName == Binary_p_Ekjh)
			|| (modelName == Gaussian_HD_p_AkjBkQkDk)
			|| (modelName == Gaussian_HD_p_AkBkQkDk)
			|| (modelName == Gaussian_HD_p_AkjBkQkD)
			|| (modelName == Gaussian_HD_p_AjBkQkD)
			|| (modelName == Gaussian_HD_p_AkjBQkD)
			|| (modelName == Gaussian_HD_p_AjBQkD)
			|| (modelName == Gaussian_HD_p_AkBkQkD)
			|| (modelName == Gaussian_HD_p_AkBQkD)
			|| modelName == Heterogeneous_p_E_L_B ||
			modelName == Heterogeneous_p_E_L_Bk ||
			modelName == Heterogeneous_p_E_Lk_B ||
			modelName == Heterogeneous_p_E_Lk_Bk ||
			modelName == Heterogeneous_p_Ej_L_B ||
			modelName == Heterogeneous_p_Ej_L_Bk ||
			modelName == Heterogeneous_p_Ej_Lk_B ||
			modelName == Heterogeneous_p_Ej_Lk_Bk ||
			modelName == Heterogeneous_p_Ek_L_B ||
			modelName == Heterogeneous_p_Ek_L_Bk ||
			modelName == Heterogeneous_p_Ek_Lk_B ||
			modelName == Heterogeneous_p_Ek_Lk_Bk ||
			modelName == Heterogeneous_p_Ekj_L_B ||
			modelName == Heterogeneous_p_Ekj_L_Bk ||
			modelName == Heterogeneous_p_Ekj_Lk_B ||
			modelName == Heterogeneous_p_Ekj_Lk_Bk ||
			modelName == Heterogeneous_p_Ekjh_L_B ||
			modelName == Heterogeneous_p_Ekjh_L_Bk ||
			modelName == Heterogeneous_p_Ekjh_Lk_B ||
			modelName == Heterogeneous_p_Ekjh_Lk_Bk) {
		res = false;
	}
	return res;
}

void editSimpleTab(double * tab, int64_t n, std::string sep, std::string before, std::ostream & flux) {
	int64_t i;
	flux << before;
	for (i = 0; i < n; i++) {
		flux << tab[i] << sep;
	}
	flux << endl ;
}

void editSimpleTab(int64_t    * tab, int64_t n, std::ostream & flux ) {
	int64_t i;
	for (i = 0; i < n; i++)
		flux << tab[i] << endl ;
}

//---------------
// Move on file fi until what is reached
// after using that function fi is just after the first time what appears
//---------------
void moveUntilReach(std::ifstream & fi, std::string  what) {
	std::string keyWord = "";
	ConvertBigtoLowString(what);
	// init reading at the beginning of file //
	fi.clear();
	fi.seekg(0, ios::beg);
	// read until finding *what* we are looking for
	do {
		fi >> keyWord ;
		ConvertBigtoLowString(keyWord);
	}

		// while( !fi.eof() && strcmp(keyWord,what)!=0 ) ;
	while ( !fi.eof() && (keyWord.compare(what) != 0) ) ;
	//  delete[] keyWord;
}

//-------------------
// read nbNbCluster file names (ex : titi ; toto;tutu)
void readTabFileName(std::ifstream & fi, int64_t nbNbCluster, std::string* tabFileName, std::string& keyWord) {
	int64_t k = 0;

	std::string c = "";
	std::string c1 = "";
	std::string tmp = "";
	std::string strBeforePv = "";
	std::string strAfterPv = "";

	fi >> c;
	// on ne convertit pas tout le nom en minuscules sinon il y a des erreurs dans le nom des fichiers

	while (!isKeyword(c) && (!fi.eof())) {

		if (c.compare(";") == 0) {
			k++;
			fi >> c;
		}
		else {
			if (c.find_first_of(';') == 0) { // si c commence par ;
				k++;
				strAfterPv = c.substr(1, c.length());
			}
			else {
				strAfterPv = c;
			}
			while ((strAfterPv.find_first_of(';') != std::string::npos)) { // ; est dans la chaine de caracteres
				tmp = strAfterPv;
				strBeforePv = tmp.substr(0, tmp.find_first_of(';'));
				strAfterPv = tmp.substr(tmp.find_first_of(';') + 1, tmp.length());


				if (tabFileName[k].length() == 0) {
					tabFileName[k] = strBeforePv;
				}
				else {
					tabFileName[k].append(" ");
					tabFileName[k].append(strBeforePv);

				}
				k++;
			}

			if (tabFileName[k].length() == 0) {
				tabFileName[k] = strAfterPv;
			}
			else {
				tabFileName[k].append(" ");
				tabFileName[k].append(strBeforePv);
			}
			fi >> c;
		}
	}

	if ( (k != nbNbCluster - 1) || (tabFileName[nbNbCluster - 1].compare("") == 0) || (tabFileName[nbNbCluster - 1].compare(" ") == 0))
		THROW(InputException, wrongPartitionFileName);

	keyWord = c;
}

void initToZero(double* tab, int64_t n) {
	double * p_tab = tab;
	int64_t i;
	for (i = 0 ; i < n ; i++, p_tab++) {
		*p_tab = 0.0;
	}
}

inline void echange(double * tab, int64_t i1, int64_t i2) {
	double tmp = tab[i1];
	tab[i1]           = tab[i2];
	tab[i2]           = tmp ;
}

inline void echange(int64_t * tab, int64_t i1, int64_t i2) {
	int64_t tmp   = tab[i1];
	tab[i1]         = tab[i2];
	tab[i2]         = tmp ;
}

void selectionSortWithOrder(double * tabRandom, int64_t * tabOrder, int64_t left, int64_t right) {
	int64_t i, j;
	int64_t min;

	for (i = left; i < right; i++) {
		min = i;
		for (j = i + 1; j <= right; j++)
			if (tabRandom[j] < tabRandom[min])
				min = j;
		echange(tabRandom, min, i);
		echange(tabOrder, min, i);
	}
}

int64_t partition(double * tabRandom, int64_t * tabOrder, int64_t left, int64_t right) {

	double val = tabRandom[left];
	int64_t lm    = left - 1;
	int64_t rm    = right + 1;
	for (; ; ) {
		do
			rm--;
		while (tabRandom[rm] > val);

		do
			lm++;
		while ( tabRandom[lm] < val);

		if (lm < rm) {
			echange(tabRandom, rm, lm);
			echange(tabOrder , rm, lm);
		}
		else
			return rm;
	}
}

void quickSortWithOrder(double * tabRandom, int64_t * tabOrder, int64_t left, int64_t right) {
	if (left < (right - SMALL_ENOUGH_TO_USE_SELECTION_SORT)) {
		int64_t split_pt = partition(tabRandom, tabOrder, left, right);
		quickSortWithOrder(tabRandom, tabOrder, left      , split_pt);
		quickSortWithOrder(tabRandom, tabOrder, split_pt + 1, right);
	}
	else selectionSortWithOrder(tabRandom, tabOrder, left, right);
}

//-------------------
//generateRandomIndex [TODO: suboptimal method near when array contain a lot of "false"]
//-------------------
int64_t generateRandomIndex(bool * tabIndividualCanBeUsedForInitRandom, double * weight, double totalWeight) {
	double rndWeight, sumWeight;
	int64_t idxSample;

	/* Generate a random integer between 0 and _nbSample-1 */
	bool IdxSampleCanBeUsed = false;  // idxSample can be used
	while (!IdxSampleCanBeUsed) {
		// get index of sample with weight //
		rndWeight = (int64_t) (totalWeight * rnd() + 1);
		sumWeight = 0.0;
		idxSample = -1;
		while (sumWeight < rndWeight) {
			idxSample++;
			sumWeight += weight[idxSample];
		}
		//cout<<"index tire au hasard :"<<idxSample<<endl;
		IdxSampleCanBeUsed = tabIndividualCanBeUsedForInitRandom[idxSample];
	}
	// on indique que cet individu ne pourra pas �tre tir� au  hasard pour une autre classe
	tabIndividualCanBeUsedForInitRandom[idxSample] = false;
	//cout<<"choisi"<<endl;
	return idxSample;
}

void inputCriterion(std::ifstream & fi, CriterionName & criterionName) {
	std::string a = "";
	fi >> a;
	if (a.compare("BIC") == 0) {
		criterionName = BIC;
	}
	else if (a.compare("CV") == 0) {
		criterionName = CV;
	}
	else if (a.compare("ICL") == 0) {
		criterionName = ICL;
	}
	else if (a.compare("NEC") == 0) {
		criterionName = NEC;
	}
	else if (a.compare("DCV") == 0) {
		criterionName = DCV;
	}
	else {
		THROW(InputException, wrongCriterionName);
	}
}

void inputCVinitBlocks(std::ifstream & fi, CVinitBlocks cVinitBlocks) {
	std::string a = "";
	fi >> a;
	if (a.compare("CV_RANDOM") == 0) {
		cVinitBlocks = CV_RANDOM;
	}
	else if (a.compare("DIAG") == 0) {
		cVinitBlocks = CV_DIAG;
	}
	else {
		THROW(InputException, wrongCVinitType);
	}
}

void inputDCVinitBlocks(std::ifstream & fi, DCVinitBlocks dCVinitBlocks) {
	std::string a = "";
	fi >> a;
	if (a.compare("DCV_RANDOM") == 0) {
		dCVinitBlocks = DCV_RANDOM;
	}
	else if (a.compare("DIAG") == 0) {
		dCVinitBlocks = DCV_DIAG;
	}
	else {
		THROW(InputException, wrongDCVinitType);
	}
}

}
