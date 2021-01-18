/***************************************************************************
                             SRC/mixmod/Utilities/exceptions/NumericException.h  description
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
/** @file XEMNumericException.h
 *  @brief Exception class for handling Numeric errors.
 *  @author Parmeet Bhatia
 **/

#ifndef XEM_NUMERICEXCEPTION_H
#define XEM_NUMERICEXCEPTION_H

#include "mixmod/Utilities/exceptions/Exception.h"

namespace XEM {

class NumericException : public Exception {

public:
	
	NumericException(std::string, int, NumericError) throw ();
	NumericException(NumericError) throw ();
	Exception * clone() throw ();
	virtual const char* what() const throw ();
	virtual bool operator==(const Exception&) const throw ();
	virtual void run(std::ostream & flux = std::cout) const throw ();

	virtual ~NumericException() throw () {
	}

	static std::map<NumericError, const char*> create_map()
	{
		std::map<NumericError, const char*> m;
		
		m.insert(std::make_pair(int64_t_max_error, "int64_t limits reached"));
		m.insert(std::make_pair(CEM_INIT_error, "Error : error in CEM_INIT initialization : in all tries occurs an error"));
		m.insert(std::make_pair(SEM_MAX_error, "Error : error in SEM_MAX initialization : in all tries occurs an error"));
		m.insert(std::make_pair(SMALL_EM_error, "Error : error in SMALL_EM initialization : in all tries occurs an error"));
		m.insert(std::make_pair(tabNkNotInteger, "tabNk is not an integer"));
		m.insert(std::make_pair(sumFiNullAndfkTPrimNull, "SumF[i]=0 in multinomial case"));
		m.insert(std::make_pair(sumFiNullInMultinomialCase, "SumF[i]=0 in multinomial case. Please report this bug to Mixmod team !"));
		m.insert(std::make_pair(nonPositiveDefiniteMatrix,"Variance matrix is non positive definite"));
		m.insert(std::make_pair(nullDeterminant, "Matrix determinant null"));
		m.insert(std::make_pair(randomProblem, "Problem in Random function call"));
		m.insert(std::make_pair(nullLikelihood, "Null likelihood"));
		m.insert(std::make_pair(noProbability, "Probability (cf. NormPdf function) not between 0 and 1"));
		m.insert(std::make_pair(pbNEC, "NEC criterion, the difference between log-likelihood for k and 1 cluster is null"));
		m.insert(std::make_pair(nullNk, "Empty mixture component"));
		m.insert(std::make_pair(numericError, "A numeric error occured"));
		m.insert(std::make_pair(errorSigmaConditionNumber, "error in variance matrix condition number"));
		m.insert(std::make_pair(minDeterminantSigmaValueError, "determinant of matrix too small"));
		m.insert(std::make_pair(minDeterminantWValueError, "determinant of matrix too small"));
		m.insert(std::make_pair(minDeterminantDiagWkValueError, "determinant of matrix too small"));
		m.insert(std::make_pair(minDeterminantDiagWValueError, "determinant of matrix too small"));
		m.insert(std::make_pair(minDeterminantBValueError, "determinant of matrix too small"));
		m.insert(std::make_pair(minDeterminantRValueError, "determinant of matrix too small"));
		m.insert(std::make_pair(minDeterminantWkValueError, "determinant of matrix too small"));
		m.insert(std::make_pair(minDeterminantShapeValueError, "determinant of matrix too small"));
		m.insert(std::make_pair(minDeterminantDiagQtmpValueError, "determinant of matrix too small"));

		return m;
	}

	static std::map<NumericError, const char*> mapErrorMsg;

protected:
	
	NumericError _errorType;
};

}

#endif /* XEMNUMERICEXCEPTION_H_ */
