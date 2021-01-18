/***************************************************************************
                             SRC/mixmod/DiscriminantAnalysis/Predict/PredictStrategy.h  description
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
#ifndef XEMPredictStrategy_H
#define XEMPredictStrategy_H

namespace XEM {

// class pre-declaration
class Model;
class Parameter;

/** 
 \class XEMPredictStrategy
 @author F. Langrognet - R Lebret
		@date 2012
		@brief XEMPredictStrategy class
 */

class  PredictStrategy {

public:

	/// Default constructor
	PredictStrategy(Parameter * classificationRule);

	/// Constructor
	PredictStrategy(const PredictStrategy & strategy);

	/// Destructor
	~PredictStrategy();

	/// Run method
	void run(Model * model);

	// verify method
	bool verify();

	// get pointer to the algorithm
	Parameter * getClassificationRule() const;

private:
	
	// classification rule
	Parameter * _classificationRule;
};

// get pointer to the algorithm
inline Parameter * PredictStrategy::getClassificationRule() const {
	return _classificationRule;
}

}

#endif
