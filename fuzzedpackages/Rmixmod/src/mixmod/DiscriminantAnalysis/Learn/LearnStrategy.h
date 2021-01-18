/***************************************************************************
                             SRC/mixmod/DiscriminantAnalysis/Learn/LearnStrategy.h  description
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
#ifndef XEMLearnStrategy_H
#define XEMLearnStrategy_H

namespace XEM {

// pre-declaration
class Model;

/** 
 \class XEMLearnStrategy
 Main class for Learn Strategy (1rst step of discriminant analysis)
 @author F. Langrognet - R Lebret
		@date 2012
		@brief XEMLearnStrategy class
 */
class  LearnStrategy {

public:

	/// Default constructor
	LearnStrategy();

	/// Constructor
	LearnStrategy(const LearnStrategy & strategy);

	/// Destructor
	~LearnStrategy();

	/// Run method
	void run(Model * model);

	// verify method
	bool verify();

};

}

#endif
