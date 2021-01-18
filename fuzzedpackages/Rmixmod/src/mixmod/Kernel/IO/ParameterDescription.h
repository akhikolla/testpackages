/***************************************************************************
                             SRC/mixmod/Kernel/IO/ParameterDescription.h  description
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

#ifndef XEMPARAMETERDESCRIPTION_H
#define XEMPARAMETERDESCRIPTION_H

#include "mixmod/Utilities/Util.h"

namespace XEM {

// pre-declaration
class Model;
class ModelOutput;
class Parameter;

/** 
 \class XEMParameterDescription
 @author F. Langrognet
		@date 2011
		@brief XEMParameterDescription class
 */
class ParameterDescription {

public:
	/// Default constructor
	ParameterDescription();

	/// Constructor
	ParameterDescription(Model* iEstimation);

	/// Constructor
	ParameterDescription(ModelOutput* iEstimation);

	// constructor for Binary
	/// Constructor
	ParameterDescription(
			int64_t nbCluster, 
			int64_t nbVariable, 
			std::vector< int64_t > nbFactor, 
			FormatNumeric::FormatNumericFile format, 
			std::string filename, 
			std::string infoName, 
			ModelName & modelName);

	// constructor for Gaussian
	/// Constructor
	ParameterDescription(
			int64_t nbCluster, 
			int64_t nbVariable, 
			FormatNumeric::FormatNumericFile format, 
			std::string filename, 
			std::string infoName, 
			ModelName & modelName);

	// constructor for Composite
	/// Constructor
	ParameterDescription(
			int64_t nbCluster,
			int64_t nbVariable_binary,
			int64_t nbVariable_gaussian,
			std::vector< int64_t > nbFactor,
			FormatNumeric::FormatNumericFile format,
			std::string filename,
			std::string infoName,
			ModelName & modelName);

	ParameterDescription(Parameter * iparam);

	//constructor for binary data
	ParameterDescription(
			int64_t nbCluster,
			int64_t nbVariable,
			ModelName& modelName,
			double * proportions,
			double ** centers,
			double *** scatters,
			std::vector< int64_t > nbFactor);

	//constructor for Gaussian data
	ParameterDescription(
			int64_t nbCluster,
			int64_t nbVariable,
			ModelName& modelName,
			double * proportions,
			double ** means,
			double *** variances);

	//constructor for Heterogeneous
	ParameterDescription(
			int64_t nbCluster, 
			int64_t nbBinaryVariable, 
			int64_t nbGaussianVariable, 
			ModelName& modelName, 
			double * proportions, 
			double ** centers, 
			double *** scatters, 
			double ** means, 
			double *** variances, 
			std::vector< int64_t > nbFactor);

	/// Destructor
	~ParameterDescription();

	/// Comparison operator
	bool operator==(ParameterDescription & paramDescription) const;

	/// getParameter
	Parameter * getParameter();

	///getInfoName
	std::string getInfoName();

	///getPbDimension
	int64_t getNbVariable();

	///getNbCluster
	int64_t getNbCluster();

	///getFormat
	FormatNumeric::FormatNumericFile getFormat();

	///getFilename
	std::string getFilename();

	///getModelType
	ModelType * getModelType();

	///getTabNbModality
	std::vector<int64_t> & getTabNbFactor();

	void saveNumericValues(std::string fileName = "");

private:

	std::string _infoName;

	int64_t _nbVariable;

	int64_t _nbCluster;

	FormatNumeric::FormatNumericFile _format; //format of  numeric file

	std::string _filename;

	std::vector<int64_t> _nbFactor;

	ModelType * _modelType;

	Parameter * _parameter;
};

inline Parameter * ParameterDescription::getParameter() {
	if (_parameter) {
		return _parameter;
	}
	else {
		THROW(OtherException, nullPointerError);
	}
}

inline int64_t ParameterDescription::getNbCluster() {
	return _nbCluster;
}

inline std::string ParameterDescription::getInfoName() {
	return _infoName;
}

inline int64_t ParameterDescription::getNbVariable() {
	return _nbVariable;
}

inline FormatNumeric::FormatNumericFile ParameterDescription::getFormat() {
	return _format;
}

inline std::string ParameterDescription::getFilename() {
	return _filename;
}

inline ModelType * ParameterDescription::getModelType() {
	return _modelType;
}

inline std::vector<int64_t> & ParameterDescription::getTabNbFactor() {
	return _nbFactor;
}

}

#endif // XEMDATADESCRIPTION_H
