/***************************************************************************
                             SRC/mixmod/Kernel/Parameter/CompositeParameter.h  description
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
#ifndef XEMCOMPOSITEPARAMETER_H_
#define XEMCOMPOSITEPARAMETER_H_
/**@file XEMCompositeParameter.h
 * @brief Composite parameter class for heterogeneous clustering.
 * @author Parmeet Bhatia
 */

#include "mixmod/Kernel/Parameter/Parameter.h"
#include "mixmod/Kernel/Parameter/GaussianEDDAParameter.h"

namespace XEM {

class GaussianParameter;
class BinaryParameter;

class CompositeParameter : public Parameter {

public:
	
	//default constructor
	CompositeParameter();
	//copy constructor
	CompositeParameter(const CompositeParameter *);
	//Actual constructor
	CompositeParameter(Model*, ModelType*, int64_t*);
	CompositeParameter(const Parameter * igaussian, const Parameter * ibinary,
			ModelType * imodelType);
	void InstantiateBinaryandGaussianParamter(ModelType*, int64_t*);
	virtual void setModel(Model * iModel);
	virtual void reset();
	virtual CompositeParameter * clone() const;
	/** Return Product of Pdf's of Gaussian and Binary parameters.*/
	virtual double getPdf(int64_t iSample, int64_t kCluster) const;
	virtual void getAllPdf(double ** tabFik, double * tabProportion) const;
	/** Return Product of Pdf's of Gaussian and Binary parameters.*/
	virtual double getPdf(Sample * x, int64_t kCluster) const;
	virtual int64_t getFreeParameter() const;

	virtual double getLogLikelihoodOne() const;
	virtual void initUSER(Parameter * iParam);
	virtual void initForInitRANDOM();
	virtual void updateForInitRANDOMorUSER_PARTITION(Sample ** tabSampleForInit, bool * tabClusterToInitialze);

	virtual void initForInitUSER_PARTITION(int64_t & nbInitializedCluster, bool * tabNotInitializedCluster, Partition * initPartition) {
		THROW(OtherException, FunctionNotYetImplemented);
	}
	virtual void MStep();
	virtual void edit();
	virtual void edit(std::ofstream & oFile, bool text = false);

	virtual void input(std::ifstream & fi);

	virtual void recopy(Parameter * otherParameter);
	virtual void updateForCV(Model * originalModel, CVBlock & CVBlock);
	/**typecast overloading for BinaryData*/
	operator BinaryParameter*();
	/**typecast overloading for XEGaussianData*/
	operator GaussianParameter*();

	virtual BinaryParameter* getBinaryParameter() {
		return (BinaryParameter*) _parameterComponent[0];
	}

	virtual GaussianParameter* getGaussianParameter() {
		return (GaussianParameter*) _parameterComponent[1];
	}
	virtual ~CompositeParameter();
	
protected:
	
	vector<Parameter*> _parameterComponent;
	vector<ModelType*> _parameterModelType;
};

inline CompositeParameter::operator BinaryParameter *() {
	return (BinaryParameter*) _parameterComponent[0];
}

inline CompositeParameter::operator GaussianParameter *() {
	return (GaussianParameter*) _parameterComponent[1];
}

}

#endif /* XEMCOMPOSITEPARAMETER_H_ */
