/***************************************************************************
                             SRC/mixmod/Kernel/Parameter/CompositeParameter.cpp  description
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
#include "mixmod/Kernel/Model/Model.h"
#include "mixmod/Kernel/Parameter/CompositeParameter.h"
#include "mixmod/Kernel/IO/GaussianSample.h"
#include "mixmod/Kernel/IO/BinarySample.h"
#include "mixmod/Kernel/Parameter/BinaryParameter.h"
#include "mixmod/Kernel/Parameter/GaussianParameter.h"
#include "mixmod/Kernel/Parameter/BinaryEParameter.h"
#include "mixmod/Kernel/Parameter/BinaryEkParameter.h"
#include "mixmod/Kernel/Parameter/BinaryEjParameter.h"
#include "mixmod/Kernel/Parameter/BinaryEkjParameter.h"
#include "mixmod/Kernel/Parameter/BinaryEkjhParameter.h"
#include "mixmod/Kernel/Parameter/GaussianDiagParameter.h"
#include "mixmod/Kernel/Model/ModelType.h"

namespace XEM {

CompositeParameter::CompositeParameter() {
	// TODO Auto-generated constructor stub
}

CompositeParameter::CompositeParameter(const CompositeParameter * param)
: Parameter(param->getModel(), param->getModelType())
{
	_parameterComponent.resize(2);
	_parameterModelType.resize(2);
	_parameterComponent[0] =
			(const_cast<CompositeParameter*> (param)->getBinaryParameter())->clone();
	_parameterComponent[1] =
			((const_cast<CompositeParameter*> (param)->getGaussianParameter()))->clone();
	_parameterModelType[0] = new ModelType(_parameterComponent[0]->getModelType()->_nameModel);
	_parameterComponent[0]->setModelType(_parameterModelType[0]);
	_parameterModelType[1] = new ModelType(_parameterComponent[1]->getModelType()->_nameModel);
	_parameterComponent[1]->setModelType(_parameterModelType[1]);
}

CompositeParameter::CompositeParameter(const Parameter * igaussian, const Parameter * ibinary,
		ModelType * imodelType) : Parameter(ibinary->getNbCluster(),
		ibinary->getPbDimension() + igaussian->getPbDimension(), imodelType)
{
	_parameterComponent.resize(2);
	_parameterModelType.resize(2);
	_parameterComponent[0] =
			(const_cast<Parameter*> (ibinary)->getBinaryParameter())->clone();
	_parameterComponent[1] =
			((const_cast<Parameter*> (igaussian)->getGaussianParameter()))->clone();
	_parameterModelType[0] = new ModelType(_parameterComponent[0]->getModelType()->_nameModel);
	_parameterComponent[0]->setModelType(_parameterModelType[0]);
	_parameterModelType[1] = new ModelType(_parameterComponent[1]->getModelType()->_nameModel);
	_parameterComponent[1]->setModelType(_parameterModelType[1]);
}

CompositeParameter::CompositeParameter(Model* iModel, ModelType* iModelType,
		int64_t * tabNbModality) : Parameter(iModel, iModelType)
{
	_parameterComponent.resize(2);
	_parameterModelType.resize(2);
	//Instantiate Binary and Gaussian models
	InstantiateBinaryandGaussianParamter(iModelType, tabNbModality);
}

CompositeParameter::~CompositeParameter() {
	for (unsigned int i = 0; i < _parameterComponent.size(); ++i) {
		if (_parameterComponent[i]) delete _parameterComponent[i];
		if (_parameterModelType[i]) delete _parameterModelType[i];
	}
}

double CompositeParameter::getLogLikelihoodOne() const
{
	// TODO: adjust likelihoods combination
	return _parameterComponent[0]->getLogLikelihoodOne()
		+ _parameterComponent[1]->getLogLikelihoodOne();
}

void CompositeParameter::InstantiateBinaryandGaussianParamter(
		ModelType* modeltype, int64_t * tabNbModality)
{
	ModelName modelname = modeltype->getModelName();
	switch (modelname) {
	case Heterogeneous_p_E_L_B:
	{
		_parameterModelType[0] = new ModelType(Binary_p_E);
		_parameterModelType[1] = new ModelType(Gaussian_p_L_B);
		_parameterComponent[0] =
				new BinaryEParameter(_model, _parameterModelType[0], tabNbModality);
		_parameterComponent[1] = new GaussianDiagParameter(_model, _parameterModelType[1]);
		break;
	}
	case Heterogeneous_p_E_Lk_B:
	{
		_parameterModelType[0] = new ModelType(Binary_p_E);
		_parameterModelType[1] = new ModelType(Gaussian_p_Lk_B);
		_parameterComponent[0] =
				new BinaryEParameter(_model, _parameterModelType[0], tabNbModality);
		_parameterComponent[1] = new GaussianDiagParameter(_model, _parameterModelType[1]);
		break;
	}
	case Heterogeneous_p_E_L_Bk:
	{
		_parameterModelType[0] = new ModelType(Binary_p_E);
		_parameterModelType[1] = new ModelType(Gaussian_p_L_Bk);
		_parameterComponent[0] =
				new BinaryEParameter(_model, _parameterModelType[0], tabNbModality);
		_parameterComponent[1] = new GaussianDiagParameter(_model, _parameterModelType[1]);
		break;
	}
	case Heterogeneous_p_E_Lk_Bk:
	{
		_parameterModelType[0] = new ModelType(Binary_p_E);
		_parameterModelType[1] = new ModelType(Gaussian_p_Lk_Bk);
		_parameterComponent[0] =
				new BinaryEParameter(_model, _parameterModelType[0], tabNbModality);
		_parameterComponent[1] = new GaussianDiagParameter(_model, _parameterModelType[1]);
		break;
	}
	case Heterogeneous_p_Ek_L_B:
	{
		_parameterModelType[0] = new ModelType(Binary_p_E);
		_parameterModelType[1] = new ModelType(Gaussian_p_L_B);
		_parameterComponent[0] =
				new BinaryEkParameter(_model, _parameterModelType[0], tabNbModality);
		_parameterComponent[1] = new GaussianDiagParameter(_model, _parameterModelType[1]);
		break;
	}
	case Heterogeneous_p_Ek_Lk_B:
	{
		_parameterModelType[0] = new ModelType(Binary_p_E);
		_parameterModelType[1] = new ModelType(Gaussian_p_Lk_B);
		_parameterComponent[0] =
				new BinaryEkParameter(_model, _parameterModelType[0], tabNbModality);
		_parameterComponent[1] = new GaussianDiagParameter(_model, _parameterModelType[1]);
		break;
	}
	case Heterogeneous_p_Ek_L_Bk:
	{
		_parameterModelType[0] = new ModelType(Binary_p_E);
		_parameterModelType[1] = new ModelType(Gaussian_p_L_Bk);
		_parameterComponent[0] =
				new BinaryEkParameter(_model, _parameterModelType[0], tabNbModality);
		_parameterComponent[1] = new GaussianDiagParameter(_model, _parameterModelType[1]);
		break;
	}
	case Heterogeneous_p_Ek_Lk_Bk:
	{
		_parameterModelType[0] = new ModelType(Binary_p_E);
		_parameterModelType[1] = new ModelType(Gaussian_p_Lk_Bk);
		_parameterComponent[0] =
				new BinaryEkParameter(_model, _parameterModelType[0], tabNbModality);
		_parameterComponent[1] = new GaussianDiagParameter(_model, _parameterModelType[1]);
		break;
	}
	case Heterogeneous_p_Ej_L_B:
	{
		_parameterModelType[0] = new ModelType(Binary_p_E);
		_parameterModelType[1] = new ModelType(Gaussian_p_L_B);
		_parameterComponent[0] =
				new BinaryEjParameter(_model, _parameterModelType[0], tabNbModality);
		_parameterComponent[1] = new GaussianDiagParameter(_model, _parameterModelType[1]);
		break;
	}
	case Heterogeneous_p_Ej_Lk_B:
	{
		_parameterModelType[0] = new ModelType(Binary_p_E);
		_parameterModelType[1] = new ModelType(Gaussian_p_Lk_B);
		_parameterComponent[0] =
				new BinaryEjParameter(_model, _parameterModelType[0], tabNbModality);
		_parameterComponent[1] = new GaussianDiagParameter(_model, _parameterModelType[1]);
		break;
	}
	case Heterogeneous_p_Ej_L_Bk:
	{
		_parameterModelType[0] = new ModelType(Binary_p_E);
		_parameterModelType[1] = new ModelType(Gaussian_p_L_Bk);
		_parameterComponent[0] =
				new BinaryEjParameter(_model, _parameterModelType[0], tabNbModality);
		_parameterComponent[1] = new GaussianDiagParameter(_model, _parameterModelType[1]);
		break;
	}
	case Heterogeneous_p_Ej_Lk_Bk:
	{
		_parameterModelType[0] = new ModelType(Binary_p_E);
		_parameterModelType[1] = new ModelType(Gaussian_p_Lk_Bk);
		_parameterComponent[0] =
				new BinaryEjParameter(_model, _parameterModelType[0], tabNbModality);
		_parameterComponent[1] = new GaussianDiagParameter(_model, _parameterModelType[1]);
		break;
	}
	case Heterogeneous_p_Ekj_L_B:
	{
		_parameterModelType[0] = new ModelType(Binary_p_E);
		_parameterModelType[1] = new ModelType(Gaussian_p_L_B);
		_parameterComponent[0] =
				new BinaryEkjParameter(_model, _parameterModelType[0], tabNbModality);
		_parameterComponent[1] = new GaussianDiagParameter(_model, _parameterModelType[1]);
		break;
	}
	case Heterogeneous_p_Ekj_Lk_B:
	{
		_parameterModelType[0] = new ModelType(Binary_p_E);
		_parameterModelType[1] = new ModelType(Gaussian_p_Lk_B);
		_parameterComponent[0] =
				new BinaryEkjParameter(_model, _parameterModelType[0], tabNbModality);
		_parameterComponent[1] = new GaussianDiagParameter(_model, _parameterModelType[1]);
		break;
	}
	case Heterogeneous_p_Ekj_L_Bk:
	{
		_parameterModelType[0] = new ModelType(Binary_p_E);
		_parameterModelType[1] = new ModelType(Gaussian_p_L_Bk);
		_parameterComponent[0] =
				new BinaryEkjParameter(_model, _parameterModelType[0], tabNbModality);
		_parameterComponent[1] = new GaussianDiagParameter(_model, _parameterModelType[1]);
		break;
	}
	case Heterogeneous_p_Ekj_Lk_Bk:
	{
		_parameterModelType[0] = new ModelType(Binary_p_E);
		_parameterModelType[1] = new ModelType(Gaussian_p_Lk_Bk);
		_parameterComponent[0] =
				new BinaryEkjParameter(_model, _parameterModelType[0], tabNbModality);
		_parameterComponent[1] = new GaussianDiagParameter(_model, _parameterModelType[1]);
		break;
	}
	case Heterogeneous_p_Ekjh_L_B:
	{
		_parameterModelType[0] = new ModelType(Binary_p_E);
		_parameterModelType[1] = new ModelType(Gaussian_p_L_B);
		_parameterComponent[0] =
				new BinaryEkjhParameter(_model, _parameterModelType[0], tabNbModality);
		_parameterComponent[1] = new GaussianDiagParameter(_model, _parameterModelType[1]);
		break;
	}
	case Heterogeneous_p_Ekjh_Lk_B:
	{
		_parameterModelType[0] = new ModelType(Binary_p_E);
		_parameterModelType[1] = new ModelType(Gaussian_p_Lk_B);
		_parameterComponent[0] =
				new BinaryEkjhParameter(_model, _parameterModelType[0], tabNbModality);
		_parameterComponent[1] = new GaussianDiagParameter(_model, _parameterModelType[1]);
		break;
	}
	case Heterogeneous_p_Ekjh_L_Bk:
	{
		_parameterModelType[0] = new ModelType(Binary_p_E);
		_parameterModelType[1] = new ModelType(Gaussian_p_L_Bk);
		_parameterComponent[0] =
				new BinaryEkjhParameter(_model, _parameterModelType[0], tabNbModality);
		_parameterComponent[1] = new GaussianDiagParameter(_model, _parameterModelType[1]);
		break;
	}
	case Heterogeneous_p_Ekjh_Lk_Bk:
	{
		_parameterModelType[0] = new ModelType(Binary_p_E);
		_parameterModelType[1] = new ModelType(Gaussian_p_Lk_Bk);
		_parameterComponent[0] =
				new BinaryEkjhParameter(_model, _parameterModelType[0], tabNbModality);
		_parameterComponent[1] = new GaussianDiagParameter(_model, _parameterModelType[1]);
		break;
	}
	case Heterogeneous_pk_E_L_B:
	{
		_parameterModelType[0] = new ModelType(Binary_pk_E);
		_parameterModelType[1] = new ModelType(Gaussian_pk_L_B);
		_parameterComponent[0] =
				new BinaryEParameter(_model, _parameterModelType[0], tabNbModality);
		_parameterComponent[1] = new GaussianDiagParameter(_model, _parameterModelType[1]);
		break;
	}
	case Heterogeneous_pk_E_Lk_B:
	{
		_parameterModelType[0] = new ModelType(Binary_pk_E);
		_parameterModelType[1] = new ModelType(Gaussian_pk_Lk_B);
		_parameterComponent[0] =
				new BinaryEParameter(_model, _parameterModelType[0], tabNbModality);
		_parameterComponent[1] = new GaussianDiagParameter(_model, _parameterModelType[1]);
		break;
	}
	case Heterogeneous_pk_E_L_Bk:
	{
		_parameterModelType[0] = new ModelType(Binary_pk_E);
		_parameterModelType[1] = new ModelType(Gaussian_pk_L_Bk);
		_parameterComponent[0] =
				new BinaryEParameter(_model, _parameterModelType[0], tabNbModality);
		_parameterComponent[1] = new GaussianDiagParameter(_model, _parameterModelType[1]);
		break;
	}
	case Heterogeneous_pk_E_Lk_Bk:
	{
		_parameterModelType[0] = new ModelType(Binary_pk_E);
		_parameterModelType[1] = new ModelType(Gaussian_pk_Lk_Bk);
		_parameterComponent[0] =
				new BinaryEParameter(_model, _parameterModelType[0], tabNbModality);
		_parameterComponent[1] = new GaussianDiagParameter(_model, _parameterModelType[1]);
		break;
	}
	case Heterogeneous_pk_Ek_L_B:
	{
		_parameterModelType[0] = new ModelType(Binary_pk_E);
		_parameterModelType[1] = new ModelType(Gaussian_pk_L_B);
		_parameterComponent[0] =
				new BinaryEkParameter(_model, _parameterModelType[0], tabNbModality);
		_parameterComponent[1] = new GaussianDiagParameter(_model, _parameterModelType[1]);
		break;
	}
	case Heterogeneous_pk_Ek_Lk_B:
	{
		_parameterModelType[0] = new ModelType(Binary_pk_E);
		_parameterModelType[1] = new ModelType(Gaussian_pk_Lk_B);
		_parameterComponent[0] =
				new BinaryEkParameter(_model, _parameterModelType[0], tabNbModality);
		_parameterComponent[1] = new GaussianDiagParameter(_model, _parameterModelType[1]);
		break;
	}
	case Heterogeneous_pk_Ek_L_Bk:
	{
		_parameterModelType[0] = new ModelType(Binary_pk_E);
		_parameterModelType[1] = new ModelType(Gaussian_pk_L_Bk);
		_parameterComponent[0] =
				new BinaryEkParameter(_model, _parameterModelType[0], tabNbModality);
		_parameterComponent[1] = new GaussianDiagParameter(_model, _parameterModelType[1]);
		break;
	}
	case Heterogeneous_pk_Ek_Lk_Bk:
	{
		_parameterModelType[0] = new ModelType(Binary_pk_E);
		_parameterModelType[1] = new ModelType(Gaussian_pk_Lk_Bk);
		_parameterComponent[0] = new
				BinaryEkParameter(_model, _parameterModelType[0], tabNbModality);
		_parameterComponent[1] = new GaussianDiagParameter(_model, _parameterModelType[1]);
		break;
	}
	case Heterogeneous_pk_Ej_L_B:
	{
		_parameterModelType[0] = new ModelType(Binary_pk_E);
		_parameterModelType[1] = new ModelType(Gaussian_pk_L_B);
		_parameterComponent[0] =
				new BinaryEjParameter(_model, _parameterModelType[0], tabNbModality);
		_parameterComponent[1] = new GaussianDiagParameter(_model, _parameterModelType[1]);
		break;
	}
	case Heterogeneous_pk_Ej_Lk_B:
	{
		_parameterModelType[0] = new ModelType(Binary_pk_E);
		_parameterModelType[1] = new ModelType(Gaussian_pk_Lk_B);
		_parameterComponent[0] =
				new BinaryEjParameter(_model, _parameterModelType[0], tabNbModality);
		_parameterComponent[1] = new GaussianDiagParameter(_model, _parameterModelType[1]);
		break;
	}
	case Heterogeneous_pk_Ej_L_Bk:
	{
		_parameterModelType[0] = new ModelType(Binary_pk_E);
		_parameterModelType[1] = new ModelType(Gaussian_pk_L_Bk);
		_parameterComponent[0] =
				new BinaryEjParameter(_model, _parameterModelType[0], tabNbModality);
		_parameterComponent[1] = new GaussianDiagParameter(_model, _parameterModelType[1]);
		break;
	}
	case Heterogeneous_pk_Ej_Lk_Bk:
	{
		_parameterModelType[0] = new ModelType(Binary_pk_E);
		_parameterModelType[1] = new ModelType(Gaussian_pk_Lk_Bk);
		_parameterComponent[0] =
				new BinaryEjParameter(_model, _parameterModelType[0], tabNbModality);
		_parameterComponent[1] = new GaussianDiagParameter(_model, _parameterModelType[1]);
		break;
	}
	case Heterogeneous_pk_Ekj_L_B:
	{
		_parameterModelType[0] = new ModelType(Binary_pk_E);
		_parameterModelType[1] = new ModelType(Gaussian_pk_L_B);
		_parameterComponent[0] =
				new BinaryEkjParameter(_model, _parameterModelType[0], tabNbModality);
		_parameterComponent[1] = new GaussianDiagParameter(_model, _parameterModelType[1]);
		break;
	}
	case Heterogeneous_pk_Ekj_Lk_B:
	{
		_parameterModelType[0] = new ModelType(Binary_pk_E);
		_parameterModelType[1] = new ModelType(Gaussian_pk_Lk_B);
		_parameterComponent[0] =
				new BinaryEkjParameter(_model, _parameterModelType[0], tabNbModality);
		_parameterComponent[1] = new GaussianDiagParameter(_model, _parameterModelType[1]);
		break;
	}
	case Heterogeneous_pk_Ekj_L_Bk:
	{
		_parameterModelType[0] = new ModelType(Binary_pk_E);
		_parameterModelType[1] = new ModelType(Gaussian_pk_L_Bk);
		_parameterComponent[0] =
				new BinaryEkjParameter(_model, _parameterModelType[0], tabNbModality);
		_parameterComponent[1] = new GaussianDiagParameter(_model, _parameterModelType[1]);
		break;
	}
	case Heterogeneous_pk_Ekj_Lk_Bk:
	{
		_parameterModelType[0] = new ModelType(Binary_pk_E);
		_parameterModelType[1] = new ModelType(Gaussian_pk_Lk_Bk);
		_parameterComponent[0] =
				new BinaryEkjParameter(_model, _parameterModelType[0], tabNbModality);
		_parameterComponent[1] = new GaussianDiagParameter(_model, _parameterModelType[1]);
		break;
	}
	case Heterogeneous_pk_Ekjh_L_B:
	{
		_parameterModelType[0] = new ModelType(Binary_pk_E);
		_parameterModelType[1] = new ModelType(Gaussian_pk_L_B);
		_parameterComponent[0] =
				new BinaryEkjhParameter(_model, _parameterModelType[0], tabNbModality);
		_parameterComponent[1] = new GaussianDiagParameter(_model, _parameterModelType[1]);
		break;
	}
	case Heterogeneous_pk_Ekjh_Lk_B:
	{
		_parameterModelType[0] = new ModelType(Binary_pk_E);
		_parameterModelType[1] = new ModelType(Gaussian_pk_Lk_B);
		_parameterComponent[0] =
				new BinaryEkjhParameter(_model, _parameterModelType[0], tabNbModality);
		_parameterComponent[1] = new GaussianDiagParameter(_model, _parameterModelType[1]);
		break;
	}
	case Heterogeneous_pk_Ekjh_L_Bk:
	{
		_parameterModelType[0] = new ModelType(Binary_pk_E);
		_parameterModelType[1] = new ModelType(Gaussian_pk_L_Bk);
		_parameterComponent[0] =
				new BinaryEkjhParameter(_model, _parameterModelType[0], tabNbModality);
		_parameterComponent[1] = new GaussianDiagParameter(_model, _parameterModelType[1]);
		break;
	}
	case Heterogeneous_pk_Ekjh_Lk_Bk:
	{
		_parameterModelType[0] = new ModelType(Binary_pk_E);
		_parameterModelType[1] = new ModelType(Gaussian_pk_Lk_Bk);
		_parameterComponent[0] =
				new BinaryEkjhParameter(_model, _parameterModelType[0], tabNbModality);
		_parameterComponent[1] = new GaussianDiagParameter(_model, _parameterModelType[1]);
		break;
	}
	default:
		THROW(InputException, wrongModelName);
	}
}

void CompositeParameter::getAllPdf(double** tabFik, double* tabProportion)const {
	int64_t nbSample = _model->getNbSample();
	int64_t i;
	int64_t k;

	for (i = 0; i < nbSample; i++) {
		for (k = 0; k < _nbCluster; k++) {
			tabFik[i][k] = getPdf(i, k) * tabProportion[k];
		}
	}
}

double CompositeParameter::getPdf(int64_t iSample, int64_t Kcluster) const {
	return _parameterComponent[1]->getPdf(iSample, Kcluster)
			* _parameterComponent[0]->getPdf(iSample, Kcluster);
}

double CompositeParameter::getPdf(Sample* x, int64_t kcluster) const {
	Sample * GSample_ = x->getGaussianSample();
	assert(GSample_ != NULL);
	Sample * BSample_ = x->getBinarySample();
	assert(BSample_ != NULL);
	return _parameterComponent[1]->getPdf(GSample_, kcluster)
			* (_parameterComponent[0]->getPdf(BSample_, kcluster));
}

void CompositeParameter::MStep() {
	//Tab proportion Update
	computeTabProportion();
	//Binary M-step
	_parameterComponent[0]->MStep();
	//Gaussian M-step
	_parameterComponent[1]->MStep();

#ifndef NDEBUG
	double * binarytabprop_ = _parameterComponent[0]->getTabProportion();
	double * gaussiantabprop_ = _parameterComponent[1]->getTabProportion();
	for (int i = 0; i < _nbCluster; ++i) {
		assert(binarytabprop_[i] == gaussiantabprop_[i]);
	}
#endif
}

CompositeParameter * CompositeParameter::clone() const {
	CompositeParameter * newparam = new CompositeParameter(this);
	return newparam;
}

int64_t CompositeParameter::getFreeParameter() const {
	int64_t freeparam = _parameterComponent[1]->getFreeParameter()
			+ _parameterComponent[0]->getFreeParameter()-(_nbCluster - 1);
	return freeparam;
}

void CompositeParameter::initForInitRANDOM() {
	_parameterComponent[0]->initForInitRANDOM();
	_parameterComponent[1]->initForInitRANDOM();
}

void CompositeParameter::updateForInitRANDOMorUSER_PARTITION(
		Sample**tabSampleForInit, bool*tabClusterToInitialize)
{
	_parameterComponent[0]->updateForInitRANDOMorUSER_PARTITION(
			tabSampleForInit, tabClusterToInitialize);
	_parameterComponent[1]->updateForInitRANDOMorUSER_PARTITION(
			tabSampleForInit, tabClusterToInitialize);
}

void CompositeParameter::recopy(Parameter* other) {
	_parameterComponent[0]->recopy(other);
	_parameterComponent[1]->recopy(other);
}

void CompositeParameter::reset() {
	_parameterComponent[0]->reset();
	_parameterComponent[1]->reset();
}

void CompositeParameter::updateForCV(Model* originalModel, CVBlock&CVBlock) {
	_parameterComponent[0]->updateForCV(originalModel, CVBlock);
	_parameterComponent[1]->updateForCV(originalModel, CVBlock);
}

void CompositeParameter::initUSER(Parameter * iParam) {
	double * iTabProportion = iParam->getTabProportion();
	for (int k = 0; k < _nbCluster; k++) {
		// proportion (no respecting model type)
		if (!hasFreeProportion(_modelType->_nameModel)) {
			_tabProportion[k] = 1.0 / _nbCluster;
		}
		else {
			_tabProportion[k] = iTabProportion[k];
		}
	}
	_parameterComponent[0]->initUSER(iParam);
	_parameterComponent[1]->initUSER(iParam);
}

void CompositeParameter::edit() {
	cout << "Binary Parameters\n";
	cout << "**********************************************************\n";
	_parameterComponent[0]->edit();
	cout << "\nGaussian Parameters\n";
	cout << "**********************************************************\n";
	_parameterComponent[1]->edit();
}

void CompositeParameter::input(std::ifstream & fi)
{
	_parameterComponent[0]->input(fi);
	_parameterComponent[1]->input(fi);
}

void CompositeParameter::edit(std::ofstream & oFile, bool text) {
	_parameterComponent[0]->edit(oFile, text);
	_parameterComponent[1]->edit(oFile, text);
}

void CompositeParameter::setModel(Model * iModel) {
	_model = iModel;
	_parameterComponent[0]->setModel(iModel);
	_parameterComponent[1]->setModel(iModel);
}

}
