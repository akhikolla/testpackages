/***************************************************************************
                             SRC/mixmod/Clustering/ClusteringStrategyInit.h  description
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
#ifndef XEMClusteringStrategyInit_H
#define XEMClusteringStrategyInit_H

#include "mixmod/Utilities/Util.h"

namespace XEM {

// pre-declaration
class Data;
class Parameter;
class Partition;
class Model;

/**
  @brief Base class for StrategyInitType(s)
  @author F Langrognet 
 */

class ClusteringStrategyInit {

public:

	/// Default constructor
	ClusteringStrategyInit();

	/// Copy constructor
	ClusteringStrategyInit(const ClusteringStrategyInit & strategyInit);

	/// Destructor
	virtual ~ClusteringStrategyInit();

	/// getStrategyInitName
	const StrategyInitName & getStrategyInitName() const;

	/// setStrategyInitName
	void setStrategyInitName(StrategyInitName initName);

	/// getNbTry
	const int64_t getNbTry() const;

	/// setNbTry
	void setNbTry(int64_t nbTry);

	/// getNbIteration
	const int64_t getNbIteration() const;

	/// set NbIteration
	void setNbIteration(int64_t nbIteration);

	/// getEpsilon
	const double getEpsilon() const;

	/// setEpsilon
	void setEpsilon(double epsilon);

	/// getStopName
	const AlgoStopName getStopName() const;

	/// setStopName
	void setStopName(AlgoStopName stopName);

	/* Parameter */
	//-----------//
	/// getNbInitParameter
	const int64_t & getNbInitParameter() const;

	/// getTabInitParameter
	const Parameter ** getTabInitParameter() const;

	/// getTabInitParameter
	Parameter * getInitParameter(int64_t index) const;

	/// setInitParam
	void setInitParam(std::string & paramFileName, int64_t position);

	/// setTabInitParam
	void setTabInitParameter(Parameter ** tabInitParameter, int64_t nbInitParameter);

	/* Partition */
	//-----------//
	///  getNbPartition
	const int64_t & getNbPartition() const;

	/// getTabPartition
	const Partition ** getTabPartition() const;

	/// getTabPartition
	Partition * getPartition(int64_t index) const;

	///set Init Partition
	void setPartition(Partition * part, int64_t position);

	///set Init Partition
	void setPartition(std::string & paramFileName, int64_t position);

	/// setTabPartition
	void setTabPartition(Partition ** tabPartition, int64_t nbPartition);

	/* Initialization */
	//----------------//

	/// Initialization by EM of the parameters of the model
	void initSMALL_EM(Model*& model);

	/// Initialization by CEM of the parameters of the model
	void initCEM_INIT(Model*& model);

	/// Initialization by SEM of the parameters of the model
	void initSEM_MAX(Model*& model);

	/* Input / Output */
	//----------------//
	// input
	void input(std::ifstream & fi, Data *& data, int64_t nbNbCluster,
			int64_t * tabNbCluster, ModelType * modelType, bool & alreadyRead);

	// verification
	bool  verify() const;

	// print out strategy initialization
	friend std::ostream & operator << (std::ostream & fo, ClusteringStrategyInit & strategyInit);

private:

	/// Initialization strategy type
	StrategyInitName _strategyInitName;

	/// nbTry
	int64_t _nbTry;

	/// stopName (for smallEm & sem_max)
	AlgoStopName _stopName;

	/// nbIteration
	int64_t _nbIteration;

	/// epsilon
	double _epsilon;

	/* USER */
	/// number of InitParameter
	int64_t _nbInitParameter;

	/// Init Parameters to initialize strategy
	Parameter ** _tabInitParameter;

	/* USER_PARTITION */
	/// number of Partition
	int64_t _nbPartition;

	/// Labels to initialize strategy
	Partition ** _tabPartition;

	bool _deleteTabParameter ;

	void oneRunOfSmallEM(Model*& model, double & logLikelihood);
};

inline  const StrategyInitName & ClusteringStrategyInit::getStrategyInitName() const {
	return _strategyInitName;
}

inline const AlgoStopName ClusteringStrategyInit::getStopName() const {
	return _stopName;
}

inline const int64_t & ClusteringStrategyInit::getNbInitParameter() const {
	return _nbInitParameter;
}

inline const int64_t & ClusteringStrategyInit::getNbPartition() const {
	return _nbPartition;
}

inline const Parameter ** ClusteringStrategyInit::getTabInitParameter() const {
	return const_cast<const Parameter**> (_tabInitParameter);
}

inline Parameter * ClusteringStrategyInit::getInitParameter(int64_t index) const {
	return _tabInitParameter[index];
}

inline const Partition** ClusteringStrategyInit::getTabPartition() const {
	return  const_cast<const Partition**> (_tabPartition);
}

inline Partition* ClusteringStrategyInit::getPartition(int64_t index) const {
	return _tabPartition[index];
}

inline const int64_t ClusteringStrategyInit::getNbTry() const {
	return _nbTry;
}

inline const int64_t ClusteringStrategyInit::getNbIteration() const {
	return _nbIteration;
}

inline const double ClusteringStrategyInit::getEpsilon() const {
	return _epsilon;
}

}

#endif
