/***************************************************************************
                             SRC/mixmod/Kernel/IO/BinaryData.h  description
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
#ifndef XEMBINARYDATA_H
#define XEMBINARYDATA_H

#include "mixmod/Kernel/IO/Data.h"

namespace XEM {

// pre-declaration
class Partition;
class BinarySample;

/**
  @brief Base class for Binary Data
  @author F Langrognet 
 */

class BinaryData : public Data {

public:

	/// Default constructor
	BinaryData();

	/// Constructor
	BinaryData(const BinaryData & iData);

	/// Constructor
	BinaryData(int64_t nbSample, int64_t pbDimension, 
			const std::string & dataFileName, int64_t * tabNbModality);

	/// Constructor
	BinaryData(int64_t nbSample, int64_t pbDimension, std::vector<int64_t> nbModality);

	/// Constructor
	BinaryData(int64_t nbSample, int64_t pbDimension, 
			std::vector<int64_t> nbModality, int64_t ** matrix);

	/// Constructor for dataReduce
	BinaryData(int64_t nbSample, int64_t pbDimension, int64_t * tabNbModality, 
			double weightTotal, Sample **& matrix, double * weight);


	/// Constructor (used in DCV context)
	BinaryData(int64_t nbSample, int64_t pbDimension, Data * originalData, CVBlock & block);

	/// Desctructor
	virtual ~BinaryData();



	virtual void input(const DataDescription & dataDescription);

	/** @brief copy
		@return A copy of data
	 */
	virtual Data * clone() const;

	/**  @brief Copy
		 @return A copy data matrix
	 */
	virtual Sample ** cloneMatrix();

	/** @brief Read data from binary data file
		@param fi Binary Data file to read
	 */
	virtual void input(std::ifstream & fi);

	/** @brief Write binary data in output file
		@param fo Output file to write into
	 */
	virtual void output(std::ostream & fo);


	virtual bool verify()const;

	/** @brief Get matrix of data Sample
		@return A vector of XEMSample
	 */
	Sample ** getDataMatrix() const;

	/** @brief Get matrix of data Sample
		@param idxSample Index of sample to get values
		@return A vector of XEMSample
	 */
	int64_t * getDataTabValue(int64_t idxSample) const;

	/** @brief Get tab modality
		@return A vector of number of modality
	 */
	int64_t * getTabNbModality() const;

	Data * reduceData(std::vector<int64_t> & correspondcenceOriginDataToReduceData, 
			Partition * knownPartition, Partition * initPartition, 
			Partition *& oKnownPartition, Partition *& oInitPartition);

protected:

	/// Array of modality for each dimension
	int64_t * _tabNbModality;
	Data * _reducedData;
};

//---------------
// inline methods
//---------------

inline Sample ** BinaryData::getDataMatrix() const {
	return _matrix;
}

inline int64_t * BinaryData::getTabNbModality() const {
	return _tabNbModality;
}

}

#endif
