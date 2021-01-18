/***************************************************************************
                             SRC/mixmod/Kernel/Algo/Algo.cpp  description
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
#include "mixmod/Kernel/Algo/Algo.h"
#include "mixmod/Kernel/Algo/EMAlgo.h"

namespace XEM {

//-----------
//Constructor
//-----------
Algo::Algo() {
	_indexIteration = 0;
	_algoStopName = defaultAlgoStopName;
	_epsilon = defaultEpsilon;
	_nbIteration = defaultNbIteration;
	_xml_old = 0;
	_xml = 0;
}

// copy constructor
Algo::Algo(const Algo & algo) {
	_indexIteration = algo._indexIteration;
	_algoStopName = algo._algoStopName;
	_epsilon = algo._epsilon;
	_nbIteration = algo._nbIteration;
	_xml_old = algo._xml_old;
	_xml = algo._xml;
}

Algo::Algo(AlgoStopName algoStopName, double espsilon, int64_t nbIteration) {
	_indexIteration = 1;
	_algoStopName = algoStopName;
	setEpsilon(espsilon);
	setNbIteration(nbIteration);
	_xml_old = 0;
	_xml = 0;
}

//----------
//Destructor
//----------
Algo::~Algo() {
}

void Algo::setEpsilon(double epsilon) {
	if (epsilon < minEpsilon) {
		THROW(InputException, epsilonTooSmall);
	}
	else if (epsilon > maxEpsilon) {
		THROW(InputException, epsilonTooLarge);
	}
	else {
		_epsilon = epsilon;
	}
}

//--------------
// continueAgain
//--------------
/* Stopping rule for algorithm : continueAgain */
bool Algo::continueAgain() {
	//cout<<"Algo::continueAgain"<<endl;
	bool result, res1, res2;
	double diff;
	ofstream progressFile;

	if (_indexIteration == 1) {
		result = true;
		return result;
	}
	else {

		if (_indexIteration > maxNbIteration) {
			result = false;
			return result;
		}
		else {
			switch (_algoStopName) {
				case NBITERATION:
					if (MASSICCC == 10) {
						progressFile.open ("progress.json");
						progressFile << "{ \"Progress\" : " << ((double)_indexIteration - 1.0)/(double)_nbIteration * 100.0 << "}";
						progressFile.close();
					}
					result = (_indexIteration <= _nbIteration);
					break;

				case EPSILON:
					if (MASSICCC == 10) {
						progressFile.open ("progress.json");
						progressFile << "{ \"Progress\" : " << ((double)_indexIteration - 1.0)/1000.0 * 100.0 << "}";
						progressFile.close();
					}
					if (_indexIteration <= 3) {
						result = true;
					}
					else {
						diff = fabs(_xml - _xml_old);
						result = (diff >= _epsilon);
					}
					if (!result) {
						if (MASSICCC == 10) {
							progressFile.open ("progress.json");
							progressFile << "{ \"Progress\" : 100 }";
							progressFile.close();
						}
						break;
					}

				case NBITERATION_EPSILON:
					if (MASSICCC == 10) {
						progressFile.open ("progress.json");
						progressFile << "{ \"Progress\" : " << ((double)_indexIteration - 1.0)/(double)_nbIteration * 100.0 << "}";
						progressFile.close();
					}
					res1 = (_indexIteration <= _nbIteration);
					if (_indexIteration <= 3) {
						res2 = true;
					}
					else {
						diff = fabs(_xml - _xml_old);
						res2 = (diff >= _epsilon);
					};
					result = (res1 && res2);
					if (!result) {
						if (MASSICCC == 10) {
							progressFile.open ("progress.json");
							progressFile << "{ \"Progress\" : 100 }";
							progressFile.close();
						}
						break;
					}

			default: result = (_indexIteration <= _nbIteration);
			}
			return result;
		}
	}
}

//-----------
// ostream <<
//-----------
std::ostream & operator <<(std::ostream & fo, Algo & algo) {
	AlgoName algoName = algo.getAlgoName();
	fo << "\t  Type : " << AlgoNameToString(algoName);

	fo << "\t  Stopping rule : ";

	switch (algo._algoStopName) {

	case NBITERATION:
		fo << "NBITERATION" << endl;
		fo << "\t  Number of iterations : " << algo._nbIteration << endl;
		break;
	case NBITERATION_EPSILON:
		fo << "NBITERATION_EPSILON" << endl;
		fo << "\t  Number of iterations : " << algo._nbIteration << endl;
		fo << "\t  Set tolerance (xml criterion) : " << algo._epsilon << endl;
		break;
	case EPSILON:
		fo << "EPSILON" << endl;
		fo << "\t  Set tolerance (xml criterion) : " << algo._epsilon << endl;
		break;
	case NO_STOP_NAME:
		break;
	default:
		break;
	}
	return fo;
}

void Algo::edit(std::ostream & out) {
	AlgoName algoName = getAlgoName();
	out << "\t  Type : " << AlgoNameToString(algoName) << endl;

	out << "\t  Stopping rule : ";

	switch (_algoStopName) {

	case NBITERATION:
		out << "NBITERATION" << endl;
		out << "\t  Number of iterations : " << _nbIteration << endl;
		break;
	case NBITERATION_EPSILON:
		out << "NBITERATION_EPSILON" << endl;
		out << "\t  Number of iterations : " << _nbIteration << endl;
		out << "\t  Set tolerance (xml criterion) : " << _epsilon << endl;
		break;
	case EPSILON:
		out << "EPSILON" << endl;
		out << "\t  Set tolerance (xml criterion) : " << _epsilon << endl;
		break;
	case NO_STOP_NAME:
		break;
	default:
		break;
	}
}

//----------------
// others functions
//
Algo * createDefaultClusteringAlgo() {
	if (defaultClusteringAlgoName != EM) {
		THROW(OtherException, internalMixmodError);
	}
	Algo * algo = new EMAlgo();
	return algo;
}

}
