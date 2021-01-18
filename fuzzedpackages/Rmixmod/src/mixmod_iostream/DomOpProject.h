/***************************************************************************
                             SRC/MIXMOD_IOSTREAM/XEMDomClusteringProject.h  description
    copyright            : (C) MIXMOD Team - 2001-2011
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

#ifndef XEM_DOMOPPROJECT_H
#define XEM_DOMOPPROJECT_H

#include "mixmod_iostream/DomProject.h"

namespace XEM {
  class ClusteringMain;
///use to create .mixmod file in Clustering case
  class DomOpProject : public DomProject {

public : 
    
	///constructor by default
	DomOpProject();    

	///destructor
	virtual ~DomOpProject();

	///constructor by initialization
	DomOpProject(xmlpp::Element *root);

	///fill the xmlpp::Document to create the .mixmod file from a ClusteringInput and ClusteringOutput
	//void writeClustering(string & s, ClusteringMain * cMain);    
	void writeMixmodXml(string & s, ClusteringMain * cMain);
	void writeMixmodXml(string & s, LearnMain * lMain);    
    void writeMixmodXml(string& s, PredictMain * pMain);
	///read a XML file and fill ClusteringInput
	//void readClustering(ClusteringInput * cInput);
	void readXmlFillIn(ClusteringInput * cInput);
	//void readXmlFillIn(LearnInput * cInput);        
    template<class T>
    void readXmlFillIn(T  *cInput);
    PredictInput * readXmlPredictInput();
    //template<class ClusteringInput*>
    //void readXmlFillIn(ClusteringInput *cInput);          
	///read a XML file and fill ClusteringOutput
	//void readClustering(ClusteringOutput * cOutput);
    //void readXmlFillOut(ClusteringOutput * cOutput);
    //void readXmlFillOut(LearnOutput * cOutput);
    template<typename T, typename U>
      void readXmlFillOut(T  * cOutput, Input *inp);
    //template<typename T>    
    //void readXmlFillOut(T * cOutput);    
};

} //end namespace

#endif // XEM_DOMOPPROJECT_H
