/***************************************************************************
@ F. Rousset 2005-2007

francois.rousset@umontpellier.fr

This file is part of Genepop'007
This software is a computer program whose purpose is to perform statistical analyses.

This software is governed by the CeCILL license under French law and
abiding by the rules of distribution of free software.  You can  use,
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info".

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability.

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or
data to be ensured and,  more generally, to use and operate it in the
same conditions as regards security.

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.

 ***************************************************************************/
#ifndef BOOTSTRAP_H
#define BOOTSTRAP_H

#include <fstream>
#include <string>
#include <vector>
#include <limits>



class CABCbootstrap {
        size_t ABCloc;
  size_t nb_units;
        std::string testLegend;
	    std::vector<double> delta;
	    double bidullevel,z,ahat;
        double cancelland(double Pvalue);
        double (*estimFnPtr)(std::vector<double> d);
        double testPoint;
        friend double cancellandWrapper(double unidirPvalue);
        double seuil_inf,seuil_sup;
public:
        int bootstrapOverLoci(double (*estimatingFnPtr)(std::vector<double> d),std::string legend,std::string boot0utfile,bool clear_screen=true);
        double Pvalue(double testPt,bool unidir,bool verbose);
        double t0,tinf,tsup,testPointPvalue;
        CABCbootstrap(size_t units);
        ~CABCbootstrap() {}
};

std::vector<double> bisection_search(double (*func)(double d),double x1,double x2,bool verbose=true);

#endif
