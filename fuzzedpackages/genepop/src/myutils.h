/***************************************************************************
@ F. Rousset 2005-2006

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
#ifndef MYUTILS_H
#define MYUTILS_H

// proc�dures d'affichages de vecteurs et de temps �coul�

#include <vector>
#include <iterator>



template <typename Type>
int ostream_vector(const std::vector<Type> &v, std::ostream& stream){
std::ostream_iterator<Type> locout(stream," ");
copy(v.begin(),v.end(),locout);
//stream << endl;
return 0;
}

template <typename Type>
int ofstream_vector(std::vector<Type> &v, std::ofstream& stream){
std::ostream_iterator<Type> locout(stream," ");
copy(v.begin(),v.end(),locout);
//stream << endl;
return 0;
}

template <typename Type>
int ofstream_vector_mll(std::vector<Type> &v,Type mll, std::ofstream& stream){
std::ostream_iterator<Type> locout(stream," ");
copy(v.begin(),v.end(),locout);
stream<< mll << std::endl;
return 0;
}

template <typename Type1,typename Type2>
int ofstream_vector_mll(std::vector<Type1> &v,std::vector<Type2> &mll, std::ofstream& stream){
std::ostream_iterator<Type1> locout1(stream," ");
std::ostream_iterator<Type2> locout2(stream," ");
copy(v.begin(),v.end(),locout1);
copy(mll.begin(),mll.end(),locout2);
stream<<std::endl;
return 0;
}


template <typename Type>
int ostream_vec_vector(std::vector<std::vector<Type> > &v, std::ostream& stream){
std::ostream_iterator<Type> locout(stream," ");
    // Without the keyword typename, the ANSI rules require that
    // the qualified identifer T::A is not the name of a type,
    // because it is dependent upon a template paramter.
    //
    // The comment shows where the ANSI rules require the keyword typename.
///////    typedef /*typename*/ T::A type;
for(typename std::vector<std::vector<Type> >::iterator it=v.begin();it<v.end();it++) {
	copy((*it).begin(),(*it).end(),locout);
	stream << std::endl;
}
return 0;
}

/*
Pas simplement la class declaration mais aussi definition de functions
=> pas sa place dans header file
#include <ctime>


class ShowTime{
protected:
string description;
clock_t avant,mtnt;

public:
ShowTime(string description) {
	this->description=description;
	avant=clock();
}
~ShowTime() {
	mtnt=clock();
	cout<<float(mtnt- avant)/CLOCKS_PER_SEC<<"s "<<description;cout.flush();
}
};
*/

#endif //MYUTILS_H
