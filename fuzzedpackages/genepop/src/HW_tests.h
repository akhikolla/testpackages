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
#ifndef HW_TESTS_H
#define HW_TESTS_H

// comme j'ai defini hardy1 pour passer les valeurs de variables, j'ai
// tr�s peu besoin d'extern

const int lowSwitchNbr=1000;

int hardymin();
int hardy1(bool defbool,bool prbool, bool globbool,bool hwbool,std::string& hw_outfile);
void ecriture_sample_HW(std::string& hw_outfile);
int lecture_fich_PL(bool testBool,std::vector<int>& effallnbr);
int enumeration_test(int nn, double uobs);
int calcul_proba(size_t nn);
void traitement_result_fichiers(std::vector<std::string>& markName,
                                std::vector<std::string>& exacName,
                                std::vector<int>& effallnbr,
                                std::string& hw_outfile);
int enum_test_et_affich(std::vector<std::string>& exacName);
int HW_Pvalues_estimate();
void alonzy();
int HWfile_info();
int global_U_initialize(std::vector<std::vector<bool> >& indic,std::vector<double>& Uloc,std::vector<double>& Upop);
int HW_Pvalues_compile(std::vector<std::vector<bool> >& indic,std::vector<double>& Uloc,std::vector<double>& Upop,std::string& hw_outfile);
int HW_Pvalues_chains(std::vector<std::string>& markName);
void fic_lect();
void ecriture_result(std::string& hw_outfile);
void delete_proba();
int Genclean_HW();

void initializeHWtests();
void cleanHWtests();

#endif
