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
#ifndef F_EST_H
#define F_EST_H
extern std::string statname,_logdist;
extern double mindist,maxdist;
extern bool _a_stat,_e_stat;
extern bool geoDistFromGeoFile;
extern std::string geoDistFile;
//extern std::vector<double>diversities;

template <typename Type>
struct    MStype{size_t bbordel;
    				 Type mmsp;
    				 Type mmsi;
    				 Type mmsg;
    				 Type nnc;
};

#ifdef FLOAT_CALC
typedef float MSreal;
#else
typedef double MSreal;
#endif //float_calc

namespace datamatrix {
    extern size_t	nb_sam_migf;
    extern std::vector<std::vector<long double> >data;
}  //extern in namespace


namespace NS_F_est {
  extern std::size_t global_pop_it, global_pop2_it;
    extern size_t nb_sam;
    extern size_t nb_locus;
    extern size_t nb_pair_sam;
    extern bool _first_of_repl;
    extern MStype<MSreal> *houlaMS,*houlaMSit;
    extern double MSp2P,MSg2Pw,MSi2P,MSi2Pw;
}  //extern in namespace


namespace NS_FFF_slmt {
    extern double SSiTot,SSgTot,MSg2P;
}

//prototype fn aillaurs
int controle_choix();
//prototypes fn locales
int ident_size();
int genotip2();
void lecture_floc();
void calc_sfreqs_Nc(int indic, std::vector<double>& tailleMoy,
                    std::vector<unsigned long int>& scgsp,
                    std::vector<std::vector<double> >& sfreqs,
                    size_t& Ntot);
void calculSSetMS(std::vector<double>& tailleMoy,
                  std::vector<unsigned long int>& scgsp,
                  std::vector<std::vector<double> >& sfreqs,
                  size_t& Ntot);
void lecturePaires(void);
int main2x2(std::vector<bool>& ploidBool);
int pairwMS(std::vector<bool>& ploidBool);
void idxinf(std::vector<std::vector<long double> >& locdata,
            std::vector<std::vector<long double> >& indx);
void idxsup(std::vector<std::vector<long double> >& locdata,
            std::vector<std::vector<long double> >& indx);
void mantelTest(bool clearscreen, bool rankBool);
std::vector<double> calcwritecorw();
std::vector<double> isoldeproc(const char nom_fich_mig[]);
void readGGFile(const char nom_fich_mig[]);
void writedat(std::vector<std::vector<long double> >m,char nom_fich_mig[]);
void writepma();
void conversion(void);
void MS_for_isolde();
void lecture_popi_popj();
int ecriture_pop_tot(std::ofstream& f_mig);
int set_options(bool perf,bool indiv,bool identbool);
int set_phylipmatrix(bool pmat);
int set_ptrs();
int delete_ptrs();
//int choix_stat();
int create_matrices(const char nom_fich_mig[]);
void lecture_Xi_Xj_pmoy();
//int set_astat_dim(bool astat,bool logdist);
void tabFtotabM(std::vector<std::vector<unsigned long int> > * tabF);
void FisParPop(bool identitybool,size_t iLoc,std::ofstream &fic_out,
               std::vector<double>& MSgTotHez,
               std::vector<int>& NLocHez,
               std::vector<int>& NLoc,
               std::vector<int>& NtotLoc,
               std::vector<double>& SSiTotLoc,
               std::vector<double>& SSgTotLoc);
void delete_tabM_tabCode();
void hierFstat(bool identitybool,int Indic,size_t iLoc,std::ofstream &fic_out, std::vector<std::vector<double> > * FFF);
void conversionFst();
void conversionGeo();
void writeGraOnly(const char nom_fich_mig[]);

void initializeFest();
void cleanFest();


#endif
