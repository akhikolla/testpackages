/***************************************************************************
© R. Leblois 2001-2004
© F. Rousset 2005--

francois.rousset@umontpellier.fr

This file is part of Genepop'007
This software is a computer program whose purpose is to perform statistical analyses.

This software is governed by the CeCILL license under French law and
abiding by the rules of distribution of free software.  You can  use,
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt".

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
#ifndef H_GENEPOP
#define H_GENEPOP
#include <iostream>
#include <iomanip>
#include <cstring>
#include <cstdlib>
#include <list>
#include <string>
#include <map>
#include <vector>
#include <fstream>
#include <cctype>
#include <cstdio>

extern std::ostream cout_abyss;
#ifdef COMPATIBILITYRCPP
#define noR_cout cout_abyss // cout_abyss to be initialized so that the stream contents are ignored. 
#define RnoR_cerr Rcpp::Rcerr 
#define Ronly_cout Rcpp::Rcout
#else
#define noR_cout std::cout
#define RnoR_cerr cerr 
#define Ronly_cout cout_abyss
#endif

extern bool inputCheckBool,cinGetOnError;
extern const std::string fichierIn;

void _gotoxy(int x,int y);
int wherey();
int wherex();
void effacer_ecran();
int cmp_nocase(const std::string& s, const std::string& s2);
void rtrim(std::string *s);
size_t minAllele(ssize_t genotype,char coding);
size_t maxAllele(ssize_t genotype,char coding);

// CBR: previously known as CFichier_genepop::set_eof_check_EOLtype(bool set_eof /*=true*/)
int set_eof_check_EOLtype(std::string EOLFileName, bool set_eof);
// CBR: previously part of CFichier_genepop::parseFile()
void set_UNIX_EOLtype(std::string EOLFileName);

void initializegenepop();
void cleangenepop();

// ------------------------------------------------------------------------------------------
// ------------ d['e]clarations ----------------------------------------------------------------
// ------------------------------------------------------------------------------------------
// CLocus stocke les effectifs alleliques !
// CLocusGP enfait une map !
// CFichier_genepop a vector<CLocusGP *>loci !
//    et vector<CPopulation *>pops qui a aussi vector<CLocusGP *>loci !


// classe de base d'un locus
class CLocus {
	public:
		std::string locusName; // nom du locus
		size_t identif; // n° d'ordre du locus
		int alleleMax; // taille de l'all[`e]le le plus grand (NOM d'allele max)
		int galleleMax; // taille de l'all[`e]le le plus grand (NOM d'allele max)
		CLocus(size_t identifiant, std::string locName); // constructeur
		virtual ~CLocus();
};


// classe de stockage d'un all[`e]le
class CAllele {
	public:
		int _identif; // identifiant de l'all[`e]le // used as index of vector
		unsigned long int _effective; // effectif de cet all[`e]le
		CAllele(int identif, unsigned long int effective); // constructeur
		void incEffective(); // incr['e]mente l'effectif
		unsigned long int getEffective(); // renvoie l'effectif
};

// classe d['e]riv['e] d'un locus pour le stockage des all[`e]les
class CLocusGP : public CLocus {
	private:
		int galleleExists(int num); // renvoie vrai si l'all[`e]le <num> existe
// mpa of integer values to pointers to CAlleles... no pointers created here
	public:
		int alleleExists(int num); // renvoie vrai si l'all[`e]le <num> existe
		std::map<int, CAllele *>alleles; // stockage des all[`e]les des genos complets (table [`a] acc[`e]s rapide !)
		std::map<int, CAllele *>galleles; // stockage des all[`e]les des genos (in)complets (table [`a] acc[`e]s rapide !)
		std::map<int, CAllele *>::const_iterator iter; // it['e]rateur sur les cl['e]s (de galleles [`a] priori)
		CLocusGP(size_t identifiant, std::string locName); // constructeur
		~CLocusGP();
		size_t getNumber(); // renvoie le nombre d'all[`e]les dans les genos complets
		unsigned long int getEffective(int num); // renvoie l'effectif de l'all[`e]le <num> in full genotypes
		unsigned long int getgEffective(int num); // renvoie l'effectif de l'(g)all[`e]le <num> if there are full diploid genotypes
		unsigned long int getggEffective(int num); //(migraine) renvoie l'effectif de l'(g)all[`e]le <num> even if only in incomplete genotypes
		void resetIterator(); // positionne iter au debut de alleles
		void resetgIterator(); // positionne iter au debut de galleles
		int getNext(); // va chercher les types alleliques presents dans alleles (renvoie -1 [`a] la fin)
		int getgNext(); // va chercher les types alleliques presents dans galleles (renvoie -1 [`a] la fin)
		size_t declareAllele(int num); // cr['e]e ou incr['e]mente l'all[`e]le <num> ; renvoie le nombre courant d'all[`e]les
		size_t declaregAllele(int num); // cr['e]e ou incr['e]mente l'(g)all[`e]le <num> ; renvoie le nombre courant d'(g)all[`e]les
		int AlleleIdentif(int num);
};

// classe d['e]riv['e] d'un locus pour le stockage des matrices g['e]notypiques
// PAS UTILISEE POUR L'INSTANT
//class Cmatrice_geno : public CLocus{
//private:
//      Cmatrice_geno(int pop, int loc);
//      ~Cmatrice_geno();
//}; //class matrice_geno;


// classe de stockage d'un typage
class CTypage {
	public:
		bool valid;
		int gene1; //
		int gene2; //
		int minAllele(); // renvoie le plus petit des deux all[`e]les
		int maxAllele(); // renvoie le plus grand des deux all[`e]les
};

// classe de stockage d'un individu
class CIndividual {
	private:
		std::string _indName; // nom de l'individu
		std::vector<CTypage>_typages; // vecteur de typages, dans l'ordre des locus
		std::vector<CTypage>_gtypages; // vecteur de typages, dans l'ordre des locus
	public:
		CIndividual(std::string indName, size_t nbLoc); // constructeur
		~CIndividual(); // destructeur
		void addTypage(int allele1, int allele2); // ajoute un typage diploide
		void addgTypage(int allele1, int allele2); // ajoute un typage diploide event. incomplet
		void addTypage(int allele1); // ajoute un typage haploide
		void addEmptyTypage(); // ajoute un typage vide
		bool isValid(size_t iLoc); // renvoie la validit['e] du typage au locus <iLoc>
		int getTypage(size_t iLoc,int all); // renvoie le typage...
		int getMinAllele(size_t iLoc); // renvoie le plus petit des deux all[`e]les au locus <iLoc>
		int getMaxAllele(size_t iLoc); // renvoie le plus grand des deux all[`e]les au locus <iLoc>
		std::string indName(); // renvoie le nom de l'individu
};

// classe de stockage d'une population
class CPopulation {
	public:
		std::vector<CLocusGP *>loci; // vecteur de stockage des locus et des all[`e]les compt['e]s pour cette pop
		std::vector<CIndividual *>inds; // vecteur de stockage des individus
		std::string popName(); // renvoie le nom de la pop (= nom du dernier individu)
		CPopulation(std::vector<CLocusGP *>lociSource); // constructeur
		~CPopulation(); // destructeur
		size_t addIndividual(std::string indName); // ajoute un individu ; renvoie le nombre courant d'individus
		CIndividual *lastIndividual(); // renvoie le dernier individu de la pop (courant lors du parsing)
};




// classe de cr['e]ation et de comptage de g['e]notypes
class CGenotypes { //complete genotypes only!! (NB null allele is info, >00)
	private:
		std::map<size_t, unsigned long int>genotypes; // compteur de tous les g['e]notypes existants (type,effectifs)
		bool genotypeExists(size_t genotype); // renvoie vrai si le g['e]notype <genotype> existe
		std::map<size_t, unsigned long int>::const_iterator iter; // it['e]rateur sur les cl['e]s
		size_t _sum;
	public:
		size_t declareGenotype(size_t genotype); // cr['e]e ou incr['e]mente le g['e]notype <genotype> ; renvoie le nb de g['e]notypes
		void clear();
		unsigned long int getEffective(ssize_t genotype); // renvoie l'effectif du g['e]notype <genotype>
		unsigned long int getEffective(size_t genotype); // renvoie l'effectif du g['e]notype <genotype>
		void fillGenotypes(size_t iLoc, CPopulation *pop,char coding); // remplit le truc avec les machins fournis [`a] la construction
		void resetIterator(); // reset l'it['e]rateur "iter" de "genotypes" [`a] begin()
		ssize_t getNext(); // renvoie la cl['e] suivant, -1 si end() atteint
		size_t getSum(); // renvoie le nombre total de g['e]notypes
		size_t getNumber(); // nombre de g['e]notypes diff['e]rents
    size_t getMaxIdx(); //largest genotype
		void printValues(CGenotypes *model, std::string *output, unsigned effectiveWidth, int foldLines); // imprime les effectifs des g['e]notypes
};

class CGenobilocus {
private:
		std::map<size_t, CGenotypes>::const_iterator iter; // it['e]rateur sur les cl['e]s
public:
        std::map<size_t,CGenotypes> mapmap; // (genotype locus1, CGeno=(genotype locus2,effectif 1x2))
        CGenotypes marginal;
		void resetIterator(); // reset l'it['e]rateur "iter" de "genotypes" [`a] begin()
		ssize_t getNext(); // renvoie la cl['e] suivant, -1 si end() atteint
	    char coding1,coding2; //nb de chiffres par genotype
        std::vector<std::vector<unsigned long int> > tabule();
        size_t getMinDim();
        size_t getcolNbr();
        size_t getligNbr();
        CGenobilocus() {coding1=0;coding2=0;marginal.clear();}
};

// classe de stockage des donn['e]es d'un fichier genepop
class CFichier_genepop {
	public:
    std::vector<char>coding; //nb de chiffres par genotype
		std::string fileName; // nom du fichier lu
		std::string fileTitle; // titre du fichier (1[`e]re ligne)
		std::vector<CLocusGP *>loci; // vecteur de locus et des all[`e]les compt['e]s pour l'ensemble du fichier
		std::vector<CPopulation *>pops; // vecteur des populations
		std::vector<std::vector<double> >coord_pop;
		CFichier_genepop(std::string fn); // constructeur minimaliste
		~CFichier_genepop(); // destructeur
        int checkName(); // checks that the file (or file +".txt") can be opened. CAN change fileName to add ".txt"
		size_t addPopulation(); // ajoute une population ; renvoie le nombre courant de populations
		CPopulation *lastPop(); // renvoie la derni[`e]re population (courante lors du parsing)
		int parseFile(); // lit et stocke le fichier genepop <filename> ; renvoie True si succ[`e]s
		CGenobilocus read_bilocus(size_t pop,size_t loc1,size_t loc2); // lecture pair loci pour une pop
		int affiche_nb_alleles(); // affichage ['e]cran nb alleles/locus
		size_t addLocus(std::string locName); // ajoute un locus au vecteur de locus ; renvoie le nombre courant de locus
		int createFichierIN();//cree le fichier.IN a partir des donn['e]es lu dans le fichier genepop
		int extract_coord_pop(); //remplit coord_pops
		int computeCheckWriteDistMat(const char file[]); // ecrit une demi matrice
};



template <typename Type>
void enligne(std::vector<Type> &v, std::ostream& stream, int coding) {
//.fill(' ') ne semble pas marcher avec stringstream (un autre character marche...)
//chercher de l'info sur stringstream
//colCount = output->length(); et donc une chaine etait necess pour compter le nb de colonne
// ici on vire la limite de 255 colonnes... la simplif est drastique
	for (typename std::vector<Type>::iterator ii=v.begin();ii<v.end();ii++) stream<<std::left<<std::setw(coding)<<(*ii);
}


extern CFichier_genepop *fichier_genepop;


#endif
