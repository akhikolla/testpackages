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
#ifndef CT_TESTS_H
#define CT_TESTS_H
#include "GenepopS.h"
#include "genepop.h"

//ceux la pour le template below...
namespace NS_GG { // for genic-genotypic tests
extern   std::vector<std::vector<int> >types; // potentiellement virable mais laborieux
//extern   std::vector<std::vector<int> >atable;
//extern   int allmax;
}

class Cctable {
public:
   Cctable(); //empty
   Cctable(std::vector<std::vector<unsigned long int> > table);
   virtual ~Cctable();
   int print(std::ostream& ost);
   bool purgeZeros(bool genobool);
   bool verifInfo();
   size_t maxCellCount();
   std::vector<double>Proba_test();
   std::vector<double>G_test();
   std::vector<double>GG_test();
   void switchSP();
   void switchSP_GG();
//en genic-genotypique la magouille rapide sur 4 termes de Gobs est incorrecte qunad on la transpose � 8 termes:
// par ex un homoz va mener � deux fois le m�me terme du G
   double calc_Gobs();
   double calc_GGobs(); // calcule aussi table allelique
   void cumul(double& privfreq,long int& nbpriv,std::vector<double>& ssize); //pour private alleles
   double calc_alleleNbr_trend(); // for EricI
   double calc_GG_geneDiv_trend(); // for EricI
   double calc_G_geneDiv_trend(std::ostream& ost);//ost OPTIONAL // for CecileA
   double calc_G_geneDiv_trend();
   size_t get_nb_lig() {return nb_lig;}
   void filltypesGG(CGenotypes* allGenotypes,char coding);
   void cleanVar();
private:
   double calc_GG(); //
  std::vector<std::vector<int> >atable;
  std::vector<std::vector<unsigned long int> > ctable;//the table itself
  std::vector<std::vector<int> > typesGG;
   size_t nb_lig,nb_col; //dimensions de la table
   size_t total;
   std::vector<unsigned long int>ligmarg,colmarg;
   std::vector<std::vector<double> > expected; //so called expected values (denominators)
   size_t lig1,lig2,col1,col2;
   void choix();
};


//------------------------------------------------------------------------------
// fonction template overloaded. Voir int crunchLocTable(...,vector<vector<double> > * tabF) in CT_tests.cpp
// l'autre tabF sert a stocker les resultats des tests (double) dans les options 1 et 3,
// ce tabF set � stocker la table geno (unsigned long int) dans l'option 5.2
template <typename Type> // et ceci est pour pouvoir detourner l'output ailleurs que dans le fichier de sortie
int crunchLocTable(int statindic,size_t iLoc,Type& fichier_out,
                   std::vector<std::vector<unsigned long int> > * tabF) {
  using namespace NS_GG;
  //version option 5.2, 6.1-4, (8.1 ??)
  // on doit pouvoir utiliser �a pour zapper les fichier P_L_...
  // tabF sert a stocker la table geno dans l'option 5.2
  // pas reussi � faire un template propre pour cette option � cause de _warning_ sur copy de testresult (c'est le bordel si ou veut les d�clarer Type)
  bool genobool=false;
  if (statindic==2) genobool=true;
  else if (statindic!=1) noR_cout<<"Incorrect statindic in cruncLocTable()";
  //	int nb_sam=fichier_genepop->pops.size();
  char coding=fichier_genepop->coding[iLoc];
  int allelecoding=2+std::max(coding/2,coding%4);   // 2 3 4 6 => 4 5 4 5
  std::vector<unsigned long int>dummyvec;
  int ligmarg;
  ssize_t genotype;
  int allele;
  std::vector<std::vector<int> >table; //table pop X x genotypes counts et types alleliques
  std::vector<std::vector<unsigned long int> >toutable; //table pop X x genotypes counts et types alleliques
  CGenotypes *allGenotypes=NULL; // TOUS les gentoypes du locus courant dans le fichier
  std::vector<CGenotypes *>genopops; // gentoypes pour chaque population au locus courant
  std::vector<CGenotypes *>::iterator pG; // iterateurs sur les effectifs
  std::vector<CPopulation *>::iterator p,pp; // CPopulation contient les popName()
  std::vector<CGenotypes *>::iterator p1,p2; // CGenotypes contient toute l'info sauf les popName()
  // construire structure pour tests genotypiques versus genic
  if (genobool) { //genotypes
    allGenotypes = new CGenotypes;
    // instanciation des structures de stockage des gentoypes
    genopops.resize(fichier_genepop->pops.size());
    // purge de la structure de gentoypes globale
    allGenotypes->clear();
    // it�rations sur les populations pour le locus courant
    pG = genopops.begin(); // initialisation de l'it�rateur sur les populations
    for(p = fichier_genepop->pops.begin(); p != fichier_genepop->pops.end(); p++) {
      allGenotypes->fillGenotypes(iLoc, *p,coding); // remplissage de la structure de gentoypes cumul�s total sample
      (*pG) = new CGenotypes; // instanciation de la structure de stockage de gentoypes de la pop courante pour le locus courant
      (*pG)->clear(); // nettoyage
      (*pG)->fillGenotypes(iLoc, *p,coding);
      pG++;
    }
    fichier_out << "Pop       Genotypes:" << std::endl;
  } else fichier_out << "Pop       Alleles:" << std::endl;
  fichier_out << "          ";
  for(int tiret=0; tiret<24; tiret ++){fichier_out << "----------";} // 240 -
  fichier_out << std::endl<<"          ";
  // ecriture labels colonnes tables: genotypic else genic
  if (genobool) {
    std::vector<std::vector<int> >typesloc(2);
    typesloc[0].resize(0);typesloc[1].resize(0);   //vide le contenu
    allGenotypes->resetIterator();
    while ((genotype = allGenotypes->getNext()) >= 0) { //sur les genos complets seulement
      typesloc[0].push_back(int(minAllele(genotype,coding)));
      typesloc[1].push_back(int(maxAllele(genotype,coding)));
    }
    enligne(typesloc[0],fichier_out,allelecoding);    // def dans genepop.cpp. Ecriture destypes alleliques
    fichier_out<<std::endl;
    fichier_out<<"          ";
    enligne(typesloc[1],fichier_out,allelecoding);
  } else {
    fichier_genepop->loci[iLoc]->resetIterator();
    std::vector<int>dummyintvec(0);
    while ((allele = fichier_genepop->loci[iLoc]->getNext()) >= 0)   dummyintvec.push_back(allele);
    enligne(dummyintvec,fichier_out,allelecoding);
  }
  fichier_out<<"  Total"<<std::endl<<std::endl;
  // valuation et ecriture des lignes de la toutable  (toute, pairbool ou non)
  toutable.resize(0);
  if (genobool) {
    pG = genopops.begin();
    for(std::vector<CPopulation *>::iterator ii=fichier_genepop->pops.begin();ii<fichier_genepop->pops.end();ii++) {
      allGenotypes->resetIterator();
      dummyvec.resize(0);
      ligmarg=0;
      while ((genotype = allGenotypes->getNext()) >= 0) {
        dummyvec.push_back((*pG)->getEffective(genotype));
        ligmarg+=dummyvec.back();
      }
      fichier_out<<std::setw(10)<<(*ii)->popName().substr(0,9);
      enligne(dummyvec,fichier_out,allelecoding);
      fichier_out<<"   "<<ligmarg<<std::endl;
      toutable.push_back(dummyvec);
      pG++;
    }
    fichier_out<<std::endl;
    fichier_out<<"Total:    ";
    allGenotypes->resetIterator();
    dummyvec.resize(0);
    ligmarg=0;
    while ((genotype = allGenotypes->getNext()) >= 0) {
      dummyvec.push_back(allGenotypes->getEffective(genotype));
      ligmarg+=dummyvec.back();
    }
    enligne(dummyvec,fichier_out,allelecoding);
    fichier_out<<"   "<<ligmarg<<std::endl;
  } else { //NOT genobool
    for(std::vector<CPopulation *>::iterator ii=fichier_genepop->pops.begin();ii<fichier_genepop->pops.end();ii++) {
      fichier_genepop->loci[iLoc]->resetIterator();
      dummyvec.resize(0);
      ligmarg=0;
      while ((allele = fichier_genepop->loci[iLoc]->getNext()) >= 0) {
        dummyvec.push_back((*ii)->loci[iLoc]->getEffective(allele));
        ligmarg+=dummyvec.back();
      }
      fichier_out<<std::setw(10)<<(*ii)->popName().substr(0,9);
      enligne(dummyvec,fichier_out,allelecoding);
      fichier_out<<"   "<<ligmarg<<std::endl;
      toutable.push_back(dummyvec);
    }
    fichier_out<<std::endl;
    fichier_out<<"Total:    ";
    fichier_genepop->loci[iLoc]->resetIterator();
    dummyvec.resize(0);
    ligmarg=0;
    while ((allele = fichier_genepop->loci[iLoc]->getNext()) >= 0) {
      dummyvec.push_back(fichier_genepop->loci[iLoc]->getEffective(allele));
      ligmarg+=dummyvec.back();
    }
    enligne(dummyvec,fichier_out,allelecoding);
    fichier_out<<"   "<<ligmarg<<std::endl;
  } // fin genobool or not
  tabF->resize(toutable.size()+2); //+1 suffirait pour haploide mais forcerait � bidouiller les procedures suivant crunchloctable
  if (genobool) { // ecriture finale dans tabF
    size_t tabFit=0;
    for (std::vector<std::vector<unsigned long int> >::iterator ii=toutable.begin();ii<toutable.end();ii++) {
      (*tabF)[tabFit].resize((*ii).size()); // moui... toutes les lignes ont la m�me lgr, important pour tabFtotabM...
      copy((*ii).begin(),(*ii).end(),(*tabF)[tabFit].begin());
      tabFit++;
    }
    (*tabF)[tabFit].resize(0);
    (*tabF)[tabFit+1].resize(0);
    allGenotypes->resetIterator();
    while ((genotype = allGenotypes->getNext()) >= 0) { //sur les genos complets seulement
      (*tabF)[tabFit].push_back(minAllele(genotype,fichier_genepop->coding[iLoc]));
      (*tabF)[tabFit+1].push_back(maxAllele(genotype,fichier_genepop->coding[iLoc]));
    }
    for(pG = genopops.begin(); pG != genopops.end(); pG++) delete (*pG);
    delete allGenotypes;
  } else { //haploide
    size_t tabFit=0;
    for (std::vector<std::vector<unsigned long int> >::iterator ii=toutable.begin();ii<toutable.end();ii++) {
      (*tabF)[tabFit].resize((*ii).size());
      copy((*ii).begin(),(*ii).end(),(*tabF)[tabFit].begin());
      tabFit++;
    }
    (*tabF)[tabFit].resize(0);
    (*tabF)[tabFit+1].resize(0);
    fichier_genepop->loci[iLoc]->resetIterator();
    while ((allele = fichier_genepop->loci[iLoc]->getNext()) >= 0) {
      (*tabF)[tabFit].push_back(size_t(allele)); // >0
      (*tabF)[tabFit+1].push_back(0); //cf commentaire ci dessus sur taille tabF (du coup ca devient indic haploidie)
    }
  }
  return 0;
}

extern void (Cctable::*switchFnPtr)();

int crunchLocTable(int statindic,size_t iLoc,std::ostream& fichier_out,std::vector<std::vector<double> > * tabF);
int struc();
void initializeCTtests();
void cleanCTtests();

#endif
