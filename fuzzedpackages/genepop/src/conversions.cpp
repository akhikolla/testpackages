/***************************************************************************
@ F. Rousset 2005-2006

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
#include <sstream>
#include "GenepopS.h" // seulement pour Zegenepopsound...
#include "genepop.h"
#include "tools.h"

#include <libgen.h>
#ifdef COMPATIBILITYRCPP
#include <Rcpp.h>
#endif
#include "conversions.h"
using namespace std;


//conversion to their input file formats: options 7.1--7.4-------------------------------------

int conversion(char Indic){
  /** Indic: 1 => Fstat; 2 =>  Biosys L; 3 => Biosys N; 4 => linkdos; 5 => diploidisation of haploid data
  6=> relabelling alleles; 7,8=> one pop per indiv; 9 => haploid sampling from diploid data **/
  //    char d[MAX_LINE_SIZE];
  //	char toto[MAX_LINE_SIZE + 1];
  string buf,typst;
  stringstream stst(stringstream::in | stringstream::out);
  string outName,numName;
  string::size_type pos,lpos;
  size_t nb_sam=fichier_genepop->pops.size();
  size_t nb_locus=fichier_genepop->loci.size();
  string corrgp_file=fichier_genepop->fileName;
  std::size_t iLoc, pop_ctr;
  int maxAll=0, //sur les loci
    all,typ;
  vector<CIndividual *>::iterator ind_ptr;
  vector<vector<int> >renum;
  for (vector<CLocusGP *>::iterator ii=fichier_genepop->loci.begin();ii<fichier_genepop->loci.end();ii++) {
    if((all=(*ii)->galleleMax)>maxAll) maxAll = all;
  }
  size_t randAllele=0;
  size_t lhap=0;
  
  
  noR_cout<<endl;
  
  if(nb_locus < 1){
#ifdef COMPATIBILITYRCPP
    // Rcpp::Rcerr<<"At least one locus is required. Check your input file"<<endl;
#else
    cerr<<"At least one locus is required. Check your input file"<<endl;
    if (cinGetOnError) cin.get();
#endif
    
    genepop_exit(0, "At least one locus is required. Check your input file");
  } // ELSE
  // ouverture fichier entr�e
  ifstream inFile(corrgp_file.c_str(), ios::in);
  if(!inFile.is_open()){
#ifdef COMPATIBILITYRCPP
    // Rcpp::Rcerr<<"Cannot open "<<corrgp_file.c_str()<<endl;
#else
    cerr<<"Cannot open "<<corrgp_file.c_str()<<endl;
    if (cinGetOnError) cin.get();
#endif
    
    genepop_exit(-1, "Cannot open file");
  }
  // ouverture fichier sortie
  outName=corrgp_file; // passer w a wp
  std::string pathFile = corrgp_file.substr(0, corrgp_file.find_last_of("/\\")+1); // add by jimmy
  std::string nameFile = corrgp_file.substr(corrgp_file.find_last_of("/\\")+1); // add by jimmy
  if (Indic==1) outName+=".DAT";
  else if (Indic==2 || Indic==3) outName+=".BIO";
  else if (Indic==4) outName+=".LKD";
  else if (Indic==5) {outName=pathFile+"D"+nameFile;}
  else if (Indic==9) {outName=pathFile+"H"+nameFile;}
  else if (Indic==6) {outName=pathFile+"N"+nameFile;}
  else if (Indic==7 || Indic==8 ) {outName=pathFile+"I"+nameFile;}
  ofstream wpOut(outName.c_str(), ios::out);
  if(!wpOut.is_open()){
#ifdef COMPATIBILITYRCPP
    // Rcpp::Rcerr<<"Cannot open file "<<outName<<endl;
#else
    cerr<<"Cannot open file "<<outName<<endl;
    if (cinGetOnError) cin.get();
#endif
    
    genepop_exit(0, "Cannot open file ");
  }
  // ecriture prambule fichier sortie
  if (Indic==1) {
    wpOut<<nb_sam<<" "<<nb_locus<<" "<<maxAll<<" 2"<<endl;
    for (vector<CLocusGP *>::iterator ii=fichier_genepop->loci.begin();ii<fichier_genepop->loci.end();ii++) {
      wpOut<<(*ii)->locusName<<endl;
    }
  } else if (Indic==2 || Indic==3) {
    wpOut<<fichier_genepop->fileTitle<<endl;
    wpOut<<"NOTU="<<nb_sam<<",NLOC="<<nb_locus<<",NALL="<<maxAll<<",CRT;\n";
    size_t maxLengthNomLocus=0;
    for (vector<CLocusGP *>::iterator ii=fichier_genepop->loci.begin();ii<fichier_genepop->loci.end();ii++) {
      if((*ii)->locusName.length()>maxLengthNomLocus) maxLengthNomLocus=(*ii)->locusName.length();
    }
    wpOut<<"("<<nb_locus<<"(1X,A"<<maxLengthNomLocus<<"/))\n";
    for (vector<CLocusGP *>::iterator ii=fichier_genepop->loci.begin();ii<fichier_genepop->loci.end();ii++) {
      wpOut<<" "<<(*ii)->locusName<<endl;
    }
  } else if (Indic==4) {
    wpOut<<nb_sam<<" "<<nb_locus<<endl;
    for (vector<CPopulation *>::iterator ii=fichier_genepop->pops.begin();ii<fichier_genepop->pops.end();ii++) {
      wpOut<<(*ii)->popName()<<endl;
      wpOut<<(*ii)->inds.size()<<endl;
    }
    for (vector<CLocusGP *>::iterator ii=fichier_genepop->loci.begin();ii<fichier_genepop->loci.end();ii++) {
      wpOut<<(*ii)->locusName<<" "<<(*ii)->galleleMax<<endl;
    }
    wpOut<<left;
    
  } else if (Indic==5) {
    wpOut<<"Diploidisation of "<<corrgp_file.c_str()<<" ("<<fichier_genepop->fileTitle<<")"<<endl;
    for (vector<CLocusGP *>::iterator ii=fichier_genepop->loci.begin();ii<fichier_genepop->loci.end();ii++) {
      wpOut<<(*ii)->locusName<<endl;
    }
  } else if (Indic==9) {
    wpOut<<"Random sampling of haploid genotypes from diploid ones in "<<corrgp_file.c_str()<<" ("<<fichier_genepop->fileTitle<<")"<<endl;
    for (vector<CLocusGP *>::iterator ii=fichier_genepop->loci.begin();ii<fichier_genepop->loci.end();ii++) {
      wpOut<<(*ii)->locusName<<endl;
    }
  } else if (Indic==7 || Indic==8) {
    wpOut<<"'Individualisation' of 'pop' data for "<<corrgp_file.c_str()<<" ("<<fichier_genepop->fileTitle<<")"<<endl;
    for (vector<CLocusGP *>::iterator ii=fichier_genepop->loci.begin();ii<fichier_genepop->loci.end();ii++) {
      wpOut<<(*ii)->locusName<<endl;
    }
  } else if (Indic==6) {
    wpOut<<fichier_genepop->fileTitle<<endl;
    for (vector<CLocusGP *>::iterator ii=fichier_genepop->loci.begin();ii<fichier_genepop->loci.end();ii++) {
      wpOut<<(*ii)->locusName<<endl;
    }
    // remplit tableau des renumerotations d'alleles et ecriture dans fichier .num
    numName=corrgp_file+".NUM";
    ofstream numfile(numName.c_str(), ios::out);
    if(!numfile.is_open()){
      RnoR_cerr<<"Cannot create file "<<numName<<endl;
      if (cinGetOnError) cin.get();
      genepop_exit(0, "Cannot open file ");
    }
    numfile<<"Old File: "<<corrgp_file<<"\nNew file: "<<outName<<"\n\n("<<fichier_genepop->fileTitle;
    numfile<<")\n\nAlleles have been renumbered\n(For each locus: first line = old names; second line = new names)\n";
    renum.resize(nb_locus);
    // pour un acces plus rapide que les maps, je fais des vectors... ce n'est pas tr�s �l�gant vu leur dimension
    for (iLoc=0;iLoc<nb_locus;iLoc++) {
      numfile<<"\nLocus: "<<fichier_genepop->loci[iLoc]->locusName<<endl;
      numfile<<" "<<fichier_genepop->loci[iLoc]->galleles.size()<<" alleles:\n    "<<left;
      renum[iLoc].resize(size_t(fichier_genepop->loci[iLoc]->galleleMax+1));
      renum[iLoc][0]=0; // commode de mettre l'allele nul dans la table pour les calculs sur 'typ' ci dessous
      int localit=1;
      int allelecoding=1+std::max(fichier_genepop->coding[iLoc]/2,fichier_genepop->coding[iLoc]%4);   // 2 3 4 6 => 3 4 3 4
      for(map<int,CAllele*>::iterator ii=fichier_genepop->loci[iLoc]->galleles.begin();
          ii!=fichier_genepop->loci[iLoc]->galleles.end();ii++) {
        numfile<<setw(allelecoding)<<ii->first;
        renum[iLoc][size_t(ii->first)]=localit;
        localit++;
      }
      numfile<<"\n    ";
      for(map<int,CAllele*>::iterator ii=fichier_genepop->loci[iLoc]->galleles.begin();
          ii!=fichier_genepop->loci[iLoc]->galleles.end();ii++) {
        numfile<<setw(allelecoding)<<renum[iLoc][size_t(ii->first)];
      }
      numfile<<"\n    ";
    }
    numfile<<"\nNormal ending.";
    numfile.close();
  }
  if (Indic==2) wpOut<<"\nSTEP DATA:\nDATYP=1,ALPHA;\n(A4,1X,"<<nb_locus<<"(2A1,1X))\n";
  else if (Indic==3) wpOut<<"\nSTEP DATA:\nDATYP=1;\n(A4,1X,"<<nb_locus<<"(2I2,1X))\n";
  //	inFile.get(toto, MAX_LINE_SIZE); inFile.ignore(1); // skip file title
  getline(inFile,buf); // skip file title
  if (inFile.eof()) {
    RnoR_cerr << "The file '" << corrgp_file.c_str() << "' exists but does not even contain a title line" << endl;
#ifdef COMPATIBILITYRCPP
#else
    if (cinGetOnError) cin.get();
#endif
    
    genepop_exit(-1, "The file does not even contain a title line");
  }
  do { // next skips loc names
    //		inFile.get(toto, MAX_LINE_SIZE); buf = toto; inFile.ignore(1);
    getline(inFile,buf);
    if (inFile.eof()) {
#ifdef COMPATIBILITYRCPP
      // Rcpp::Rcerr << "(!) From conversion(): The file '" << corrgp_file.c_str() << "' exists but does not contain data" << endl;
#else
      cerr << "(!) From conversion(): The file '" << corrgp_file.c_str() << "' exists but does not contain data" << endl;
      if (cinGetOnError) cin.get();
#endif
      
      genepop_exit(-1, "(!) From conversion(): The file exists but does not contain data");
    }
    lpos =min(buf.find('\r'), buf.length());
    buf=buf.substr(0, lpos); // remove a possible final \r
    while ((buf[0] == ' ') || (buf[0] == '\t')) {buf.erase(0, 1);} // vire les blancs initiaux
#ifdef TESTBLANK
    if (buf.length()==0) goto next_line;
#else
    if (buf.length()==0) {
#ifdef COMPATIBILITYRCPP
      // Rcpp::Rcerr << "The file '" << corrgp_file.c_str() << "' contains a blank line. Check this file." << endl; //FER 08/2006
#else
      cerr << "The file '" << corrgp_file.c_str() << "' contains a blank line. Check this file." << endl; //FER 08/2006
      if (cinGetOnError) cin.get();
#endif
      
      genepop_exit(-1, "The file contains a blank line. Check this file.");
    }
#endif
    stst.clear();
    stst.str(buf);
    stst>>typst; // lit un mot dans la stst
    if (cmp_nocase(typst, "POP") != 0) { // pas un "pop" (mais un nom de locus peut commencer par pop)
      pos = buf.find(',');
      if (pos == string::npos) {}
      else { while (pos!=string::npos)
      {buf.erase(0, pos+1); pos = buf.find(',');}
      }
    }
  } while (cmp_nocase(typst, "POP") != 0);
  // a ce point on a trouv� la premiere pop...
  pop_ctr=1;
  if (Indic==2 || Indic==3) wpOut<<"POP # "<<pop_ctr<<" ("<<fichier_genepop->pops[pop_ctr-1]->popName()<<")\n";
  else if (Indic==5 || Indic==6 || Indic==9) {wpOut<<"pop\n"; ind_ptr=fichier_genepop->pops[pop_ctr-1]->inds.begin();}
  else if (Indic==8) {ind_ptr=fichier_genepop->pops[pop_ctr-1]->inds.begin();}
  do {
    //		inFile.get(toto, MAX_LINE_SIZE); buf = toto; inFile.ignore(1);
    getline(inFile,buf);
    lpos =min(buf.find('\r'), buf.length());
    buf=buf.substr(0, lpos);
    while ((buf[0] == ' ') || (buf[0] == '\t')) {buf.erase(0, 1);} //shiny adds an initial blank !!
    pos = buf.find(',');  // rajout FER pour verif qu'on n'est pas sur la ligne d'un genotype
    if (cmp_nocase(buf.substr(0,3), "POP") == 0 && pos == string::npos) { // begin subsequent pop
      pop_ctr++;
      if (Indic==2 || Indic==3) wpOut<<"NEXT\nPOP # "<<pop_ctr<<" ("<<fichier_genepop->pops[pop_ctr-1]->popName()<<")\n";
      else if (Indic==5 || Indic==6 || Indic==9 ) {
        wpOut<<"pop\n";
        ind_ptr=fichier_genepop->pops[pop_ctr-1]->inds.begin();
      } else if (Indic==8) {ind_ptr=fichier_genepop->pops[pop_ctr-1]->inds.begin();}
    } else { // individu, pas pop
      if (Indic==1) wpOut<<" "<<pop_ctr<<"  ";
      else if (Indic==2 || Indic==3) wpOut<<"     ";
      else if (Indic==5 || Indic==6 || Indic==9) {
        // if the "pop"s are incorrectly (re-)read here, the ind_ptr goes beyon range =>crash
        if((*ind_ptr)->indName().size()>0) wpOut<<" "<<(*ind_ptr)->indName()<<", "; else wpOut<<"     , ";
      } else if (Indic==7) {
        wpOut<<"Pop"<<endl; //each individual in a pop
        wpOut<<" "<<fichier_genepop->pops[pop_ctr-1]->inds.back()->indName()<<", "; // coordinates of POP
      } else if (Indic==8) {
        wpOut<<"Pop"<<endl; //each individual in a pop
        wpOut<<" "<<(*ind_ptr)->indName()<<", "; // coordinates of INDIV
      }
      buf.erase(0, pos+1);
      iLoc=0;
      do {
        while ((buf[0] == ' ') || (buf[0] == '\t')) {buf.erase(0, 1);} // vire les blancs
        while (buf.length()==0) { // si on est la c'est qu'on a pas tous les loci sur une ligne
          //                    inFile.get(toto, MAX_LINE_SIZE); buf = toto; inFile.ignore(1);
          getline(inFile,buf);
          lpos =min(buf.find('\r'), buf.length());
          buf=buf.substr(0, lpos);
          while ((buf[0] == ' ') || (buf[0] == '\t')) {buf.erase(0, 1);} // vire les blancs
        } //� ce point on doit avoir qq chose
        lpos =min(buf.find(' '), min(buf.find('\t'), buf.length()));
        if (Indic==1 || Indic==3) wpOut<<buf.substr(0, lpos).c_str()<<" ";
        else if (Indic==2) {
          typ=atoi(buf.substr(0, lpos).c_str());
          //cout<<atoi(buf.substr(0, lpos).c_str())<<" "<<typ<<" "<<typ/100<<" "<<char(64+typ/100)<<typ%100<<" "<<char(64+typ%100)<<" "<<char(64+(typ%100));getchar();
          if (typ/100==0) wpOut<<" "; else wpOut<<char(64+typ/100);
          if (typ%100==0) wpOut<<" "; else wpOut<<char(64+typ%100);
          wpOut<<" ";
        } else if (Indic==4) {
          typ=atoi(buf.substr(0, lpos).c_str());
          if (lpos==4) {
            wpOut<<setw(2)<<typ/100<<" "<<setw(2)<<typ%100<<" ";
          } else if (lpos==6) {
            wpOut<<setw(3)<<typ/1000<<" "<<setw(3)<<typ%1000<<" ";
          } else if (lpos==2) {
            wpOut<<setw(2)<<typ<<" ";
          } else if (lpos==3) {
            wpOut<<setw(3)<<typ<<" ";
          }
        } else if (Indic==5) {
          if (lpos>3) wpOut<<buf.substr(0, lpos).c_str()<<" ";  // already diploid
          else wpOut<<buf.substr(0, lpos).c_str()<<buf.substr(0, lpos).c_str()<<" "; //haploid to diploid
        } else if (Indic==9) {
          if (lpos<4) wpOut<<buf.substr(0, lpos).c_str()<<" ";  // already haploid
          else {
            lhap=lpos/2;
            randAllele=size_t(alea.randInt(1));
            wpOut<<buf.substr(lhap*randAllele, lhap).c_str()<<" "; //random haploid from diploid
          }
        } else if (Indic==7 || Indic==8) {
          wpOut<<buf.substr(0, lpos).c_str()<<" "; // no specific conversion in this case
        } else if (Indic==6) {
          typ=atoi(buf.substr(0, lpos).c_str());
          if (lpos==4) {
            typ=100*renum[iLoc][size_t(typ/100)]+renum[iLoc][typ%100];
            if (typ<10) wpOut<<"000";
            else if (typ<100) wpOut<<"00";
            else if (typ<1000) wpOut<<"0";
            wpOut<<typ<<" ";
          } else if (lpos==6) {
            typ=1000*renum[iLoc][size_t(typ/1000)]+renum[iLoc][typ%1000];
            if (typ<10) wpOut<<"00000";
            else if (typ<100) wpOut<<"0000";
            else if (typ<1000) wpOut<<"000";
            else if (typ<10000) wpOut<<"00";
            else if (typ<100000) wpOut<<"0";
            wpOut<<typ<<" ";
          } else if (lpos==2) {
            typ=renum[iLoc][size_t(typ)];
            if (typ<10) wpOut<<"0";
            wpOut<<typ<<" ";
          } else if (lpos==3) {
            typ=renum[iLoc][size_t(typ)];
            if (typ<10) wpOut<<"00";
            else if (typ<100) wpOut<<"0";
            wpOut<<typ<<" ";
          }
        }
        buf.erase(0,lpos);
        // r�cup�ration des typages
        iLoc ++;
      } while (iLoc<nb_locus);
      if (Indic==5 || Indic==6 || Indic==8 || Indic==9) ind_ptr++;
      wpOut<<endl;
      //cout<<endl;
    }
  } while (! inFile.eof());
  if (Indic==2 || Indic==3) wpOut<<"NEXT\nEND;\n";
  inFile.close();
  wpOut.close();
  
  noR_cout<<"Normal ending."<<endl;
  noR_cout<<"The file "<<outName<<" is ready"<<endl;
  
  if (!perf) ZeGenepopSound();
  if (pauseGP) {
    noR_cout<<"(Return) to continue"<<endl; getchar();
  }
  return 0;
} // main
