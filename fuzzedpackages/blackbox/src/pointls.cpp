#include <fstream>
#ifdef NO_R_CONSOLE
#include <iostream>
#endif
#include "pointls.h"
#define R_NO_REMAP
#include "R.h" // pour Rprintf

Cpointls::Cpointls(Cpointls& ptls) {
    this->xy.resize(0);
    for (std::vector<std::vector<ioType> >::iterator ii=ptls.xy.begin();ii!=ptls.xy.end();ii++)
        this->xy.push_back((*ii));
    std::stringstream ost(std::stringstream::in | std::stringstream::out);
    ost<<"Copy of "<<ptls.pointlsFileName;
    ost>>pointlsFileName;
    ost.clear();
}

Cpointls::Cpointls(double* xyarray,int *nrow,int *ncol) { // better code? :
  fittedparamnbr=(*ncol)-1;
  this->xy.resize(*nrow);
  for (int ii=0;ii<(*nrow);ii++) {
    this->xy[ii].resize(0);
    for (int jj=0;jj<(*ncol);jj++) {
      this->xy[ii].push_back(xyarray[ii*(*ncol)+jj]);
#ifdef NO_R_CONSOLE
#else
//      Rprintf("%d %d %f \n",ii,jj,this->xy[ii].back());
#endif
    }
  }
  std::stringstream ost(std::stringstream::in | std::stringstream::out);
  ost>>pointlsFileName;
  ost.clear();
}


int Cpointls::read_pointls(std::string filename) {
	std::stringstream strstr(std::stringstream::in | std::stringstream::out);
    std::string buf;
    std::vector<ioType> tempv;
    ioType templd;
    char chr;
    pointlsFileName=filename;

	std::ifstream inFile(filename.c_str(), std::ios::in);
  	if (! inFile.is_open()) {
#ifdef NO_R_CONSOLE
  		std::cerr << "Unable to open file '" << filename.c_str() << "'" << std::endl;
  		if (batchDebug) std::cin.get();
        inFile.clear();
  		exit(-1);
#else
  		throw Rcpp::exception("Unable to open file.");
  		//Rf_error("Unable to open file... I exit");
#endif
  	}
  	while (!inFile.eof()) {
        std::getline(inFile,buf);
        strstr.str(buf);
        chr=strstr.peek();
        if (isdigit(chr)) {//ligne commence par nombre
            tempv.resize(0);
            while (!strstr.eof()) {
//truc bizarre si pas le test suivant: blancs en fin de ligne c'est "lu" mais il ne reste que l'ancienne valeur dans templd
            	if (isdigit(strstr.peek())) {
                   strstr>>templd;
//std::cout<<"!!"<<templd;getchar();
            	   tempv.push_back(templd);
                } else strstr.get();
            }
            xy.push_back(tempv);
        }
        strstr.clear();
    }
    fittedparamnbr=xy[0].size()-1;
    return 0;
}

int Cpointls::minuslogLTologL() { //tr?s rustique mais...
    //from positive to negative
    for (std::vector<std::vector<ioType> >::iterator ii=xy.begin();ii!=xy.end();ii++) (*ii).back()*=-1.;
return 0;
}

int Cpointls::select_columns(int fittedp) { // primitif : garde les fittedp 1e cols et la derni?re
    for (std::vector<std::vector<ioType> >::iterator ii=xy.begin();ii!=xy.end();ii++) {
        ii->erase(ii->begin()+fittedp,ii->end()-1);
    }
return 0;
}

int Cpointls::selectTop() { // suppose derniere col est y
   ioType extremum=xy[0].back();
   for (unsigned int ii=1;ii<xy.size();ii++) {
           if (xy[ii].back()>extremum) extremum=xy[ii].back();    //ll
   }
#ifdef NO_R_CONSOLE
std::cout<<extremum<<std::endl;
#endif
   ioType threshold=extremum+10.;  //mll
   std::vector<std::vector<ioType> >::iterator ii=xy.begin();
   while (ii<xy.end()) {
//ostream_vector((*ii),std::cout);getchar();
//  std::cout<<int(ii-xy.begin());
           if (ii->back()>threshold) ii=xy.erase(ii); else ii++;   //mll
   }
#ifdef NO_R_CONSOLE
   std::cout<<"selectTop() selected "<<xy.size()<<" points"<<std::flush;
#endif

return 0;
}

