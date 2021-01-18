#include <Rcpp.h>
using namespace Rcpp;


void erasefirstspace (std::vector< std::string> &vecstr)
{
  for (std::string &stri : vecstr)
  {if (stri.size()>0)
    stri.erase(stri.begin(),stri.begin()+1);}
}

// [[Rcpp::export]]
Rcpp::DataFrame transactiontoBitmax( std::vector<std::string> transac, char deli)
{ int nb_individus = transac.size();

  std::vector <std::string> listename;
  std::unordered_map<std::string,int> maptr;
  for (std::string s :  transac)
  { std::string q ="";
    for (char c :s)
      { if (c!= deli) {q+= c;}
      else {maptr[q]; q.erase(q.begin(),q.end());}
      }
    maptr[q];
  }

  int nb_var = maptr.size();

  std::unordered_map<std::string,int>::iterator itmap;
  std::unordered_map<std::string,int>::iterator itendmap = maptr.end();
  int a = 0;
  for (itmap = maptr.begin();itmap != itendmap;itmap++)
  {itmap->second = a++;
    listename.push_back(itmap->first);}

  std::vector<std::vector<int>> mytab (nb_var, std::vector<int>(nb_individus));
  int i = 0;
  for (std::string s :  transac)
  { std::string q ="";
    for (char c :s)
      { if (c!= deli) {q+= c;}
      else {mytab[maptr[q]][i]=1; q.erase(q.begin(),q.end());}
      }
    mytab[maptr[q]][i++]=1;
  }

  Rcpp::DataFrame Dataout (mytab);
  Dataout.attr("names") = listename;

  return Dataout;
}


class regroupe_function;

struct freq
{ public:
  freq (std::string name, int nbf) : nom (name), support(nbf){};
  std::string nom;
  int support;
  std::vector<int>* indyces = new std::vector<int> (support);
  freq * fils = NULL;
  freq * frere = NULL;
  void (*parcours)(freq&,short * tab);
  ~freq() {};
};

class regroupe_function
{ public :
  void (*racine)(freq&, short*);
  void (*feuille)(freq&, short * );
  void (*frere)(freq&, short * );
  void (*frerefils)(freq&, short * );
  void (*fils)(freq&,short * );

};


int nbfreq = 0;
int Sup = 0;
regroupe_function repertoire;
int nbind = 0;
std::string curname = "root";


void extract_and_erase_set (freq&alpha,  std::vector<std::string> &namevalue,std::vector<int>& supvalue, std::vector<float> &relativesupvalue, int &sit)
{
  namevalue[sit]=alpha.nom;
  supvalue[sit]=alpha.support;
  relativesupvalue[sit]=float(alpha.support)/float(nbind);
  sit++;
  if (alpha.fils != NULL)
  {extract_and_erase_set(*alpha.fils,namevalue,supvalue,relativesupvalue,sit);
    delete alpha.fils;}
  if (alpha.frere != NULL)
  {extract_and_erase_set(*alpha.frere,namevalue,supvalue,relativesupvalue,sit);
    delete alpha.frere;}
  delete alpha.indyces;
}

void tri_tableau (std::vector<short*> &dataframept, std::vector<int> mytab, int nb_elem, std::vector<std::string> &listenoms)
{
  int i,j,min,imin,temp;
  short* tempt;
  std::string tempstr;
  for (i = 0; i<nb_elem-1;i++)
  {imin = i; min=mytab[i];
  for (j=i+1;j<nb_elem;j++) {if (mytab[j]<min) {min = mytab[j];imin=j;}}
  temp=mytab[imin];mytab[imin]=mytab[i];mytab[i]=temp;
  tempt = dataframept[imin]; dataframept[imin] = dataframept[i];dataframept[i]=tempt;
  tempstr = listenoms[imin]; listenoms[imin] = listenoms[i];listenoms[i]=tempstr;
  }
}

std::vector<short*> init_prefixtree (std::vector<std::vector<short>> &mydataframe,std::vector<std::string>& listenoms, int support )
{ int n = mydataframe.size();
  std::vector<int> somme_vec (n);
  std::vector<short*> mydataframept (n);
  int tokeep = 0;
  for (int h = 0; h <n; h++ )
  {int a = std::accumulate(mydataframe[h].begin(),mydataframe[h].end(),0);
    somme_vec[h]=a;
    mydataframept[h]=&mydataframe[h][0];
    if (a > support ){tokeep++;}}
  tri_tableau(mydataframept, somme_vec,n,listenoms);
  tokeep = n-tokeep;

  listenoms.erase(listenoms.begin(),listenoms.begin()+tokeep);
  listenoms.shrink_to_fit();
  mydataframept.erase(mydataframept.begin(),mydataframept.begin()+tokeep);
  mydataframept.shrink_to_fit();
  return mydataframept;
}


void root (freq& alpha, short * binaryvec)
{ int sum = 0;
  nbfreq++;
  for (int r=0; r<nbind;r++)
  {sum+= binaryvec[r];}
  freq *nouv = new freq (curname,sum);
  curname = ' ' +curname;
  std::vector<int>::iterator hs = nouv->indyces->begin();
  for (int r = 0; r< nbind; r++)
  { if(binaryvec[r])  *hs++= r;
  }
  nouv->parcours = repertoire.frere;
  nouv->frere = alpha.fils;
  alpha.fils->parcours(*alpha.fils,binaryvec);
  alpha.fils = nouv;
}

void leaf (freq& alpha,short * binaryvec)
{  int sum = 0;
  std::vector<int>::iterator it;
  std::vector<int>::iterator end = alpha.indyces->end();
  for (it=alpha.indyces->begin();it!=end;it++)
  { sum+= binaryvec[*it];
  }
  if (sum > Sup)
  { nbfreq++;
    std::string name = alpha.nom+curname;
    alpha.fils = new freq (name,sum);
    std::vector<int>::iterator hs = alpha.fils->indyces->begin();
    for (it=alpha.indyces->begin();it!=end;it++)
    { if(binaryvec[*it]) *hs++=*it;
    }
    alpha.parcours = repertoire.fils;
    alpha.fils->parcours = repertoire.feuille;
  }
}


void bro (freq& alpha, short * binaryvec)
{ int sum = 0;
  std::vector<int>::iterator it;
  std::vector<int>::iterator end = alpha.indyces->end();
  for (it=alpha.indyces->begin();it!=end;it++)
  { sum+= binaryvec[*it];
  }
  if (sum > Sup)
  { nbfreq++;
    std::string name = alpha.nom+curname;
    alpha.fils  = new freq (name,sum);
    std::vector<int>::iterator hs = alpha.fils->indyces->begin();
    for (it=alpha.indyces->begin();it!=end;it++)
    { if(binaryvec[*it]) *hs++=*it;
    }

    alpha.parcours = repertoire.frerefils;
    alpha.fils->parcours = repertoire.feuille;
  }
  alpha.frere->parcours(*alpha.frere,binaryvec);
}



void broson (freq& alpha, short * binaryvec)
{  int sum = 0;
  std::vector<int>::iterator it;
  std::vector<int>::iterator end = alpha.indyces->end();
  for (it=alpha.indyces->begin();it!=end;it++)
  { sum+= binaryvec[*it];
  }

  if (sum > Sup)
  { nbfreq++;
    std::string name = alpha.nom+curname;
    freq *nouv = new freq (name,sum);
    std::vector<int>::iterator hs = nouv->indyces->begin();
    for (it=alpha.indyces->begin();it!=end;it++)
    { if(binaryvec[*it]) *hs++=*it;
    }
    nouv->parcours = repertoire.frere;
    nouv->frere = alpha.fils;
    alpha.fils->parcours(*alpha.fils,binaryvec);
    alpha.fils = nouv;
  }
  alpha.frere->parcours(*alpha.frere,binaryvec);
}


void son (freq& alpha, short * binaryvec)
{  int sum = 0;
  std::vector<int>::iterator it;
  std::vector<int>::iterator end = alpha.indyces->end();
  for (it=alpha.indyces->begin();it!=end;it++)
  { sum+= binaryvec[*it];
  }

  if (sum > Sup)
  { nbfreq++;
    std::string name = alpha.nom+curname;
    freq *nouv = new freq (name,sum);
    std::vector<int>::iterator hs = nouv->indyces->begin();
    for (it=alpha.indyces->begin();it!=end;it++)
    { if(binaryvec[*it]) *hs++=*it;
    }

    nouv->parcours = repertoire.frere;
    nouv->frere = alpha.fils;
    alpha.fils->parcours(*alpha.fils,binaryvec);
    alpha.fils = nouv;
  }

}




// [[Rcpp::export]]

DataFrame prefrecset (std::vector<std::vector<short>> Bitmax, std::vector<std::string> varnames, float relativeSup)
{ nbind= Bitmax[0].size();
  Sup = nbind *relativeSup -1 ;
  std::vector<short*> datapt = init_prefixtree(Bitmax, varnames, Sup);
  Rcpp::Rcout <<"the Supportvalue is " << Sup+1  << std::endl;
  if(datapt.size() ==0) {
    std::vector<std::string> outempty;
    outempty.push_back("empty");
    Rcpp::Rcout << "Not a single frequent itemSet, please try with a lower relative minSup value. " << std::endl;
    return Rcpp::DataFrame::create(Rcpp::Named("no itemSet") = outempty);
  }

  short * binarynewvec = datapt[0];
  nbfreq = 1;
  freq racine ("",0);

  int suminit = 0;
  for (int r=0; r<nbind;r++)
  {suminit += binarynewvec[r];}
  freq* init = new freq (varnames[0], suminit);
  std::vector<int>::iterator it = init->indyces->begin();
  for (int h=0; h < nbind;h++)
  {if (binarynewvec[h]) *it++=h;}

  nbfreq++;
  repertoire.racine=&root;
  repertoire.frere = &bro;
  repertoire.feuille = &leaf;
  repertoire.frerefils = &broson;
  repertoire.fils = &son;
  racine.parcours = repertoire.racine;
  racine.fils = init;
  racine.fils->parcours =repertoire.feuille;
  int nbvar=datapt.size();
  Rcpp::Rcout << " Start PrefRec with  " << nbvar << " frequents variables "<<  std::endl;
  std::clock_t stimer, etimer;
  stimer = std::clock();
  for (int h =1; h<nbvar ;h++)
  {
    curname = varnames[h];
    binarynewvec = datapt[h];
    racine.parcours(racine, binarynewvec);

  }
  etimer = std::clock();
  double diffe = double(etimer-stimer)/double(CLOCKS_PER_SEC);
  Rcpp::Rcout <<" PrefRec end in " << diffe <<" sec " << std::endl;
  Rcpp::Rcout << " Number of frequent set:  " << nbfreq <<  std::endl;
  int sizeout = nbfreq+1;
  int sit = 0;
  std::vector<int> supportvalue (sizeout);
  std::vector<float> relativesupportvalue(sizeout);
  std::vector<std::string> namevalue (sizeout);
  extract_and_erase_set(racine,namevalue,supportvalue,relativesupportvalue,sit);

  return Rcpp::DataFrame::create(Rcpp::Named("id") = namevalue,
                                 Rcpp::Named("support") = supportvalue,
                                 Rcpp::Named("relative_support")=relativesupportvalue);

}



struct frek {
public :
  frek (): unit_name(), support() {};
  frek ( std::vector<std::string> un, int s) : unit_name(un), support(s) {};
  std::vector<std::string> unit_name;
  int support;
  ~ frek(){};
};

struct rules {
  rules (std::string a, std::string c, float cf): antecedant (a), consequent(c), confiance (cf){};
  std::string antecedant;
  std::string consequent;
  float confiance;
};

std::unordered_map<std::string,frek> Mapfrek;
std::list <rules> Ruleslistes;
float conf =0;

void extraction_rules ( std::vector<std::string>& antec, std::vector<std::string>& conseq , std::vector<float>& rulesvalues)
{
  std::list<rules>::iterator itrulesend = Ruleslistes.end();
  std::list<rules>::iterator itrules;
  int i=0;
  for (itrules = Ruleslistes.begin(); itrules!= itrulesend; itrules++,i++)
  {antec[i]=itrules->antecedant;
    conseq[i]=itrules->consequent;
    rulesvalues[i]=itrules->confiance;}

}

std::vector<std::string> inter_vec (std::vector<std::string> &b, std::vector<std::string> &a)
{
  std::vector<std::string> out;
  std::vector<std::string>::iterator it;
  for (std::string go:a)
  { it =std::find(b.begin(),b.end(),go);
    if (it!=b.end()) {out.push_back(go);} }
  return out;
}

std::string creaantecedant (std::vector<std::string> vecant, int a)
{ vecant.erase(vecant.begin()+a);
  std::string outn = "";
  for (std::string str:vecant)
  {outn+= str;}
  return outn;
}

void  rules_test(std::vector<std::vector<std::string>> candid, double sup_consequent , float Minconf, std::string consequent)
{ std::vector<std::vector<std::string>> transi_candid;
  size_t h;
  for ( h=0; h < candid.size()-1;h++)
  { transi_candid.erase(transi_candid.begin(),transi_candid.end());
    size_t g;
    for( g=h+1; g<candid.size();g++)
    {
      std::vector<std::string> neo = inter_vec(candid[h], candid[g]);

      if (neo.size() >0)
      {
        std::string antecedent = "";

        for (std::string astr : neo)
        {antecedent += astr;}


        float nouvconf;
        nouvconf= sup_consequent/Mapfrek[antecedent].support;
        if ( nouvconf> Minconf )     {
          transi_candid.push_back(neo);

          Ruleslistes.emplace_front(rules(antecedent,consequent,nouvconf));}}
      if (transi_candid.size()>1) {rules_test(transi_candid, sup_consequent, Minconf, consequent);}
    }}
};



void Gen_rules (frek& candid, std::string consequent, double sup_consequent, float Minconf )
{ int vt = candid.unit_name.size();
  std::vector < std::vector < std::string>> tab_unt_name;
  for (int i = 0; i <vt;i++)
  { std::string antecedent = creaantecedant(candid.unit_name,i);
    float nouvconf;
    nouvconf= sup_consequent / Mapfrek[antecedent].support;
    if (nouvconf > Minconf) {
      Ruleslistes.emplace_front(rules(antecedent,consequent,nouvconf));
      tab_unt_name.push_back(Mapfrek[antecedent].unit_name);
    }
  }
  if (tab_unt_name.size() > 1) {
    rules_test(tab_unt_name,sup_consequent,Minconf,consequent);}

};

void rootr (freq& alpha, short * binaryvec)
{   int sum = 0;
  nbfreq++;
  for (int r=0; r<nbind;r++)
  {sum+= binaryvec[r];}
  freq *nouv = new freq (curname,sum);
  std::vector<int>::iterator hs = nouv->indyces->begin();

  for (int r = 0; r< nbind; r++)
  { if(binaryvec[r])  *hs++= r;
  }
  std::vector<std::string> tr;
  tr.push_back(curname);
  Mapfrek[curname] = frek(tr,sum);
  nouv->parcours = repertoire.frere;
  nouv->frere = alpha.fils;

  alpha.fils->parcours(*alpha.fils,binaryvec);
  alpha.fils = nouv;
}

void leafr (freq& alpha,short * binaryvec)
{
  int sum = 0;
  std::vector<int>::iterator it;
  std::vector<int>::iterator end = alpha.indyces->end();

  for (it=alpha.indyces->begin();it!=end;it++)
  { sum+= binaryvec[*it];
  }

  if (sum > Sup)
  { nbfreq++;

    std::string name = alpha.nom+curname;
    freq *nouv = new freq (name,sum);
    std::vector<int>::iterator hs = nouv->indyces->begin();
    for (it=alpha.indyces->begin();it!=end;it++)
    { if(binaryvec[*it]) *hs++=*it;
    }
    std::vector<std::string> tr (Mapfrek[alpha.nom].unit_name);
    tr.push_back(curname);
    frek trfk (tr,sum);
    Mapfrek[name] = trfk;
    Gen_rules(trfk, name, sum,conf);
    alpha.fils = nouv;
    alpha.parcours = repertoire.fils;
    alpha.fils->parcours = repertoire.feuille;
  }
}


void bror (freq& alpha, short * binaryvec)
{   int sum = 0;
  std::vector<int>::iterator it;
  std::vector<int>::iterator end = alpha.indyces->end();
  for (it=alpha.indyces->begin();it!=end;it++)
  { sum+= binaryvec[*it];
  }

  if (sum > Sup)
  { nbfreq++;
    std::string name = alpha.nom+curname;
    freq *nouv = new freq (name,sum);
    std::vector<int>::iterator hs = nouv->indyces->begin();

    for (it=alpha.indyces->begin();it!=end;it++)
    { if(binaryvec[*it]) *hs++=*it;
    }

    std::vector<std::string> tr = Mapfrek[alpha.nom].unit_name;
    tr.push_back(curname);
    frek trfk (tr,sum);
    Mapfrek[name] = trfk;
    Gen_rules(trfk, name, sum,conf);
    alpha.fils = nouv;
    alpha.parcours = repertoire.frerefils;
    alpha.fils->parcours = repertoire.feuille;
  }
  alpha.frere->parcours(*alpha.frere,binaryvec);
}



void brosonr (freq& alpha, short * binaryvec)
{ int sum = 0;
  std::vector<int>::iterator it;
  std::vector<int>::iterator end = alpha.indyces->end();
  for (it=alpha.indyces->begin();it!=end;it++)
  { sum+= binaryvec[*it];
  }

  if (sum > Sup)
  { nbfreq++;
    std::string name = alpha.nom+curname;
    freq *nouv = new freq (name,sum);
    std::vector<int>::iterator hs = nouv->indyces->begin();
    for (it=alpha.indyces->begin();it!=end;it++)
    { if(binaryvec[*it]) *hs++=*it;
    }
    std::vector<std::string> tr = Mapfrek[alpha.nom].unit_name;
    tr.push_back(curname);
    frek trfk (tr,sum);
    Mapfrek[name] = trfk;
    Gen_rules(trfk, name, sum,conf);
    nouv->parcours = repertoire.frere;
    nouv->frere = alpha.fils;
    alpha.fils->parcours(*alpha.fils,binaryvec);
    alpha.fils = nouv;
  }
  alpha.frere->parcours(*alpha.frere,binaryvec);
}


void sonr (freq& alpha, short * binaryvec)
{   int sum = 0;
  std::vector<int>::iterator it;
  std::vector<int>::iterator end = alpha.indyces->end();
  for (it=alpha.indyces->begin();it!=end;it++)
  { sum+= binaryvec[*it];
  }

  if (sum > Sup)
  { nbfreq++;
    std::string name = alpha.nom+curname;
    freq *nouv = new freq (name,sum);
    std::vector<int>::iterator hs = nouv->indyces->begin();
    for (it=alpha.indyces->begin();it!=end;it++)
    { if(binaryvec[*it]) *hs++=*it;
    }

    std::vector<std::string> tr = Mapfrek[alpha.nom].unit_name;
    tr.push_back(curname);
    frek trfk (tr,sum);
    Mapfrek[name] = trfk;
    Gen_rules(trfk, name, sum,conf);
    nouv->parcours = repertoire.frere;
    nouv->frere = alpha.fils;
    alpha.fils->parcours(*alpha.fils,binaryvec);
    alpha.fils = nouv;
  }

}




// [[Rcpp::export]]
List prefrecrules (std::vector<std::vector<short>> Bitmax, std::vector<std::string> varnames,float relativeSup,  float Minconf)
{ nbfreq = 1;
  conf =  float(Minconf);
  Ruleslistes.erase(Ruleslistes.begin(),Ruleslistes.end());
  Mapfrek.erase(Mapfrek.begin(),Mapfrek.end());
  nbind= Bitmax[0].size();
  Sup = nbind *relativeSup -1 ;
  std::vector<short*> datapt = init_prefixtree(Bitmax, varnames, Sup);
  Rcpp::Rcout <<"the Supportvalue is " << Sup+1 << " and the conf value is " << conf << std::endl;

  if(datapt.size() ==0) {
    std::vector<std::string> outempty;
    outempty.push_back("empty");
    Rcpp::Rcout << "Not a single frequent itemSet, please try with a lower relative minSup value. " << std::endl;
    DataFrame freqset = Rcpp::DataFrame::create(Rcpp::Named("no itemset") = outempty);
    DataFrame rulesset = Rcpp::DataFrame::create(Rcpp::Named("no rules") = outempty);
    return Rcpp::List::create( Rcpp::Named("frequent_itemset")= freqset,
                               Rcpp::Named("confident_rules")=rulesset);

  }

  short * binarynewvec = datapt[0];
  freq racine ("",0);
  int suminit = 0;
  for (int r=0; r<nbind;r++)
  {suminit += binarynewvec[r];}
  std::string init_name = ' ' +varnames[0];
  freq* init = new freq (init_name, suminit);
  std::vector<int>::iterator it = init->indyces->begin();
  for (int h=0; h < nbind;h++)
  {if (binarynewvec[h]) *it++=h;}
  nbfreq++;
  std::vector<std::string> tr;
  tr.push_back(init_name);
  Mapfrek[init_name] = frek(tr,suminit);

  repertoire.racine=&rootr;
  repertoire.frere = &bror;
  repertoire.feuille = &leafr;
  repertoire.frerefils = &brosonr;
  repertoire.fils = &sonr;
  racine.parcours = repertoire.racine;
  racine.fils = init;
  racine.fils->parcours =repertoire.feuille;

  int nbvar=datapt.size();
  Rcpp::Rcout << " Start PrefRec with  " << nbvar << " frequents variables "<<  std::endl;
  std::clock_t stimer, etimer;
  stimer = std::clock();
  for (int h =1; h<nbvar ;h++)
  {
    curname = ' '+varnames[h];
    binarynewvec = datapt[h];
    racine.parcours(racine, binarynewvec);

  }
  etimer = std::clock();
  double diffe = double(etimer-stimer)/double(CLOCKS_PER_SEC);
  Rcpp::Rcout <<" prefrecRules ends in  " << diffe <<" secs " << std::endl;
  Rcpp::Rcout << " Number of frequent set : " << nbfreq <<"number of confidents rules " << Ruleslistes.size() <<   std::endl;
  int sizeout = nbfreq+1;
  int sit = 0;
  std::vector<int> supportvalue (sizeout);
  std::vector<float> relativesupportvalue(sizeout);
  std::vector<std::string> namevalue (sizeout);

  extract_and_erase_set(racine,namevalue,supportvalue,relativesupportvalue,sit);

  std::vector<std::string> vecconsequent (Ruleslistes.size());
  std::vector<std::string> vecantecedant (Ruleslistes.size());
  std::vector<float> veconf (Ruleslistes.size());

  extraction_rules(vecantecedant,vecconsequent,veconf);

  erasefirstspace(namevalue);
  erasefirstspace(vecantecedant);
  erasefirstspace(vecconsequent);

  DataFrame freqset = Rcpp::DataFrame::create(Rcpp::Named("id") = namevalue,
                                              Rcpp::Named("support") = supportvalue,
                                              Rcpp::Named("relative_support")=relativesupportvalue);

  DataFrame rulesset = Rcpp::DataFrame::create(Rcpp::Named("antecedant") = vecantecedant,
                                               Rcpp::Named("consequent") = vecconsequent,
                                               Rcpp::Named("confiance") = veconf);


  return Rcpp::List::create( Rcpp::Named("frequent_itemset")= freqset,
                             Rcpp::Named("confident_rules")=rulesset);
}



