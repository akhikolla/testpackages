/**@file RankCluster.cpp
 * * * * * * * * ** @brief Implementation of methods of the class @c RankCluster
 * * * * * * * * **/
 
 #include "RankCluster.h"
 
 #include <cstdlib>
 #include <algorithm>
 #include <string>
 #include <limits>
 #include <cmath>
 #include <map>
 #include <iostream>
 #include <ctime>
 #include <Rmath.h>
 
 using namespace std;
 using namespace Eigen;
 //using namespace Rcpp;
 
 //constructor
 RankCluster::RankCluster()
 {}
 
 //constructor
 RankCluster::RankCluster(std::vector<std::vector<int> > const& X,int g, vector<int> const& m, SEMparameters const& param)
                          : m_(m),n_(X.size()),d_(m.size()),g_(g),
                            data_(d_,vector<PartialRank> (n_)),
                            z_(n_),
                            mu_(d_,vector<vector<int> > (g_)),
                            p_(d_,vector<double>(g_)),
                            proportion_(g),
                            parameter_(param),
                            partial_(false),
                            dataOk_(true),
                            indexPb_(m.size())
 {
   //  try
   //  {
   // convert data in the good notation and create information about missing and partial
   conversion2data(X);
   //  }
   //  catch(string const& chaine)
   //    {dataOk_=false;}
   
 }
 
 //constructor
 RankCluster::RankCluster(vector<vector<int> > const& X, vector<int> const& m, SEMparameters const& param,
                          vector<double> const& proportion, vector<vector<double> > const& p, vector<vector<vector<int> > > const& mu)
                          : m_(m), n_(X.size()), d_(m.size()), g_(proportion.size()),
                            data_(d_,vector<PartialRank> (n_)),
                            z_(n_), mu_(mu), p_(p),
                            proportion_(proportion),
                            parameter_(param),
                            partial_(false),
                            dataOk_(true),
                            indexPb_(m.size())
 {
  //  try
  //  {
  // convert data in the good notation and create information about missing and partial
  conversion2data(X);
  //  }
  //  catch(string const& chaine)
  //    {dataOk_=false;}
 }
 
 //copy constructor
 RankCluster::RankCluster(RankCluster& rankClusterObject)
                          : m_(rankClusterObject.m()),
                            n_(rankClusterObject.n()),
                            d_(rankClusterObject.d()),
                            g_(rankClusterObject.g()),
                            data_(rankClusterObject.data()),
                            mu_(rankClusterObject.mu()),
                            p_(rankClusterObject.p()),
                            proportion_(rankClusterObject.proportion()),
                            parameter_(rankClusterObject.parameter()),
                            partial_(rankClusterObject.partial()),
                            dataOk_(rankClusterObject.dataOk())
 {
   //nothing to do
 }
 
 
 //destructor
 RankCluster::~RankCluster()
 {
   // nothing
 }
 
 void RankCluster::readRankingRank(vector<vector<int> > const& X, int const& dim, int const& j, vector<int> const& indM)
 {
   //initialization
   int indiceElement = 0;
   data_[dim][j].isNotFull = false;
   
   //multi dim rank temporary
   vector<vector<int> > temp(m_[dim]+1);
   
   for(int i = indM[dim]; i < indM[dim+1]; i++)
   {
     temp[X[j][i]].push_back(indiceElement+1);
     indiceElement++;
   }
   
   //vector containing index of partial element
   vector<int> partialIndex;
   
   int skip = 0;
   //index 0 is for missing, we don't manage in this loop
   for(int i = 1; i < (int) temp.size(); i++)
   {
     if(skip)
     {
       if(temp[i].size()!=0)
       {
         dataOk_ = false;
         indexPb_[dim].push_back(j+1);
         //throw string("Problem with data.");
       }
       
       skip--;
     }
     else
     {
       //tied case
       if(temp[i].size() > 1)
       {
         data_[dim][j].isNotFull = true;
         partial_ = true;
         skip = temp[i].size() - 1;
         data_[dim][j].missingData.push_back(temp[i]);
         
         vector<int> missingIndex(temp[i].size());
         
        for(int ii = 0; ii < (int) temp[i].size(); ii++)
          missingIndex[ii] = i + ii - 1;
         
         data_[dim][j].missingIndex.push_back(missingIndex);
       }
       else
       {
         //normal case
         if(temp[i].size()==1)
          data_[dim][j].rank[i-1] = temp[i][0];
         else//temp[i].size=0//partial case
          partialIndex.push_back(i-1);
       }
     }
   }
   
   //problem with the data : index of 0 et element at missing position don't match
   if(temp[0].size() != partialIndex.size())
   {
     dataOk_ = false;
     indexPb_[dim].push_back(j+1);
     //throw string("Problem with data.");
   }
   
   //add partial
   if(temp[0].size() !=0 )
   {
     data_[dim][j].isNotFull = true;
     partial_ = true;
     data_[dim][j].missingData.push_back(temp[0]);
     data_[dim][j].missingIndex.push_back(partialIndex);
   }
 }
 
 void RankCluster::conversion2data(vector<vector<int> > const& X)
 {
   //size of a row of X
   vector<int> indM(d_+1,0);
   for(int i = 0; i < d_; i++)
   indM[i+1] = indM[i] + m_[i];
   
   //resize data
   for(int i = 0; i < d_; i++)
    for(int j = 0; j < n_; j++)
      data_[i][j].rank.resize(m_[i]);
   
   //begin the read of the data row by row
   for(int j = 0; j < n_; j++)
   {
     //dim by dim
     for(int dim(0);dim<d_;dim++)
     {
       //read rank j of dim dim
       readRankingRank(X,dim,j,indM);
     }
   }
   
 }
 
 
 void RankCluster::initialization()
 {
   double alea;
   
   //zik initialization with multinomial of equal proba
   if(g_ != 1)
   {
     for(int i = 0; i < n_; i++)
     {
       alea=runif(0.,1.);
       for(int j = 0; j < g_; j++)
       if((alea>(double) j/g_) & (alea<(double) (j+1)/g_))
       {
         z_[i]=j;
         break;
       }
     }
   }
   else
   {
     for(int i(0); i < n_; i++)
     z_[i]=0;
   }
   
   
   //mu & p  initialization
   for(int k = 0; k < d_; k++)
   {
     for(int i = 0; i < g_; i++)
     {
       //initialization of p_ with double between 0.5 and 1
       alea=(double) runif(0.5,1.);
       p_[k][i]=alea;
       //initialization of mu_ with alea rank of size m_
       mu_[k][i].resize(m_[k]);
       for(int j = 0; j < m_[k]; j++)
        mu_[k][i][j]=j+1;
       random_shuffle(mu_[k][i].begin(),mu_[k][i].end(),randWrapper);
     }
   }
   
   //proportion initialization
   for(int i = 0; i < n_; i++)
   proportion_[z_[i]]++;
   
   for(int i = 0; i < g_; i++)
   proportion_[i]/=(double) n_;
   
   
   //partial data and order of presentation initialization
   for(int dim = 0; dim < d_; dim++)
   {
     vector<int> rankTemp(m_[dim]);
     for(int i = 0; i < m_[dim]; i++)
      rankTemp[i]=i+1;
     for(int ind = 0; ind < n_; ind++)
     {
       //initialization of y
       random_shuffle(rankTemp.begin(),rankTemp.end(),randWrapper);
       data_[dim][ind].y=rankTemp;
       
       if(data_[dim][ind].isNotFull)
       {
         for(int ii = 0; ii < (int) data_[dim][ind].missingIndex.size(); ii++)
         {
           //initialization of Partial Rank
           vector<int> rankTemp2(data_[dim][ind].missingIndex[ii]);
           random_shuffle(rankTemp2.begin(),rankTemp2.end(),randWrapper);
           
           for(int iii = 0; iii < (int) data_[dim][ind].missingData[ii].size(); iii++)
           data_[dim][ind].rank[rankTemp2[iii]] = data_[dim][ind].missingData[ii][iii];
         }
       }
     }
   }
   
   
   indexPartialData_=vector<vector<int> > (d_);
   for(int dim = 0; dim < d_; dim++)
   {
     for(int ind = 0; ind < n_; ind++)
     {
       if(data_[dim][ind].isNotFull)
       indexPartialData_[dim].push_back(ind);
     }
   }
   
   vector<vector<vector<int> > > donneesPartiel(d_);
   for(int dim = 0; dim < d_; dim++)
   for(vector<int>::iterator it = indexPartialData_[dim].begin(); it != indexPartialData_[dim].end(); it++)
   donneesPartiel[dim].push_back(data_[dim][*it].rank);
   
   
   //sauvegarde initialisation
   output_.initialPartialRank = donneesPartiel;
   output_.initialP = p_;
   output_.initialZ = z_;
   output_.initialMu = mu_;
   output_.initialProportion = proportion_;
 }
 
 
 void RankCluster::SEstep()
 {
   //simulation of order of presentation for each dimension
   for(int dim = 0; dim < d_; dim++)
   gibbsY(dim);
   
   //simulation of z
   zSimulation();
   
   //simulation of partial rank for each dimension
   for(int dim = 0; dim < d_; dim++)
   gibbsX(dim);
 }
 
 
 void RankCluster::gibbsY(int indexDim)
 {
   double p1(0),p2(0),alea(0);
   set<int>::iterator itset;
   
   //rank 1 2..m
   vector<int> yTemp(m_[indexDim]);
   for(int j = 0; j < m_[indexDim]; j++)
   yTemp[j] = j+1;
   
   for(int ind = 0; ind < n_; ind++)
   {
     //Gibbs sampling
     vector<int> y(m_[indexDim]),y2(m_[indexDim]),y1(m_[indexDim]);
     
     //initialization of p1 and y1
     y=yTemp;
     random_shuffle(y.begin(),y.end(),randWrapper);//simulation of alea rank
     y1=y;
     p1=probaCond(data_[indexDim][ind].rank,y1,mu_[indexDim][z_[ind]],p_[indexDim][z_[ind]]);
     
     //start of iteration
     for(int iter(0);iter<parameter_.nGibbsSE[indexDim];iter++)
     {
       for(int k(0);k<m_[indexDim]-1;k++)
       {
         //new y to test (old y with inversion of 2 adjacents elements)
         y2=y;
         y2[k]=y[k+1];
         y2[k+1]=y[k];
         
         //compute the probability of accept the changement of y
         p2=probaCond(data_[indexDim][ind].rank,y2,mu_[indexDim][z_[ind]],p_[indexDim][z_[ind]]);
         
         alea=(double) runif(0,(p1+p2));
         
         if(alea<p2)//changement acceptation
         {
           y=y2;
           p1=p2;
           y1=y;
         }
         else
         y=y1;
       }
     }
     data_[indexDim][ind].y=y;
   }
   
 }
 
 
 void RankCluster::zSimulation()
 {
   
   if(g_ != 1)
   {
     double alea(0),sumTik(0);
     vector<double> lim(g_+1,0),tik(g_);
     
     for(int ind(0);ind<n_;ind++)
     {
       //computation of the probability to belong to each cluster
       for(int k(0);k<g_;k++)
       tik[k]=1;
       
       sumTik=0;
       for(int k(0);k<g_;k++)
       {
         for(int l(0);l<d_;l++)
         tik[k]*=probaCond(data_[l][ind].rank,data_[l][ind].y,mu_[l][k],p_[l][k]);
         tik[k]*=proportion_[k];
         sumTik+=tik[k];
       }
       
       for(int i(0);i<g_;i++)
       tik[i]/=sumTik;
       
       //z follow a multinomial law of parameter tik
       for(int k(1);k<g_+1;k++)
       lim[k]=lim[k-1]+tik[k-1];
       
       alea=(double) runif(0.,1.);
       
       for(int j(0);j<g_;j++)
       if((lim[j]<=alea)&&(alea<=lim[j+1]))
       {
         z_[ind]=j;
         break;
       }
     }    
   }
   else
   {
     for(int i=0; i < (int) z_.size(); i++)
     z_[i]=0;
   }
   
   
   
 }
 
 void RankCluster::gibbsX(int indexDim)
 {
   double p1(0),p2(0),alea(0);
   
   for(int ind = 0; ind < n_; ind++)
   {
     if(data_[indexDim][ind].isNotFull)
     {
       //Algorithme de Gibbs
       vector<int> x(m_[indexDim]),x2(m_[indexDim]),x1(m_[indexDim]);
       
       //initialisation de mu et p pour Gibbs
       x = data_[indexDim][ind].rank;
       x1 = x;
       p1 = probaCond(x1,data_[indexDim][ind].y,mu_[indexDim][z_[ind]],p_[indexDim][z_[ind]]);
       
       for(int iter = 0; iter < parameter_.nGibbsSE[indexDim]; iter++)
       {
         for(int ii = 0; ii < (int) data_[indexDim][ind].missingIndex.size(); ii++)
         {
           for(int i = 0; i < (int) data_[indexDim][ind].missingIndex[ii].size()-1; i++)
           {
             //nouveau x à tester, ancien x auquel on inverse 2 éléments partiels
             x2 = x;
             x2[data_[indexDim][ind].missingIndex[ii][i]] = x[data_[indexDim][ind].missingIndex[ii][i+1]];
             x2[data_[indexDim][ind].missingIndex[ii][i+1]] = x[data_[indexDim][ind].missingIndex[ii][i]];
             
             p2 = probaCond(x2,data_[indexDim][ind].y,mu_[indexDim][z_[ind]],p_[indexDim][z_[ind]]);
             
             alea = (double) runif(0.,p1+p2);
             
             if(alea < p2)//acceptation du changement
             {
               x = x2;
               p1 = p2;
               x1 = x;
             }
             else
             x = x1;
           }
         }
         
       }
       data_[indexDim][ind].rank = x;
     }
   }
 }
 
 
 void RankCluster::Mstep()
 {
   //MAJ proportion
   for(int k(0);k<g_;k++)
   {
     proportion_[k]=0;
     for(int ind(0);ind<n_;ind++)
     {
       if(z_[ind]==k)
       proportion_[k]++;
     }
     proportion_[k]/=(double) n_;
     if(proportion_[k]==0)
     throw string("NON CONVERGENCE DE L'ALGORITHME : a proportion is equal to 0");
   }
   
   //simulation of mu for ezach dim and cluster
   for(int dim(0);dim<d_;dim++)
   {
     for(int numCl(0);numCl<g_;numCl++)
     simuM(dim,numCl);
   }
   
 }
 
 
 void RankCluster::simuM(int indexDim,int indCl)
 {
   long double p1(0),p2(0),alea(0),lnp1(0),lnp2(0),lnp1Plusp2(0);
   double s1(0),lim(-numeric_limits<double>::max());
   vector<int> tabFact(m_[indexDim]),comp(2);
   
   tabFact=tab_factorial(m_[indexDim]);
   
   vector<vector<int> > MU(parameter_.nGibbsM[indexDim]);
   vector<double> P(parameter_.nGibbsM[indexDim],0),L(parameter_.nGibbsM[indexDim],0);
   int indMu,compteur(0);
   int tailleComposante(0);
   
   
   for(int i(0);i<n_;i++)
   if(z_[i]==indCl)
   tailleComposante++;
   
   vector<double> G(tailleComposante),A_G(tailleComposante);
   
   
   map<int,vector<double> > muTeste;
   map<int,vector<double> >::iterator itmapmu;
   
   //Algorithme de Gibbs
   vector<int> mu(m_[indexDim]),mu2(m_[indexDim]),mu1(m_[indexDim]);
   
   //initialization of mu  and p
   mu=mu_[indexDim][indCl];
   mu1=mu;
   lnp1=0;
   
   //initial proba
   for(int ind(0);ind<n_;ind++)
   {
     if(z_[ind]==indCl)
     {
       p1=probaCond(data_[indexDim][ind].rank,data_[indexDim][ind].y,mu1,p_[indexDim][indCl]);
       lnp1+=log(p1);
     }
   }
   
   for(int iter(0);iter<parameter_.nGibbsM[indexDim];iter++)
   {
     //new mu
     for(int k(0);k<m_[indexDim]-1;k++)
     {
       //new mu to test
       mu2=mu;
       mu2[k]=mu[k+1];
       mu2[k+1]=mu[k];
       
       lnp2=0;
       
       //new proba
       for(int ind(0);ind<n_;ind++)
       {
         if(z_[ind]==indCl)
         {
           p2=probaCond(data_[indexDim][ind].rank,data_[indexDim][ind].y,mu2,p_[indexDim][indCl]);
           lnp2+=log(p2);
         }
       }
       
       //p1+p2
       //ln(p1+p2)~ln(
       long double diffln;
       if(lnp1>lnp2)
       {
         diffln=lnp2-lnp1;
         lnp1Plusp2=lnp1;
       }
       else
       {
         diffln=lnp1-lnp2;
         lnp1Plusp2=lnp2;
       }
       
       //taylor development of order for estiamting lnp1Plusp2
       for(int ordre(1);ordre<6;ordre++)
       {
         //* (long double) std::pow((int) -1,(int) ordre-1)
         int sign=1;
         if(ordre%2==1)
         sign=-1;
         lnp1Plusp2+=(long double) sign/ordre*exp(diffln*ordre);
       }
       
       alea=(long double) runif(0.,1.);
       
       // acceptaion of change or not
       if(alea<exp(lnp2-lnp1Plusp2))//accept the changement
       {
         mu=mu2;
         lnp1=lnp2;
         mu1=mu;
       }
       else
       mu=mu1;
       
     }//fin parcours mu
     MU[iter]=mu;
     
     //MAJ p
     indMu=rank2index(MU[iter],tabFact);
     itmapmu=muTeste.find(indMu);
     
     if(itmapmu==muTeste.end())//if we have already tested this mu, we do not redo computation
     {
       s1=0;
       compteur=0;
       //computation of G and A-G (number of good and bad comparison) for loglikelihood
       double somG(0),somA_G(0);
       for(int ind(0);ind<n_;ind++)
       {
         if(z_[ind]==indCl)
         {
           comp=comparaison(data_[indexDim][ind].rank,data_[indexDim][ind].y,MU[iter]);
           G[compteur]=comp[1];
           A_G[compteur]=comp[0]-comp[1];
           s1+=comp[0];
           somG+=G[compteur];
           somA_G+=A_G[compteur];
           compteur++;
         }
       }
       
       P[iter]=somG/s1;
       
       if((P[iter]!=0) & (P[iter]!=1))
       {
         //L[iter]+=(G[i]*log(P[iter])+A_G[i]*log(1-P[iter]));
         L[iter]=somG*log(P[iter])+somA_G*log(1-P[iter]);
       }
       else
       {
         if((P[iter]==0) & (somG==0))
         L[iter]=0;
         else
         {
           if((P[iter]==1) & (somA_G==0))
           L[iter]=0;
           else
           L[iter]=lim;
         }
       }
       /*vector<double> stock(3);
       stock[0]=somG;
       stock[1]=somA_G;
       stock[2]=s1;*/
       vector<double> stock(2);
       stock[0]=P[iter];
       stock[1]=L[iter];
       muTeste[indMu]=stock;//MAJ des mu test�
     }
     else
     {
       L[iter]=(itmapmu->second)[1];
       P[iter]=(itmapmu->second)[0];
     }
     
     p_[indexDim][indCl]=P[iter];
   }
   
   
   int indice(max_element(L.begin(),L.end())-L.begin());
   p_[indexDim][indCl]=P[indice];
   mu_[indexDim][indCl]=MU[indice];
   
 }
 
 
 
 typedef struct ListeMu ListeMu;
 void RankCluster::likelihood(vector<vector<vector<vector<int> > > > &listeMu,vector<vector<vector<double> > > &resP,vector<vector<double> > &resProp)
 {
   //we put the same mu together and make the mean of their parameters
   //double t1,t2,tL(0);
   
   //t1=clock();
   struct ListeMu;
   struct ListeMu
   {
     double compteur;//number of same mu
     std::vector<std::vector<std::vector<int> > > rangComplet;
     std::vector<std::vector<double> > p;
     std::vector<double> prop;
     ListeMu* suivant;
   };
   
   bool continuer(true),egaliteRang;
   
   
   ListeMu*headMu=new ListeMu;
   ListeMu*currMu=headMu;
   ListeMu*next=0;
   currMu->compteur=1;
   currMu->suivant=0;
   currMu->rangComplet=listeMu[0];
   currMu->p=resP[0];
   currMu->prop=resProp[0];
   
   int nbMu(1);
   
   for (int j(1);j<parameter_.maxIt-parameter_.burnAlgo;j++)//we see all the mu
   {
     continuer=true;
     currMu=headMu;
     while(continuer)
     {
       egaliteRang=true;
       //look if the j-th mu is the same that the current mu
       for(int J(0);J<d_;J++)
       for(int k(0);k<g_;k++)
       {
         for(int i(0);i<m_[J];i++)
         if(currMu->rangComplet[J][k][i]!=listeMu[j][J][k][i])
         {
           egaliteRang=false;
           break;
         }
       }
       
       
       if(egaliteRang)
       {
         //same mu
         currMu->compteur++;
         //we sum the proportion and p
         for(int compt1(0);compt1<g_;compt1++)
         {
           currMu->prop[compt1]+=resProp[j][compt1];
           for(int compt2(0);compt2<d_;compt2++)
           currMu->p[compt2][compt1]+=resP[j][compt2][compt1];
         }
         continuer=false;//no need to check other mu of the struct
       }
       else
       {
         //not the same mu
         if(currMu->suivant==0)
         { //if the current mu is the last, we add thz j-th mu in the struct
         nbMu++;
         continuer=false;
         next=new ListeMu;
         currMu->suivant=next;
         currMu=next;
         currMu->compteur=1;
         currMu->suivant=0;
         currMu->rangComplet=listeMu[j];
         currMu->prop=resProp[j];
         currMu->p=resP[j];
         }
         else
         currMu=currMu->suivant;//we will test the next mu
       }
     }
   }
   
   //t2=clock();
   //cout<<"Temps regroupement mu: "<<(double) (t2-t1)/CLOCKS_PER_SEC<<"s"<<endl;
   int compteur(0);
   
   //calcul logvraisemblance
   //if(parameter_.detail)
   //cout<<"Number of reference rank which must compute the log-likelihood: "<<nbMu<<endl;
   
   double Llast(-numeric_limits<double>::max()), L;
   
   
   vector<vector<vector<int> > > Y(d_,vector<vector<int> > (n_)),xPartialTemp(output_.initialPartialRank);
   vector<vector<vector<double> > > scoreTemp(output_.initialPartialRank.size());
   for(int ii = 0; ii < (int) scoreTemp.size(); ii++)
   {
     scoreTemp[ii].resize(output_.initialPartialRank[ii].size());
     for(int iii = 0; iii < (int) scoreTemp[ii].size(); iii++)
     {
       scoreTemp[ii][iii].resize(output_.initialPartialRank[ii][iii].size());
     }
   }
   
   //Now, we have the list of all the different Mu
   currMu=headMu;
   ArrayXXd tik(n_,g_);
   ArrayXXd probabilities(n_,g_);//estimate probability for an individual to belong to each cluster
   
   //for each mu, we will compute the associate log likelihood
   while(currMu->suivant!=0)
   {
     //if(parameter_.detail)
     //cout<<"*";
     
     //mean of the parameter
     for(int compt1(0);compt1<g_;compt1++)
     {
       currMu->prop[compt1]/=currMu->compteur;
       for(int compt2(0);compt2<d_;compt2++)
       currMu->p[compt2][compt1]/=currMu->compteur;
     }
     
     //compute the log likelihood
     //t1=clock();
     L=computeLikelihood(currMu->rangComplet,currMu->p,currMu->prop,tik,Y,xPartialTemp,probabilities,scoreTemp);
     //t2=clock();
     //tL+=t2-t1;
     
     if(L>Llast)
     {
       //the current mu has a better loglikelihood, we save the parameter
       Llast = L;
       mu_ = currMu->rangComplet;
       p_ = currMu->p;
       proportion_ = currMu->prop;
       output_.tik = tik;
       output_.L = L;
       output_.probabilities = probabilities;
       output_.partialRankScore = scoreTemp;
       for(int dim = 0; dim < d_; dim++)
       {
         for(int ind = 0; ind < n_; ind++)
         data_[dim][ind].y=Y[dim][ind];
         
         int compteur(0);
         for(vector<int>::iterator it = indexPartialData_[dim].begin(); it != indexPartialData_[dim].end(); it++)
         {
           data_[dim][*it].rank=xPartialTemp[dim][compteur];
           compteur++;
         }
       }
     }
     
     next=currMu->suivant;
     delete currMu;//delete the mu
     currMu=next;
     compteur++;
   }
   
   //the last mu of the struct
   //mean of the parameter
   for(int compt1 = 0; compt1 < g_; compt1++)
   {
     currMu->prop[compt1] /= currMu->compteur;
     for(int compt2 = 0; compt2 < d_; compt2++)
     currMu->p[compt2][compt1] /= currMu->compteur;
   }
   
   //if(parameter_.detail)
   //cout<<"*"<<endl;
   
   //compute log likelihood
   //t1=clock();
   L=computeLikelihood(currMu->rangComplet,currMu->p,currMu->prop,tik,Y,xPartialTemp,probabilities,scoreTemp);
   
   //t2=clock();
   //tL+=t2-t1;
   compteur++;
   
   if(L > Llast)
   {
     //the current mu has a better loglikelihood, we save the parameter
     Llast=L;
     mu_=currMu->rangComplet;
     p_=currMu->p;
     proportion_=currMu->prop;
     output_.tik=tik;
     output_.L=L;
     output_.probabilities=probabilities;
     output_.partialRankScore = scoreTemp;
     vector<int> compteurPartiel(d_,0);
     
     for(int dim = 0; dim < d_; dim++)
     {
       for(int ind = 0; ind < n_; ind++)
       data_[dim][ind].y=Y[dim][ind];
       
       int compteur(0);
       for(vector<int>::iterator it = indexPartialData_[dim].begin(); it != indexPartialData_[dim].end(); it++)
       {
         data_[dim][*it].rank=xPartialTemp[dim][compteur];
         compteur++;
       }
     }
   }
   delete currMu;//delete the last mu
   
   
   //if(parameter_.detail)
   //cout<<"Computing time for log-likelihood approximation: "<<(double) tL/CLOCKS_PER_SEC<<"s ("<<(double) tL/CLOCKS_PER_SEC/compteur<<"s per mu)."<<endl;
   
 }
 
 
 //LL gibbs
 double RankCluster::computeLikelihood(vector<vector<vector<int> > > const& mu,vector<vector<double> > const& p,
 vector<double> const& proportion,ArrayXXd &tik,vector<vector<vector<int> > > &Y,vector<vector<vector<int> > > &xTemp, ArrayXXd &probabilities,
 vector<vector<vector<double> > > &score)
 {
   long double p1(0),p2(0),p1x(0),p2x(0),alea(0),l(0),li(0);
   double div((double) (parameter_.nGibbsL-parameter_.burnL));
   vector<int> compteur(d_,0);
   vector<int> x1,x2;
   
   //objet pour stockage de calcul pour éviter répétition
   ArrayXXd proba1(d_,g_),proba2(d_,g_);
   ArrayXd proba1X(g_),proba2X(g_);
   
   //store proportion
   ArrayXd prop(g_);
   for(int i = 0; i < g_; i++)
   prop(i)=proportion[i];
   ArrayXXd propb(1,g_);
   propb=prop.transpose();
   
   vector<vector<int> > yTemp(d_),x(d_);
   vector<vector<vector<int> > > scoreCount(d_);
   for(int j = 0; j < d_; j++)
   {
     //génération rang 1 2 3 ..m par dimension
     yTemp[j].resize(m_[j]);
     for(int i = 0; i < m_[j]; i++)
     yTemp[j][i]=i+1;
     
     //initialize score
     scoreCount[j].resize(m_[j],vector<int> (m_[j],0));
   }
   
   //start estimation
   for(int ind = 0; ind < n_; ind++)
   {
     //cout<<"ind "<<ind<<endl;
     vector<vector<int> > y(d_),y2(d_),y1(d_);
     li=0;
     //algorithme de Gibbs pour simuler yi
     
     //initialisation de y et p pour Gibbs
     y=yTemp;
     for(int j = 0; j < d_; j++)
     {
       for(int jj = 0; jj < m_[j]; jj++)
       for(int jjj = 0; jjj < m_[j]; jjj++)
       scoreCount[j][jj][jjj]=0;
       
       random_shuffle(y[j].begin(),y[j].end(),randWrapper);//permutation de 1 2 3 ..m
       x[j]=data_[j][ind].rank;
     }
     
     y1 = y;
     
     for(int k = 0; k < g_; k++)
     {
       tik(ind,k) = 0;
       for(int j = 0; j < d_; j++)
       proba1(j,k) = probaCond(x[j],y1[j],mu[j][k],p[j][k]);
     }
     
     p1 = (long double) (propb*proba1.colwise().prod()).sum();
     proba2 = proba1;
     
     //start gibbs for sample ind
     for(int iter = 0; iter< parameter_.nGibbsL; iter++)
     {
       
       /*simulation des y*/
       for(int J = 0; J < d_; J++)
       {
         for(int K = 0; K < m_[J]-1; K++)
         {
           //"état" à tester (inversion de 2 éléments adjacents)
           y2 = y;
           y2[J][K] = y[J][K+1];
           y2[J][K+1] = y[J][K];
           
           for(int k = 0; k < g_; k++)//tester un stockage des proba calcul� pour �viter r�p�tition de calculs dans la boucle
           proba2(J,k) = probaCond(x[J],y2[J],mu[J][k],p[J][k]);
           
           
           p2 = (long double) (propb*proba2.colwise().prod()).sum();
           
           alea = (long double) runif(0.,p1+p2);//unif(0,p1+p2)(double) 
           
           if(alea < p2)//acceptation du changement de y
           {
             y[J] = y2[J];
             p1 = p2;
             proba1.row(J) = proba2.row(J);
             y1[J] = y[J];
           }
           else//on ne modifie pas y
           {
             y[J] = y1[J];//rajout J
             proba2.row(J) = proba1.row(J);
           }
         }
       }
       //y_i est mis à jour
       
       
       
       /*simulation des x_i^j qui sont partiels*/
       for(int J = 0; J < d_; J++)
       {
         if(data_[J][ind].isNotFull)//simulation de xi si partiel
         {
           x1 = x[J];
           proba1X = proba1.row(J);
           p1x = (proba1X * prop).sum();
           for(int kk = 0; kk < (int) (data_[J][ind].missingIndex).size()-1; kk++)
           {
             for(int k = 0; k < (int) (data_[J][ind].missingIndex[kk]).size()-1; k++)//Gibbs sur les x
             {
               //nouveau x à tester
               x2 = x[J];
               x2[data_[J][ind].missingIndex[kk][k]] = x[J][data_[J][ind].missingIndex[kk][k+1]];
               x2[data_[J][ind].missingIndex[kk][k+1]] = x[J][data_[J][ind].missingIndex[kk][k]];
               
               for(int l = 0; l < g_; l++)
               proba2X(l) = probaCond(x2,y[J],mu[J][l],p[J][l]);
               
               p2x = (proba2X * prop).sum();
               
               alea = (double) runif(0.,p1x+p2x);
               
               if(alea < p2)//acceptation du changement
               {
                 x[J] = x2;
                 p1x = p2x;
                 proba1X = proba2X;
                 x1 = x[J];
               }
               else
               x[J] = x1;
             }
           }
           
           proba1.row(J) = proba1X;
         }
         
       }
       
       if(iter >= parameter_.burnL)
       {
         ArrayXd calculInter(g_);
         for(int cl = 0; cl < g_; cl++)
         {
           calculInter(cl)=1;
           for(int dim = 0; dim < d_; dim++)
           calculInter(cl) *= proba1(dim,cl);
           
           probabilities(ind,cl) = calculInter(cl);
           calculInter(cl) *= propb(cl);
         }
         
         double den(calculInter.sum());
         tik.row(ind) += (calculInter/den);
         
         li+=(long double) 1/den;
         
         //compute score partial rank
         for(int dim = 0; dim < d_; dim++)
         {
           if(data_[dim][ind].isNotFull)
           {
             for(int indElem = 0; indElem < m_[dim]; indElem++)
             scoreCount[dim][indElem][x[dim][indElem]-1]++;
           }
         }
       }//end if not burn
       
       
     }//fin du gibbs pour l'individu ind
     
     
     //l -= log(li*div);
     l -= log(li);
     
     tik.row(ind) /= div;
     probabilities.row(ind) /= div;
     
     //sauvegarde des nouveau y et x
     for(int j(0);j<d_;j++)
     {
       Y[j][ind]=y[j];
       if(data_[j][ind].isNotFull)
       {
         xTemp[j][compteur[j]]=x[j];
         for(int elem = 0; elem < m_[j]; elem++)
         score[j][compteur[j]][elem] = ( (double) (scoreCount[j][elem][x[j][elem]-1]) / (double) (parameter_.nGibbsL - parameter_.burnL) );
         compteur[j]++;
       }
       
     }
     
   }//fin boucle sur n
   
   l += n_ * log(div);
   return l;
 }
 
 
 //compute the final partition
 void RankCluster::computePartition()
 {
   if(g_>1)
   {//calcul partition
   double max;
   for(int ind(0);ind<n_;ind++)
   {
     max=output_.tik(ind,0);
     z_[ind]=0;
     for(int k(1);k<g_;k++)
     {
       if(output_.tik(ind,k)>max)
       {
         max=output_.tik(ind,k);
         z_[ind]=k;
       }
     }
     //cout<<z_[ind]<<" "<< output_.tik(ind,z_[ind])<<endl;
   }
   }
 }
 
 //return probability to belong to the final cluster
 ArrayXd RankCluster::probability() const
 {
   ArrayXd probability(n_);
   for(int ind(0);ind<n_;ind++)
   probability(ind)=output_.probabilities(ind,z_[ind]);
   
   
   return probability;
 }
 
 
 
 void RankCluster::computeDistance(vector<vector<double> > const& resProp,vector<vector<vector<double> > > const& resP,
 vector<vector<vector<vector<int> > > > const& resMu,vector<vector<int> > const& resZ,
 vector<vector<vector<vector<int> > > > const& resDonneesPartiel)
 {
   int const iterTotal(parameter_.maxIt-parameter_.burnAlgo);
   
   //initialization of container
   output_.distProp=vector<vector<double> >(iterTotal,vector<double> (g_));
   output_.distP=vector<vector<vector<double> > >(iterTotal,vector<vector<double> > (d_,vector<double> (g_)));
   output_.distMu=vector<vector<vector<int> > >(iterTotal,vector<vector<int> > (d_,vector<int> (g_)));
   output_.distZ=vector<double>(iterTotal);
   
   //compute the distance between the output parameters and parameters from each iteration
   for(int i(0);i<iterTotal;i++)
   {
     //distance between partition
     output_.distZ[i]=indiceRand(z_,resZ[i]);
     
     for(int cl(0);cl<g_;cl++)
     {
       //distance between proportion
       output_.distProp[i][cl]=pow(resProp[i][cl]-proportion_[cl],2);
       
       for(int dim(0);dim<d_;dim++)
       {
         //distance between p
         output_.distP[i][dim][cl]=pow(resP[i][dim][cl]-p_[dim][cl],2);
         
         //distance between mu
         output_.distMu[i][dim][cl]=distanceKendall(mu_[dim][cl],resMu[i][dim][cl]);
       }
     }
   }
   
   //distance between partial rank
   vector<vector<vector<int> > > distRangPartiel (iterTotal,vector<vector<int> > (d_));
   if(partial_)
   {
     for(int i(0);i<iterTotal;i++)
     {
       for(int dim(0);dim<d_;dim++)
       {
         int compteur(0);
         //for(int k(0);k<resDonneesPartiel[i][dim].size();k++)
         for(vector<int>::iterator it=indexPartialData_[dim].begin();it!=indexPartialData_[dim].end();it++)
         {
           distRangPartiel[i][dim].push_back(distanceKendall(data_[dim][*it].rank,resDonneesPartiel[i][dim][compteur]));
           compteur++;
         }
       }
     }
   }
   
   //changement de format
   vector<int> compteurElemPartiel(d_,0);
   output_.distPartialRank=vector<vector<vector<int> > > (resDonneesPartiel.size());
   vector<int> rangTemp(d_);
   
   for(int iter(0); iter < (int) distRangPartiel.size(); iter++)
   {
     for(int dim(0);dim<d_;dim++)
     compteurElemPartiel[dim]=0;
     
     for(int ind(0);ind<n_;ind++)
     {
       for(int dim(0);dim<d_;dim++)
       {
         if(data_[dim][ind].isNotFull)
         {
           rangTemp[dim]=distRangPartiel[iter][dim][compteurElemPartiel[dim]];
           compteurElemPartiel[dim]++;
         }
         else
         rangTemp[dim]=0;
       }
       output_.distPartialRank[iter].push_back(rangTemp);
     }
   }
 }
 
 void RankCluster::run()
 {
   convergence_=false;
   int nbTry(0);
   while(!convergence_ && nbTry<parameter_.maxTry)
   {
     try
     {
       //double t0,t1,t2,t3,tM(0),tSE(0);
       
       //if(parameter_.detail)
       //{
       //cout<<"##########################################################"<<endl;
       //cout<<"#  SEM-Gibbs Algorithm for multivariate partial ranking  #"<<endl;
       //cout<<"##########################################################"<<endl;
       //}
       //t0=clock();
       initialization();
       //t1=clock();
       
       //if(parameter_.detail)
       //cout<<"Initialization: "<<(double) (t1-t0)/CLOCKS_PER_SEC<<"s."<<endl;
       
       
       //objects for storing the estimated parameters at each iteration
       vector<int> indrang(g_);
       vector<vector<vector<double> > > resP(parameter_.maxIt-parameter_.burnAlgo,vector<vector<double> >(d_,vector<double> (g_)));
       vector<vector<double> > resProp(parameter_.maxIt-parameter_.burnAlgo,(vector<double> (g_)));
       vector<vector<int> > resZ(parameter_.maxIt-parameter_.burnAlgo,vector<int> (n_));
       vector<vector<vector<vector<int> > > > resMu(parameter_.maxIt-parameter_.burnAlgo,mu_);
       vector<vector<vector<vector<int> > > > resDonneesPartiel(parameter_.maxIt-parameter_.burnAlgo,output_.initialPartialRank);
       
       
       //algorithm
       //if(parameter_.detail)
       for(int iter(0);iter<parameter_.maxIt;iter++)
       {
         //if(parameter_.detail)
         //cout<<"*";
         
         //t2=clock();
         SEstep();
         //t3=clock();
         //tSE+=t3-t2;
         
         //t2=clock();
         Mstep();
         //t3=clock();
         //tM+=t3-t2;
         
         //we store the estimated parameters
         if(iter>=parameter_.burnAlgo)
         {
           
           for(int l(0);l<d_;l++)
           {
             for(int k(0);k<g_;k++)
             {
               if(p_[l][k]<0.5)
               {
                 p_[l][k]=1-p_[l][k];
                 inverseRang(mu_[l][k]);
               }
             }
           }
           
           for(int k(0);k<g_;k++)
           indrang[k]=rank2index(mu_[0][k],tab_factorial(m_[0]));
           
           //the first cluster must be the cluster with the more little index of mu
           tri_insertionMulti(indrang,proportion_,p_,mu_,z_,g_,d_,n_);//tri selon les mu pour que 2 3=3 2
           
           //store parameters
           resP[iter-parameter_.burnAlgo]=p_;
           resProp[iter-parameter_.burnAlgo]=proportion_;
           resMu[iter-parameter_.burnAlgo]=mu_;
           resZ[iter-parameter_.burnAlgo]=z_;
           
           for(int dim(0);dim<d_;dim++)
           {
             int compteur(0);
             for(vector<int>::iterator it=indexPartialData_[dim].begin();it!=indexPartialData_[dim].end();it++)
             {
               resDonneesPartiel[iter-parameter_.burnAlgo][dim][compteur]=data_[dim][*it].rank;
               compteur++;
             }
           }
           
         }//end storage
       }//end SEM
       
       //if(parameter_.detail)
       //cout<<endl<<endl<<"Loglikelihood estimation"<<endl;
       //t2=clock();
       //if(parameter_.detail)
       //{
       //cout<<"Computing time for SE step: "<<(double) tSE/CLOCKS_PER_SEC<<"s ( "<<(double) tSE/CLOCKS_PER_SEC/parameter_.maxIt<<"s per step)."<<endl;
       //cout<<"Computing time for M step: "<<(double) tM/CLOCKS_PER_SEC<<"s ( "<<(double) tM/CLOCKS_PER_SEC/parameter_.maxIt<<"s per step )."<<endl;
       //}
       
       //compute loglikelihood and choice of the best parameters
       likelihood(resMu,resP,resProp);
       //t3=clock();
       
       //compute the partition associated to the best result
       computePartition();
       
       //compute distance between estimated parameters at each iteration
       computeDistance(resProp,resP,resMu,resZ,resDonneesPartiel);
       
       //compute criterion
       output_.bic = BIC(output_.L,n_,2*g_*d_+g_-1);
       
       output_.icl = output_.bic;
       
       output_.entropy = ArrayXd(n_);
       for(int i = 0; i < n_; i++)
       {
         output_.entropy(i)=0;
         for(int j = 0; j < g_; j++)
         {
           if(output_.tik(i,j)!=0)
           output_.entropy(i)-=output_.tik(i,j)*std::log(output_.tik(i,j));
         }
         output_.icl += 2*output_.entropy(i);
       }
       
       //if(parameter_.detail)
       //{
       //cout<<"Total computing time : "<<(double) (t3-t0)/CLOCKS_PER_SEC<<"s"<<endl;
       //}
       
       
       //if the p<0.5, we invert the associate mu and put p=1-p
       for(int j = 0; j < d_; j++)
       for(int k = 0; k < g_; k++)
       {
         if(p_[j][k]<0.5)
         {
           p_[j][k] = 1-p_[j][k];
           inverseRang(mu_[j][k]);
         }
       }
       
       convergence_=true;
       
     }//end try
     catch(string const& chaine)
     {convergence_=false;}
     nbTry++;
   }
 }
 
 
 //LL gibbs
 void RankCluster::estimateCriterion(double &L,double &bic,double &icl)
 {
   
   /*initialisation partial rank and order of presentation*/
   //partial data and order of presentation initialization
   for(int dim(0);dim<d_;dim++)
   {
     vector<int> rankTemp(m_[dim]);
     for(int i(0);i<m_[dim];i++)
     rankTemp[i]=i+1;
     for(int ind(0);ind<n_;ind++)
     {
       //initialization of y
       random_shuffle(rankTemp.begin(),rankTemp.end(),randWrapper);
       data_[dim][ind].y=rankTemp;
       
       if(data_[dim][ind].isNotFull)
       {
         for(int i = 0; i < (int) data_[dim][ind].missingData.size(); i++)
         {
           //initialization of Partial Rank
           vector<int> rankTemp2(data_[dim][ind].missingIndex[i]);
           random_shuffle(rankTemp2.begin(),rankTemp2.end(),randWrapper);
           
           for(int ii = 0; ii < (int) data_[dim][ind].missingData[i].size(); ii++)
           data_[dim][ind].rank[rankTemp2[ii]] = data_[dim][ind].missingData[i][ii];
         }
         
       }
     }
   }
   
   /*log likelihood computation*/
   ArrayXXd tik(n_,g_);
   long double p1(0),p2(0),p1x(0),p2x(0),alea(0),li(0);
   double div((double) (parameter_.nGibbsL-parameter_.burnL));
   vector<int> compteur(d_,0);
   vector<int> x1,x2;
   
   //objet pour stockage de calcul pour éviter répétition
   ArrayXXd proba1(d_,g_),proba2(d_,g_);
   ArrayXd proba1X(g_),proba2X(g_);
   
   ArrayXd prop(g_);
   for(int i(0);i<g_;i++)
   prop(i)=proportion_[i];
   ArrayXXd propb(1,g_);
   propb=prop.transpose();
   
   //génération rang 1 2 3 ..m par dimension
   vector<vector<int> > yTemp(d_),x(d_);
   for(int j(0);j<d_;j++)
   {
     yTemp[j].resize(m_[j]);
     for(int i(0);i<m_[j];i++)
     yTemp[j][i]=i+1;
   }
   
   vector<double> logL(parameter_.nGibbsL-parameter_.burnL,0);
   
   //simulation de y multi dimensionnel
   for(int ind(0);ind<n_;ind++)
   {
     
     vector<vector<int> > y(d_),y2(d_),y1(d_);
     li=0;
     //algorithme de Gibbs pour simuler yi
     
     //initialisation de y et p pour Gibbs
     y=yTemp;
     for(int j(0);j<d_;j++)
     {
       random_shuffle(y[j].begin(),y[j].end(),randWrapper);//permutation de 1 2 3 ..m
       x[j]=data_[j][ind].rank;
     }
     
     y1=y;
     
     for(int k(0);k<g_;k++)
     {
       tik(ind,k)=0;
       for(int j(0);j<d_;j++)
       proba1(j,k)=probaCond(x[j],y1[j],mu_[j][k],p_[j][k]);
     }
     
     p1=(long double) (propb*proba1.colwise().prod()).sum();
     proba2=proba1;
     
     for(int iter(0);iter<parameter_.nGibbsL;iter++)
     {
       
       /*simulation des y*/
       for(int J(0);J<d_;J++)
       {
         for(int K(0);K<m_[J]-1;K++)
         {
           //"état" à tester (inversion de 2 éléments adjacents)
           y2=y;
           y2[J][K]=y[J][K+1];
           y2[J][K+1]=y[J][K];
           
           for(int k(0);k<g_;k++)//tester un stockage des proba calcul� pour �viter r�p�tition de calculs dans la boucle
           proba2(J,k)=probaCond(x[J],y2[J],mu_[J][k],p_[J][k]);
           
           
           p2=(long double) (propb*proba2.colwise().prod()).sum();
           
           alea=(long double) runif(0.,p1+p2);//unif(0,p1+p2)
           
           if(alea<p2)//accept changement
           {
             y[J]=y2[J];
             p1=p2;
             proba1.row(J)=proba2.row(J);
             y1[J]=y[J];
           }
           else//do not change y
           {
             y[J]=y1[J];//rajout J
             proba2.row(J)=proba1.row(J);
           }
         }
       }
       //y_i is updated
       
       
       /*simulation of partial rank with a gibbs sampler*/
       for(int J = 0; J < d_; J++)
       {
         if(data_[J][ind].isNotFull)//simulation of xi if it is a partial rank
         {
           x1 = x[J];
           proba1X = proba1.row(J);
           p1x = (proba1X * prop).sum();
           
           for(int kk = 0; kk < (int) (data_[J][ind].missingIndex).size()-1; kk++)//Gibbs sur les x
           {
             for(int k = 0; k < (int) (data_[J][ind].missingIndex[kk]).size()-1; k++)
             {
               //new x to test
               x2 = x[J];
               x2[data_[J][ind].missingIndex[kk][k]] = x[J][data_[J][ind].missingIndex[kk][k+1]];
               x2[data_[J][ind].missingIndex[kk][k+1]] = x[J][data_[J][ind].missingIndex[kk][k]];
               
               for(int l = 0; l < g_; l++)
               proba2X(l) = probaCond(x2,y[J],mu_[J][l],p_[J][l]);
               
               p2x = (proba2X * prop).sum();
               
               alea = (double) runif(0.,p1x+p2x);
               
               if(alea < p2)//we accept the changement
               {
                 x[J] = x2;
                 p1x = p2x;
                 proba1X = proba2X;
                 x1 = x[J];
               }
               else
               x[J] = x1;
             }
           }
           proba1.row(J)=proba1X;
         }
       }
       
       if(iter >= parameter_.burnL)
       {
         ArrayXd calculInter(g_);
         for(int cl = 0; cl < g_; cl++)
         {
           calculInter(cl) = 1;
           for(int dim = 0; dim < d_; dim++)
           calculInter(cl) *= proba1(dim,cl);
           calculInter(cl) *= propb(cl);
         }
         
         double den = calculInter.sum();
         tik.row(ind) += (calculInter/den);
         
         li += (long double) 1/den;
       }
       
     }//end gibbs sampling for sample ind
     
     //L -= log(li*div);
     L -= log(li);
     
     tik.row(ind) /= div;
     
   }//end loop on sample
   
   L += (double) n_ * log(div);
   
   output_.L = L;
   
   output_.bic = BIC(output_.L,n_,2*g_*d_+g_-1);
   
   output_.icl = output_.bic;
   
   ArrayXd entropy(n_);
   for(int i = 0; i < n_; i++)
   {
     entropy(i)=0;
     for(int j = 0; j < g_; j++)
     {
       if(tik(i,j)!=0)
       entropy(i) -= tik(i,j)*std::log(tik(i,j));
     }
     output_.icl += 2*entropy(i);
   }
   
   bic = output_.bic;
   icl = output_.icl;
 }
 
