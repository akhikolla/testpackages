#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <set>
#include <Rmath.h>

#include "test.h"

using namespace std;

//fonction pour simuler un N-échantillon d'ISR(mu,p)
void simulMixtureISR(vector<vector<int> > &simul,int const& n,int const& m,vector<vector<int> > const& mu,
		vector<double> const& p,vector<double> const& prop)
{
    vector<int> s(m),rgTemp(m);
    int l,classe(0);
    double correct,alea(0);
    bool compar,avance;

    vector<double> limite(prop.size()+1,0);
    for(int i(1); i < (int) limite.size(); i++)
    	limite[i]=limite[i-1]+prop[i-1];

    for(int i(0);i<m;i++)
        rgTemp[i]=i+1;

    for(int i(0);i<n;i++)
    {
    	//tirage aléatoire de la classe
    	alea=(double) runif(0.,1.);
        for(int j(0); j < (int) prop.size(); j++)
        {
            if((alea>limite[j]) & (alea<limite[j+1]))
            {
            	classe=j;
                break;
            }
        }

        //simulation d'un rang aléatoire: permutation du vecteur 1 2..m
        s=rgTemp;
        random_shuffle(s.begin(),s.end(),randWrapper);

        simul[i][0]=s[0];
        for (int j(1);j<m;j++)
        {
            l=0;
            avance=true;
            while (avance && l<j)
            {
                correct = (double) runif(0.,1.);
                compar=(positionRank(mu[classe],s[j])<positionRank(mu[classe],simul[i][l]));
                if ((compar && correct<p[classe])||(!compar && correct>p[classe]))
                {
                    for(int k(j-1);k>=l;k--)
                        simul[i][k+1]=simul[i][k];

                    simul[i][l]=s[j];
                    avance=false;
                }
                else
                    l++;
            }
            if (l==j)
                simul[i][l]=s[j];
        }
    }

}



double khi2(vector<vector<int> > const& data,vector<double> const& p,vector<double> const& prop,vector<vector<int> > const& mu,int const& nBoot)
{
	int const g(prop.size()),m(data[0].size()),n(data.size());
	int factM(factorial(m));
	double dkhi2(0),pvalue(0),mult((double) 2/factM),prob(0);

	//************** cacul des effectifs theoriques
	vector<double> effTheo(factM,0);
	vector<vector<int> > x(factM,vector<int> (m));

	//generation of the index of y
	vector<int> tabFact(tab_factorial(m)),listeY;
	listeY=listeSigma(m,tabFact);


	//generation of all the rank for the dim
	for(int i(0);i<factM;i++)
		x[i]=index2rank(i+1,m);

	for(int ind(0);ind<factM;ind++)
	{
		for(int cl(0);cl<g;cl++)
		{
			prob=0;
			for(vector<int>::iterator it=listeY.begin();it!=listeY.end();it++)
				prob+=probaCond(x[ind],x[*it-1],mu[cl],p[cl]);
			prob*=mult;
			prob*=prop[cl];
			effTheo[ind]+=prob;
		}
		effTheo[ind]=(double) effTheo[ind]*n;

	}

	//************** cacul des effectifs empiriques
	vector<double> effEmp(factM,0);
	vector<vector<int> > freq;
	int indexData;
	for(int i(0);i<n;i++)
	{
		indexData=rank2index(data[i],tabFact);
		effEmp[indexData-1]++;
	}

	//************** distance du khi2
	for(int i(0);i<factM;i++)
		dkhi2+=pow(effTheo[i]-effEmp[i],2);


	//************** evaluation de la loi de dKhi2 sous H0: modele = ISR
	vector<vector<int> > simulation(n,vector<int> (m,0));
	int dkhi2Sim(0);
	for(int iter(0);iter<nBoot;iter++)
	{
		vector<double> effSim(factM,0);
		simulMixtureISR(simulation,n,m,mu,p,prop);
		for(int i(0);i<n;i++)
		{
			indexData=rank2index(simulation[i],tabFact);
			effSim[indexData-1]++;
		}

		dkhi2Sim=0;
		for(int i(0);i<factM;i++)
			dkhi2Sim+=pow(effTheo[i]-effSim[i],2);

		if(dkhi2Sim>dkhi2)
			pvalue++;
	}

	pvalue/=nBoot;

	return pvalue;
}


double khi2partial(vector<Rank > &data,vector<double> const& p,vector<double> const& prop,vector<vector<int> > const& mu,int const& nBoot)
{
	int const g(prop.size()),m(data[0].rank.size()),n(data.size());
	int factM(factorial(m));
	double dkhi2(0),pvalue(0),mult((double) 2/factM),prob(0);

	//************** cacul des effectifs theoriques
	vector<double> effTheo(factM,0);
	vector<vector<int> > x(factM,vector<int> (m));

	//generation of the index of y
	vector<int> tabFact(tab_factorial(m)),listeY;
	listeY=listeSigma(m,tabFact);

	//generation of all the rank for the dim
	for(int i(0);i<factM;i++)
		x[i]=index2rank(i+1,m);


	for(int ind(0);ind<factM;ind++)
	{
		for(int cl(0);cl<g;cl++)
		{
			prob=0;
			for(vector<int>::iterator it=listeY.begin();it!=listeY.end();it++)
				prob+=probaCond(x[ind],x[*it-1],mu[cl],p[cl]);
			prob*=mult;
			prob*=prop[cl];
			effTheo[ind]+=prob;
		}
		effTheo[ind]=(double) effTheo[ind]*n;
//cout<<x[ind]<<"     "<<effTheo[ind]<<endl;
	}

	//************** cacul des effectifs empiriques
	vector<double> effEmp(factM,0);
	vector<vector<int> > freq;
	int indexData;

	//pour les rang non partiels
	//doit être fait avant car utilisé pour les rangs partiels
	for(int i(0);i<n;i++)
	{
		if(!data[i].isPartial)
		{
			indexData=rank2index(data[i].rank,tabFact);
			effEmp[indexData-1]++;
		}
	}
//cout<<effEmp<<endl;
	//pour les rangs partiels
	vector<double> effEmpPartiel(factM,0);

	for(int i(0);i<n;i++)
	{
		if(data[i].isPartial)
		{
			int nbMiss=data[i].missingIndex.size();
			int indexPerm(0);
			double somme(0);
			vector<double> effTemp(factorial(nbMiss));
			vector<int> indexTemp(factorial(nbMiss));
			do
			{
				int compteur(0);
				for(set<int>::iterator it=data[i].missingNumber.begin(); it!=data[i].missingNumber.end(); it++)
				{
					data[i].rank[data[i].missingIndex[compteur]]=*it;
					compteur++;
				}

				indexData=rank2index(data[i].rank,tabFact);
				indexTemp[indexPerm]=indexData;
				effTemp[indexPerm]=effEmp[indexData-1];
				somme+=effEmp[indexData-1];
				indexPerm++;
			}
			while (next_permutation(data[i].missingIndex.begin(),data[i].missingIndex.end()) );

			if(somme!=0)//si un des rangs possibles apparait: effemp+= freq d'apparition
			{
				for(int j(0); j< (int) indexTemp.size(); j++)
					effEmpPartiel[indexTemp[j]-1]=effTemp[j]/somme;
			}
			else//sinon effemp +=1/nb rang possible
			{	double div=(double) 1/indexTemp.size();
				for(int j(0);j < (int) indexTemp.size(); j++)
					effEmpPartiel[indexTemp[j]-1]+=div;
			}
		}

	}

	for(int i(0);i<factM;i++)
		effEmp[i] += effEmpPartiel[i];

	//************** distance du khi2
	for(int i(0);i<factM;i++)
		dkhi2 += pow(effTheo[i]-effEmp[i],2);

//cout<<effTheo<<endl;
//cout<<effEmp<<endl;

	//************** evaluation de la loi de dKhi2 sous H0: modele = ISR
	vector<vector<int> > simulation(n,vector<int> (m,0));
	double dkhi2Sim(0);

	for(int iter(0);iter<nBoot;iter++)
	{
		vector<double> effSim(factM,0);
		simulMixtureISR(simulation,n,m,mu,p,prop);
		for(int i(0); i < n; i++)
		{
		  //on réintroduit des 0 dans la simulation au même endroit
			if(data[i].isPartial)
			{
			  //nbMiss : nombre d'element manquant
				int nbMiss(data[i].missingIndex.size()),indexPerm(0);
				double somme(0);

				vector<double> effTemp(factorial(nbMiss));

				vector<int> indexTemp(factorial(nbMiss));
				vector<int> missingNumber(nbMiss);

				vector<int> index(nbMiss);
				for(int j = 0; j< nbMiss; j++)
				  index[j]=j;

				for(int j(0);j<nbMiss;j++)//on prends les numéro manquant
					missingNumber[j]=simulation[i][data[i].missingIndex[j]];
//				cout<<simulation[i]<<endl;
//				cout<<missingNumber<<endl;
				do
				{
					for(int compteur=0; compteur < nbMiss; compteur++)
						simulation[i][data[i].missingIndex[compteur]]=missingNumber[index[compteur]];
//					cout<<simulation[i]<<endl;
					indexData=rank2index(simulation[i],tabFact);
					indexTemp[indexPerm]=indexData;
					effTemp[indexPerm]=effTheo[indexData-1];
					somme+=effTheo[indexData-1];
					indexPerm++;
				}
				while (next_permutation(index.begin(),index.end()) );

				//cout<<indexTemp<<endl;


				for(int j(0);j < (int) indexTemp.size(); j++)
				  effSim[indexTemp[j]-1] += effTemp[j]/somme;



			}
			else//cas normal
			{
				indexData=rank2index(simulation[i],tabFact);
				effSim[indexData-1]++;
			}

		}
//cout<<effSim<<endl;
//cout<<effTheo<<endl;
		dkhi2Sim=0;
		for(int i(0); i < factM; i++)
			dkhi2Sim += pow(effTheo[i]-effSim[i],2);
//cout<<dkhi2Sim<<"  "<<dkhi2<<endl;
		if(dkhi2Sim>dkhi2)
			pvalue++;

	}
//cout<<endl<<dkhi2<<endl;
//cout<<pvalue<<endl;
//cout<<effTheo<<endl;
//cout<<effEmp<<endl;
	pvalue/=nBoot;

	return pvalue;
}



//---------------------------- divergence de kullback
void updateD(double &divKL,vector<int> &index, vector<vector<vector<double> > > const& p1,vector<vector<vector<double> > >  const& p2,int const& d,int const& g,
		vector<double> const& proportion1,vector<double> const& proportion2)
{
	double p1b(0),p2b(0);
	for(int k(0);k<g;k++)
	{
		int compteur(0);
		double p1Temp(1),p2Temp(1);
		for(vector<int>::iterator it=index.begin();it!=index.end();it++)
		{
			p1Temp*=p1[compteur][k][*it];
			p2Temp*=p2[compteur][k][*it];
			compteur++;
		}
		p1b+=p1Temp*proportion1[k];
		p2b+=p2Temp*proportion2[k];
	}
	divKL+=p1b*log(p1b/p2b);
}

void updateIndex(vector<int> &index,int i,vector<int> const& factm,bool &stop)
{
	if(i<0)
		stop=true;
	else
	{
		if(index[i]<factm[i]-1)
			index[i]++;
		else
		{
			updateIndex(index,i-1,factm,stop);
			index[i]=0;
		}
	}
}

void computePQ(vector<vector<vector<double> > > &p, vector<vector<vector<double> > > &q,vector<vector<vector<int> > > const& mu1,
		vector<vector<vector<int> > > const& mu2,vector<vector<double> > const& p1,vector<vector<double> > const& p2,vector<int> const& m,int d, int g)
{
	bool isDimDiff(true);
	int n=factorial(m[0]);
	double mult(1);
	vector<int> listeY;
	vector<vector<int> > x(n,vector<int> (m[0]));
	for(int dim(0);dim<d;dim++)
	{
		if(isDimDiff)
		{
			n=factorial(m[dim]);
			x=vector<vector<int> > (n,vector<int> (m[0]));
			//generation of the index of y
			vector<int> tabFact(tab_factorial(m[dim]));
			listeY=listeSigma(m[dim],tabFact);

			//generation of all the rank for the dim
			for(int i(0);i<n;i++)
				x[i]=index2rank(i+1,m[dim]);

			mult=(double) 2/factorial(m[dim]);
		}

		if(dim!=d-1)
		{
			if(m[dim]==m[dim+1])
				isDimDiff=false;
		}


		for(int cl(0);cl<g;cl++)
		{
			for(int ind(0);ind<n;ind++)
			{
				p[dim][cl][ind]=0;
				q[dim][cl][ind]=0;
				for(vector<int>::iterator it=listeY.begin();it!=listeY.end();it++)
				{
					p[dim][cl][ind]+=probaCond(x[ind],x[*it],mu1[dim][cl],p1[dim][cl]);
					q[dim][cl][ind]+=probaCond(x[ind],x[*it],mu2[dim][cl],p2[dim][cl]);
				}

				p[dim][cl][ind]*=mult;
				q[dim][cl][ind]*=mult;
			}
		}
	}
}


double divKL(vector<int> const& m,vector<vector<vector<int> > > const& mu1,vector<vector<vector<int> > > const& mu2,
		vector<vector<double> > const& p1, vector<vector<double> >  const& p2,vector<double> const& proportion1,vector<double> const& proportion2)
{
	double divKL(0);
	int const d=m.size();
	int const g=proportion1.size();
	vector<int> factm(d);
	for(int i(0);i<d;i++)
		factm[i]=factorial(m[i]);
	//dans p et q on stocke les probas pour tt les rangs possibles de ttes les dims
	vector<vector<vector<double> > > p(d,vector<vector<double> >(g));
	for(int i(0);i<d;i++)
		for(int j(0);j<g;j++)
			p[i][j].resize(factorial(m[i]));

	vector<vector<vector<double> > > q(p);

	computePQ(p,q,mu1,mu2,p1,p2,m,d,g);

	vector<int> index(d,0);
	int i(d-1);
	bool stop(false);

	while(!stop)
	{
		updateIndex(index,i,factm,stop);
		updateD(divKL,index,p,q,d,g,proportion1,proportion2);
	}

	return divKL;
}
