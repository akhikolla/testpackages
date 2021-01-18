/*
 * cplfunctionvec.hpp
 *
 *  Created on: 16 avr. 2013
 *      Author: robin
 */

#ifndef CPLFUNCTIONVEC_HPP_
#define CPLFUNCTIONVEC_HPP_


class cplfunctionvec {

  private:
  std::vector<cplfunction> MycplfunctionList_;

  public:
  // Destructor
  ~cplfunctionvec(){
    MycplfunctionList_.clear();
  };

  //Constructors
  cplfunctionvec() : MycplfunctionList_(){};
  cplfunctionvec(int i) : MycplfunctionList_(i){};

  //Wrapper to base functions
  std::vector<cplfunction>::iterator begin(){return(MycplfunctionList_.begin());};
  std::vector<cplfunction>::iterator end(){return(MycplfunctionList_.end());};
  std::vector<cplfunction>::reverse_iterator rbegin(){return(MycplfunctionList_.rbegin());};
  void vec_set( int i,cplfunction value ) { MycplfunctionList_.at(i) = value; };
  cplfunction vec_get( int i) { return(MycplfunctionList_.at(i)); };
  int size(){ return(MycplfunctionList_.size()); };
  
  void push_back(cplfunction func){MycplfunctionList_.push_back(func);};

  // serialized push
  void SerialPush_1Breaks_Functions(Rcpp::NumericVector S1, Rcpp::NumericVector B1)
  {
	  int length=S1.size();
	  Rcpp::NumericVector Slopes(1);
	  Rcpp::NumericVector BreakPoints(1);
	  for (int compteur=0; compteur<length; compteur++){
		Slopes[0]=S1[compteur];
		BreakPoints[0]=B1[compteur];
		//vectorofcplfunctions_.push_back(cplfunction(Slopes,BreakPoints,0));
		MycplfunctionList_.push_back(cplfunction(Slopes,BreakPoints,0.));
	  }
  }


  void SerialPush_2Breaks_Functions(Rcpp::NumericVector S1,Rcpp::NumericVector S2, Rcpp::NumericVector B1,Rcpp::NumericVector B2)
  {
	  int length=S1.size();
	  Rcpp::NumericVector Slopes(2);
	  Rcpp::NumericVector BreakPoints(2);
	  for (int compteur=0; compteur<length; compteur++){
		Slopes[0]=S1[compteur];Slopes[1]=S2[compteur];
		BreakPoints[0]=B1[compteur];BreakPoints[1]=B2[compteur];
		//vectorofcplfunctions_.push_back(cplfunction(Slopes,BreakPoints,0));
		MycplfunctionList_.push_back(cplfunction(Slopes,BreakPoints,0));
	  }
  };

  void SerialPush_nBreaks_Functions(Rcpp::NumericMatrix S, Rcpp::NumericMatrix B)
  {
	  int length=S.nrow(),nbbreak=S.ncol();
	  Rcpp::NumericVector Slopes(nbbreak);
	  Rcpp::NumericVector BreakPoints(nbbreak);
	  for (int compteur=0; compteur<length; compteur++){
		Slopes=S(compteur,_);BreakPoints=B(compteur,_);
		//vectorofcplfunctions_.push_back(cplfunction(Slopes,BreakPoints,0));
		MycplfunctionList_.push_back(cplfunction(Slopes,BreakPoints,0));
	  }
  };


  void SerialPush_Store_Functions(double gamma1, double gamma2, Rcpp::NumericVector Pflex, Rcpp::NumericVector Prod)
  {// gamma1 rendement stockage, gamma2 rendement d?stockage
	  int length=Pflex.size();

	  Rcpp::NumericVector Slopes2(3);
	  Rcpp::NumericVector BreakPoints2(3);
	  Rcpp::NumericVector Slopes1(2);
	  Rcpp::NumericVector BreakPoints1(2);
	  std::vector<cplfunction> res2;
	  cplfunction tmp1,tmp2,res;
	  double Prodtmp;
	  for (int compteur=0; compteur<length; compteur++){
		  //calcul de f1=(x*1(x>0)+gamma2 * x*1(x>0)-P/gamma1)_+
		  // + f2=(x*1(x>0)+gamma2 * x*1(x>0)-P/gamma1)_+
		  Prodtmp=Prod[compteur]+Pflex[compteur];
		if (Prod[compteur]<0)
		{// deux breakpoints
			Slopes2[0]=0;Slopes2[1]=gamma2;Slopes2[2]=1;
			BreakPoints2[0]=-std::numeric_limits<double>::infinity();
			BreakPoints2[1]=Prod[compteur]/gamma1; BreakPoints2[2]=0;
			res2.push_back(cplfunction(Slopes2,BreakPoints2,0));

		}else
		{
			Slopes1[0]=0;Slopes1[1]=1;
				BreakPoints1[0]=-std::numeric_limits<double>::infinity();
				BreakPoints1[1]=Prod[compteur]/gamma1;
				res2.push_back(cplfunction(Slopes1,BreakPoints1,0));
		}


		if (Prodtmp<0)
		{// deux breakpoints
			Slopes2[0]=0;Slopes2[1]=gamma2;Slopes2[2]=1;
			BreakPoints2[0]=-std::numeric_limits<double>::infinity();
			BreakPoints2[1]=Prodtmp/gamma1; BreakPoints2[2]=0;
			//(MycplfunctionList_.rbegin())->Sumf(cplfunction(Slopes2,BreakPoints2,0))
			MycplfunctionList_.push_back(cplfunction(Slopes2,BreakPoints2,0));
      
		}else
		{
			Slopes1[0]=0;Slopes1[1]=1;
				BreakPoints1[0]=-std::numeric_limits<double>::infinity();
				BreakPoints1[1]=Prodtmp/gamma1;
				MycplfunctionList_.push_back(cplfunction(Slopes1,BreakPoints1,0));
		}

		// res=Suml(tmp1,tmp2);
    //cplfunction tmp3=Suml(tmp2,tmp1);
    //MycplfunctionList_.push_back(tmp1);
		(MycplfunctionList_.rbegin())->Sumf(*(res2.rbegin()));
		//itend=MycplfunctionList_.rbegin();
		//itend->Sumf(tmp2);
		//vectorofcplfunctions_.push_back(cplfunction(Slopes,BreakPoints,0));
		//MycplfunctionList_.push_back(Suml(tmp1,tmp2));
	  }
  };

  void SerialPenalize(Rcpp::NumericVector alpha,
		  Rcpp::NumericVector inf,Rcpp::NumericVector sup)
  {
	  int length=MycplfunctionList_.size();
	  //assert(alpha.size()==length);
	  Rcpp::NumericVector Slopes(2);
	  Rcpp::NumericVector BreakPoints(2);
	  std::vector<cplfunction> f;
    double zero=0;
	  cplfunction tmp1,tmp2,tmp3;
	  for (int compteur=0; compteur<length; compteur++){
		Slopes[0]=alpha[compteur];Slopes[1]=std::numeric_limits<double>::infinity();
		BreakPoints[0]=inf[compteur];BreakPoints[1]=sup[compteur];
		tmp1=MycplfunctionList_[compteur];
		//f.push_back(cplfunction(Slopes,BreakPoints,zero));
		//vectorofcplfunctions_.push_back(cplfunction(Slopes,BreakPoints,0));
		vec_set(compteur,Suml(tmp1,cplfunction(Slopes,BreakPoints,zero)));
	  }
  };


  //Optim problem solving

Rcpp::List OptimMargInt(NumericVector Pmoins,NumericVector Pplus,NumericVector Cmoins,NumericVector Cplus)
{
	//cplfunctionvec Couts =*Coutsptr;
	int length=Pmoins.size();
	int compteur=0;
	std::vector<double> xEtoile(length);
	std::vector<cplfunction> f;
	cplfunction tmpfunc,tmpfunc2,tmpfunc3;
	std::vector<cplfunction>::iterator itend,it = MycplfunctionList_.begin();
	tmpfunc2=*it;
	tmpfunc2.Squeeze(Pmoins[compteur],Pplus[compteur]);
	f.push_back(tmpfunc2);
	compteur++; ++it;
	itend=MycplfunctionList_.end();
	while ( it!=itend)
	{
		tmpfunc=*it;
		cplfunction tmpfunc2= *(f.rbegin());
		if (tmpfunc.is_an_infinite_line())
		{
			tmpfunc2.Squeeze(Cmoins[compteur-1],Cplus[compteur-1]);
			tmpfunc2.EpiSum_Withline(Pmoins[compteur],Pplus[compteur],tmpfunc.FirstSlopeVal_);
			f.push_back(tmpfunc2);
		}else
		{
			tmpfunc.Squeeze(Pmoins[compteur],Pplus[compteur]);
			tmpfunc.Legendre();
			tmpfunc2.Squeeze(Cmoins[compteur-1],Cplus[compteur-1]);
			tmpfunc2.Legendre();
			if (tmpfunc.Breakpoints_.size()>tmpfunc2.Breakpoints_.size())
			{
				tmpfunc.Sumf(tmpfunc2);
				tmpfunc.Legendre();
				f.push_back(tmpfunc);
			}else
			{
				tmpfunc2.Sumf(tmpfunc);
				tmpfunc2.Legendre();
				f.push_back(tmpfunc2);
			}
		}
		compteur++; ++it;
	}

	std::vector<cplfunction>::reverse_iterator itr,itf,itfrend;
	itr = MycplfunctionList_.rbegin();
	itf= f.rbegin();
	compteur=length-1;
	tmpfunc= *(itf);  ++itf;
	tmpfunc.Squeeze(Cmoins[compteur],Cplus[compteur]);
	xEtoile[compteur]=tmpfunc.Argmin();
	double z=xEtoile[compteur];

	itfrend=f.rend();
	while(itf!= itfrend)
	{
		--compteur;
		tmpfunc=*itr; ++itr;
		tmpfunc2=*itf; ++itf;
		tmpfunc.Squeeze(Pmoins[compteur+1],Pplus[compteur+1]);
		tmpfunc2.Squeeze(Cmoins[compteur],Cplus[compteur]);
		tmpfunc2.Swap(z);
		if (tmpfunc.Breakpoints_.size()>tmpfunc2.Breakpoints_.size())
		{
			tmpfunc.Sumf(tmpfunc2);
			xEtoile[compteur]=tmpfunc.Argmin();
		}else
		{
			tmpfunc2.Sumf(tmpfunc);
			xEtoile[compteur]=tmpfunc2.Argmin();
		}
		z=z-xEtoile[compteur];
		xEtoile[compteur]=z;
	}
	double tmpval,tmpval1=0;
	for (int i=0;i<length;i++)
	{
		tmpval=xEtoile[i];
		xEtoile[i]=xEtoile[i]-tmpval1;
		tmpval1=tmpval;
	}
	return Rcpp::List::create(Rcpp::Named("xEtoile") = Rcpp::wrap(xEtoile));

	};

	Rcpp::NumericVector OptimPriceMarket_(NumericVector Pplus, double Conso)
	{
	int length=MycplfunctionList_.size();
	int compteur=0;
	std::vector<double> xEtoile(length);


	std::vector<cplfunction> f;
	cplfunction tmpfunc,tmpfunc2,tmpfunc3;
	std::vector<cplfunction>::iterator itend,it = MycplfunctionList_.begin();
	tmpfunc2=*it;
	//Rcout<<"Compteur1 "<<compteur<<endl;
	tmpfunc2.Squeeze((double)0.,Pplus[compteur]);
	f.push_back(tmpfunc2);
	compteur++; ++it;

	itend=MycplfunctionList_.end();
	while ( it!=itend)
	{
			tmpfunc=*it;
			cplfunction tmpfunc2= *(f.rbegin());
			if (tmpfunc.is_an_infinite_line())
			{
				tmpfunc2.EpiSum_Withline((double)0.,Pplus[compteur],tmpfunc.FirstSlopeVal_);
				f.push_back(tmpfunc2);
			}else
			{
				tmpfunc.Squeeze((double)0.,Pplus[compteur]);
				tmpfunc.Legendre();
				tmpfunc2.Legendre();
				if (tmpfunc.Breakpoints_.size()>tmpfunc2.Breakpoints_.size())
				{
					tmpfunc.Sumf(tmpfunc2);
					tmpfunc.Legendre();
					f.push_back(tmpfunc);
				}else
				{
					tmpfunc2.Sumf(tmpfunc);
					tmpfunc2.Legendre();
					f.push_back(tmpfunc2);
				}
			}

			compteur++; ++it;
  	}
  	std::vector<cplfunction>::reverse_iterator itr,itf,itfrend;
     itr = MycplfunctionList_.rbegin();
     itf= f.rbegin();
     compteur=length-1;
     tmpfunc= *(itf);  ++itf;
   	xEtoile[compteur]=Conso;
   	double z=xEtoile[compteur];
	itfrend=f.rend();
	while(itf!= itfrend)
	{
			--compteur;
			tmpfunc=*itr; ++itr;
			tmpfunc2=*itf; ++itf;
			tmpfunc.Squeeze((double)0.,Pplus[compteur+1]);
			tmpfunc2.Swap(z);
			if (tmpfunc.Breakpoints_.size()>tmpfunc2.Breakpoints_.size())
			{
				tmpfunc.Sumf(tmpfunc2);
				xEtoile[compteur]=tmpfunc.Argmin();
			}else
			{
				tmpfunc2.Sumf(tmpfunc);
				xEtoile[compteur]=tmpfunc2.Argmin();
			}
			z=z-xEtoile[compteur];
			xEtoile[compteur]=z;
		}

  	 double tmpval,tmpval1=0;
  	 for (int i=0;i<length;i++){
  		 tmpval=xEtoile[i];
  		 xEtoile[i]=xEtoile[i]-tmpval1;
  		 tmpval1=tmpval;
  	 }
     return(Rcpp::wrap(xEtoile));
	};

};
Rcpp::List OptimPriceStorage_(NumericVector Prices,NumericVector Pmoins,NumericVector Pplus,NumericVector Cmoins,NumericVector Cplus)
{
	int length=Pmoins.size();
	NumericVector::iterator Pmoins_it=Pmoins.begin(), Cmoins_it=Cmoins.begin();
	NumericVector::iterator Pplus_it=Pplus.begin(),Cplus_it=Cplus.begin();
	NumericVector::iterator Prices_it=Prices.begin(),Prices_itend = Prices.end();
	std::vector<cplfunction> f;
	cplfunction tmpfunc,tmpfunc3;

	cplfunction tmpfunc2=cplfunction(*Pmoins_it,*Pplus_it,0,*Prices_it,numeric_limits<double>::infinity());
	f.push_back(tmpfunc2);
	++Pplus_it;++Pmoins_it; ++Prices_it;

	while ( Prices_it!=Prices_itend)
	{
		//cplfunction tmpfunc2= *(f.rbegin());
		tmpfunc2.Squeeze(*Cmoins_it,*Cplus_it);
		tmpfunc2.EpiSum_Withline(*Pmoins_it,*Pplus_it,*Prices_it);
		f.push_back(tmpfunc2);
		++Prices_it;++Cmoins_it;++Cplus_it;++Pmoins_it; ++Pplus_it;
	}

	std::vector<cplfunction>::reverse_iterator itr,itf,itfend;
	itfend=f.rend();
	itf= f.rbegin();
	int compteur=length-1;
	itf->Squeeze(Cmoins[compteur],Cplus[compteur]);
	//tmpfunc.Squeeze(Pmoins[length-1],Pplus[length-1]);
	NumericVector xEtoile(length);
	xEtoile[compteur]=itf->Argmin();  ++itf;
	double z=xEtoile[compteur];
	while(itf!= itfend)
	{
		--compteur;
		cplfunction tmpfunc=cplfunction(Pmoins[compteur+1],Pplus[compteur+1],0,Prices[compteur+1],numeric_limits<double>::infinity());
		itf->Squeeze(Cmoins[compteur],Cplus[compteur]);
		itf->Swap(z);
		itf->Sumf(tmpfunc);
		xEtoile[compteur]=itf->Argmin();
		++itf;
		z=z-xEtoile[compteur];
		xEtoile[compteur]=z;
	}
	double tmpval,tmpval1=0;
	for (int i=0;i<length;i++)
	{
		tmpval=xEtoile[i];
		xEtoile[i]=xEtoile[i]-tmpval1;
		tmpval1=tmpval;
	}
	return Rcpp::List::create(
			Rcpp::Named("xEtoile") = Rcpp::wrap(xEtoile));
	}





double evalf_(NumericVector BreakPoints, NumericVector Prices,double x)
{
	if (abs(x)<=50)
	{
		return(x*Prices[1]);
	}else if (x>0)
	{
		return 50*Prices[1]+(x-50)*Prices[2];
	}else
	{
		return -50*Prices[1]+(x+50)*Prices[0];
	}
}


Rcpp::NumericVector SerialOptimPriceStorage(NumericMatrix Prices,NumericMatrix BreakPoints,NumericVector Pmoins,NumericVector Pplus,NumericVector Cmoins,NumericVector Cplus)
{
	  cplfunctionvec f(Prices.nrow());
	  for (int compteur=0; compteur<Prices.nrow(); compteur++){
		//vectorofcplfunctions_.push_back(cplfunction(Slopes,BreakPoints,0));
		  f.vec_set(compteur,cplfunction(Prices(compteur,_),BreakPoints(compteur,_),0));
	  }

	  int ncases=Pmoins.size();
	  NumericVector benefit(ncases);

	  for (int compteur=0; compteur<ncases; compteur++){
		//vectorofcplfunctions_.push_back(cplfunction(Slopes,BreakPoints,0));
		NumericVector Pmoinstmp(Prices.nrow(),Pmoins[compteur]);
		NumericVector Cmoinstmp(Prices.nrow(),Cmoins[compteur]);
		NumericVector PPlustmp(Prices.nrow(),Pplus[compteur]);
		NumericVector Cplustmp(Prices.nrow(),Cplus[compteur]);
		NumericVector res=(f.OptimMargInt(Pmoinstmp,PPlustmp,Cmoinstmp,Cplustmp))["xEtoile"];
		benefit[compteur]=0;
		 for (int compteur2=0; compteur2<Prices.nrow(); compteur2++){
				//vectorofcplfunctions_.push_back(cplfunction(Slopes,BreakPoints,0));
			 benefit[compteur]=benefit[compteur]+evalf_(BreakPoints(compteur2,_),Prices(compteur2,_),res[compteur2]);
		}
	  }
	return benefit;
}





Rcpp::NumericVector OptimPriceMarket_l_(NumericVector Prices,NumericVector Pplus, double Conso)
{
int length=Prices.size();
NumericVector::iterator Pplus_it=Pplus.begin(),Prices_itend = Prices.end(),Prices_it=Prices.begin();
std::vector<cplfunction> f;
cplfunction tmpfunc,tmpfunc3;

cplfunction tmpfunc2=cplfunction((double)0.,*Pplus_it,0,*Prices_it,numeric_limits<double>::infinity());
f.push_back(tmpfunc2);
++Pplus_it; ++Prices_it;
while ( Prices_it!=Prices_itend)
{
	tmpfunc2.EpiSum_Withline((double)0.,*Pplus_it,*Prices_it);
	f.push_back(tmpfunc2);
	++Prices_it; ++Pplus_it;
}
std::vector<cplfunction>::reverse_iterator itr,itf,itfend;
itfend=f.rend();
itf= f.rbegin();
int compteur=length-1;
NumericVector xEtoile(length);
xEtoile[compteur]=Conso;  ++itf;
double z=xEtoile[compteur];
while(itf!= itfend)
{
	--compteur;
	cplfunction tmpfunc=cplfunction((double)0.,Pplus[compteur+1],0,Prices[compteur+1],numeric_limits<double>::infinity());
	itf->Swap(z);
	itf->Sumf(tmpfunc);
	xEtoile[compteur]=itf->Argmin();
	++itf;
	z=z-xEtoile[compteur];
	xEtoile[compteur]=z;
}
double tmpval,tmpval1=0;
for (int i=0;i<length;i++)
{
	tmpval=xEtoile[i];
	xEtoile[i]=xEtoile[i]-tmpval1;
	tmpval1=tmpval;
}
 return(Rcpp::wrap(xEtoile));
};


Rcpp::NumericMatrix OptimPriceMarket_l(NumericMatrix OffresPrix,NumericMatrix Availability,NumericVector Conso)
{
	int nbpasTemps=OffresPrix.nrow(),nbProd=OffresPrix.ncol();
	Rcpp::NumericMatrix Power(nbpasTemps,nbProd);

	for (int compteurt=0; compteurt<nbpasTemps; compteurt++){
		  Power(compteurt,_)=OptimPriceMarket_l_(OffresPrix(compteurt,_),Availability(compteurt,_),Conso[compteurt]);
	}

	return Power;
}



#endif /* CPLFUNCTIONVEC_HPP_ */
