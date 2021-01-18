/*
 * cplqfunctionvec.hpp
 *
 *  Created on: 2 mai 2013
 *      Author: robin
 */

#ifndef CPLQFUNCTIONVEC_HPP_
#define CPLQFUNCTIONVEC_HPP_


class cpqfunctionvec {

  private:
  std::vector<cpqfunction> MycpqfunctionList_;

  public:
  // Destructor
  ~cpqfunctionvec(){
    MycpqfunctionList_.clear();
  };

  //Constructors
  cpqfunctionvec() : MycpqfunctionList_(){};
  cpqfunctionvec(int i) : MycpqfunctionList_(i){};

  //Wrapper to base functions
  std::vector<cpqfunction>::iterator begin(){return(MycpqfunctionList_.begin());};
  std::vector<cpqfunction>::iterator end(){return(MycpqfunctionList_.end());};
  std::vector<cpqfunction>::reverse_iterator rbegin(){return(MycpqfunctionList_.rbegin());};
  void vec_set( int i,cpqfunction value ) { MycpqfunctionList_.at(i) = value; };
  cpqfunction vec_get( int i) { return(MycpqfunctionList_.at(i)); };
  int size(){ return(MycpqfunctionList_.size()); };
  void push_back(cpqfunction func){MycpqfunctionList_.push_back(func);};

  // serialized push
  void SerialPush_1Breaks_Functions(Rcpp::NumericVector S0, Rcpp::NumericVector S1,Rcpp::NumericVector B1)
  {
	  int length=S1.size();
	  Rcpp::NumericVector Slopes0(1),Slopes1(1);
	  Rcpp::NumericVector BreakPoints(2);
	  for (int compteur=0; compteur<length; compteur++){
		Slopes0[0]=S0[compteur];Slopes1[0]=S1[compteur];
		BreakPoints[0]=B1[compteur];
		BreakPoints[1]=numeric_limits<double>::infinity();
		//vectorofcpqfunctions_.push_back(cpqfunction(Slopes,BreakPoints,0));
		MycpqfunctionList_.push_back(cpqfunction(Slopes0,Slopes1,BreakPoints,0.0));
	  }
  }

  void SerialPush_0Breaks_Functions(Rcpp::NumericVector S0, Rcpp::NumericVector S1)
  {
	  int length=S1.size();
	  Rcpp::NumericVector Slopes0(1),Slopes1(1);
	  Rcpp::NumericVector BreakPoints(2);
	  for (int compteur=0; compteur<length; compteur++){
		Slopes0[0]=S0[compteur];Slopes1[0]=S1[compteur];
		BreakPoints[0]=-numeric_limits<double>::infinity();
		BreakPoints[1]=numeric_limits<double>::infinity();
		//vectorofcpqfunctions_.push_back(cpqfunction(Slopes,BreakPoints,0));
		MycpqfunctionList_.push_back(cpqfunction(Slopes0,Slopes1,BreakPoints,0.0));
	  }
  }

  void SerialPenalize(Rcpp::NumericVector alpha,
		  Rcpp::NumericVector inf,Rcpp::NumericVector sup)
  {
	  int length=MycpqfunctionList_.size();
	  //assert(alpha.size()==length);
	  Rcpp::NumericVector Slopes(2);
	  Rcpp::NumericVector BreakPoints(2);
	  std::vector<cpqfunction> f;
    double zero=0;
	  cpqfunction tmp1,tmp2,tmp3;
	  for (int compteur=0; compteur<length; compteur++){
		Slopes[0]=alpha[compteur];Slopes[1]=std::numeric_limits<double>::infinity();
		BreakPoints[0]=inf[compteur];BreakPoints[1]=sup[compteur];
		tmp1=MycpqfunctionList_[compteur];
		f.push_back(cpqfunction(Slopes,Slopes,BreakPoints,zero));
		//vectorofcplfunctions_.push_back(cplfunction(Slopes,BreakPoints,0));
		vec_set(compteur,Sumq(tmp1,f[compteur]));
	  }
  };



  //Optim problem solving
  Rcpp::List OptimMargInt(NumericVector Pmoins,NumericVector Pplus,NumericVector Cmoins,NumericVector Cplus){
      //cpqfunctionvec Couts =*Coutsptr;
      int length=Pmoins.size();
      int compteur=0;
      std::vector<double> xEtoile(length);
      std::vector<cpqfunction> f;

      cpqfunction tmpfunc,tmpfunc2,tmpfunc3;
      std::vector<cpqfunction>::iterator it = MycpqfunctionList_.begin();

      //cpqfunctionvec::iterator itprec = Couts.begin();
      tmpfunc=*it;
    	tmpfunc.Squeeze(Pmoins[compteur],Pplus[compteur]);
    	 //tmpfunc.Squeeze(Cmoins[0],Cplus[0]);
    	f.push_back(tmpfunc);
      compteur++; ++it;
    	while ( it!=MycpqfunctionList_.end())
    	{
			tmpfunc=*it;
			cpqfunction tmpfunc2= *(f.rbegin());
			tmpfunc.Squeeze(Pmoins[compteur],Pplus[compteur]);
			tmpfunc.Etoile();
			tmpfunc2.Squeeze(Cmoins[compteur-1],Cplus[compteur-1]);
			tmpfunc2.Etoile();
			cpqfunction tmpfunc3 = Sumq(tmpfunc,tmpfunc2);
			tmpfunc3.Etoile();
			f.push_back(tmpfunc3);
			compteur++; ++it;
    	}
       std::vector<cpqfunction>::reverse_iterator itr,itf;
       itr = MycpqfunctionList_.rbegin();
       itf= f.rbegin();
       compteur=length-1;
       tmpfunc= *(itf);  ++itf;
    	 tmpfunc.Squeeze(Cmoins[compteur],Cplus[compteur]);
    	 //tmpfunc.Squeeze(Pmoins[length-1],Pplus[length-1]);
    	 xEtoile[compteur]=tmpfunc.Argmin();
    	 double z=xEtoile[compteur];
		while(itf!= f.rend())
		{
			--compteur;
			tmpfunc=*itr; ++itr;
			tmpfunc2=*itf; ++itf;
			tmpfunc.Squeeze(Pmoins[compteur+1],Pplus[compteur+1]);
			tmpfunc2.Squeeze(Cmoins[compteur],Cplus[compteur]);
			cpqfunction tmpfunc3=InfConfFunctq(tmpfunc,tmpfunc2,z);
			xEtoile[compteur]=tmpfunc3.Argmin();
			z=z-xEtoile[compteur];
			xEtoile[compteur]=z;
		}
    	 double tmpval,tmpval1=0;
    	 for (int i=0;i<length;i++){
    		 tmpval=xEtoile[i];
    		 xEtoile[i]=xEtoile[i]-tmpval1;
    		 tmpval1=tmpval;
    	 }
       return Rcpp::List::create(
    		Rcpp::Named("xEtoile") = Rcpp::wrap(xEtoile));
     };

  Rcpp::NumericVector OptimPriceMarket_(NumericVector Pplus, double Conso)
  {
  	int length=MycpqfunctionList_.size();
  	int compteur=0;
    std::vector<double> xEtoile(length);

  	std::vector<cpqfunction> f;
    cpqfunction tmpfunc,tmpfunc2,tmpfunc3;
    std::vector<cpqfunction>::iterator it = MycpqfunctionList_.begin();
    tmpfunc2=*it;
    tmpfunc2.Squeeze((double)0.,Pplus[compteur]);
    f.push_back(tmpfunc2);
    compteur++; ++it;
    	while ( it!=MycpqfunctionList_.end())
    	{

  			tmpfunc=*it;
  			cpqfunction tmpfunc2= *(f.rbegin());
  			tmpfunc.Squeeze((double)0,Pplus[compteur]);
  			tmpfunc.Etoile();
  			tmpfunc2.Etoile();
  			cpqfunction tmpfunc3 = Sumq(tmpfunc,tmpfunc2);
  			tmpfunc3.Etoile();
  			f.push_back(tmpfunc3);
  			compteur++; ++it;
    	}
       std::vector<cpqfunction>::reverse_iterator itr,itf;
       itr = MycpqfunctionList_.rbegin();
       itf= f.rbegin();
       compteur=length-1;
       tmpfunc= *(itf);  ++itf;
    	xEtoile[compteur]=Conso;
    	double z=xEtoile[compteur];
  		while(itf!= f.rend())
  		{
  			--compteur;
  			tmpfunc=*itr; ++itr;
  			tmpfunc2=*itf; ++itf;
  			tmpfunc.Squeeze((double)0.,Pplus[compteur+1]);
  			cpqfunction tmpfunc3=InfConfFunctq(tmpfunc,tmpfunc2,z);
  			xEtoile[compteur]=tmpfunc3.Argmin();
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

Rcpp::NumericVector OptimPriceMarket_q_(NumericVector Pricesa,NumericVector Pricesb,NumericVector Pplus, double Conso)
{
	int length=Pricesa.size();
	NumericVector::iterator Pplus_it=Pplus.begin(),Prices_itend = Pricesa.end();
	NumericVector::iterator Pricesb_it=Pricesb.begin(),Pricesa_it=Pricesa.begin();
	Rcpp::NumericVector Slopes0(1), Slopes1(1),BreakPoints(2);
	std::vector<cpqfunction> f;
	cpqfunction tmpfunc,tmpfunc3;

	Slopes1[0]=*Pricesb_it; Slopes0[0]=*Pricesa_it;
	BreakPoints[0]=0; BreakPoints[1]=*Pplus_it;
	cpqfunction tmpfunc2=cpqfunction(Slopes0,Slopes1,BreakPoints,0.0);

	f.push_back(tmpfunc2);
	++Pplus_it; ++Pricesa_it; ++Pricesb_it;
	while ( Pricesa_it!=Prices_itend)
	{
		Slopes1[0]=*Pricesb_it; Slopes0[0]=*Pricesa_it;
		BreakPoints[0]=0;  BreakPoints[1]=*Pplus_it;
		cpqfunction tmpfunc=cpqfunction(Slopes0,Slopes1,BreakPoints,0.0);
		tmpfunc.Etoile();
		tmpfunc2.Etoile();
		cpqfunction tmpfunc3 = Sumq(tmpfunc,tmpfunc2);
		tmpfunc3.Etoile();
		f.push_back(tmpfunc3);
		++Pricesa_it;++Pricesb_it; ++Pplus_it;
	}
	std::vector<cpqfunction>::reverse_iterator itr,itf,itfend;
	itfend=f.rend();
	itf= f.rbegin();
	int compteur=length-1;
	NumericVector xEtoile(length);
	xEtoile[compteur]=Conso;  ++itf;
	double z=xEtoile[compteur];
	while(itf!= itfend)
	{
		--compteur;
		Slopes1[0]=Pricesb[compteur+1]; Slopes0[0]=Pricesa[compteur+1];
		BreakPoints[0]=0;  BreakPoints[1]=Pplus[compteur+1];
		cpqfunction tmpfunc=cpqfunction(Slopes0,Slopes1,BreakPoints,0.0);
		tmpfunc2=*(itf); ++itf;
		cpqfunction tmpfunc3=InfConfFunctq(tmpfunc,tmpfunc2,z);
		xEtoile[compteur]=tmpfunc3.Argmin();
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


Rcpp::NumericMatrix OptimPriceMarket_q(NumericMatrix OffresPrixa,NumericMatrix OffresPrixb,NumericMatrix Availability,NumericVector Conso)
{
	int nbpasTemps=OffresPrixa.nrow(),nbProd=OffresPrixa.ncol();
	Rcpp::NumericMatrix Power(nbpasTemps,nbProd);

	for (int compteurt=0; compteurt<nbpasTemps; compteurt++){
		  Power(compteurt,_)=OptimPriceMarket_q_(OffresPrixa(compteurt,_),OffresPrixb(compteurt,_),Availability(compteurt,_),Conso[compteurt]);
	}

	return Power;
}




#endif /* CPLQFUNCTIONVEC_HPP_ */
