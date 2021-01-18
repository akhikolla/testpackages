/*
 * cpqfunction.hpp
 *
 *  Created on: 16 avr. 2013
 *      Author: robin
 */

#ifndef CPQFUNCTION_HPP_
#define CPQFUNCTION_HPP_

class  cpqfunction {
    /* private fields */

    public:
    map<double,pair<double,double> > Breakpoints_; // breakpoints, where we have a polynom here polynom is 1/2 ax^2+bx+c a=.first b=.second
    double FirstBreakVal_; // firstbreakval

    ~cpqfunction(){
      Breakpoints_.clear();
    };

    cpqfunction() : Breakpoints_(), FirstBreakVal_(0){}

    cpqfunction(int NbCoefficients, pair<double,double> * Coefficients, double * BreakPoints,double FirstBreakVal) {
		  create_cpqfunction(NbCoefficients,Coefficients,BreakPoints,FirstBreakVal);
	  };

    void create_cpqfunction(int NbCoefficients, pair<double,double> * Coefficients, double * BreakPoints,double FirstBreakVal) {
  	  for (int i=0; i<NbCoefficients; i++){
  	    Breakpoints_[BreakPoints[i]]=Coefficients[i];
  	  }
  	  Breakpoints_[BreakPoints[NbCoefficients]]=pair<double,double>(numeric_limits<double>::infinity(),numeric_limits<double>::infinity());
  	  FirstBreakVal_= FirstBreakVal;
    };

    cpqfunction(pair<double,double> * Coefficients,double val){
  	  // This function constructs a simple quadratic function, no breaks.
  	  int NbCoefficients=1;
  	  double BreakPoints[1]={-numeric_limits<double>::infinity()};
      create_cpqfunction(NbCoefficients,Coefficients,BreakPoints,val);
    }

    cpqfunction(Rcpp::NumericVector Slopes0,Rcpp::NumericVector Slopes1, Rcpp::NumericVector BreakPoints,double FirstBreakVal){
  	int NbSlopes=  Slopes1.size();
		if (NbSlopes+1==BreakPoints.size()){
				if (isincreasing(BreakPoints)){
					for (int i=0; i<NbSlopes; i++){
						pair<double,double> coeffs=Slopes2Coeffs(Slopes0[i],Slopes1[i]);
						  if (
								  Slopes0[i]<=Slopes1[i] &&
								  	  (i==0 ||
								  	   getSlope(coeffs,BreakPoints[i])>= getSlope(Slopes2Coeffs(Slopes0[i-1],Slopes1[i-1]),BreakPoints[i])-1e-7)
							 )
						  {
							Breakpoints_[BreakPoints[i]]=Slopes2Coeffs(Slopes0[i],Slopes1[i]);
						  }else{
								Rcout<<"getSlope(coeffs,BreakPoints[i])"<<getSlope(coeffs,BreakPoints[i])<<endl;
								Rcout<<"getSlope(Slopes2Coeffs(Slopes0[i-1],Slopes1[i-1]),BreakPoints[i])"<<getSlope(Slopes2Coeffs(Slopes0[i-1],Slopes1[i-1]),BreakPoints[i])<<endl;

							Rprintf( "Error: non increasing Slopes" ) ;
									throw nonincreasingslopes() ;
						  }
					  }
					  Breakpoints_[BreakPoints[NbSlopes]]=pair<double,double>(numeric_limits<double>::infinity(),numeric_limits<double>::infinity());
					  FirstBreakVal_= FirstBreakVal;
				}else{
					Rprintf( "Error: non increasing breakpoints" ) ;
					throw nonincreasingbreakpoints() ;
				}
		}else{
			Rprintf( "Error: number of Slopes must be number of breaks -1 " ) ;
			throw nonincreasingslopes() ;
		}
	}

    cpqfunction(cpqfunction const & x) : Breakpoints_(x.Breakpoints_), FirstBreakVal_(x.FirstBreakVal_) {
	  }

    cpqfunction* clone() const {
        return new cpqfunction(*this) ;
    }

    cpqfunction(double uniquebreak,pair<double,double> * Coefficients,double val){
	   // This function constructs a simple quadratic function bounded from below breaks.
	   int NbCoefficients=1;
	   double BreakPoints [1]={uniquebreak};
	   create_cpqfunction(NbCoefficients,Coefficients,BreakPoints,val);
   };

    cpqfunction(double * twobreaks,pair<double,double> Coefficient, double val){
	    // Simple quadratic function bounded from above and below
	    int NbCoefficients=1;
	    pair<double,double> Coefficients [2];
	    Coefficients[0]=Coefficient; Coefficients[1]=pair<double,double>(numeric_limits<double>::infinity(),numeric_limits<double>::infinity());
	    create_cpqfunction(NbCoefficients,Coefficients,twobreaks,val);
    };
//    cpqfunction(Rcpp::NumericVector TwoBreaks,Rcpp::NumericVector Slopes)
//    {
//    	int NbCoefficients=1;
//    	pair<double,double> Coefficients [2];
//    	Coefficients[0]=Slopes2Coeffs(Slopes[0],Slopes[1]);
//    	Coefficients[1]=pair<double,double>(numeric_limits<double>::infinity(),numeric_limits<double>::infinity());
//
//    	create_cpqfunction(NbCoefficients,Coefficients,twobreaks,val);
//    }
    //stack exceptions, implement std::exception
    class emptyfunc : public std::exception {
     public:
      const char * what() { return "empty function"; }
    };

    class nonincreasingslopes : public std::exception {
     public:
      const char * what() { return "non increasing slopes"; }
    };

    class nonincreasingbreakpoints : public std::exception {
     public:
      const char * what() { return "non increasing breakpoints"; }
    };

    Rcpp::List get_BreakPoints(){
      std::vector<double> Breakpoints;
    	std::vector<double> Slopes1;
      std::vector<double> Slopes0;
	 	 	map<double,pair<double,double> >::iterator it=Breakpoints_.begin();
	 	 	int nbSlopes=0,compteur=0;
	 	 	while(it != Breakpoints_.end()) {it++; nbSlopes++;}
	 	 	nbSlopes--;
  			it=Breakpoints_.begin();
  			compteur=0;
  			while(it != Breakpoints_.end()) {
  				Breakpoints.push_back( it->first );
  				if (compteur != (nbSlopes+1)){
  					Slopes0.push_back( it->second.second );
            Slopes1.push_back( it->second.first+it->second.second);
  				}
  				it++; compteur++;
  			}

  			return Rcpp::List::create(
				Rcpp::Named("Breakpoints") = Rcpp::wrap(Breakpoints),
				Rcpp::Named("Slopes0") = Rcpp::wrap(Slopes0),
        Rcpp::Named("Slopes1") = Rcpp::wrap(Slopes1));
  	}


    void print()
    {
		std::vector<double> Breakpoints;
		std::vector<double> Slopes1;
		std::vector<double> Slopes0;
		map<double,pair<double,double> >::iterator it=Breakpoints_.begin();
		map<double,pair<double,double> >::iterator itplus1=Breakpoints_.begin(); ++itplus1;
		while(itplus1 != Breakpoints_.end())
		{
			Rcout<<"|"<<it->first<<"|";
			Rcout<<"__("<< getSlope(it->second,0);
			Rcout<<","<< getSlope(it->second,1)<<")__";
			++it; ++itplus1;
		}
		//last point
		Rcout<<"|"<<it->first<<"|";
		Rcout<<"__("<<it->second.first;
		Rcout<<","<< it->second.second<<")__";
		Rcout<<endl;

		it=Breakpoints_.begin();
		itplus1=Breakpoints_.begin(); ++itplus1;

		while(itplus1 != Breakpoints_.end())
		{
			Rcout<<"|"<<it->first<<"|";
			Rcout<<"__("<< it->second.first;
			Rcout<<","<< it->second.second<<")__";
			++it; ++itplus1;
		}
		//last point
		Rcout<<"|"<<it->first<<"|";
		Rcout<<"__("<<it->second.first;
		Rcout<<","<< it->second.second<<")__";
		Rcout<<endl;

  	}


    cpqfunction & operator = (cpqfunction const & s) {
    /* Cleanup current data */
    if(this != &s) {
     Breakpoints_.clear();
     /* copy needed data, call copy constructor
      * not efficient but will call copy constructor
      * */
     Breakpoints_=s.Breakpoints_;
     FirstBreakVal_=s.FirstBreakVal_;
    }

    return *this;
   }

    double evalf(double x)
    {
    	map<double, pair<double,double> >::iterator it = Breakpoints_.begin();
    	map<double, pair<double,double> >::reverse_iterator rev_it = Breakpoints_.rbegin();
    	if (it->first>x || rev_it->first<x)
    	{
    		return numeric_limits<double>::infinity();
    	}
    	else
    	{
    		  std::map<double, pair<double,double> >::const_reverse_iterator
    		  last_element_not_greater_than(Breakpoints_.upper_bound(x));

    		  if (Breakpoints_.rend() == last_element_not_greater_than) {
    		    return -1;
    		  }
    		double diff=FirstBreakVal_-getVal(Breakpoints_.begin()->second,Breakpoints_.begin()->first);
    		return getVal(last_element_not_greater_than->second,x)+diff;
    	}
    };
    void AddSimple(double const & breakpoint,pair<double,double> const & left,pair<double,double> const & right,double const & val){
     map<double, pair<double,double> >::iterator i = Breakpoints_.begin();

	   if ((left.first==right.first)&&(left.second==right.second))
	   {//in this case the break is not a real break...
		   FirstBreakVal_=FirstBreakVal_+getVal(left,i->first)-getVal(left,breakpoint)+val;
		   while(i != Breakpoints_.end()) {
			   (*i).second.first=i->second.first+left.first;
			   (*i).second.second=i->second.second+left.second;
		   	   ++i;
		   }

	   }else{
           if (breakpoint<=(*Breakpoints_.begin()).first){
               //BreakPoint is out of the domain, on the left
        	   FirstBreakVal_=FirstBreakVal_+getVal(right,i->first)-getVal(right,breakpoint)+val;
    		   while(i != Breakpoints_.end()) {
    			   (*i).second.first=i->second.first+right.first;
    			   (*i).second.second=i->second.second+right.second;
    		   	   ++i;
    		   }
           }
           else{
        	   FirstBreakVal_=FirstBreakVal_+getVal(left,i->first)-getVal(left,breakpoint)+val;
               if (breakpoint>=(*Breakpoints_.rbegin()).first){
          		   while(i != Breakpoints_.end()) {
        			   (*i).second.first=i->second.first+left.first;
        			   (*i).second.second=i->second.second+left.second;
          			   ++i;
          		   }
               }else{/*here the new breakpoint is inside the domain of this and
               the rightslope and left Coefficients are different*/
        		   map<double, pair<double,double> >::iterator it,ittmp;
        		   unsigned int initialsize=Breakpoints_.size();
        		   //insert the new breakpoint

          		   it=Breakpoints_.insert(pair<double, pair<double,double> > (breakpoint, pair<double,double>(0.0,0.0))).first;
          		   it--; ittmp=it; it++;
          		   if (Breakpoints_.size()!=initialsize){
          			(*it).second = (*ittmp).second;
          		   }

          		   map<double, pair<double,double> >::iterator i = Breakpoints_.begin();
          		   while(i != it) {
        			   (*i).second.first=i->second.first+left.first;
        			   (*i).second.second=i->second.second+left.second;
          			   ++i;
          		   }
          		   while(i != Breakpoints_.end()) {
        			   (*i).second.first=i->second.first+right.first;
        			   (*i).second.second=i->second.second+right.second;
          			   ++i;
          		   }
	           }
	   	   }
	   }

   };

	bool push_quad_if_not_linear(map<double,pair<double,double> >::iterator it)
	{// returns true if was not linear (and pushed)
		// the polynom is 1/2 ax^2 +bx+c
		bool res=false;
		double a=(it->second).first;
		if (a!=0)
		{
			res=true;
			Breakpoints_[getSlope(it->second,it->first)]=pair<double,double>(1/a,-(it->second).second/a);
		}
		return(res);
	}

	bool push_linear_if_diff(map<double,pair<double,double> >::iterator it,
			map<double,pair<double,double> >::iterator itplus1)
	{// return true if dif was true and linear was pushed
		bool res=false;
		double x0=itplus1->first;
		double slopeplus1=getSlope(itplus1->second,x0);
		double slope=getSlope(it->second,x0);
		if ((slopeplus1-slope)!=0)
		{
			res=true;
			Breakpoints_[slope]=pair<double,double>(0.,x0);
		}
		return (res);
	}


	bool is_a_point()
	{
		return ((Breakpoints_.size()==1)&&
				(Breakpoints_.begin()->second).first==numeric_limits<double>::infinity());
	}

	bool is_an_infinite_line()
	{
		return ((Breakpoints_.size()==1)&&
				(Breakpoints_.begin()->second).first==0);
	}

    void Etoile()
    {
    	//Rcout << __FUNCTION__ << endl;
    	//this->print();
		cpqfunction tmp(*this);
		Breakpoints_.clear();
		bool done=false;
		bool is_last_linear=false;
		if (tmp.is_a_point())
		{// a point get transformed into a line
			Breakpoints_[-numeric_limits<double>::infinity()]=pair<double,double>(0.,Breakpoints_.begin()->first);
			done=true;
		}

		if (tmp.is_an_infinite_line())
		{// an infinite line get transformed into a point
			Breakpoints_[((tmp.Breakpoints_.begin())->second).second]=pair<double,double>(0.,numeric_limits<double>::infinity());
			done=true;
		}

		if (done)
		{
			// do nothing
		}
		else
		{
			map<double,pair<double,double> >::iterator it=tmp.Breakpoints_.begin();
			map<double,pair<double,double> >::iterator itplus1=tmp.Breakpoints_.begin(); ++itplus1;
			bool is_last_infinite=false;

			if ((it->first!=-numeric_limits<double>::infinity()))
			{// first break is not infinity : this will give a slope and a first break with inf
				Breakpoints_[-numeric_limits<double>::infinity()]=pair<double,double>(0.,it->first);
				is_last_linear=true;
			}else
			{// nothing here for now

			}

			while (itplus1!=tmp.Breakpoints_.end())
			{
				is_last_linear=!push_quad_if_not_linear(it) && is_last_linear;
				if (itplus1->first!=numeric_limits<double>::infinity())
				{
					is_last_linear=push_linear_if_diff(it,itplus1);
				}
				++it; ++itplus1;
			}

			if (is_last_linear)
			{
				Breakpoints_[getSlope(it->second,0)]=pair<double,double>(numeric_limits<double>::infinity(),numeric_limits<double>::infinity());
			}else
			{
				Breakpoints_[numeric_limits<double>::infinity()]=pair<double,double>(numeric_limits<double>::infinity(),numeric_limits<double>::infinity());
			}
		}



		//treating the first non null breakval
		map<double,pair<double,double> >::iterator it=tmp.Breakpoints_.begin();
	    map<double,pair<double,double> >::iterator itobj=Breakpoints_.begin();
		if (itobj->first==-numeric_limits<double>::infinity()){itobj++;}
		if (it->first==-numeric_limits<double>::infinity()){it++;}

		pair<double,double> Coef=pair<double,double>(itobj->second.first,itobj->second.second);

		double c_init=tmp.FirstBreakVal_-getVal(it->second,it->first);
		FirstBreakVal_= getVal(Coef,itobj->first)-c_init+(it->second.second)*(it->second.second)/(2*it->second.first);


		//Rcout << __FUNCTION__<< " out" <<endl;
		//this->print();

    }


    bool eq(cpqfunction  const & cpqfunction1){
	   if (FirstBreakVal_!=cpqfunction1.FirstBreakVal_)
	   {
		   return(false);
	   }
	   if (Breakpoints_.size()!=cpqfunction1.Breakpoints_.size())
	   {
		   return(false);
	   }else{
		   map<double, pair<double,double> > mybreak=Breakpoints_;
  		   map<double, pair<double,double> >::iterator it = Breakpoints_.begin();
  				map<double, pair<double,double> >::const_iterator    it2=cpqfunction1.Breakpoints_.begin();
  		   while(it != Breakpoints_.end())
  		   {
  			   if (	it->first==it2->first&&
  					   it->second.first==it2->second.first&&
  					   it->second.second==it2->second.second)
  			   {
  				 ++it;++it2;
  			   }else
  			   {
  				   return(false);
  			   }
  		   }
  		   return(true);
	   }
   }

    double Argmin(){
	  // Rcout << __FUNCTION__ << endl;
	   //this->print();
	   double res;
	   cpqfunction tmp(*this);
	   int NbCoefficients=tmp.Breakpoints_.size()-1;
	   if (NbCoefficients<2){
		   if (NbCoefficients==1){
		       if (getSlope(tmp.Breakpoints_.rbegin()->second,tmp.Breakpoints_.rbegin()->first)<=0){
		    	  res =tmp.Breakpoints_.rbegin()->first;
		       }else{
		    	   if (getSlope(tmp.Breakpoints_.begin()->second,tmp.Breakpoints_.begin()->first)>0){
		    		   res =tmp.Breakpoints_.begin()->first;
		    	   }else{// here f'(0)<0 and f'(1)>0, in particular f'(0)!=f'(1).
		    		   res=getXetoile(tmp.Breakpoints_.begin()->second);
		    	   }
		       }
		   }else{
			   if (NbCoefficients==0){
				   res =tmp.Breakpoints_.begin()->first;
			   }else{
				   Rcout<<"NbCoefficients="<<NbCoefficients<<endl;
				   throw emptyfunc();
			   }
		   }
	   }else{
	       if (getSlope(tmp.Breakpoints_.begin()->second,tmp.Breakpoints_.begin()->first)>0){
	    	   res =tmp.Breakpoints_.begin()->first;
	       }else{
      		map<double, pair<double,double> >::iterator i = tmp.Breakpoints_.begin();
			if (i->second.first==i->second.second){
				res=i->first;
			}else{
				res=getXetoile(i->second);
			}

      		++i;
      		   while(i != tmp.Breakpoints_.end()) {
      			   if (res>i->first){
      				   res=i->first;
      			   }
         			if (getSlope(i->second,i->first)>0){
         				break;
         			}else{
						if (i->second.first==i->second.second){
							res=i->first;
						}else{
							res=getXetoile(i->second);
						}
         			}
      			++i;
      		   }
	       }
	   }
	  // Rcout<<"res="<<res<<endl;
	   return(res);
   }

    void Squeeze(double leftBreak,double rightBreak){
  	  //Rcout << __FUNCTION__ << "("<<leftBreak<<","<<rightBreak<<")"<<endl;
	  //this->print();
	   cpqfunction tmp(*this);

	   if (tmp.Breakpoints_.size()<1 ||leftBreak>=rightBreak ||tmp.Breakpoints_.begin()->first>=rightBreak ||tmp.Breakpoints_.rbegin()->first<=leftBreak){
		   if (tmp.Breakpoints_.begin()->first==rightBreak||tmp.Breakpoints_.begin()->first==-numeric_limits<double>::infinity()){
			   Breakpoints_.clear();
			   Breakpoints_[tmp.Breakpoints_.begin()->first]=pair<double,double>(numeric_limits<double>::infinity(),numeric_limits<double>::infinity());
		   }else{
			   if (tmp.Breakpoints_.rbegin()->first==leftBreak||tmp.Breakpoints_.rbegin()->first==numeric_limits<double>::infinity()){
				   Breakpoints_.clear();
				   Breakpoints_[tmp.Breakpoints_.rbegin()->first]=pair<double,double>(numeric_limits<double>::infinity(),numeric_limits<double>::infinity());
			   }else{
				   Rcout<<"tmp.Breakpoints_.size()"<<tmp.Breakpoints_.size()<<endl;
				   Rcout<<"rightBreak"<<rightBreak<<endl;
				   Rcout<<"leftBreak"<<leftBreak<<endl;
				   Rcout<<"tmp.Breakpoints_.begin()->first"<<tmp.Breakpoints_.begin()->first<<endl;
				   Rcout<<"tmp.Breakpoints_.rbegin()->first"<<tmp.Breakpoints_.rbegin()->first<<endl;
				   Rcout<<"empty function as a result of Squeeze"<<endl;
				   throw emptyfunc();
			   }
		   }
	   }else{
		   if (tmp.Breakpoints_.size()==1){
			   Breakpoints_.clear();
			   Breakpoints_[tmp.Breakpoints_.rbegin()->first]=pair<double,double>(numeric_limits<double>::infinity(),numeric_limits<double>::infinity());
		   }else{
	   		   map<double, pair<double,double> >::iterator itleft,itright,itb;
			   unsigned int initialsize=tmp.Breakpoints_.size();

			   //insert the new breakpoint
			   if (tmp.Breakpoints_.begin()->first<leftBreak){

				pair<map<double, pair<double,double> >::iterator,bool> breakinsertion=Breakpoints_.insert(pair<double, pair<double,double> > (leftBreak, pair<double,double>(0.0,1.0)));
				itleft=breakinsertion.first;
			   if (breakinsertion.second){
				   --itleft; pair<double,double> u=itleft->second; ++itleft;
				  //<<(*ittmp).second<<"left B"<<leftBreak<<endl;
				   (*itleft).second = u;
			   }
			   itb=Breakpoints_.begin();
			   Breakpoints_.erase(itb,itleft);
			   }

			   if (tmp.Breakpoints_.rbegin()->first>rightBreak){
			   initialsize=Breakpoints_.size();
			   itright=Breakpoints_.insert(pair<double, pair<double,double> > (rightBreak, pair<double,double>(0.0,1.0))).first;
			   itright++;
			   itb=Breakpoints_.end();
			   if (itright!=itb) Breakpoints_.erase(itright,itb);
			   map<double, pair<double,double> >::reverse_iterator irev=Breakpoints_.rbegin();
			   irev->second=pair<double,double>(numeric_limits<double>::infinity(),numeric_limits<double>::infinity());
			   }
		   }
	   }
	 //  Rcout<<"out : "<<endl;
	  // this->print();
   }

    void Sumf(cpqfunction  &  cpqfunction1){
  	   //Rcout << __FUNCTION__ <<endl;
  	   //this->print();
  	   //cpqfunction1.print();
  	  cpqfunction tmp(*this),tmp1=cpqfunction1;

  	  (*this).Squeeze(tmp1.Breakpoints_.begin()->first , tmp1.Breakpoints_.rbegin()->first);

      if (tmp1.Breakpoints_.size()<=2){
        if (tmp1.Breakpoints_.size()==1){
  			  if (tmp1.Breakpoints_.begin()->first!=Breakpoints_.begin()->first){
  				  Rcout<<"in Sumf"<<endl;
  				  throw emptyfunc();
  			  }
        }else{
          //Rcout << tmp1.Breakpoints_.rbegin()->first;
        	(*this).AddSimple(tmp1.Breakpoints_.begin()->first,
        					tmp1.Breakpoints_.begin()->second,
        					tmp1.Breakpoints_.begin()->second,
        					tmp1.FirstBreakVal_);
  			/*  map<double,pair<double,double> >::iterator it=Breakpoints_.begin();
  			FirstBreakVal_=FirstBreakVal_+tmp1.FirstBreakVal_;
  			  double a,b;
  			  while (it != Breakpoints_.end()){
  				  a=it->second.first; b=it->second.second;
  				  (*it).second.first=a+tmp1.Breakpoints_.begin()->second.first;
  				  (*it).second.second=b+tmp1.Breakpoints_.begin()->second.second;
  				  ++it;
  			  }*/
  		  }
      }else{
  		  double a,b;
  		  map<double,pair<double,double> >::iterator it=tmp1.Breakpoints_.begin();
  		  ++it;
  		  map<double, pair<double,double> >::iterator itplus=it,itplus2;
  		  pair<double,double> ab,zero=pair<double,double>(0.0,0.0);
  		  ++it;itplus2=it;
  		  it=tmp1.Breakpoints_.begin();
  		  //double const & breakpoint,pair<double,double> const & left,pair<double,double> const & right,double const & val
  	    (*this).AddSimple(itplus->first,it->second,itplus->second,tmp1.evalf(itplus->first));
  	    ++itplus;++it;++itplus2;
  	    while (itplus2!=tmp1.Breakpoints_.end()){
  	      a=itplus->second.first-it->second.first;
  	      b=itplus->second.second-it->second.second;
  	      ab=pair<double,double>(a,b);
  	      (*this).AddSimple(itplus->first,zero,ab,0.0);
          ++itplus;++it;++itplus2;
  	    }
      }
  	  //  Rcout<<"out :"<<endl;
  	   // this->print();
    }

    void Swap(double y) {
     //Rcout << __FUNCTION__ << " " << y << endl;
	 //this->print();
	   if(Breakpoints_.size() < 1)
		   throw emptyfunc();
	   map<double,pair<double,double> >::reverse_iterator rit;
	   cpqfunction tmp(*this);
	   Breakpoints_.clear();
	   rit = tmp.Breakpoints_.rbegin();
	   double last_first = rit->first;
	   ++rit;
	   while(rit != tmp.Breakpoints_.rend()){
		   double a=(rit->second.first);
		   double b=(rit->second.second);
		   Breakpoints_[y-last_first] = pair<double,double>(a,-a*y-b);
		   last_first = rit->first;
		   ++rit;
	   }
	   Breakpoints_[y-last_first] = pair<double,double>(numeric_limits<double>::infinity(),numeric_limits<double>::infinity());
	   //Rcout<<"out:"<<endl;
	   //this->print();
	   //return(*this);
   }

};// end of class cpqfunction definition



#endif /* CPQFUNCTION_HPP_ */
