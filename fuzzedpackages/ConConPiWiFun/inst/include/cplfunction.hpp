/*
 * cplfunction.hpp
 *
 *  Created on: 27 avr. 2013
 *      Author: robin
 */

#ifndef cplfunction_HPP_
#define cplfunction_HPP_


/*
 * cplfunction.hpp
 *
 *  Created on: 16 avr. 2013
 *      Author: robin
 */



class cplfunction {
	// this class implements the convex continuous piecewise functions with a map
	// this allows nlog(n) sum of two such function
	// FirstSlopeVal_ is the absolute slope associated to the first breakpoint
	// Breakpoints_ is
	//
	//

    public:
    std::map<double,double> Breakpoints_; // breakpoints
    double FirstBreakVal_; // firstbreakval
    double FirstSlopeVal_ ;

    ~cplfunction(){
      Breakpoints_.clear();
    };

    cplfunction()
    	: Breakpoints_(),FirstBreakVal_(0),
    	  FirstSlopeVal_(-std::numeric_limits<double>::infinity()){};

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

    cplfunction(Rcpp::NumericVector Slopes, Rcpp::NumericVector BreakPoints,double FirstBreakVal){
  		int NbSlopes=  Slopes.size();
  		if (NbSlopes==BreakPoints.size()){
  			if (isincreasing(Slopes)){
  				if (isincreasing(BreakPoints)){
  					FirstSlopeVal_=Slopes[0];
  					Breakpoints_[BreakPoints[0]]=0;
  					for (int i=1; i<NbSlopes; i++)
  					{
  						Breakpoints_[BreakPoints[i]]=Slopes[i]-Slopes[i-1];
  					}
  					FirstBreakVal_= FirstBreakVal;
  				}else{
  					Rprintf( "Error: non increasing breakpoints" ) ;
  					throw nonincreasingbreakpoints() ;
  				}
  			}else{
  				Rprintf( "Error: non increasing Slopes" ) ;
  				throw nonincreasingslopes() ;
  			}
  		}else{
  			Rprintf( "Error: number of Slopes must be number of breaks+1  " ) ;
  			throw nonincreasingslopes() ;
  		}
	  };

	  cplfunction(cplfunction const & x) :
		  Breakpoints_(x.Breakpoints_),FirstBreakVal_(x.FirstBreakVal_),
		  FirstSlopeVal_(x.FirstSlopeVal_){};

    cplfunction* clone() const {
        return new cplfunction(*this) ;
    };

    cplfunction(double uniquebreak,double val)//a single point has 1 break and infinite FirstSlopeVal_
    	:Breakpoints_(),FirstBreakVal_(val),FirstSlopeVal_(std::numeric_limits<double>::infinity())
    {
    	Breakpoints_[uniquebreak]=0;
    };

    //a half line is just 1 point and a slope
    cplfunction(double uniquebreak,double val,double Slope1)
    	:Breakpoints_(),
    	 FirstBreakVal_(val),FirstSlopeVal_(Slope1)
    {
    	Breakpoints_[uniquebreak]=0;
    };


    //a "V" function is 2 points and 2 slopes the first point being infinity
    cplfunction(double uniquebreak,double val,double Slope1, double Slope2)
		:Breakpoints_(),
		 FirstBreakVal_(val),FirstSlopeVal_(Slope1)
    {
    	Breakpoints_[-std::numeric_limits<double>::infinity()]=0;
    	Breakpoints_[uniquebreak]=Slope2;
    };
    cplfunction(double breakleft,double breakright,double val,double Slope1,double Slope2)
    	:Breakpoints_(),
    	 FirstBreakVal_(val),FirstSlopeVal_(Slope1)
    {
    	Breakpoints_[breakleft]=0;
    	Breakpoints_[breakright]=Slope2-Slope1;
    };
    
    
    Rcpp::List get_BreakPoints()
    {
    	int nbBreaks=Breakpoints_.size();
    	std::vector<double> Breakpoints(nbBreaks);
  		std::vector<double> Slopes(nbBreaks);
  		int compteur=0;
  		std::map<double,double>::iterator Breakpoints_it=Breakpoints_.begin();
  		while(Breakpoints_it != Breakpoints_.end())
  		{
  			Breakpoints[compteur]=Breakpoints_it->first;
  			if (compteur==0)
  			{
  				Slopes[compteur]=FirstSlopeVal_;
  			}else
  			{
				Slopes[compteur]=Slopes[compteur-1]+Breakpoints_it->second;
  			}
  			Breakpoints_it++; compteur++;
  		}

		return Rcpp::List::create(
			Rcpp::Named("Breakpoints") = Breakpoints,
			Rcpp::Named("Slopes") = Slopes);
  	}

    void print()
     {
     	int nbBreaks=Breakpoints_.size();
     	std::vector<double> Breakpoints(nbBreaks);
   		std::vector<double> Slopes(nbBreaks);
   		int compteur=0;
   		std::map<double,double>::iterator Breakpoints_it=Breakpoints_.begin();
   		if (Breakpoints_it->second!=0)
   		{
   			Rcout<<"Warning first Slope diff non null =  "<< Breakpoints_it->second <<", ";
   		}

   		while(Breakpoints_it != Breakpoints_.end())
   		{
   			Rcout<<"|"<<Breakpoints_it->first<<"|";
   			if (compteur==0)
   			{
   				Slopes[compteur]=FirstSlopeVal_;
   				Rcout<<"__"<<FirstSlopeVal_<<"__";
   			}else
   			{
   				Slopes[compteur]=Slopes[compteur-1]+Breakpoints_it->second;
   				Rcout<<"__"<<Slopes[compteur]<<"__";
   			}
   			Breakpoints_it++; compteur++;
   		}
   		Rcpp::Rcout<<std::endl;
   	}



/*
    cplfunction(double * twobreaks,double slope, double val){
 	   int NbSlopes=1;
 	   double Slopes [2]={slope, numeric_limits<double>::infinity()};
 	   create_cplfunction(NbSlopes,Slopes,twobreaks,val);
    }
*/
    /* cplfunction(simplefunction sfunc){
 	   int NbSlopes=2;
 	   double Slopes [2]={sfunc.leftslope_, sfunc.rightslope_};
 	   double BreakPoints [3]={-numeric_limits<double>::infinity(),sfunc.breakpoint_,numeric_limits<double>::infinity()};
 	   create_cplfunction(NbSlopes,Slopes,BreakPoints,sfunc.val_);
    };*/

    cplfunction & operator = (cplfunction & s) {
     /* Cleanup current data */
     if(this != &s) {
      Breakpoints_.clear();
      /* copy needed data, call copy constructor
       * not efficient but will call copy constructor
       * */
      Breakpoints_=s.Breakpoints_;
      FirstBreakVal_=s.FirstBreakVal_;
      FirstSlopeVal_=s.FirstSlopeVal_;
     }
     return *this;
    }

    void AddSimple(double leftslope, double rightslope, double val, double breakpoint)
    {
    	//Rcout << __FUNCTION__ << "("<<leftslope<< ","<<rightslope<<","<<breakpoint<<")"<<"in "<<endl;
    	//this->print();
    	std::map<double, double>::iterator i = Breakpoints_.begin();
    	FirstBreakVal_=FirstBreakVal_+val;

		if (breakpoint<=Breakpoints_.begin()->first)
		{//BreakPoint is out of the domain, on the left
			  FirstSlopeVal_=FirstSlopeVal_+rightslope;
		}else
		{
		  if (breakpoint>=Breakpoints_.rbegin()->first && Breakpoints_.rbegin()->second== std::numeric_limits<double>::infinity()){
			  FirstSlopeVal_=FirstSlopeVal_+leftslope;
		  }else
		  {
			/*here the new breakpoint is inside the domain of this*/

			FirstSlopeVal_=FirstSlopeVal_+leftslope;
			double diff=rightslope-leftslope;
			std::pair<std::map<double, double>::iterator,bool> tmp_insert=Breakpoints_.insert(pair<double, double> (breakpoint, diff));
			if (!tmp_insert.second)
			{//insert the new breakpoint if it does not exist and if it exists increment :
				double tmpval=tmp_insert.first->second;
				(*tmp_insert.first).second=tmpval+rightslope-leftslope;
			}
		  }
		}
		//Rcout << __FUNCTION__ << "out "<<endl;
		//this->print();
    }

    bool eq(cplfunction  const & cplfunction1){
 	   if (FirstBreakVal_!=cplfunction1.FirstBreakVal_){
 		   return(false);
 	   }
 	   if (FirstSlopeVal_!=cplfunction1.FirstSlopeVal_){
 		  return(false);
 	   }
 	   if (Breakpoints_.size()!=cplfunction1.Breakpoints_.size()){
 		   return(false);
 	   }else{
   		   std::map<double, double>::iterator i = Breakpoints_.begin();
				std::map<double, double>::const_iterator	   i2=cplfunction1.Breakpoints_.begin();
   		   while(i != Breakpoints_.end()) {
   			   if (i->first==i2->first&&i->second==i2->second){
   				 ++i;++i2;
   			   }else{
   				   return(false);
   			   }
   		   }
   		   return(true);
 	   }
    }


//    double evalf(double x)
//	{
//		map<double, double >::iterator it = Breakpoints_.begin();
//		map<double, double >::reverse_iterator rev_it = Breakpoints_.rbegin();
//		double res,precBreak=it->first;
//		if (it->first!=numeric_limits<double>::infinity())
//		{
//			res=FirstBreakVal_;
//		}else{
//			res=-numeric_limits<double>::infinity();
//		}
//
//		double precSlope=FirstSlopeVal_;
//		if (is_an_infinite_line() || is_a_point()||Breakpoints_.size()<2)
//		{
//			return(FirstBreakVal_+FirstSlopeVal_*x);
//		}
//		else if (it->first==numeric_limits<double>::infinity() && it->first )
//		{
//
//		}
//		else if (it->first>x || (rev_it->first<x && !is_last_infinity()) )
//		{
//			return numeric_limits<double>::infinity();
//		}
//		else
//		{
//			  while (it!=Breakpoints_.end() && it->first<x)
//			  {
//				if (res!=-numeric_limits<double>::infinity())
//				{
//				  res=res+(it->second-precBreak)*precSlope;
//				  precSlope=precSlope+it->second;
//				  precBreak=it->first;
//				  ++it;
//				}else
//				{
//
//				  res=res+(it->second-precBreak)*precSlope;
//				  precSlope=precSlope+it->second;
//				  precBreak=it->first;
//				  ++it;
//				}
//
//			  }
//
//			  if (it==Breakpoints_.end())
//			  {
//				  return res;
//			  }else
//			  {
//				  return res;
//			  }
//
//		}
//	};

    bool is_last_infinity()
    {
    	if (Breakpoints_.size()==1)
    	{
    		return(FirstSlopeVal_!=std::numeric_limits<double>::infinity());
    	}else
    	{
    		return(((Breakpoints_.rbegin()->second!=std::numeric_limits<double>::infinity()) &&
    				FirstSlopeVal_!=std::numeric_limits<double>::infinity()));
    	}
    }

    bool is_a_point()
    {
    	return(FirstSlopeVal_==std::numeric_limits<double>::infinity());
    }

    bool is_an_infinite_line()
    {
    	return(FirstSlopeVal_!=std::numeric_limits<double>::infinity() &&
    			Breakpoints_.size()==1 &&
    			(Breakpoints_.begin())->first==-std::numeric_limits<double>::infinity());
    }

    // Etoile is replaced by "Legendre"
    void Etoile()
    {
    //	Rcout << __FUNCTION__<< " in" <<endl;
    //	this->print();
		cplfunction tmp(*this);
		Breakpoints_.clear();
		bool done=false;

		if (tmp.is_a_point())
		{// a point get transformed into a line
			FirstBreakVal_=-tmp.FirstBreakVal_;
			FirstSlopeVal_=(tmp.Breakpoints_.begin())->first;
			Breakpoints_[-std::numeric_limits<double>::infinity()]=0;
			done=true;
		}

		if (tmp.is_an_infinite_line())
		{// an infinite line get transformed into a point
			FirstBreakVal_=-tmp.FirstBreakVal_;
			Breakpoints_[(tmp.Breakpoints_.begin())->second+FirstSlopeVal_]=0;
			FirstSlopeVal_=std::numeric_limits<double>::infinity();
			done=true;
		}
		if (tmp.Breakpoints_.size()==1 && !done)
		{// only one breakpoint not a line not a point : this is a half line
		 // gives a half line with 2 breaks
			FirstBreakVal_=-tmp.FirstBreakVal_;
			FirstSlopeVal_=(tmp.Breakpoints_.begin())->first;
			Breakpoints_[-std::numeric_limits<double>::infinity()]=0;
			Breakpoints_[tmp.FirstSlopeVal_]=std::numeric_limits<double>::infinity();
			done=true;
		}
		if (done)
		{
			// do nothing
		}
		else
		{
			std::map<double,double>::iterator it=tmp.Breakpoints_.begin();
			std::map<double,double>::iterator itplus1=tmp.Breakpoints_.begin(); ++itplus1;
			FirstBreakVal_=-tmp.FirstBreakVal_;
			double NewBreak,NewSlope,PastBreakPrec,NewBreakPrec;
			if ((it->first!=-std::numeric_limits<double>::infinity()))
			{// first break is not infinity : this will give a slope and a first break with inf
				FirstSlopeVal_=it->first;
				Breakpoints_[-std::numeric_limits<double>::infinity()]=0;
			}else
			{
				FirstSlopeVal_=itplus1->first;
			}
			PastBreakPrec=FirstSlopeVal_;
			NewBreakPrec=tmp.FirstSlopeVal_;

			while (itplus1!=tmp.Breakpoints_.end())
			{
				NewSlope=itplus1->first-PastBreakPrec;
				PastBreakPrec=itplus1->first;
				NewBreak=it->second+NewBreakPrec;
				NewBreakPrec=NewBreak;
				Breakpoints_[NewBreak]=NewSlope;
				++it; ++itplus1;
			}
			if (tmp.is_last_infinity())
			{
				NewBreak=it->second+NewBreakPrec;
				Breakpoints_[NewBreak]=std::numeric_limits<double>::infinity();
			}
		}
	//	Rcout << __FUNCTION__<< " out" <<endl;
	//	this->print();
    };

    void EpiSum_Withline(double lowerbound,double upperbound,double Slope)
    {
    	//Rcout << __FUNCTION__<< " in" <<endl;
    	//this->print();
    	bool done=false;
		if (is_a_point())
		{// a point get transformed into a line
			double uniquebreak=Breakpoints_.begin()->first;
			Breakpoints_.clear();
			FirstSlopeVal_=Slope;
			Breakpoints_[lowerbound+uniquebreak]=0;
			Breakpoints_[upperbound+uniquebreak]=std::numeric_limits<double>::infinity();
			done=true;
		}

		if (is_an_infinite_line())
		{// an infinite line get transformed into a point a point remains the same and back it is the line
			// ... no change
			done=true;
		}
		if (Breakpoints_.size()==1 && !done)
		{// only one breakpoint not a line not a point : this is a half line
		 // gives a half line with 2 breaks
			double uniquebreak=Breakpoints_.begin()->first;
			double uniqueslope=Breakpoints_.begin()->second+FirstSlopeVal_;
			Breakpoints_.clear();
			if (Slope<uniqueslope)
			{// f*+g* has two breaks and two slopes
				FirstSlopeVal_=Slope;
				Breakpoints_[lowerbound+uniquebreak]=0;
				Breakpoints_[upperbound+uniquebreak]=uniqueslope-Slope;
			}else
			{
				FirstSlopeVal_=uniqueslope;
				Breakpoints_[lowerbound+uniquebreak]=0;
			}

			done=true;
		}
		if (done)
		{
			// do nothing
		}
		else
		{// there are 2 breaks or more.

	    	std::map<double,double>::iterator it=Breakpoints_.begin();
	    	std::map<double,double>::iterator itplus1=Breakpoints_.begin();++itplus1;
	    	std::map<double,double>::iterator itend=Breakpoints_.end();
	    	double CurrentSlope=FirstSlopeVal_;
	    	double newbreakval,newbreakslope,x;
	    	bool newbreak;
	    	if (CurrentSlope<Slope)
	    	{
				while (itplus1!=itend && (CurrentSlope+itplus1->second)<Slope )
				{// for the breakpoint that have slope < Slope just shift their key "on the left" by "lowerbound"
					CurrentSlope=CurrentSlope+itplus1->second;
					x=it->first;
					const_cast<double&>(it->first) = x+lowerbound;
					++it;++itplus1;
				}

		    	if (itplus1==itend)
		    	{
		    		newbreak=false;
		    	}else
		    	{//CurrentSlope<Slope <=CurrentSlope+itplus1->second
					x=it->first;
					const_cast<double&>(it->first) = x+lowerbound;
					if ((CurrentSlope+itplus1->second)==Slope)
					{
						newbreak=false;
					}else
					{
						newbreak=true;
						newbreakval=itplus1->first+lowerbound;
						newbreakslope=Slope-CurrentSlope;
						if (itplus1->second!=std::numeric_limits<double>::infinity())
						{
							itplus1->second=itplus1->second-newbreakslope;
						}

					}

					while (itplus1!=itend)
					{
						x=itplus1->first;
						const_cast<double&>(itplus1->first) = x+upperbound;
						++itplus1;
					}
		    	}
	    	}else
	    	{// Slope were lower or equal to the first Slope
	    		//Slope <=CurrentSlope
	    		if (it->first==-std::numeric_limits<double>::infinity())
	    		{
	    			newbreak=false;
					while (it!=itend)
					{
						x=it->first;
						const_cast<double&>(it->first) = x+upperbound;
						++it;
					}
	    		}else
	    		{//it->first>-numeric_limits<double>::infinity()
	    			if (CurrentSlope==Slope)
	    			{
	    				newbreak=false;
						x=it->first;
						const_cast<double&>(it->first) = x+lowerbound;
	    				while (itplus1!=itend)
	    				{
	    					x=itplus1->first;
	    					const_cast<double&>(itplus1->first) = x+upperbound;
	    					++itplus1;
	    				}
	    			}else
	    			{//CurrentSlope>Slope and it->first>-numeric_limits<double>::infinity()
						newbreak=true;
						newbreakval=it->first+lowerbound;
						it->second=it->second+FirstSlopeVal_-Slope; // from zero to true val
						FirstSlopeVal_=Slope;
						newbreakslope=0;

						while (it!=itend)
						{
							x=it->first;
							const_cast<double&>(it->first) = x+upperbound;
							++it;
						}
	    			}
	    		}
	    	}


	    	if (newbreak)
	    	{
	    		Breakpoints_[newbreakval]=newbreakslope;
	    	}

		}
		//Rcout << __FUNCTION__<< " out" <<endl;
		//this->print();
    }

    double Argmin(){
 	 // Rcout << __FUNCTION__ << endl;
 	 // this->print();
 	   double res;
 	   double precslope=FirstSlopeVal_;
 	   int NbSlopes=Breakpoints_.size();

 	   if (is_a_point()){res=Breakpoints_.begin()->first;}
 	   if (is_an_infinite_line()){
 		   if (FirstSlopeVal_==0){res=0;}
 		   else if (FirstSlopeVal_<0){res=std::numeric_limits<double>::infinity();}
 		   else if (FirstSlopeVal_>0){res=-std::numeric_limits<double>::infinity();}
 	   }

	   if (precslope>0){
		   res =Breakpoints_.begin()->first;
	   }else{
		 std::map<double, double>::iterator i = Breakpoints_.begin();
		 std::map<double, double>::iterator iend = Breakpoints_.end();
		++i;
		  while(i != iend) {
				res=i->first;
				precslope=precslope+i->second;
				if (precslope>0){ break;}
				++i;
				if(precslope<0 && is_last_infinity())
				{
					res=std::numeric_limits<double>::infinity();
				}
		  }
	   }
	//  Rcout<<"res="<<res<<endl;
	   return(res);

 	 //  return(res);
    };

    void Squeeze(double leftBreak,double rightBreak)
    {
     	  // Rcout << __FUNCTION__ << "("<<leftBreak<< ","<<rightBreak<<")"<<" in "<<endl;
     	//   Rcout<< "left : "<<leftBreak<<", right : "<<rightBreak<<endl;
     	 //  Rcout<< "this left "<< Breakpoints_.begin()->first << "this right : "<<Breakpoints_.rbegin()->first<<endl;
     	  //Rcout<<  "this right slope : "<<Breakpoints_.rbegin()->second<<endl;

     	  // this->print();
     	  // Rcout<<"FirstSlopeVal_ : "<<FirstSlopeVal_<<endl;
    	// test for empty interval or empty intersection of function and interval
    	bool done=false;
		if (  (leftBreak>rightBreak) ||
				(Breakpoints_.rbegin()->first<leftBreak && !is_last_infinity() ) ||
				(Breakpoints_.begin()->first>rightBreak && !is_last_infinity()) )
		{
			double arg_ = 1e-7;
			if ((leftBreak-Breakpoints_.rbegin()->first<arg_))
			{
				double tmpbreak=Breakpoints_.rbegin()->first;
				Breakpoints_.clear();
				Breakpoints_[tmpbreak]=std::numeric_limits<double>::infinity();
				done=true;
			}else if (Breakpoints_.begin()->first-rightBreak<arg_){
				double tmpbreak=Breakpoints_.begin()->first;
				Breakpoints_.clear();
				Breakpoints_[tmpbreak]=std::numeric_limits<double>::infinity();
				done=true;
			}else
			{
				if (leftBreak>=rightBreak){Rcout<<"leftBreak>=rightBreak"<<std::endl;}
				Rcpp::Rcout<<"Empty function thrown in Squeeze"<<std::endl;
				Rcpp::Rcout << __FUNCTION__ << "("<<leftBreak<< ","<<rightBreak<<")"<<" in "<<std::endl;
				this->print();
				Rcpp::Rcout<< "Breakpoints_.begin()->first-rightBreak"<<Breakpoints_.begin()->first-rightBreak<<std::endl;
				throw emptyfunc();
			}

		}
		else if (is_a_point()||done)
		{

		}else
		{
			/// taking care of left break
			std::map<double, double>::iterator it=Breakpoints_.begin();
			std::map<double, double>::iterator itend=Breakpoints_.end();
			if (it->first<leftBreak)
			{// something will have to be cut on the left
				while (it!=itend && it->first<leftBreak )
				{
					FirstSlopeVal_=FirstSlopeVal_+it->second;
					++it;
				}
				if (it==itend)
				{
					Breakpoints_.erase(Breakpoints_.begin(),it);
					Breakpoints_.insert(pair<double, double> (leftBreak, 0.0));
				}else
				{
					if (it!=Breakpoints_.begin()){
						Breakpoints_.erase(Breakpoints_.begin(),it);
					}

					if (Breakpoints_.begin()->first==leftBreak){
						FirstSlopeVal_=FirstSlopeVal_+Breakpoints_.begin()->second;
						Breakpoints_.begin()->second=0.0;
					}else{
						Breakpoints_.insert(std::pair<double, double> (leftBreak, 0.0));
					}
				}
			}

			// taking care of right break
			if (is_last_infinity())
			{
				if (rightBreak!=std::numeric_limits<double>::infinity())
				{// something has to be cut on the right
					std::pair<std::map<double, double>::iterator,bool> tmp_insert=Breakpoints_.insert(std::pair<double, double> (rightBreak, numeric_limits<double>::infinity()));
					std::map<double, double>::iterator tmp_insert_it=tmp_insert.first;
					if (!tmp_insert.second)
					{
						(tmp_insert_it)->second=std::numeric_limits<double>::infinity();
					}

					++tmp_insert_it;
					if (tmp_insert_it!=Breakpoints_.end())
					{
						Breakpoints_.erase(tmp_insert_it,Breakpoints_.end());
					}
					if (Breakpoints_.size()==1)
					{
						FirstSlopeVal_=std::numeric_limits<double>::infinity();
					}
				}
			}else
			{
				if (rightBreak!=std::numeric_limits<double>::infinity() && Breakpoints_.rbegin()->first>rightBreak )
				{// something has to be cut on the right
					pair<map<double, double>::iterator,bool> tmp_insert=Breakpoints_.insert(pair<double, double> (rightBreak, numeric_limits<double>::infinity()));

					map<double, double>::iterator tmp_insert_it=tmp_insert.first;
					if (!tmp_insert.second)
					{
						(tmp_insert_it)->second=std::numeric_limits<double>::infinity();
					}

					++tmp_insert_it;
					if (tmp_insert_it!=Breakpoints_.end())
					{
						Breakpoints_.erase(tmp_insert_it,Breakpoints_.end());
					}
					if (Breakpoints_.size()==1)
					{
						Breakpoints_.begin()->second=0;
						FirstSlopeVal_=std::numeric_limits<double>::infinity();
					}
				}
			}
		}

		//Rcout << __FUNCTION__ << "("<<leftBreak<< ","<<rightBreak<<")"<<" out "<<endl;
		//this->print();
		//Rcout << "FirstSlopeVal_: "<<FirstSlopeVal_<<endl;
     };


    std::vector<double> Squeeze2(double leftBreak,double rightBreak)
        {
         	   //Rcout << __FUNCTION__ << "("<<leftBreak<< ","<<rightBreak<<")"<<" in "<<endl;
         	  // Rcout<< "left : "<<leftBreak<<", right : "<<rightBreak<<endl;
         	  // Rcout<< "this left "<< Breakpoints_.begin()->first << "this right : "<<Breakpoints_.rbegin()->first<<endl;
         	  //Rcout<<  "this right slope : "<<Breakpoints_.rbegin()->second<<endl;

         	   //this->print();
         	   //Rcout<<"FirstSlopeVal_ : "<<FirstSlopeVal_<<endl;
        	// test for empty interval or empty intersection of function and interval
    	std::vector<double> res(2);
    	res[0]=0; res[1]=0;
    		if (  (leftBreak>rightBreak) ||
    				(Breakpoints_.rbegin()->first<leftBreak && !is_last_infinity() ) ||
    				(Breakpoints_.begin()->first>rightBreak && !is_last_infinity()) )
    		{
    			if (leftBreak>=rightBreak){Rcout<<"leftBreak>=rightBreak"<<std::endl;}
    			Rcout<<"Empty function thrown in Squeeze"<<std::endl;
    			Rcout << __FUNCTION__ << "("<<leftBreak<< ","<<rightBreak<<")"<<" in "<<std::endl;
    			this->print();
    			throw emptyfunc();
    		}
    		else if (is_a_point())
    		{

    		}else
    		{
    			/// taking care of left break
    			map<double, double>::iterator it=Breakpoints_.begin();
    			map<double, double>::iterator itend=Breakpoints_.end();
    			if (it->first<leftBreak)
    			{// something will have to be cut on the left
    				while (it!=itend && it->first<leftBreak )
    				{
    					FirstSlopeVal_=FirstSlopeVal_+it->second;
    					++it;
    				}
    				if (it==itend)
    				{
    					Breakpoints_.erase(Breakpoints_.begin(),it);
    					Breakpoints_.insert(pair<double, double> (leftBreak, 0.0));
    				}else
    				{
    					if (it!=Breakpoints_.begin()){
    						Breakpoints_.erase(Breakpoints_.begin(),it);
    					}

    					if (Breakpoints_.begin()->first==leftBreak){
    						FirstSlopeVal_=FirstSlopeVal_+Breakpoints_.begin()->second;
    						Breakpoints_.begin()->second=0.0;
    					}else{
    						Breakpoints_.insert(pair<double, double> (leftBreak, 0.0));
    					}
    				}
    			}

    			// taking care of right break
    			if (is_last_infinity())
    			{
    				if (rightBreak!=numeric_limits<double>::infinity())
    				{// something has to be cut on the right
    					pair<map<double, double>::iterator,bool> tmp_insert=Breakpoints_.insert(pair<double, double> (rightBreak, numeric_limits<double>::infinity()));
    					map<double, double>::iterator tmp_insert_it=tmp_insert.first;
    					if (!tmp_insert.second)
    					{
    						(tmp_insert_it)->second=numeric_limits<double>::infinity();
    					}

    					++tmp_insert_it;
    					if (tmp_insert_it!=Breakpoints_.end())
    					{
    						Breakpoints_.erase(tmp_insert_it,Breakpoints_.end());
    					}
    					if (Breakpoints_.size()==1)
    					{
    						FirstSlopeVal_=numeric_limits<double>::infinity();
    					}
    				}
    			}else
    			{
    				if (rightBreak!=numeric_limits<double>::infinity() && Breakpoints_.rbegin()->first>rightBreak )
    				{// something has to be cut on the right
    					pair<map<double, double>::iterator,bool> tmp_insert=Breakpoints_.insert(pair<double, double> (rightBreak, numeric_limits<double>::infinity()));

    					map<double, double>::iterator tmp_insert_it=tmp_insert.first;
    					if (!tmp_insert.second)
    					{
    						(tmp_insert_it)->second=numeric_limits<double>::infinity();
    					}

    					++tmp_insert_it;
    					if (tmp_insert_it!=Breakpoints_.end())
    					{
    						Breakpoints_.erase(tmp_insert_it,Breakpoints_.end());
    					}
    					if (Breakpoints_.size()==1)
    					{
    						FirstSlopeVal_=numeric_limits<double>::infinity();
    					}
    				}
    			}
    		}

    		double tmp=Breakpoints_.begin()->first-leftBreak;
    		if (tmp>0)
    		{
    			res[0]=tmp;
    		}
    		tmp=rightBreak-Breakpoints_.rbegin()->first;
    		if (tmp>0)
    		{
    			res[1]=tmp;
    		}
    		return(res);
			//Rcout << __FUNCTION__ << "("<<leftBreak<< ","<<rightBreak<<")"<<" out "<<endl;
    		//this->print();
    		//Rcout << "FirstSlopeVal_: "<<FirstSlopeVal_<<endl;
         };


    void Sumf(cplfunction & cplfunction1){
   // Rcout << __FUNCTION__ <<" in "<<endl;
   //	this->print();
   //	cplfunction1.print();
   	//  Rcout<<"cplfunction1.FirstSlopeVal_ : "<<cplfunction1.FirstSlopeVal_<<endl;
   	  if (cplfunction1.is_last_infinity())
   	  {
   		(*this).Squeeze(cplfunction1.Breakpoints_.begin()->first,numeric_limits<double>::infinity());
   	  }else
   	  {
  	   	(*this).Squeeze(cplfunction1.Breakpoints_.begin()->first,cplfunction1.Breakpoints_.rbegin()->first);
   	  }

   	  if (	(cplfunction1.Breakpoints_.size()==1) ||
   			  (cplfunction1.Breakpoints_.size()==2 && !cplfunction1.is_last_infinity()))
   	  {// linear function ... only one slope
   		  FirstSlopeVal_=FirstSlopeVal_+cplfunction1.FirstSlopeVal_;
   		  FirstBreakVal_=FirstBreakVal_+cplfunction1.FirstBreakVal_;
   	  }else if (is_a_point())
   	  {

   	  }else
   	  {// at least two slopes with 2 breaks or more
		map<double,double>::const_iterator it=cplfunction1.Breakpoints_.begin();
		++it;
		map<double, double>::const_iterator itplus=it;
		it=cplfunction1.Breakpoints_.begin();

		(*this).AddSimple(cplfunction1.FirstSlopeVal_,
				itplus->second+cplfunction1.FirstSlopeVal_,
				cplfunction1.FirstBreakVal_,itplus->first);

		++itplus;++it;
		while (itplus!=cplfunction1.Breakpoints_.end())
		{// enter this loop if there are more than 2 slopes (3 breaks or more)...
			(*this).AddSimple(0.0,itplus->second,0.0,itplus->first);
			++itplus;++it;
		}
      }

   	//Rcout << __FUNCTION__ <<" out "<<endl;
   	//this->print();
    };

    void Swap(double y)
    {
		//Rcout << __FUNCTION__ << " " << y << " in "<< endl;
		//Rcout << "y-0.6" << " " << y-0.6 << " in "<< endl;
		//this->print();
 	   cplfunction tmp(*this);
 	   Breakpoints_.clear();
 	   map<double,double>::reverse_iterator rit = tmp.Breakpoints_.rbegin();

 	   if (tmp.is_a_point())
 	   {
 		  Breakpoints_[y-rit->first]=0;
 		  FirstSlopeVal_=numeric_limits<double>::infinity();
 	   }else
 	   {
 	 	   map<double,double>::reverse_iterator ritplus1 = tmp.Breakpoints_.rbegin();
 	 	   ++ritplus1;
 			double LastSlopeVal=0.;
 			if (tmp.is_last_infinity())
 			{
 				Breakpoints_[-numeric_limits<double>::infinity()]=0;
 			}else
 			{
 				Breakpoints_[y-rit->first]=0; ++rit;++ritplus1;
 			}

 			while(ritplus1 != tmp.Breakpoints_.rend()){
				Breakpoints_[y-rit->first] = rit->second;
				LastSlopeVal=LastSlopeVal+(rit->second);
				++rit;++ritplus1;
 			}
 			if (rit->first!=-numeric_limits<double>::infinity())
 			{
 				Breakpoints_[y-rit->first] = numeric_limits<double>::infinity();
 			}

 			FirstSlopeVal_=-(FirstSlopeVal_+LastSlopeVal);
 	   }
 	 // Rcout << __FUNCTION__ << " " << y << " out "<< endl;
 	  // this->print();
 	   //return(*this);
    };

    void Legendre()
    {
        	//Rcout << __FUNCTION__<< " in" <<endl;
        	//this->print();
    	FirstBreakVal_=-FirstBreakVal_;
		if (is_a_point())
		{// a point get transformed into a line
			FirstSlopeVal_=(Breakpoints_.begin())->first;
			Breakpoints_.clear();
			Breakpoints_[-numeric_limits<double>::infinity()]=0;
		} else if (is_an_infinite_line())
		{
			double breakval=(Breakpoints_.begin())->second+FirstSlopeVal_;
			Breakpoints_.clear();
			Breakpoints_[breakval]=0;
			FirstSlopeVal_=numeric_limits<double>::infinity();
		} else if (Breakpoints_.size()==1)
		{// only one breakpoint not a line not a point : this is a half line
		 // gives a half line with 2 breaks
			double breakval=FirstSlopeVal_;
			FirstSlopeVal_=(Breakpoints_.begin())->first;
			Breakpoints_.clear();
			Breakpoints_[-numeric_limits<double>::infinity()]=0;
			Breakpoints_[breakval]=numeric_limits<double>::infinity();
		}else{// 2 breaks or more
			double breakval;
	    	map<double,double>::iterator it=Breakpoints_.begin();
	    	if (it->first==-numeric_limits<double>::infinity())
	    	{
	    		breakval=FirstSlopeVal_+it->second;
	    		++it;
	    		FirstSlopeVal_=FirstSlopeVal_+it->second;
	    		it->second=0;
	    		Breakpoints_.erase(Breakpoints_.begin());
	    	}else
	    	{
	    		breakval=-numeric_limits<double>::infinity();
	    	}

		    if (is_last_infinity())
		    {
	    		double lastbreak=flip_push_left(breakval);
	    		Breakpoints_[lastbreak]=numeric_limits<double>::infinity();
	    	}
	    	else
	    	{
	    		flip_push_left(breakval);
	    	}
		}

			//Rcout << __FUNCTION__<< " out" <<endl;
			//this->print();
    };


    inline double flip_push_left( double left_val)
    {
	// input Breakpoints
	// x0,		x1,				x2,				...,		xn
	// f(x_0)	f(x1)-f(x0),	f(x2)-f(x1),	...,		f(xn)-f(xn-1)
	//output Breakpoints
	// left_val, 	f(x_0),		f(x1),		...,		f(xn-1),
	// x0,			x1-x0,		x2-x1,		...,		xn-xn-1,
    //	returns f(xn)
    	std::map<double,double>::iterator itbegin = Breakpoints_.begin();
    	std::map<double,double>::iterator itend = Breakpoints_.end();
    	std::map<double, double>::iterator i=itbegin;
    	//assert(left_val< itbegin->first);


		double x = FirstSlopeVal_;
		double y = i->first;
    	FirstSlopeVal_=i->first;
		i->second = 0;
		const_cast<double&>(i->first) = left_val;
		double old_second=0;
		++i;

		while ( i != itend)
		{
			double old_second = i->second+x;
			double old_first=i->first;
			i->second = i->first-y;
			const_cast<double&>(i->first) = x;
			x = old_second;
			y = old_first;
			++i;
		}

		return x;

    }


  inline void shift_Breakpoint( std::map<double,double>::iterator it,double shift)
  {
	  // about changing the key of an element in a map
	  // http://stackoverflow.com/questions/5743545/what-is-the-fastest-way-to-change-a-key-of-an-element-inside-stdmap
    std::swap(Breakpoints_[it->first+shift], it->second);
    Breakpoints_.erase(it);
  }

 inline void shift_Breakpoint( std::map<double,double>::reverse_iterator it,double shift)
  {
    double first=it->first;
    double second=it->second;
    map<double,double>::reverse_iterator j = it ; ++j;
    Breakpoints_.erase(j.base()); // http://stackoverflow.com/questions/1830158/how-to-call-erase-with-a-reverse-iterator
    Breakpoints_[first+shift]=second;
  }


 };// end of class cplfunction definition




#endif /* cplfunction_HPP_ */
