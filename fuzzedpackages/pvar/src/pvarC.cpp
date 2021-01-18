// p-variation calculus for piecewise monotone functions
// The main function is `pvarC`.
// Author and maintainer: Vygantas Butkus <Vygantas.Butkus@gmail.com>
// Please do not hesitate to contact me in any question.


#include <numeric>     
#include <cmath> 
#include <queue>        
#include <Rcpp.h>       

using namespace Rcpp;

// ######################################################################################### //
// ############################### inner(C++) functions #################################### //
// ######################################################################################### //

// -------------------------------- definitions of types  ---------------------------------- // 

// p-variation point. An object with necessary info.
struct pvpoint {
  int id;  
  double val;
  double pvdiff;
};

typedef std::list<pvpoint> PrtList;    // the list of p-variation point
typedef PrtList::iterator it_PrtList;  // the iterator of element in list of p-variation point 

// p-variation temporary points, used in specific calculates - it saves extra value need for calculations.
struct pvtemppoint{
  it_PrtList it;
  double ev;
};



// --------------------- printing functions (used only for debugging) ---------------------- // 

//template <class T>
//void print (T& obj, std::string ext1=" ", std::string ext2 = "\n"){
//  Rcout << ext1 << obj << ext2;
//}
//
//void print (pvpoint& obj, std::string ext1=" ",  std::string ext2="\n"){
//  Rcout << ext1 << "[" << obj.id << "]" << obj.val << ext2;
//}
//
//void print (it_PrtList& obj, std::string ext1=" ",  std::string ext2="\n"){
//  print((*obj), ext1, ext2);
//}
//
//template <class T>
//void printList (T& obj, std::string ext1="", std::string ext2="\n"){
//  typename T::iterator it=obj.begin();
//  Rcout << ext1;
//  for(it; it!=obj.end(); it++){
//     print(*it, " ", ", ");
//  }  
//  Rcout << ext2;
//}
//
//               
//void printList (std::list<pvtemppoint>& obj, std::string ext1="", std::string ext2="\n"){
//  std::list<pvtemppoint>::iterator it; 
//  Rcout << ext1;
//  for(it=obj.begin(); it!=obj.end(); it++){
//     print(*(*it).it, " ", ", ");
//  }  
//  Rcout << ext2;
//}



// ---------------------------- calculation functions ----------------------------- // 

// last iterator of the list. It is iterator to `obj.end()-1`
template <class T>
typename T::iterator last(T& obj){
  typename T::iterator it = obj.end();
  --it;
  return(it);  
}


// the difference used in p-variation, i.e. the abs power of diff.
double pvar_diff(double diff, double p){  
  return std::pow(std::abs(diff), p);
}

// finds change point of the vector and put it in list
PrtList ChangePoints(const NumericVector& x){
    
  // Main principle:
  // if point pt[i] is in increasing interval and pt[i]>pt[i+1], then 
  // pt[i] is change point (and vise versa).

  int dir = 0;
  int n = x.size();  
  
  pvpoint pvp;
  pvp.id = 0;
  PrtList out (1, pvp); // the first and last points are always included by definition.

  for(int i = 1; i < n; ++i) {
    if(x[i-1]<x[i]){
      if(dir<0){
        pvp.id = i-1;
        out.push_back (pvp);
      }     
      dir = 1; 
    }
    if(x[i-1]>x[i]){
      if(dir>0){
        pvp.id = i-1;
        out.push_back (pvp);
      }     
      dir = -1; 
    }        
  }  
  
  pvp.id =n-1;
  out.push_back (pvp); // the first and last points are always included by definition.

  return(out);
}

// create `prt` list from x and partition that is already known. No checking, just getting iterators and 
PrtList CreateBasicPrt(const NumericVector& x, const NumericVector& partition){
  
  int n = partition.size();  
  int i = 0;
  PrtList out (n);  
  
  for(PrtList::iterator it=out.begin(); it!=out.end(); it++){
    (*it).id = partition[i]-1;
    i++;
  }
     
  return out;
}


// finds(updates) all necessary attributes of `prt` list. `prt` must have good ids.
void prepare_prt(const NumericVector& x, std::list<pvpoint>& prt,  const double& p){

  PrtList::iterator it1_prt, it2_prt;
  it1_prt = it2_prt = prt.begin();
  ++it2_prt;

  // getting first element
  (*it1_prt).val = x[(*it1_prt).id];
  (*it1_prt).pvdiff = 0;
  
  // getting all other elements
  for ( ; it2_prt != prt.end(); it1_prt++, it2_prt++){
    (*it2_prt).val = x[(*it2_prt).id];
    (*it2_prt).pvdiff = pvar_diff(x[(*it1_prt).id] - x[(*it2_prt).id], p);
  }    
}

// sequentially checks small intervals of length d
void CheckSmallIntervalsOnce(PrtList& prt, const double& p,  const int& d){
  
  // Main principle:
  // if |pt[i] - pt[i+ d]|^p > sum_{j={i+1}}^d   |pt[j] - pt[j-1]|^p
  // then all middle points (i.e. p[j], j=i+1,...,i+d-1) are insignificant 

  int dcount = 0;
  double csum = 0;
  double fjoinval;
  
  PrtList::iterator it1_prt, it2_prt;
  it1_prt = it2_prt = prt.begin();
  ++it2_prt;
  
  for ( ; it2_prt != prt.end(); it2_prt++){
    ++dcount; 
    csum += (*it2_prt).pvdiff ;  
    if(dcount==d){
      fjoinval = pvar_diff((*it1_prt).val - (*it2_prt).val, p);
      ++it1_prt;  // in this stage it1_prt is always significant
      if(csum < fjoinval){ // mid points are insignificant, delete all
        dcount = 0;
        csum = 0;    
        it1_prt = prt.erase(it1_prt, it2_prt);
        (*it2_prt).pvdiff = fjoinval;
      }else{ // all significant, move one step foward forward
        csum -= (*it1_prt).pvdiff ;
        --dcount;
      }      
    }
  }  
}

// checks small intervals of up till length dn. 
// After this function, all points are significant in any small interval (i.e. interval with length no greater dn).
void CheckSmallIntervals(PrtList& prt, const double& p,  const int& dn){

  // Main principle:
  // apply CheckSmallIntervalsOnce starting form d=3 (because 3 is the minimal length worth checking)
  // If there was no change, apply CheckSmallIntervalsOnce with d=d+2 (because insignificant points goes only in pears) 
  // If there was change start from d=3 again (because `prt` changed, therefore, we are not sure if all smaller intervals are good).

  int LastSize = 0;
  int CurSize = prt.size();
  int d = 3;
  
  while((LastSize!=CurSize) & (CurSize>3) & (d<=dn)){
    d = 3;
    LastSize = CurSize;
    CheckSmallIntervalsOnce(prt, p, d);
    CurSize = prt.size();
    while((LastSize==CurSize) & (CurSize>d+2) & (d<dn)){
      d = d + 2;
      LastSize = CurSize;
      CheckSmallIntervalsOnce(prt, p, d);
      CurSize = prt.size();
    }
  }
}


// merge two intervals ([a, v] and [v, b]) which are known to be good.
void Merge2GoodInt(PrtList& prt,  const double& p, it_PrtList a, it_PrtList v, it_PrtList b){
  
  // Main principle:
  // 1. Find potential points in intervals [a,v) and (v, b] 
  //    (i.e. the points that could make a new f-joint with any point form opposite interval).
  //    Those points are find using cummin and cummac starting from v. 
  //     Some points might be dropped out before actual checking, but experiment showed, that it is not worthwhile.   
  // 2. Sequentially check all possible joints. If any increase is detected, then all middle points are insignificant.

  if (a==v or v==b) return ; // nothing to calculate, exit the procedure.

  double amin, amax, bmin, bmax, ev, balance, maxbalance, jfoin, takefjoin;
  it_PrtList prt_it, prt_ait, prt_bit;
  std::list<pvtemppoint> av, vb; 
  std::list<pvtemppoint>::iterator ait, bit, tit, tait, tbit, bitstart;
  pvtemppoint pvtp;
  
  // 1. ### Find potential points 

  // --- in interval [a,v) (starting from v).  
  ev = 0;
  prt_it = v;
  amin = amax= (*v).val;
  while(prt_it!=a){
    ev += (*prt_it).pvdiff;
    --prt_it;
    if((*prt_it).val>amax){
      amax=(*prt_it).val;
      pvtp.it = prt_it;
      pvtp.ev = ev;
      av.push_back (pvtp);
    }
    if((*prt_it).val<amin){
      amin=(*prt_it).val;
      pvtp.it = prt_it;
      pvtp.ev = ev;
      av.push_back (pvtp);
    }
  }
  // printList(av, "av :");
  
  // --- in interval (v,b] (starting from v). 
  ev = 0;                 
  prt_it = v;
  bmin = bmax = (*v).val;
  while(prt_it!=b){
    ++prt_it;
    ev += (*prt_it).pvdiff;
    if((*prt_it).val>bmax){
      bmax=(*prt_it).val;
      pvtp.it = prt_it;
      pvtp.ev = ev;
      vb.push_back (pvtp);
    }
    if((*prt_it).val<bmin){
      bmin=(*prt_it).val;
      pvtp.it = prt_it;
      pvtp.ev = ev;
      vb.push_back (pvtp);
    }    
  }
  // printList(vb, "vb :");

  // 2. ### Sequentially check all possible joints: finding the best i,j \in [a, v)x(v,b] that could be joined
  takefjoin = 0;
  maxbalance = 0;
  for(ait=av.begin(); ait!=av.end(); ait++){
    for(bit=vb.begin(); bit!=vb.end(); bit++){
      // std::cout <<  (*(*ait).it).id << " - " << (*(*bit).it).id << ":\n";
      jfoin = pvar_diff( (*(*ait).it).val - (*(*bit).it).val, p );
      balance = jfoin - (*bit).ev - (*ait).ev ;
      if (balance>maxbalance){
        maxbalance = balance;
        takefjoin = jfoin;
        tait = ait;
        tbit = bit;
      }
    } 
  }  
       
  // if we found any point, join it by erasing all middle points
  if(maxbalance>0){
    // joining:
    prt_ait = (*tait).it;
    ++prt_ait;
    prt_it = prt.erase(prt_ait, (*tbit).it);     
    (*prt_it).pvdiff = takefjoin;
  }   
}

// Modifies prt to become the partition of p-variation by merging all small in [a, b]
void PvarByMerging(PrtList& prt,  const double& p, PrtList::iterator a, PrtList::iterator b, int LSI=2){
  
  // Main principle:
  // 1. find all the intervals that should be merged.
  // 2. Apply merging by pears of interval. Repeat it until all intervals are merged.

  it_PrtList it1_prt, it2_prt, v, it, it2;
  std::list<it_PrtList> IterList;
  std::list<it_PrtList>::iterator it_IterList, a_IL, v_IL, b_IL;

  // 1. ### Finding all the intervals that will be merged
  it = a;
  int count = 0;
  while(it!=b){
    if(count % LSI == 0){
      IterList.push_back (it);
    }
    ++count;
    ++it;
  }
  IterList.push_back (it);

  // ### 2. Apply merging by pears of interval until everything is merged.
  while(IterList.size()>2){
    a_IL = v_IL = b_IL = IterList.begin();
    ++v_IL; 
    ++b_IL; 
    ++b_IL; 
    while((b_IL!=IterList.end()) & (v_IL!=IterList.end())){
      Merge2GoodInt(prt, p, *a_IL, *v_IL, *b_IL);
      a_IL = IterList.erase(v_IL); 
      v_IL = b_IL = a_IL;
      ++v_IL; 
      ++b_IL; 
      ++b_IL;      
    }
  }
}



// ######################################################################################### //
// ################# functions exported to R for normal use ################################ //
// ######################################################################################### //


//' Change Points of a \code{numeric} vector
//'
//' Finds changes points (i.e. corners) in the \code{numeric} vector.
//'
//' The end points of the vector will be always included in the results. 
//'
//' @return The vector of index of change points.
//' @param x \code{numeric} vector.
//' @export
//' @examples
//' x <- rwiener(100)
//' cid <- ChangePoints(x)
//' plot(x, type="l")
//' points(time(x)[cid], x[cid], cex=0.5, col=2, pch=19)
//[[Rcpp::export("ChangePoints")]]
IntegerVector ChangePoints_fromR(const NumericVector& x){
  
  // ### input checking
  // x!=NA
  for (int i = 0; i < x.size(); ++i) {
    if(NumericVector::is_na(x[i])){
      stop("In `ChangePoints` function, `x` must not have NA values!"); 
    }
  }    
  // the length of x
  if(x.size()<2){
    if(x.size()<1){
      stop("In `ChangePoints` function, the length of `x` must be strictly positive!");      
    }
    IntegerVector out(1);
    out[0] = 1;
    return (out);
  }
    
  // ### program it self:  
  
  PrtList::iterator it_prt;
  IntegerVector::iterator it_IV;

  PrtList prt = ChangePoints(x);

  IntegerVector out(prt.size());  
  for (it_prt = prt.begin(), it_IV=out.begin(); it_prt != prt.end(); it_prt++, it_IV++){
    *it_IV = (*it_prt).id + 1;
  }  
  return(out);  
}

//' p-variation calculation (in C++)
//' 
//' An internal function(written in C++) that calculates p-variation. 
//' 
//' This is a waking horse of this packages, nonetheless, users should 
//' not call this function directly (rather use \code{\link{pvar}}).
//' 
//' @return An object of the class \code{pvar}.
//' @keywords internal
//' @inheritParams  pvar
//' @export
//[[Rcpp::export("pvarC")]]
List pvarC(const NumericVector& x, double& p, int LSI=3){
  
  // ##### Checking of possible errors. Must be done, otherwise this could crash R.
  // p!=NA
  if(NumericVector::is_na(p)){
    stop("In `pvarC` function, the value of `p` must not be NA!");
  }  
  // p>1
  if(p<=1){
    stop("In `pvarC` function, the value of `p` must be greater then 1!");
  } 
  // LSI!=NA
  if(NumericVector::is_na(LSI)){
    stop("In `pvarC` function, the value of `LSI` must not be NA!");
  }    
  // LSI = 1, 2, ..
  if(LSI<1){
    stop("In `pvarC` function, the `LSI` must be positive integer!");
  }    
  // x!=NA
  for (int i = 0; i < x.size(); ++i) {
    if(NumericVector::is_na(x[i])){
      stop("In `pvarC` function, `x` must not have NA values!"); 
    }
  }  
  // the length of x
  if(x.size()<2){
    if(x.size()<1){
      stop("In `pvarC` function, the length of `x` must be strictly positive!");      
    }
    List out = List::create(
        _["value"] = NumericVector::create(_["p-variation"]=0),
        _["x"] = x,
        _["p"] = p,
        _["partition"] = 1
    ) ;          
    out.attr("class") = "pvar";
    return out;
  }
  
  
  // ##### Program it self:
  
  PrtList::iterator it_prt;   
  PrtList prt = ChangePoints(x);
  prepare_prt(x, prt, p);  
  
  CheckSmallIntervals(prt, p, LSI);
  PvarByMerging(prt, p, prt.begin(), last(prt), LSI+1);  
    
  // output:
  double pvalue=0;;
  NumericVector partition(prt.size());  
  int i = 0;
  for(it_prt=prt.begin(); it_prt!=prt.end(); it_prt++){
    pvalue += (*it_prt).pvdiff;
    partition[i] = (*it_prt).id + 1;
    ++i;    
  }
  
  List out = List::create(
      _["value"] = NumericVector::create(_["p-variation"]=pvalue),
      _["x"] = x,
      _["p"] = p,
      _["partition"] = partition
  ) ;          
  out.attr("class") = "pvar";
  return out;

}

//' Addition of p-variation (in C++)
//' 
//' An internal function(written in C++) that merges two objects of pvar and effectively recalculates the p-variation of joined sample.
//' 
//' This is an internal function, therefore, users should 
//' not call this function directly (rather use \code{\link{AddPvar}} or \code{pv1 + pv2}).
//' 
//' @return An object of the class \code{pvar}.
//' @keywords internal
//' @inheritParams  AddPvar
//' @export
//[[Rcpp::export("AddPvarC")]]
List AddPvar(List PV1, List PV2, bool AddIfPossible=true){
    
  // ##### Checking of possible errors.  
  // p1 = p2
  NumericVector pr1 = PV1["p"];
  NumericVector pr2 = PV2["p"];
  if(pr1[0]!=pr1[0]){
    stop("In `AddPvar` function, `p` attributes in PV1 and PV2 must be equal!");
  }  
  double* pp = &pr1[0];
  double p = *pp;
  // p>1
  if(p<=1){
    stop("In `AddPvarC` function, the value of `p` must be greater then 1!");
  }   
    
  // ##### Program:
  PrtList::iterator it_prt, it2_prt, v1, v2;  
  bool Add;
  int newn;

  NumericVector x1 = PV1["x"];
  NumericVector x2 = PV2["x"];
  int n1 = x1.size();    
  int n2 = x2.size();
  
  NumericVector partition1 = PV1["partition"];
  NumericVector partition2 = PV2["partition"];  
  
  PrtList prt1 = CreateBasicPrt(x1, partition1); 
  PrtList prt2 = CreateBasicPrt(x2, partition2); 
  
  prepare_prt(x1, prt1, p);  
  prepare_prt(x2, prt2, p);  
  
  Add = AddIfPossible & (x1[n1-1]==x2[0]);
  newn = n1 + n2;
  if(Add) --newn;
  NumericVector newx(newn);

  // concatenating x:
  std::copy (x1.begin(), x1.end(), newx.begin());
  NumericVector::iterator it1 = newx.begin();
  NumericVector::iterator it2 = x2.begin();
  std::advance(it1, n1);
  if(Add) ++it2;
  std::copy (it2, x2.end(), it1);    
  
  // concatenating prt:  
  v1 = v2 = last(prt1);  //saving joining vertex, before merge  
  int addid = n1;
  if(Add){
    prt2.erase(prt2.begin());
    addid--;
  }
  if(!prt2.empty()){
    for(it_prt=prt2.begin(); it_prt!=prt2.end(); it_prt++){
      (*it_prt).id = (*it_prt).id + addid;      
    }
    prt1.splice(prt1.end(), prt2); 
    ++v2;
    (*v2).pvdiff = pvar_diff(newx[(*v2).id] - newx[(*v1).id], p); // safe way is to use :prepare_prt(newx, prt1, p); 
  }  
  
  // recalculating
  if(Add){
    Merge2GoodInt(prt1, p, prt1.begin(), v1, last(prt1));    
  }else{
    v2 = v1;
    ++v2;
    Merge2GoodInt(prt1, p, prt1.begin(), v1, v2);
    Merge2GoodInt(prt1, p, prt1.begin(), v2, last(prt1));
  }

  // removing monotonic points
  if(prt1.size()>2){
    it_prt = prt1.begin();
    ++it_prt;
    while(it_prt!=last(prt1)){
      if((*it_prt).pvdiff==0){
        it_prt = prt1.erase(it_prt);
      }else{
        it_prt++;
      }
    }
    // checking last point (it always significants, so if pvdiff=0 then the one before last is insignificant)
    if((prt1.size()>2) & ((*it_prt).pvdiff==0)){  
      --it_prt;
      it_prt = prt1.erase(it_prt);
    }
  }
  prepare_prt(newx, prt1, p); // updating the values
  
  // output
  double pvalue=0;
  NumericVector partition(prt1.size());  
  int i = 0;
  for(PrtList::iterator it_prt=prt1.begin(); it_prt!=prt1.end(); it_prt++){
    pvalue += (*it_prt).pvdiff;
    partition[i] = (*it_prt).id + 1;
    ++i;    
  }
  
  List out = List::create(
      _["value"] = NumericVector::create(_["p-variation"]=pvalue),
      _["x"] = newx,
      _["p"] = p,
      _["partition"] = partition
  ) ;          
  out.attr("class") = "pvar";
  return out;  
}

 
// ######################################################################################### //
// ############# function exported only for testing purpose ################################ //
// ######################################################################################### //


//[[Rcpp::export("CheckSmallIntervals")]]
NumericVector test_CheckSmallIntervals(const NumericVector& x, const double& p,  const int& dn){
  std::list<pvpoint> prt = ChangePoints(x);
  prepare_prt(x, prt, p);  
  CheckSmallIntervals(prt, p,  dn);
  
  // output:
  NumericVector out(prt.size());
  NumericVector::iterator nv_it;
  std::list<pvpoint>::iterator prt_it;  
  for (prt_it = prt.begin(), nv_it=out.begin(); prt_it != prt.end(); prt_it++, nv_it++){
    *nv_it = (*prt_it).id + 1;
  }  
  return(out);  
}


//[[Rcpp::export("prepare_prt")]]
List test_prepare_prt(const NumericVector& x, const double& p){
  std::list<pvpoint> prt = ChangePoints(x);
  prepare_prt(x, prt, p);  
  
  std::list<pvpoint>::iterator prt_it=prt.begin();
  
  int i = 0;
  int n = prt.size();
  IntegerVector id(n);
  LogicalVector type(n);
  NumericVector val(n);
  NumericVector pvdiff(n);
  
  for ( ; prt_it != prt.end(); prt_it++){    
    id[i] = (*prt_it).id+1;
    val[i] = (*prt_it).val;
    pvdiff[i] = (*prt_it).pvdiff;
    ++i;
  }  
  
  return List::create(
      Named("x") = val,
      Named("id") = id,
      Named("type") = type,
      Named("pvdiff") = pvdiff
  ) ;  
}
