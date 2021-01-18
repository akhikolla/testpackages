// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#define container vector         // Could be 'vector' or 'deque' (both work but there may be performance differences)
#define USE_UNORDERED_MAP true   // set to true for unordered_map; comment out to use plain stl map.

#include <Rcpp.h>
#include <cmath>

#include <string.h>
#include <iostream>
#include <unordered_map>
#include <vector>
#include <deque>
#include <utility>
#include <iterator>

using namespace std;
using namespace Rcpp; 

typedef container<signed int> mycont;  // a mycont  is a container [vector or deque] of *signed* ints.

#ifdef USE_UNORDERED_MAP
class MyVecHasher
{
public:
       size_t operator()(const mycont& vec) const
       {
              // thanks to Steffan Hooper for advice
              std::size_t seed = 0;
              for (auto& i : vec){
                  seed ^= i + 0x9e3779b9 + (seed << 6) + (seed >> 2);
              }
              return seed;
       }
};
 
typedef std::unordered_map<mycont, double , MyVecHasher > spray;
#else 
typedef map<mycont, double > spray;
#endif



spray prepare(const IntegerMatrix M, const NumericVector d){
    spray S;
    mycont v;
    signed int i,j;
    spray::iterator it;

    for(i=0; i<M.nrow() ; i++){
        if(d[i] != 0){
                v.clear();
                for(j=0; j<M.ncol(); j++){
                    v.push_back(M(i,j));
                }
                S[v] += d[i];
        }
    }  // i loop closes

    // Now remove zero entries:
    it = S.begin();
    while(it != S.end()){
        if(it->second == 0){
            it = S.erase(it); //  in C++11, erase() returns *next* iterator
        } else {
            ++it;  // else just increment the iterator
        }
    }
    return(S);
}

IntegerMatrix makeindex(const spray S){  // takes a spray, returns the matrix of indices
    const unsigned int ncol = S.begin()->first.size();
    IntegerMatrix  out(S.size(),ncol);   // index
    mycont v;
    unsigned int row=0, col=0;
    mycont::const_iterator ci;  // ci iterates through a container
    spray::const_iterator it;   // it iterates through a sparse array

    for(it=S.begin(); it != S.end(); ++it){
        v = it->first;
        col = 0;
        for(ci=v.begin() ; ci != v.end() ; ++ci){
            out(row,col++) = *ci;
        }
        row++;
    }
    return(out);
}

NumericVector makevalue(const spray S){  // takes a spray, returns data
    NumericVector  out(S.size());   // data
    unsigned int i=0;
    spray::const_iterator it;   // it iterates through a sparse array
    
    for(it=S.begin(); it != S.end(); ++it){
        out(i++) = it->second;   // initialize-and-fill is more efficient than  out.push_back(it->second) 
    }
    return(out);
}

List retval (const spray &S){  // used to return a list to R

  // In this function, returning a zero-row matrix results in a
  // segfault ('memory not mapped').  So we check for 'S' being zero
  // size and, if so, return a special Nil value.  This corresponds to
  // an empty spray object.
  
    if(S.size() == 0){
        return List::create(Named("index") = R_NilValue,
                            Named("value") = R_NilValue
                            );
    } else {
        return List::create(Named("index") = makeindex(S),
                            Named("value") = makevalue(S)
                            );
    }
}


// [[Rcpp::export]]
List spray_maker
(
  const IntegerMatrix &M, const NumericVector &d
 ){
    return retval(prepare(M,d));
}

// [[Rcpp::export]]
List spray_add
(
 const IntegerMatrix &M1, const NumericVector &d1,
 const IntegerMatrix &M2, const NumericVector &d2
 ){
     spray S1, S2;   // basic idea is S1 = S1+S2;
     spray::const_iterator it;   // it iterates through a sparse array
     mycont v;
     
    S1 = prepare(M1, d1);
    S2 = prepare(M2, d2);
    
    // the "meat" of this function, namely S1=S1+S2 (S1 += S2):
    for (it=S2.begin(); it != S2.end(); ++it){
        v = it->first;
        S1[v] += S2[v];   //S1 to be returned.  NB: S1 may have increased in size.
        if(S1[v]==0){S1.erase(v);}
    }

    return retval(S1);
}


spray prod //
(
 const spray S1, const spray S2
 ){
    spray Sout;
    spray::const_iterator it1,it2;
    mycont vsum;
    unsigned int i;

    // the "meat" of this function:  Sout=S1*S2
    for (it1=S1.begin(); it1 != S1.end(); ++it1){
        const mycont v1 = it1->first;
        const double x1 = it1->second;
        for (it2=S2.begin(); it2 != S2.end(); ++it2){
            const mycont v2 = it2->first;
            const double x2 = it2->second;
            vsum.clear();
            for(i=0; i<v1.size(); i++){
                vsum.push_back(v1[i] + v2[i]);  // meat 1: powers add
            }
            Sout[vsum] += x1*x2;                // meat 2: coefficients multiply
        }
    }
    //    for(spray::iterator it=S3.begin(); it != S3.end(); ++it){
    //    if(it->second ==0){S3.erase(it);}
    // }

    return Sout;

}    

spray unit //
(
 unsigned int n
){
    const IntegerMatrix M(1,n);
    NumericVector one(1);
    one[0] = 1;
    return prepare(M,one);
}

// [[Rcpp::export]]
List spray_mult  // multiply two sprays
(
 const IntegerMatrix &M1, const NumericVector &d1,
 const IntegerMatrix &M2, const NumericVector &d2
 ){
    return retval(prod(prepare(M1,d1),prepare(M2,d2)));
}

// [[Rcpp::export]]
List spray_overwrite // something like S1[ind(S2)] <- S2
(
 const IntegerMatrix &M1, const NumericVector &d1,
 const IntegerMatrix &M2, const NumericVector &d2
 ){
    spray S1,S2;
    mycont v;
    spray::const_iterator it;
    
    S1 = prepare(M1, d1);
    S2 = prepare(M2, d2);    

    for (it=S2.begin(); it != S2.end(); ++it){
        v = it->first;
        S1[v] = S2[v];   // the meat
    }

    return retval(S1);
}

// [[Rcpp::export]]
NumericVector spray_accessor // returns S1[]
(
 const IntegerMatrix &M, const NumericVector &d,
 const IntegerMatrix &Mindex
 ){
    spray S;
    mycont v;
    signed int i,j,k=0;
    NumericVector out(Mindex.nrow());
    
    S = prepare(M, d);

    for(i=0; i<Mindex.nrow() ; i++){
        v.clear();
        for(j=0; j<Mindex.ncol(); j++){
            v.push_back(Mindex(i,j));
        }
        out[k++] = S[v];
    }
    return out;
}


// [[Rcpp::export]]
List spray_setter // effectively S[M] <- d; return S
(
 const IntegerMatrix &M1, const NumericVector &d1,
 const IntegerMatrix &M2, const NumericVector &d2    // M2 -> index ; d2 -> value
 ){
    spray S1,S2;
    mycont v;
    signed int i,j;
    
    S1 = prepare(M1, d1);
    S2 = prepare(M2, d2);

    for(i=0; i<M2.nrow() ; i++){
        v.clear();
        for(j=0; j<M2.ncol(); j++){
            v.push_back(M2(i,j));
        }
        S1[v] = S2[v];
    }
    return retval(S1);
}


// [[Rcpp::export]]
bool spray_equality // S1 == S2
(
 const IntegerMatrix &M1, const NumericVector &d1,
 const IntegerMatrix &M2, const NumericVector &d2
 ){
    mycont v;
    spray S1, S2;  
    spray::const_iterator it;
    
    S1 = prepare(M1, d1);
    S2 = prepare(M2, d2);

    if(S1.size() != S2.size()){
        return FALSE;
    }

    for (it=S1.begin(); it != S1.end(); ++it){
        v = it->first;
        if(S1[v] != S2[v]){
            return FALSE;
        } else {
            S2.erase(v);
        }
    }
    // at this point, S1[v] == S2[v] for every index 'v' of S1;  S1\subseteq S2.
    // We need to check that every element of S2 has been accounted for:
    
    if(S2.empty()){
        return TRUE;
    } else {
        return FALSE;
    }
}

// [[Rcpp::export]]
List spray_asum_include
(
 const IntegerMatrix &M, const NumericVector &d,
 const IntegerVector &n
 ){
    spray S;
    mycont v;
    signed int i,j,k;

    for(i=0; i<M.nrow() ; i++){
        v.clear();
        for(j=0; j<M.ncol(); j++){
            v.push_back(M(i,j));
        }
        for(k=0 ;  k<n.size() ; k++){
            v[n[k]-1] = 0;    // off-by-one issue dealt with here: if n[k]=0 this means dimension 1.
        }
        S[v] += d[i];
    }
    return retval(S);
}


// [[Rcpp::export]]
List spray_asum_exclude 
(
 const IntegerMatrix &M, const NumericVector &d,
 const IntegerVector &n
 ){
    spray S;
    mycont v;
    signed int i,j;

    for(i=0; i<M.nrow() ; i++){
        v.clear();
        for(j=0; j<n.size(); j++){
            v.push_back(M(i,n[j]-1));   //off-by-one error
        }
        S[v] += d[i];
    }
    return retval(S);
}

//spray(matrix(sample(0:5,100,rep=T),ncol=5) )-> S

// [[Rcpp::export]]
List spray_deriv
(
 const IntegerMatrix &M,  const NumericVector &d,   
 const IntegerVector &n
){
    spray S;
    mycont v;
    signed int i,j,nn;

    IntegerMatrix Mout(M.nrow(),M.ncol());
    NumericVector dout(d.size());

    // create local copies of M and d:
    for(i=0 ; i<M.nrow() ; i++){
        dout[i] = d[i];
        for(j=0 ; j<M.ncol() ; j++){
            Mout(i,j) = M(i,j);
        }
    }
    
    for(i=0; i<Mout.nrow() ; i++){
        for(j=0; j<Mout.ncol() ; j++){
            nn = n[j];
            while( (nn>0) & (d[i]!=0)  ){  // while loop because it might not run at all
                dout[i] *= Mout(i,j);  // multiply d first, then decrement M (!)
                Mout(i,j)--;
                nn--;
            }
        }
        v.clear();
        for(j=0; j<Mout.ncol() ; j++){
            v.push_back(Mout(i,j));
        }                    
        S[v] += dout[i];   // increment because v is not row-unique any more
    }  // i loop closes
    return retval(S);
}

// [[Rcpp::export]]
List spray_pmax
(
 const IntegerMatrix &M1, const NumericVector &d1,
 const IntegerMatrix &M2, const NumericVector &d2 
 ){
    spray S1,S2;
    spray::const_iterator it;   // it iterates through a sparse array
    
    S1 = prepare(M1, d1);
    S2 = prepare(M2, d2);

    for (it=S1.begin(); it != S1.end(); ++it){
        const mycont v = it->first;
        if(S2[v] > S1[v]){ S1[v] = S2[v];} // S1[v] = max(S1[v],S2[v]);
        S2.erase(v); // not S2[v] = 0;  // OK because the iterator is it1 and this line modifies S2
    }
            
    for (it=S2.begin(); it != S2.end(); ++it){ //iterate through S2 keys not in S1
        const mycont v = it->first;
        if(S2[v] > 0){ S1[v] = S2[v]; }
    }

    return retval(S1);
}

// [[Rcpp::export]]
List spray_pmin
(
 const IntegerMatrix &M1, const NumericVector &d1,
 const IntegerMatrix &M2, const NumericVector &d2 
 ){
    spray S1,S2;
    spray::const_iterator it;   // it iterates through a sparse array
    
    S1 = prepare(M1, d1);
    S2 = prepare(M2, d2);

    for (it=S1.begin(); it != S1.end(); ++it){
        const mycont v = it->first;
        if(S2[v] < S1[v]){ S1[v] = S2[v]; }// S1[v] = min(S1[v],S2[v]);
        S2.erase(v);
    }
            
    for (it=S2.begin(); it != S2.end(); ++it){
        const mycont v = it->first;
        if(S2[v] < 0){S1[v] = S2[v]; } // S1[v] = min(S2[v],0);
    }

    return retval(S1);
}


// [[Rcpp::export]]
List spray_power
(
 const IntegerMatrix &M, const NumericVector &d, const NumericVector &pow
 ){
    spray out = unit(M.ncol());
    const spray S = prepare(M,d);
    unsigned int n=pow[0];

    if(n == 0){
        return retval(out);
    } else if (n==1){
        return retval(S);
    } else {  // n>1
        for( ; n>0; n--){
            out = prod(S,out);
        }
    }
    return retval(out);
}



