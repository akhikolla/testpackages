// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include <set>
#include <map>
#include <iostream>
#include <Rcpp.h>
#include <iterator>
#include <vector>
#include <assert.h>

using namespace std;
using namespace Rcpp;

typedef set<unsigned int> bracket;
typedef map<bracket, long double> hyper2;

// again it might be nice to use unsigned_map above, but this would
// need a hash function and this would require further work.

hyper2 prepareL(const List &L, const NumericVector &d){
    hyper2 H;   
    bracket b;
    unsigned int i,j;
    const unsigned n=L.size();
    Rcpp::IntegerVector iv;
    
    for(i=0; i<n ; i++){
        if(d[i] != 0){
            b.clear();
            iv = as<Rcpp::IntegerVector> (L[i]);
            for(j=0 ; j<(unsigned int) iv.size(); j++){
                b.insert(iv[j]);
            }
        H[b] += d[i];
        }
    }
    return(H);
}

List makebrackets(const hyper2 H){  // takes a hyper2, returns the brackets
    List out;
    hyper2::const_iterator ih;
    
    for(ih=H.begin(); ih != H.end(); ++ih){
        out.push_back(ih->first);
    }
    return(out);
}

NumericVector makepowers(const hyper2 H){  // takes a hyper2, returns powers
    NumericVector out(H.size());   // data
    unsigned int i=0;
    hyper2::const_iterator it;   // it iterates through a hyper2 object
    
    for(it=H.begin(); it != H.end(); ++it){
        out(i++) = it->second;   // initialize-and-fill is more efficient than  out.push_back(it->second) 
    }
    return(out);
}

List retval(const hyper2 &H){  // used to return a list to R
    
        return List::create(Named("brackets") =  makebrackets(H),
                            Named("powers")   =  makepowers(H)
                            );
}

// [[Rcpp::export]]
List identityL(const List &L, const NumericVector &p){
    const hyper2 out = prepareL(L,p);
    return retval(out);
}
 
//[[Rcpp::export]]
List addL(
          const List L1, const NumericVector p1,
          const List L2, const NumericVector p2  // p==powers
          ){
    hyper2 h1 = prepareL(L1,p1);
    hyper2 h2 = prepareL(L2,p2);
    hyper2::const_iterator it;
    bracket b;
    if(L1.size() > L2.size()){ // L1 is bigger, so iterate through L2
        for (it=h2.begin(); it != h2.end(); ++it){
            b = it->first;
            h1[b] += h2[b];  
        }
        return(retval(h1));
    } else {  // L2 is bigger
        for (it=h1.begin(); it != h1.end(); ++it){
            b = it->first;
            h2[b] += h1[b];  
        }
        return(retval(h2));
    }
}

//[[Rcpp::export]]
bool equality(  // modelled on spray_equality()
          const List L1, const NumericVector p1,
          const List L2, const NumericVector p2  // p==powers
           ){
    hyper2 h1,h2;
    hyper2::const_iterator it;
    bracket b;
   
    if(L1.size() != L2.size()){
        return false;
    }

    h1 = prepareL(L1,p1);
    h2 = prepareL(L2,p2);

    for (it=h1.begin(); it != h1.end(); ++it){
            b = it->first;
            if(h1[b] != h2[b]){
                return false;
            } else {
                h2.erase(b);
            }
        }

    if(h2.empty()){
        return true;
    } else {
        return false;
    }
}

//[[Rcpp::export]]
List accessor(
              const List L,
              const NumericVector powers,
              const List Lwanted
              ){

    const hyper2 h=prepareL(L,powers);
    hyper2 out;
    bracket b;
    const unsigned int n=Lwanted.size();
    unsigned int i,j;
    Rcpp::IntegerVector iv;

    for(i=0; i<n ; i++){
            b.clear();
            iv = as<Rcpp::IntegerVector> (Lwanted[i]);
            for(j=0 ; j<(unsigned int) iv.size(); j++){
                b.insert(iv[j]);
            }
            if(h.count(b)>0){
                out[b] = h.at(b);
            }
    }
    return retval(out);
}
    
//[[Rcpp::export]]
List overwrite(  // H1[] <- H2
              const List L1, const NumericVector powers1,
              const List L2, const NumericVector powers2
              ){

          hyper2 h1=prepareL(L1,powers1);
    const hyper2 h2=prepareL(L2,powers2);
    bracket b;
    hyper2::const_iterator it;

    for(it=h2.begin(); it != h2.end(); ++it){
        b = it->first;
        h1[b] = h2.at(b);
    }
    return retval(h1);
}

//[[Rcpp::export]]
List assigner(  // H[L] <- v
            const List L, const NumericVector p,
            const List L2,
            const NumericVector value
              ){
    hyper2 h=prepareL(L,p);
    bracket b;
    hyper2::const_iterator it;
    const unsigned int n=L2.size();
    unsigned int i,j;
    Rcpp::IntegerVector iv;

    for(i=0 ; i<n ; i++){
        b.clear();
        iv = as<Rcpp::IntegerVector> (L2[i]);
        for(j=0 ; j<(unsigned int) iv.size(); j++){
            b.insert(iv[j]);
        } 
        h[b] = value[i]; // RHS might be zero in which case this entry is deleted from h
    }
    return retval(h);
}
                     
//[[Rcpp::export]]
double evaluate(  // returns log-likelihood
                const List L,
                const NumericVector powers,
                const NumericVector probs
                  ){

    const hyper2 h = prepareL(L,powers);
    hyper2::const_iterator it;
    bracket::const_iterator ib;
    bracket b;
    double out=0;
    double bracket_total=0;
    
    for (it=h.begin(); it != h.end(); ++it){
        b = it->first;
        bracket_total = 0;
        for (ib=b.begin(); ib != b.end(); ++ib){
            bracket_total += probs[*ib-1];  //NB off-by-one error!
        }
        out += log(bracket_total) * (it->second);
    }
    return out;
}

double differentiate_single( // d(log-likelihod)/dp
                 const hyper2 h,
                 unsigned int i,   // component to diff WRT
                 const unsigned int n,   // p_1 + ...+p_n = 1
                 const NumericVector probs
                 /* probs[0] == p_1, ..., probs[n-1] = p_n.  
                    probs[n-1] = p_n = fillup value.

                    bracket b has elements in the range 1,2,...,n.
                    If, say, 3 is an element of b, then this
                    corresponds to p_3==probs[2].  In particular, if
                    n\in b, then this is the fillup.
                  */
                      ){

    hyper2::const_iterator it;
    bracket::const_iterator ib;
    bracket b;
    double out;
    double bracket_total;
    double power;
    unsigned int no_of_diff_terms; // number of p_i terms in the bracket
    unsigned int no_of_fill_terms; // number of p_n terms (=fillup) in the bracket

    i++; // off-by-one
    assert(i > 0); 
    assert(i < n);  // sic; strict

    out = 0;
    for (it=h.begin(); it != h.end(); ++it){  // standard hyper2 loop
        b = it->first;

        no_of_diff_terms = b.count(i);
        no_of_fill_terms = b.count(n);
            
        if(no_of_diff_terms == no_of_fill_terms){continue;} // bracket has same number (maybe 0!) of diff and fillups

        bracket_total = 0;
        for (ib=b.begin(); ib != b.end(); ++ib){
            bracket_total += probs[*ib-1];  //NB off-by-one error! ... p_1 == probs[0]
        }

        power = it->second;
        // the 'meat':
        out += no_of_diff_terms*power/bracket_total;
        out -= no_of_fill_terms*power/bracket_total; 
    }
    return out;
}

//[[Rcpp::export]]
List differentiate(  // returns gradient of log-likelihood
                const List L,
                const NumericVector powers,
                const NumericVector probs,
                const unsigned int n
                  ){

    unsigned int i;
    NumericVector out(n-1);  
    const hyper2 h=prepareL(L,powers);

    for(i=0; i<n-1; i++){
        out[i] = differentiate_single(h,i,n,probs);
    }
        return List::create(Named("grad_comp") =  out);
}
