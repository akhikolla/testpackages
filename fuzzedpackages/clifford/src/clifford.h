// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// [[Rcpp::depends(BH)]]


// NB: terminology.  Here "blade" means "basis blade"

#include <map>
#include <iostream>
#include <Rcpp.h>
#include <iterator>
#include <vector>
#include <assert.h>
#include <tuple>
#include <boost/dynamic_bitset.hpp>

using namespace std;
using namespace Rcpp;

typedef boost::dynamic_bitset<> blade;  
typedef map<blade, long double> clifford;
typedef std::tuple<blade, int> blade_and_sign;

clifford remove_zeros(clifford &C){
   for(auto it=C.begin() ; it != C.end() ;){
        if(it->second == 0){
            it = C.erase(it); //increments pointer
        } else {
            ++it; // increment anyway
        }
    } 
    return C;
}

clifford prepare(const List &L, const NumericVector &d, const NumericVector &m){
    clifford out;
    const size_t n=L.size();
    for(size_t i=0 ; i<n ; i++){
        if(d[i] != 0){
            Rcpp::IntegerVector iv = as<Rcpp::IntegerVector> (L[i]);
            blade b;
            b.resize(m[0]+1);  //off-by-one
            for(int j=0 ; j < iv.size(); j++){
                b[iv[j]] = 1;
            }
            out[b] += d[i];  // the meat
        } // if d[i] closes
    }  // i loop closes
    return remove_zeros(out);
}

Rcpp::IntegerVector which(const blade b){ // takes a blade, returns which(blade)
    Rcpp::IntegerVector out;
    unsigned int i;
    for(i=0 ; i<b.size() ; ++i){
        if((bool) b[i]){
            out.push_back(i); // the meat; off-by-one here
        }
    }
    return out;
}

List Rblades(const clifford C){  // takes a clifford object, returns a list of which(blades); used in retval()
    List out;
    clifford::const_iterator ic;
    for(ic=C.begin(); ic != C.end(); ++ic){
        out.push_back(which(ic->first));
    }
    return out;
}

NumericVector coeffs(const clifford C){  // takes a clifford object, returns the coefficients
    NumericVector out(C.size());
    unsigned int i=0;
    clifford::const_iterator ic;   // it iterates through a hyper2 object
    
    for(ic=C.begin(); ic != C.end(); ++ic){
        out[i] = ic->second;
        i++;
    }
    return out;
}

List retval(const clifford &C){  // used to return a list to R
        return List::create(Named("blades") =  Rblades(C),
                            Named("coeffs") =  coeffs(C)
                            );
}

clifford c_add(clifford cliff1, clifford cliff2){
    clifford::iterator ic;
    if(cliff1.size() > cliff2.size()){ // #1 is bigger, so iterate through #2
        for (ic=cliff2.begin(); ic != cliff2.end(); ++ic){
            const blade b = ic->first;
            cliff1[b] += cliff2[b];  
        }
        return remove_zeros(cliff1);
    } else {  // L2 is bigger
        for (ic=cliff1.begin(); ic != cliff1.end(); ++ic){
            const blade b = ic->first;
            cliff2[b] += cliff1[b];  
        }
        return remove_zeros(cliff2);
    }
}

blade_and_sign juxtapose(blade b1, blade b2, const signed int signature){//juxtaposes two blades, returns reduction and sign
    int sign = 1;
    blade bout;
    const size_t m = max(b1.size(),b2.size());

    b1.resize(m, false);
    b2.resize(m, false);
    bout.resize(m, false);

    for(size_t i=0 ; i<m ; ++i){
        
        if       (((bool)~b1[i]) & ((bool)~b2[i])){ bout[i] = false;  // neither
        } else if(((bool) b1[i]) & ((bool)~b2[i])){ bout[i] = true;   // just b1
        } else if(((bool)~b1[i]) & ((bool) b2[i])){ bout[i] = true;   // just b2
        } else if(((bool) b1[i]) & ((bool) b2[i])){ bout[i] = false;  // both, but...
            if(signature>0){ // ...swap sign!  
                if((signed int) i > (signed int) signature){ // NB check for off-by-one error
                    sign *= -1;
                }
            } else if(signature < 0){
                sign = 0; // exterior product, repeated index -> 0
            } else { // signature == 0
                ;  // no swaps!
            }
        }
    }
        
    for(size_t i=0 ; i<m ; ++i){
        for(size_t j=i ; j<m ; ++j){
            if((b2[i] & b1[j]) && (i<j)){
                sign *= -1;}
        }
    }
    return std::make_tuple(bout,sign);
}

clifford c_general_prod(const clifford C1, const clifford C2, const NumericVector signature, bool (*chooser)(const blade, const blade)){

    clifford out;
    clifford::const_iterator ic1,ic2;
    blade b;
    int sign;
    for(ic1=C1.begin(); ic1 != C1.end(); ++ic1){
        const blade b1 = ic1->first;
        for(ic2=C2.begin(); ic2 != C2.end(); ++ic2){
            const blade b2 = ic2->first;
            if(chooser(b1,b2)){
                tie(b, sign) = juxtapose(b1, b2, signature[0]);
                out[b] += sign*(ic1->second)*(ic2->second); // the meat
            }
        }
    }
    return remove_zeros(out);
}


bool c_equal(clifford C1, clifford C2){
    // modelled on spray_equality()
    if(C1.size() != C2.size()){
        return false;
    }

    for (clifford::const_iterator ic=C1.begin(); ic != C1.end(); ++ic){
        blade b = ic->first;
        if(C1[b] != C2[b]){
            return false;
        }
    }
    
    return true;
}

clifford c_grade(const clifford C, const NumericVector &n){
    clifford out;
    for(size_t i=0 ; i < (size_t) n.length() ; ++i){
        for(clifford::const_iterator ic=C.begin() ; ic != C.end() ; ++ic){
            const blade b = ic->first;
            if(b.count() == (size_t) n[i]){
                out[b] = ic->second;
            }
        }
    }
    return out;
}

NumericVector c_coeffs_of_blades(clifford C,
                                 const List &B,
                                 const NumericVector &m
                                 ){  
    Rcpp::NumericVector out;
    for(size_t i=0 ; i < (size_t) B.size() ; ++i){
        blade b;
        b.resize(m[0]+1);  //off-by-one; note that this code also appears in prepare()
        const IntegerVector iv = B[i];
        for(size_t j=0 ; j < (size_t) iv.size(); j++){
            b[iv[j]] = 1;
        }
        out.push_back(C[b]);
    }
    return out;
}

bool geometricproductchooser(const blade b1, const blade b2){return true;}
bool outerproductchooser    (const blade b1, const blade b2){return ((  b1 &  b2).count() == 0)                                                                  ;}
bool innerproductchooser    (const blade b1, const blade b2){return ((((b1 & ~b2).count() == 0) | ((~b1 & b2).count() == 0)) && (b1.count()>0) && (b2.count()>0));}
bool fatdotchooser          (const blade b1, const blade b2){return ((( b1 & ~b2).count() == 0) | ((~b1 & b2).count() == 0))                                     ;}
bool lefttickchooser        (const blade b1, const blade b2){return ( ( b1 & ~b2).count() == 0)                                                                  ;}
bool righttickchooser       (const blade b1, const blade b2){return ( (~b1 &  b2).count() == 0)                                                                  ;}

clifford c_geometricprod(const clifford C1, const clifford C2, const NumericVector &signature){ return c_general_prod(C1, C2, signature, &geometricproductchooser);}
clifford outerprod      (const clifford C1, const clifford C2, const NumericVector &signature){ return c_general_prod(C1, C2, signature, &outerproductchooser    );}
clifford innerprod      (const clifford C1, const clifford C2, const NumericVector &signature){ return c_general_prod(C1, C2, signature, &innerproductchooser    );}
clifford fatdotprod     (const clifford C1, const clifford C2, const NumericVector &signature){ return c_general_prod(C1, C2, signature, &fatdotchooser          );}
clifford lefttickprod   (const clifford C1, const clifford C2, const NumericVector &signature){ return c_general_prod(C1, C2, signature, &lefttickchooser        );}
clifford righttickprod  (const clifford C1, const clifford C2, const NumericVector &signature){ return c_general_prod(C1, C2, signature, &righttickchooser       );}

clifford overwrite(clifford C1, const clifford C2){  // C1[] <- C2
    clifford::const_iterator i;
    for(i=C2.begin(); i != C2.end(); ++i){
        C1[i->first] = i->second;
    }
    return C1;
}

clifford c_power(const clifford C, const NumericVector &power, const NumericVector &signature){  // p for power
    clifford out;
    unsigned int p = power[0];

    if(p<1){throw std::range_error("power cannot be <1");} 
    if(p==1){
        return C;
    } else {
        out = C; 
        for( ; p>1; p--){
            out = c_geometricprod(C,out, signature);
        }
    }
    return out;
}
