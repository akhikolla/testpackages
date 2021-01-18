#ifndef _haptools_COMMON_H
#define _haptools_COMMON_H

#include <vector>
#include <Rcpp.h>
#include <cassert>

struct haplotype {
  std::vector<int> profile;

  int count;

  friend std::ostream & operator<<(std::ostream &o, haplotype const& other) {
    std::ostringstream oss;

    if (!other.profile.empty())
    {
      std::copy(other.profile.begin(), other.profile.end() - 1, 
        std::ostream_iterator<int>(oss, ","));
      oss << other.profile.back();
      return o << "(" << oss.str() << ") x " << other.count;
    }

    return o;
  }
    
  bool operator==(const haplotype &other) const { 
    //return (profile.size() == other.profile.size() && profile == other.profile);
    return profile == other.profile;
  }

  bool operator () (const haplotype* h1, const haplotype* h2) const { 
    return *h1 == *h2;
  }
};

/*
  R:
  
  hash <- function(h) {
    hash_val <- 0
    for (i in seq_along(h)) {
      hash_val <- hash_val*31 + ifelse(h[i] < 0, -2*h[i], 2*h[i]+1)
    }
    return(hash_val)
  }
  
  hash(c(0, 0, 0))

  hash(c(0, 0, 1))
  hash(c(0, 1, 1))
  
  hash(c(0, 1, 0))
  hash(c(1, 0, 0))
  
  hash(c(1, 1, 1))
  
  max(table(apply(expand.grid((-5:5), (-5:5), (-5:5)), 1, hash)))
  #max(table(apply(expand.grid((-5:5), (-5:5), (-5:5), (-5:5)), 1, hash)))
  #max(table(apply(expand.grid((-3:3), (-3:3), (-3:3), (-3:3), (-3:3)), 1, hash)))
  
  max(table(apply(expand.grid((-10:10), (-10:10), (-10:10)), 1, hash)))
  max(table(apply(expand.grid((-20:20), (-20:20)), 1, hash)))
  max(table(apply(expand.grid((-50:50), (-50:50)), 1, hash)))
*/
struct haplotype_hash {
  /*
    Negative allele mapped to even numbers, positive to odd numbers
    Assume rather small range of alleles
  */
  
  size_t operator () (const haplotype &h) const { 
    size_t hash = 0;

    for (std::vector<int>::const_iterator it = h.profile.begin(); 
      it != h.profile.end(); ++it) {

      /*
        Negative allele mapped to even numbers, positive to odd numbers
        Assume rather small range of alleles
      */
      hash = hash*31 + (*it < 0 ? -2*(*it) : 2*(*it)+1);
    }
    
    return hash;
  }
};

#endif

