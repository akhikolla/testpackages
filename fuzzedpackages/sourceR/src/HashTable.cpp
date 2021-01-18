#include <Rcpp.h>
#include <iostream>
#include <string>
#include <map>

using namespace Rcpp;
//[[Rcpp::plugins(cpp11)]]

class HashTable {
public:
  HashTable(SEXP key, SEXP val)
  {
    CharacterVector key_(key);
    NumericVector val_(val);
    insert(key_, val_);
  }

  void insert(CharacterVector key, NumericVector val)
  {
    if (key.size() != val.size())
      throw Rcpp::exception("length(keys) must match length(val)");
    CharacterVector::iterator k;
    NumericVector::iterator v;
    for(k=key.begin(), v=val.begin();
        k!=key.end(); ++k, ++v) {
      std::string key_str = std::string(*k);
      ht_[key_str] = *v;
    }
  }

  NumericVector find(CharacterVector keys) const
  {
    NumericVector rv;
    for (auto k = keys.begin(); k != keys.end(); ++k)
    {
      auto val = ht_.find(std::string(*k));
      if (val != ht_.end())
        rv.push_back(val->second);
      else
        rv.push_back(NA_REAL);
    }
    return rv;
  }

  void erase(CharacterVector keys)
  {
    for (auto k = keys.begin(); k != keys.end(); ++k)
    {
      auto it = ht_.find(std::string(*k));
      if (it != ht_.end())
        ht_.erase(it);
      else
        throw Rcpp::exception("no such key in hashtable");
    }
  }

  NumericVector values() const
  {
    NumericVector vals;
    for (auto it = ht_.begin(); it!=ht_.end(); ++it)
      vals.push_back(it->second);
    return vals;
  }

  CharacterVector keys() const
  {
    CharacterVector keys;
    for(auto it = ht_.begin(); it!=ht_.end(); ++it)
      keys.push_back(it->first);
    return keys;
  }

  void print() const
  {
    for(auto it = ht_.begin(); it!=ht_.end(); ++it)
      Rcout << it->first << ": " << it->second << "\n";
    Rcout << std::flush;
  }

private:
  typedef std::map<std::string, double> HashTable_t;
  HashTable_t ht_;
};


/* hashtable is a lightweight wrapper around C++ std::map
 * for creating a hash table of character => numeric pairs.
 */
RCPP_MODULE(HashTable_module)
{
  class_<HashTable>("HashTable")

  .constructor<SEXP, SEXP>()

  .method("find", &HashTable::find)
  .method("insert", &HashTable::insert)
  .method("erase", &HashTable::erase)
  .method("keys", &HashTable::keys)
  .method("values", &HashTable::values)
  .method("print", &HashTable::print)
  ;
}
