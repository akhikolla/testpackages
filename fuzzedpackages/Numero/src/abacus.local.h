/* Created by Ville-Petteri Makinen
   email: ville.makinen@vipmak.net */

#ifndef abacus_local_INCLUDED
#define abacus_local_INCLUDED

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <climits>
#include <cmath>
#include <algorithm>
#include <random>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include "medusa.h"
#include "abacus.h"

using namespace std;
using namespace medusa;
using namespace abacus;

/* Encapsulate with redundant namespace in case in a collection
   of modules another module has the same class name(s) in use. */
namespace abacus_local {

  /*
   *
   */
  class BaseGaussian {
  protected:
    string method;
    mdreal center;
    mdreal offset;
    mdreal scale;
    mdreal factor;
    mdreal mu;
    mdreal sigma;
  public:
    BaseGaussian() {
      this->center = medusa::rnan();
      this->offset = center;
      this->scale = center;
      this->factor = center;
      this->mu = center;
      this->sigma = center;
    };
    ~BaseGaussian() {};
    vector<mdreal> parameters() const;
    void parameters(const vector<mdreal>&);
    mdsize transform(vector<mdreal>&) const;
    void apply(vector<mdreal>&, const mdreal) const;
  };

  /*
   *
   */
  class Gaussian : public BaseGaussian {
  private:
    vector<mdsize> qloci;
    vector<mdreal> values;
    vector<mdreal> weights;
    vector<mdreal> zscores;
  public:
    Gaussian() : BaseGaussian() {};
    Gaussian(const vector<mdreal>& prm) : BaseGaussian() {
      this->parameters(prm);
    };
    ~Gaussian() {};
    bool configure(const vector<mdreal>&, const vector<mdreal>&);
    string model() const {return method;};
    mdreal optimize(const string&);
    mdreal quality() const;
    mdreal distance(const mdreal, const mdreal, const mdreal) const;
  };
  
  /*
   *
   */
  class Approximation {
  private:
    mdreal mode;
    Gaussian positive;
    Gaussian negative;
  public:
    Approximation() {this->mode = medusa::rnan();};
    Approximation(void* ptr) {
      Approximation* p = (Approximation*)ptr;
      this->mode = p->mode;
      this->positive = p->positive;
      this->negative = p->negative;
    };
    ~Approximation() {};
    void fit(const vector<mdreal>&, const vector<mdreal>&);
    vector<mdreal> parameters() const;
    bool parameters(const vector<mdreal>&);
    mdreal transform(const mdreal) const;
  };
    
  /*
   *
   */
  class EmpiricalBuffer {
  public:
    unsigned long ndata;
    Approximation approx;
    vector<mdreal> valsorted;
    vector<mdreal> wsorted;    
    unordered_map<mdreal, mdreal> data;
  public:
    EmpiricalBuffer() {this->ndata = 0;};
    EmpiricalBuffer(void* ptr) {
      EmpiricalBuffer* p = (EmpiricalBuffer*)ptr;
      this->ndata = p->ndata;
      this->approx = p->approx;
      this->valsorted = p->valsorted;
      this->wsorted = p->wsorted;
      this->data = p->data;
    };
    ~EmpiricalBuffer() {};
    void contents(vector<mdreal>&, vector<mdreal>&);
  };
  
  /*
   *
   */
  class Array {
  private:
    mdsize ndata;
    mdsize nelem;
    mdreal rlnan;
    vector<mdreal> full;
    map<mdsize, mdreal> sparse;
  private:
    mdsize optimize();
  public:
    Array();
    ~Array();
    mdreal operator[](const mdsize) const;
    void elements(vector<Element>&, const mdsize) const;
    mdsize length();
    mdreal remove(const mdsize);
    mdsize size() const;
    bool update(const mdsize, const mdreal, const bool);
    vector<mdreal> values() const;
  };

  /*
   *
   */
  class TwowayMap {
  private:
    unordered_map<mdsize, string> rank2name;
    unordered_map<string, mdsize> name2rank;
  public:
    void erase(const mdsize);
    void erase(const string&);
    void insert(const mdsize, const string&);
    string name(const mdsize);
    mdsize rank(const string&);
  };

  /*
   *
   */
  class MatrixBuffer  {
  public:
    bool symmflag;
    mdsize nrows;
    mdsize ncols;
    mdreal rlnan;
    TwowayMap rownames;
    TwowayMap colnames;
    unordered_map<mdsize, Array> rowdata;
  public:
    MatrixBuffer() {
      this->symmflag = false;
      this->nrows = 0;
      this->ncols = 0;
      this->rlnan = medusa::rnan();
    };
    MatrixBuffer(void* ptr) {
      MatrixBuffer* p = (MatrixBuffer*)ptr;
      this->symmflag = p->symmflag;
      this->nrows = p->nrows;
      this->ncols = p->ncols;
      this->rlnan = p->rlnan;
      this->rownames = p->rownames;
      this->colnames = p->colnames;
      this->rowdata = p->rowdata;
    };
    ~MatrixBuffer() {};
    vector<Element> elements(const int, const bool);
  };

  /*
   *
   */
  class MinimizerBuffer {
  public:
    mdsize npoints;
    mdreal epsilon;
    pair<mdreal, mdreal> limits;
  public:
    MinimizerBuffer() {
      this->npoints = 0;
      this->epsilon = 0.0;
    };
    MinimizerBuffer(void* ptr) {
      MinimizerBuffer* p = (MinimizerBuffer*)ptr;
      this->npoints = p->npoints;
      this->epsilon = p->epsilon;
    };
    ~MinimizerBuffer() {};
  };
}

using namespace abacus_local;

#endif /* abacus_local_INCLUDED */
