#ifndef H_POINTLS
#define H_POINTLS

#include <vector>
#include <sstream>
#include <Rcpp.h>
#include "Krigtypes.h"

class Cpointls{
      std::string pointlsFileName;
      std::vector<std::vector<ioType> > xy;
public:
      Cpointls() {};
      Cpointls(Cpointls& ptls); //copy constructor - does not copy pointlsFileName
      Cpointls(double* xy,int* nrow,int* ncol);
      virtual ~Cpointls() {};
      int read_pointls(std::string filename);
      int select_columns(int fittedp);
      int minuslogLTologL();
      int selectTop();
      std::vector<std::vector<ioType> > getxy() {return xy;};
      std::string getname() {return pointlsFileName;};
};

#endif //H_POINTLS

