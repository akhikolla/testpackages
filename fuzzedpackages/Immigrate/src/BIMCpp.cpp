#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List BIMCpp(Function oneboostImmigrate,NumericMatrix train_xx,NumericVector train_yy,int nIter = 10, int max_iter = 10,bool removesmall = false, double sigstart=0.02,double sigend=4) {
  int num_examples = train_yy.size();
  List matrix_list = List::create();
  std::vector<double> coeff_vector(nIter,0.0);
  
  NumericVector sample_wt(num_examples);
  for ( int i = 0; i<num_examples;i++){
    sample_wt[i] = 1./num_examples;
  }
  double siggap = exp(log(sigend/sigstart)/10);
  NumericVector sigseq(nIter);
  double signow = sigstart/siggap;
  for (int i = 0; i < nIter; i++){
      signow = signow*siggap;
      sigseq[i] = signow;
        }
  double err;
  List boost_result;
  
  for(int i=0;i<nIter;i++)
  {
    boost_result = oneboostImmigrate(train_xx,train_yy,sample_wt,Named("sig") = sigseq[i], Named("removesmall") = removesmall,Named("max_iter") = max_iter);
    std::ostringstream ss;
    ss << i;
    std::string list_name =  ss.str();

    matrix_list[list_name] = boost_result["w"];
    err = boost_result["error"];

    coeff_vector[i] = boost_result["alpha"];
    if((err>=0.5)||(err==0)) 
    {
      coeff_vector.erase( coeff_vector.begin()+i+1,coeff_vector.end() );
      break;
    }
    sample_wt =  boost_result["sample_wt"]; 
  }

  List BIM_list;
  BIM_list["matrix"] = matrix_list;
  BIM_list["weights"] = coeff_vector;
  BIM_list["sample_wt"] = sample_wt;
  BIM_list["sig"] = sigseq;
  return BIM_list;
}
