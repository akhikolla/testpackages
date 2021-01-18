#include "STCXEMspatial.h"
#include "STCXEMnonspatial.h"
 
//[[Rcpp::export]]
S4  SpaTimeClusCpp(S4 input, List inputparam, NumericMatrix matT){
  S4 * input_p=&input;
  int spa = as<S4>(input_p->slot("model")).slot("spatial");
  if (spa==1){
    STCXEMspatial xem(input, inputparam, matT);
    xem.Run();
    xem.Output(input_p);
  }else{
    STCXEMnonspatial xem(input, inputparam, matT);
    xem.Run();
    xem.Output(input_p);
  }
  return input;
}     
