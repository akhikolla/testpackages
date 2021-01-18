#define STRICT_R_HEADERS
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <Rcpp.h>
using namespace Rcpp;



// [[Rcpp::export]]
void rcpp_forcelayout(DataFrame schedule, String path){
  print(schedule);
  // StringMatrix names = internal::convert_using_rfunction(schedule.names(), "as.matrix");  
  StringMatrix scheduleMatrix = internal::convert_using_rfunction(schedule, "as.matrix");
  // StringMatrix names = internal::convert_using_rfunction(rownames(schedule), "as.data.frame");
  StringVector names = rownames(scheduleMatrix);
  print(names);
  String path_local(path);
  std::ofstream myfile;

  Rcout << path_local.get_cstring() << std::endl;  
  path_local += "/extdata/DataForceLayout.js";
  
  myfile.open(path_local.get_cstring());
  myfile << "var units = {" << std::endl;
  
  int nrows = names.size();
  
  for(int i=0; i < nrows; i++)
  {
    myfile << "\"" << i+1 << "\" : \"" << names[i] << "\"," << std::endl;
  }


  myfile << "};" << std::endl;


  myfile << "var data=[" << std::endl;

  for (int i=0; i < scheduleMatrix.ncol(); i++)
  {
    myfile << "{" << std::endl;
    for (int j=0; j<nrows;j++) {
      myfile << "\"" << j+1 << "\":\"" << scheduleMatrix(j,i) << "\","  << std::endl;
    }
    myfile << "}," << std::endl;
  }
  myfile.seekp(-3, std::ios_base::cur); /**< Backspacing to take out the last comma in the object array */
  
  myfile << "];";
  myfile.close();
  
  path_local += "/../indexForceLayout.html";
  if (system(path_local.get_cstring()) != 2)
  {
  } 
}






