// lin_polylog_framed_clust.cpp

#include <Rcpp.h>
#include "OptFramedClust.h"
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List lin_polylog_framed_clust(
    std::vector<double> & X,
    int K, int frame_width,
    int first_frame, int last_frame,
    int prev_k_f, int next_k_f)
{
  struct clustering cluster = MFC(X, frame_width,  K, first_frame,
                                  last_frame, prev_k_f,  next_k_f);

  Rcpp::List result;

  result["ID"] = cluster.Frame_ID;

  result["Border"] = cluster.Borders;

  result["centers"] = cluster.centers;

  result["withinss"] = cluster.withinss;

  result["size"] = cluster.size;

  result["totss"] = cluster.totss;

  result["tot.withinss"] = cluster.ssq  ;

  result["betweenss"] = cluster.totss - cluster.ssq;

  /*
   for(int i=0;i<frame_width;i++)
   {
   cout <<endl;
   cout << i << endl;
   cout<< endl;
   for(int j=0;j<K;j++)
   {
   cout << Cluster_Border[i][j] << endl;
   }
   cout << endl;}

   */

  return(result);
}
