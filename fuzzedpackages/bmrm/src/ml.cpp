#include <Rcpp.h>
using namespace Rcpp;

//' Find first common ancestor of 2 nodes in an hclust object
//' 
//' @param hc an hclust object
//' @param a an integer vector with the first leaf node
//' @param b an integer vector with the second leaf node (same length as a)
//' @return an integer vector of the same length as a and b identifing the first common ancestors of a and b
//' @author Julien Prados
//' @export
//' @useDynLib bmrm
//' @import Rcpp
//' @examples
//'   hc <- hclust(dist(USArrests), "complete")
//'   plot(hc)
//'   A <- outer(seq_along(hc$order),seq_along(hc$order),hclust_fca,hc=hc)
//'   H <- array(hc$height[A],dim(A))
//'   image(H[hc$order,hc$order])
//'   image(A[hc$order,hc$order])
//[[Rcpp::export]]
IntegerVector hclust_fca(List hc,IntegerVector a,IntegerVector b) {
  IntegerMatrix m = hc["merge"];
  IntegerVector o = hc["order"];
  IntegerVector parents_pos(m.nrow());
  IntegerVector parents_neg(o.length());
  for(int i=0;i<m.length();i++) {
    if (m[i]<0) {
      parents_neg[-m[i]-1] = i % m.nrow();
    } else {
      parents_pos[m[i]-1] = i % m.nrow();
    }
  }

  IntegerVector root(a.length());
  for(int i=0;i<root.length();i++) {
    int rootA = parents_neg[a[i]-1];
    int rootB = parents_neg[b[i]-1];
    while (rootA!=rootB) {
      if (rootA<rootB) {
        rootA = parents_pos[rootA];
      } else {
        rootB = parents_pos[rootB];
      }
    }
    root[i] = rootA + 1;
  }
  return root;
}



