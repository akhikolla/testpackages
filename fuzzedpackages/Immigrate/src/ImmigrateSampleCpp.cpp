
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List ImmigrateSampleCpp(Function onesampleImmigrate,NumericMatrix train_xx,NumericVector train_yy,
       NumericVector sample_wt,NumericMatrix W, double epsilon = 0.01,
            double sig = 1,int max_iter = 10,bool removesmall = false){
            double p = train_xx.ncol();
            double c0 = 0;
            double c_before = c0;
            double c_after = c0+1;
            NumericMatrix w_before = W;
            NumericMatrix w_after = W;
            NumericMatrix w;
            int iter = 0;
            
            if (removesmall == false){
            while ( (fabs(c_before - c_after)>epsilon)&&(iter< max_iter)  ){
            w_before = w_after;
            c_before = c_after;
            List res = onesampleImmigrate(train_xx,train_yy,w_before,sample_wt,sig);
            NumericMatrix tmp = res["w"];
            w_after = tmp;
            c_after = res["C"];
            iter++;
            }
            }else {
            while ( (fabs(c_before - c_after)>epsilon)&&(iter< max_iter)  ){
            w_before = w_after;
            c_before = c_after;
            List res = onesampleImmigrate(train_xx,train_yy,w_before,sample_wt,sig);
            NumericMatrix tmp = res["w"];
            w_after = tmp;
            double totalsum = 0;
            for ( int i = 0; i <w_after.ncol();i++){
            for ( int j  = w_after.nrow()-1; j >=i; j-- ){
            if (w_after(i,j) < 1/p ){
            w_after(i,j) = 0;
            w_after(j,i) = 0;
            }else{
              if (i == j){totalsum += w_after(i,j)*w_after(i,j);
              }else{
                  totalsum += 2*w_after(i,j)*w_after(i,j);
                }
            }
            }
            }

            w_after = w_after/sqrt(totalsum);
            c_after = res["C"];
            iter++;
            }
            }
            List newList;
            newList["w"]=w_after;
            newList["c"]=c_after;
            newList["iter_num"]=iter;
            return newList;
            }
