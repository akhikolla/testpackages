
#include "RcppArmadillo.h"
#include "decline_header.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]


// ********** Prediction Interface **********

// [[Rcpp::export]]
NumericMatrix decline_predict_cpp(Rcpp::List decline_lst, Rcpp::List time_lst) {

   double m3_to_bbl = 6.289814;
   double m3_to_Mft3 = 35.314667 / 1000.;
   std::string input_unit = as<std::string>(decline_lst["input_unit"]);
   std::string output_unit = as<std::string>(decline_lst["output_unit"]);
   std::string fluid = as<std::string>(decline_lst["fluid"]);
   arma::vec time = as<arma::vec>(time_lst["t"]);
   std::string cls;
   if (time_lst.inherits("day")) {
      cls = "day";
   }
   if (time_lst.inherits("month")) {
      cls = "month";
   }
   if (time_lst.inherits("year")) {
      cls = "year";
   }
   int lr = time.size();
   int ll = decline_lst.length();
   int lc = 1;
   if (decline_lst.inherits("exponential") | decline_lst.inherits("harmonic") | decline_lst.inherits("hyperbolic")) {
      if (ll == 6) {
         lc = 5;
      } else {
         lc = 7;
      }
   }
   if (decline_lst.inherits("modified_hyperbolic")) {
      if (ll == 7) {
         lc = 5;
      } else {
         lc = 7;
      }
   }
   NumericMatrix results_table_(lr,lc);
   CharacterVector colname(lc);
   if (lc == 5) {
      if (output_unit == "Field") {
         if (fluid == "oil") {
            if (cls == "day") {
               colname = {"Time_(day)", "q_(bbl/day)", "Q_(bbl)", "D_(1/day)", "Beta"};
            }
            if (cls == "month") {
               colname = {"Time_(month)", "q_(bbl/month)", "Q_(bbl)", "D_(1/month)", "Beta"};
            }
            if (cls == "year") {
               colname = {"Time_(year)", "q_(bbl/year)", "Q_(bbl)", "D_(1/year)", "Beta"};
            }
         } else {
            if (cls == "day") {
               colname = {"Time_(day)", "q_(MSCF/day)", "Q_(MMSCF)", "D_(1/day)", "Beta"};
            }
            if (cls == "month") {
               colname = {"Time_(day)", "q_(MSCF/month)", "Q_(MMSCF)", "D_(1/month)", "Beta"};
            }
            if (cls == "year") {
               colname = {"Time_(day)", "q_(MSCF/year)", "Q_(MMSCF)", "D_(1/year)", "Beta"};
            }
         }
      }
      if (output_unit == "SI") {
         if (cls == "day") {
            colname = {"Time_(day)", "q_(m3/day)", "Q_(m3)", "D_(1/day)", "Beta"};
         }
         if (cls == "month") {
            colname = {"Time_(month)", "q_(m3/month)", "Q_(m3)", "D_(1/month)", "Beta"};
         }
         if (cls == "year") {
            colname = {"Time_(year)", "q_(m3/year)", "Q_(m3)", "D_(1/year)", "Beta"};
         }
      }
   } else {
      // lc == 7
      if (output_unit == "Field") {
         if (fluid == "oil") {
            if (cls == "day") {
               colname = {"Time_(day)", "q_(bbl/day)", "Q_(bbl)", "D_(1/day)", "Beta", "time_abnd_(days)", "EUR_(bbl)"};
            }
            if (cls == "month") {
               colname = {"Time_(month)", "q_(bbl/month)", "Q_(bbl)", "D_(1/month)", "Beta", "time_abnd_(months)", "EUR_(bbl)"};
            }
            if (cls == "year") {
               colname = {"Time_(year)", "q_(bbl/year)", "Q_(bbl)", "D_(1/year)", "Beta", "time_abnd_(years)", "EUR_(bbl)"};
            }
         } else {
            if (cls == "day") {
               colname = {"Time_(day)", "q_(MSCF/day)", "Q_(MMSCF)", "D_(1/day)", "Beta", "time_abnd_(days)", "EUR_(MMSCF)"};
            }
            if (cls == "month") {
               colname = {"Time_(day)", "q_(MSCF/month)", "Q_(MMSCF)", "D_(1/month)", "Beta", "time_abnd_(months)", "EUR_(MMSCF)"};
            }
            if (cls == "year") {
               colname = {"Time_(day)", "q_(MSCF/year)", "Q_(MMSCF)", "D_(1/year)", "Beta", "time_abnd_(years)", "EUR_(MMSCF)"};
            }
         }
      }
      if (output_unit == "SI") {
         if (cls == "day") {
            colname = {"Time_(day)", "q_(m3/day)", "Q_(m3)", "D_(1/day)", "Beta", "time_abnd_(days)", "EUR_(m3)"};
         }
         if (cls == "month") {
            colname = {"Time_(month)", "q_(m3/month)", "Q_(m3)", "D_(1/month)", "Beta", "time_abnd_(months)", "EUR_(m3)"};
         }
         if (cls == "year") {
            colname = {"Time_(year)", "q_(m3/year)", "Q_(m3)", "D_(1/year)", "Beta", "time_abnd_(years)", "EUR_(m3)"};
         }
      }
   }
   if(decline_lst.inherits("exponential")) {
      results_table_ = wrap(exponential(decline_lst, time));
   }
   if(decline_lst.inherits("harmonic")) {
      results_table_ = wrap(harmonic(decline_lst, time));
   }
   if(decline_lst.inherits("hyperbolic")) {
      results_table_ = wrap(hyperbolic(decline_lst, time));
   }
   if(decline_lst.inherits("modified_hyperbolic")) {
      results_table_ = wrap(modified_hyperbolic(decline_lst, time));
   }
   if (input_unit == output_unit) {
      if (output_unit == "Field") {
         if (fluid == "gas") {
            if (lc == 5) {
               for (int j = 0; j < lr; j++) {
                  results_table_(j,2) = results_table_(j,2) / 1000.;
               }
            }
            if (lc == 7) {
               for (int j = 0; j < lr; j++) {
                  results_table_(j,2) = results_table_(j,2) / 1000.;
                  results_table_(j,6) = results_table_(j,6) / 1000.;
               }
            }
         }
      }
   } else if (output_unit == "Field") {
      if (fluid == "oil") {
         if (lc == 5) {
            for (int j = 0; j < lr; j++) {
               results_table_(j,1) = results_table_(j,1) * m3_to_bbl;
               results_table_(j,2) = results_table_(j,2) * m3_to_bbl;
            }
         }
         if (lc == 7) {
            for (int j = 0; j < lr; j++) {
               results_table_(j,1) = results_table_(j,1) * m3_to_bbl;
               results_table_(j,2) = results_table_(j,2) * m3_to_bbl;
               results_table_(j,6) = results_table_(j,6) * m3_to_bbl;
            }
         }
      }
      if (fluid == "gas") {
         if (lc == 5) {
            for (int j = 0; j < lr; j++) {
               results_table_(j,1) = results_table_(j,1) * m3_to_Mft3;
               results_table_(j,2) = results_table_(j,2) * m3_to_Mft3 / 1000.;
            }
         }
         if (lc == 7) {
            for (int j = 0; j < lr; j++) {
               results_table_(j,1) = results_table_(j,1) * m3_to_Mft3;
               results_table_(j,2) = results_table_(j,2) * m3_to_Mft3 / 1000.;
               results_table_(j,6) = results_table_(j,6) * m3_to_Mft3 / 1000.;
            }
         }
      }
   } else {
      if (fluid == "oil") {
         if (lc == 5) {
            for (int j = 0; j < lr; j++) {
               results_table_(j,1) = results_table_(j,1) / m3_to_bbl;
               results_table_(j,2) = results_table_(j,2) / m3_to_bbl;
            }
         }
         if (lc == 7) {
            for (int j = 0; j < lr; j++) {
               results_table_(j,1) = results_table_(j,1) / m3_to_bbl;
               results_table_(j,2) = results_table_(j,2) / m3_to_bbl;
               results_table_(j,6) = results_table_(j,6) / m3_to_bbl;
            }
         }
      }
      if (fluid == "gas") {
         if (lc == 5) {
            for (int j = 0; j < lr; j++) {
               results_table_(j,1) = results_table_(j,1) / m3_to_Mft3;
               results_table_(j,2) = results_table_(j,2) / m3_to_Mft3;
            }
         }
         if (lc == 7) {
            for (int j = 0; j < lr; j++) {
               results_table_(j,1) = results_table_(j,1) / m3_to_Mft3;
               results_table_(j,2) = results_table_(j,2) / m3_to_Mft3;
               results_table_(j,6) = results_table_(j,6) / m3_to_Mft3;
            }
         }
      }
   }
   colnames(results_table_) = colname;
   return(results_table_);
}

