#ifndef STUMP_H
#define STUMP_H

#include <Rcpp.h>
#include <vector>
using namespace Rcpp;

class Stump
{
public:
  Stump();
  Stump(List stump_in);

  // populate data for sboost
  static void populate_data(const NumericMatrix& f, const NumericVector& o, const NumericMatrix& oi, const NumericVector& c);
  // populate data for assess
  static void populate_data(const NumericMatrix& f, const NumericVector& o);
  // populate data for predict
  static void populate_data(const NumericMatrix& f);

  void find_stump(const NumericVector& weights);
  void set_vote(double v);

  double get_vote() const;
  void update_predictions(NumericVector& predictions) const;
  void new_predictions(NumericVector& predictions) const;
  void new_predictions_integer(NumericVector& predictions) const;
  NumericVector get_contingencies(const NumericVector& predictions) const;
  List make_list() const;

private:
  static NumericMatrix features;
  static NumericVector outcomes;
  static NumericMatrix ordered_index;
  static NumericVector categorical;

  int feature;
  int direction;
  double vote;
  int is_categorical;
  double split;
  std::vector<int> positive_categories;
  std::vector<int> negative_categories;
};


#endif
