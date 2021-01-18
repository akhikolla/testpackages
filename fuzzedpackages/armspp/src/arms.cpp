#include <Rcpp.h>
#include <random>
#include <armspp>

#include "utils.h"

using namespace Rcpp;
using armspp::ARMS;
using armspp::listToDottedPair;

class FunctionWrapper {
public:
  FunctionWrapper(Function f, DottedPair arguments)
    : f_(f),
      arguments_(arguments),
      nEvaluations_(0) {}

  double operator()(double x) {
    NumericVector output;
    if (arguments_ == R_NilValue) {
      output = f_(x);
    } else {
      output = Rcpp_eval(Rcpp_lcons(
        f_,
        grow(x, arguments_)
      ));
    }
    ++nEvaluations_;
    return output[0];
  }

  int nEvaluations() const { return nEvaluations_; }

private:
  Function f_;
  DottedPair arguments_;
  int nEvaluations_;
};

// [[Rcpp::export(name = '.arms')]]
RObject arms(
  int nSamples,
  List logPdf,
  NumericVector lower,
  NumericVector upper,
  List initial,
  NumericVector convex,
  IntegerVector maxPoints,
  IntegerVector metropolis,
  NumericVector previous,
  List arguments,
  bool includeNEvaluations
) {
  std::mt19937_64 rng(static_cast<uint_fast64_t>(UINT_FAST64_MAX * R::unif_rand()));

  int maxArgumentsSize = 1;
  for (int i = 0; i < arguments.size(); ++i){
    if (::Rf_isVector(arguments[i]) && !Rf_isMatrix(arguments[i])) {
      maxArgumentsSize = std::max(maxArgumentsSize, static_cast<int>(
        static_cast<GenericVector>(
          arguments[i]
        ).size()
      ));
    }
  }

  int nEvaluations = 0;
  NumericVector samples(nSamples);
  if (
    initial.size() == 1
    && lower.size() == 1
    && upper.size() == 1
    && logPdf.size() == 1
    && convex.size() == 1
    && maxPoints.size() == 1
    && metropolis.size() == 1
    && previous.size() == 1
    && maxArgumentsSize == 1
  ) {
    // Special case: only one distribution to sample from
    NumericVector initialVector(clone(as<NumericVector>(initial[0])));

    FunctionWrapper f(logPdf[0], listToDottedPair(arguments, 0));
    ARMS<double, FunctionWrapper, NumericVector::iterator> dist(
      f,
      lower[0], upper[0], convex[0],
      initialVector.begin(), initialVector.size(),
      maxPoints[0], metropolis[0], previous[0]
    );

    for (int i = 0; i < nSamples; ++i) {
      samples[i] = dist(rng);
    }
    nEvaluations += f.nEvaluations();
  } else {
    for (int i = 0; i < nSamples; ++i) {
      NumericVector initialVector(clone(as<NumericVector>(initial[i % initial.size()])));

      FunctionWrapper f(
        logPdf[i % logPdf.size()],
        listToDottedPair(arguments, i)
      );

      ARMS<double, FunctionWrapper, NumericVector::iterator> dist(
        f,
        lower[i % lower.size()],
        upper[i % upper.size()],
        convex[i % convex.size()],
        initialVector.begin(), initialVector.size(),
        maxPoints[i % maxPoints.size()],
        metropolis[i % metropolis.size()],
        previous[i % previous.size()]
      );
      samples[i] = dist(rng);
      nEvaluations += f.nEvaluations();
    }
  }

  if (includeNEvaluations) {
    List output;
    output["n_evaluations"] = nEvaluations;
    output["samples"] = samples;
    return output;
  } else {
    return samples;
  }
}
