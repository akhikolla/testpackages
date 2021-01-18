// [[Rcpp::interfaces(r, cpp)]]
#include <Rcpp.h>

using namespace Rcpp;

//' Get the Eddington number for cycling
//'
//' Gets the \href{https://en.wikipedia.org/wiki/Arthur_Eddington#Eddington_number_for_cycling}{Eddington number for cycling}.
//' The Eddington Number for cycling, \emph{E}, is the
//' maximum number where a cyclist has ridden \emph{E} miles in \emph{E} days.
//'
//' The Eddington Number for cycling is related to computing the rank of an
//' integer partition, which is the same as computing the side length of its
//' \href{https://en.wikipedia.org/wiki/Durfee_square}{Durfee square}. Another
//' relevant application of this metric is computing the
//' \href{https://doi.org/10.1073/pnas.0507655102}{Hirsch index} for
//' publications.
//'
//' This is not to be confused with the
//' \href{https://en.wikipedia.org/wiki/Eddington_number}{Eddington Number in
//' astrophysics}, \eqn{N_{Edd}}, which represents the number of protons in the
//' observable universe.
//'
//' @param rides A vector of mileage, where each element represents a single
//'   day.
//'
//' @return An integer which is the Eddington cycling number for the
//'   data provided.
//'
//' @seealso \code{\link{E_cum}}, \code{\link{E_next}}, \code{\link{E_req}},
//'   \code{\link{E_sat}}
//'
//' @examples
//' # Randomly generate a set of 15 rides
//' rides <- rgamma(15, shape = 2, scale = 10)
//'
//' # View the rides sorted in decreasing order
//' setNames(sort(rides, decreasing = TRUE), seq_along(rides))
//'
//' # Get the Eddington number
//' E_num(rides)
//' @export
// [[Rcpp::export]]
int E_num(NumericVector &rides) {
  int n = rides.size(), E = 0, ride = 0, above = 0;
  IntegerVector H(n);

  for (int i = 0; i < n; i++) {
    ride = (int) rides[i];
    if (ride > E) {
      above++;
      if (ride < n) H[ride]++;

      if (above > E) {
        E++;
        above -= H[E];
      }
    }
  }

  return E;
}

//' Calculate the cumulative Eddington number
//'
//' This function is much like \code{\link{E_num}} except it provides
//' a cumulative Eddington number over the vector rather than a single summary
//' number.
//'
//' @inheritParams E_num
//' @seealso \code{\link{E_next}}, \code{\link{E_num}}, \code{\link{E_req}},
//'   \code{\link{E_sat}}
//' @return An integer vector the same length as \code{rides}.
//' @export
// [[Rcpp::export]]
IntegerVector E_cum(NumericVector &rides) {
  int n = rides.size(), running = 0, ride = 0, above = 0;
  IntegerVector E(n), H(n);

  for (int i = 0; i < n; i++) {
    ride = (int) rides[i];
    if (ride > running) {
      above++;
      if (ride < n) H[ride]++;

      if (above > running) {
        running++;
        above -= H[running];
      }
    }

    E[i] = running;

  }

  return E;
}

//' Get the number of rides required to increment to the next Eddington number
//'
//' Get the number of rides required to increment to the next Eddington number.
//'
//' @inheritParams E_num
//' @seealso \code{\link{E_cum}}, \code{\link{E_num}}, \code{\link{E_req}},
//'   \code{\link{E_sat}}
//' @return A named list with the current Eddington number (\code{E}) and the
//'   number of rides required to increment by one (\code{req}).
//' @export
// [[Rcpp::export]]
List E_next(NumericVector &rides) {
  int n = rides.size(), E = 0, ride = 0, above = 0;
  IntegerVector H(n + 2);

  for (int i = 0; i < n; i++) {
    ride = (int) rides[i];
    if (ride > E) {
      above++;
      H[std::min(ride, n + 1)]++;

      if (above > E) {
        E++;
        above -= H[E];
      }
    }
  }

  List out = List::create(
    _["E"] = E,
    _["req"] = E + 1 - above
  );

  out.attr("class") = "E_next";

  return out;
}
