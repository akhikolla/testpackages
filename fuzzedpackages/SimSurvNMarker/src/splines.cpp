// [[Rcpp::depends(RcppArmadillo)]]
// #define DO_CHECKS
#define ARMA_NO_DEBUG

#include "splines.h"
#include <algorithm> // lower_bound
#include <cmath> // isnan
#include <stdexcept> // invalid_argument

inline void check_splines
  (const arma::vec &boundary_knots, const arma::vec &interior_knots,
   const int order) {
#ifdef DO_CHECKS
  if(order<1, 0)
    throw std::invalid_argument("order<1");
  if(boundary_knots.size() != 2L)
    throw std::invalid_argument("boundary_knots should have length 2");
  if(interior_knots.size()>0 && boundary_knots(0)>min(interior_knots))
    throw std::invalid_argument("boundary_knots(0)>min(interior_knots)");
  if(interior_knots.size()>0 && boundary_knots(1)<max(interior_knots))
    throw std::invalid_argument("boundary_knots(1)<max(interior_knots)");
  // TODO: check if interior_knots are in ascending order?
#endif
}

inline void throw_invalid_out(
    std::string const &cl, unsigned const dim, unsigned const dim_ex){
#ifdef DO_CHECKS
  std::stringstream msg;
  msg << cl << ": invalid 'out' (dim is " << dim << "; expected "
      << dim_ex << ')';
  throw std::invalid_argument(msg.str());
#endif
}

namespace splines {

vec basisMixin::operator()(double const x, int const ders) const {
  vec out(get_n_basis());
  operator()(out, x, ders);
  return out;
}

mat basisMixin::basis(const vec &x, const int ders,
                      const double centre) const {
#ifdef DO_CHECKS
  if (ders < 0)
    throw std::invalid_argument("ders<0");
#endif
  uword const n_basis(get_n_basis()),
              n_x    (x.n_elem);
  rowvec centering =
    (std::isnan(centre) || ders > 0 ?
     zeros(n_basis) : operator()(centre, 0)).t();

  mat out(n_x, n_basis);
  vec wrk(n_basis);
  for (uword i = 0; i < n_x; i++){
    operator()(wrk, x[i], ders);
    out.row(i) = wrk.t() - centering;
  }

  return out;
}

SplineBasis::SplineBasis(const int order): order(order), knots() {
#ifdef DO_CHECKS
  if (order<1)
    throw std::invalid_argument("order<1");
#endif
}


SplineBasis::SplineBasis(const vec knots, const int order):
  order(order), knots(knots) {
#ifdef DO_CHECKS
  if (order<1)
    throw std::invalid_argument("order<1");
#endif
}

void SplineBasis::operator()(
    vec &out, double const x, const int ders) const {
  out.zeros();
#ifdef DO_CHECKS
  if(out.n_elem != SplineBasis::get_n_basis())
    throw_invalid_out(
      "splineBasis", out.n_elem, SplineBasis::get_n_basis());
#endif

  set_cursor(x);
  int io = curs - order;
  if (io < 0 || io > nknots) {
    /* Do nothing. x is already zero by default
    for (size_t j = 0; j < (size_t)order; j++) {
      out(j+io) = double(0); // R_NaN;
    }*/
  } else if (ders > 0) { /* slow method for derivatives */
    for(uword i = 0; i < (size_t)order; i++) {
      for(uword j = 0; j < (size_t)order; j++)
        a(j) = 0;
      a(i) = 1;
      out(i + io) = slow_evaluate(x, ders);
    }
  } else { /* fast method for value */
    basis_funcs(wrk, x);
    for (uword i = 0; i < wrk.n_elem; i++)
      out(i + io) = wrk(i);
  }
}

int SplineBasis::set_cursor(const double x) const {
  /* don't assume x's are sorted */
  curs = -1; /* Wall */
  boundary = 0;
  for (int i = 0; i < nknots; i++) {
    if (knots(i) >= x)
      curs = i;
    if (knots(i) > x)
      break;
  }
  if (curs > ncoef) {
    int const lastLegit = ncoef;
    if (x == knots(lastLegit)){
      boundary = 1;
      curs = lastLegit;
    }
  }
  return curs;
}

void SplineBasis::diff_table(const double x, const int ndiff) const {
  for (int i = 0; i < ndiff; i++) {
    rdel(i) = knots(curs + i) - x;
    ldel(i) = x - knots(curs - (i + 1));
  }
}

double SplineBasis::slow_evaluate(const double x, int nder) const
{
  int ti = curs,
     lpt, apt, rpt, inner,
   outer = ordm1;
  if (boundary && nder == ordm1) /* value is arbitrary */
    return 0;
  while(nder--) {  // FIXME: divides by zero
    for(inner = outer, apt = 0, lpt = ti - outer; inner--; apt++, lpt++)
      a(apt) = (double)outer * (a(apt + 1) - a(apt)) /
        (knots(lpt + outer) - knots(lpt));
    outer--;
  }
  diff_table(x, outer);
  while(outer--)
    for(apt = 0, lpt = outer, rpt = 0, inner = outer + 1;
        inner--; lpt--, rpt++, apt++)
      // FIXME: divides by zero
      a(apt) = (a(apt + 1) * ldel(lpt) + a(apt) * rdel(rpt)) /
        (rdel(rpt) + ldel(lpt));
  return a(0);
}

void SplineBasis::basis_funcs(vec &b, const double x) const {
  diff_table(x, ordm1);
  b(0) = 1;
  for (uword j = 1; j <= (uword)ordm1; j++) {
    double saved(0);
    for (size_t r = 0; r < j; r++) { // do not divide by zero
      double const den = rdel(r) + ldel(j - 1 - r);
      if(den != 0) {
        double const term = b(r)/den;
        b(r) = saved + rdel(r) * term;
        saved = ldel(j - 1 - r) * term;
      } else {
        if(r != 0 || rdel(r) != 0)
          b(r) = saved;
        saved = double(0);
      }
    }
    b(j) = saved;
  }
}

inline arma::vec get_SplineBasis_knots
(vec const &boundary_knots, vec const &interior_knots,
 int const order) {
  check_splines(boundary_knots, interior_knots, order);

  uword const nknots = interior_knots.size() + 2 * order;
  vec knots(nknots);
  for(uword i = 0; i < (uword)order; i++) {
    knots(i) = boundary_knots(0);
    knots(nknots - i - 1) = boundary_knots(1);
  }
  if (interior_knots.size() > 0)
    for(uword i = 0; i < interior_knots.size(); i++)
      knots(i + order) = interior_knots(i);

  return knots;
}

bs::bs(const vec &bk, const vec &ik, const bool inter, const int ord):
  SplineBasis(get_SplineBasis_knots(bk, ik, ord), ord),
  boundary_knots(bk), interior_knots(ik),
  intercept(inter),
  df((int)intercept + order - 1 + interior_knots.size()) {
  check_splines(boundary_knots, interior_knots, order);
}

void bs::operator()(vec &out, double const x, const int ders) const {
#ifdef DO_CHECKS
  if(out.n_elem != bs::get_n_basis())
    throw_invalid_out("bs", out.n_elem, bs::get_n_basis());
#endif
  if (x < boundary_knots(0) || x > boundary_knots(1)) {
    double const k_pivot =
        x < boundary_knots(0) ?
        0.75 * boundary_knots(0) + 0.25 * knots(order) :
        0.75 * boundary_knots(1) + 0.25 * knots(knots.n_elem - order - 2),
      delta = x - k_pivot;

    auto add_term = [&](int const d, double const f = 1){
      bs::operator()(wrks, k_pivot, d);
      out += f * wrks;
    };

    out.zeros();
    if (ders == 0) {
      add_term(0);
      add_term(1, delta);
      add_term(2, delta * delta/2.);
      add_term(3, delta * delta * delta /6.);

    } else if (ders == 1) {
      add_term(1);
      add_term(2, delta);
      add_term(3, delta * delta / 2.);

    } else if (ders == 2) {
      add_term(2);
      add_term(3, delta);

    } else if (ders == 3)
      add_term(3);

    return;
  }

  if(intercept)
    SplineBasis::operator()(out, x, ders);
  else {
    SplineBasis::operator()(wrk, x, ders);
    for(uword i = 1; i < wrk.n_elem; ++i)
      out[i - 1L] = wrk[i];
  }
}

ns::ns(const vec &boundary_knots, const vec &interior_knots,
       const bool intercept, const int order):
  bspline(boundary_knots, interior_knots, true, order),
  intercept(intercept),
  tl0(trans(bspline(boundary_knots(0), 0))),
  tl1(trans(bspline(boundary_knots(0), 1))),
  tr0(trans(bspline(boundary_knots(1), 0))),
  tr1(trans(bspline(boundary_knots(1), 1)))
  { }

void ns::operator()(vec &out, double const x, const int ders) const {
#ifdef DO_CHECKS
  if(out.n_elem != ns::get_n_basis())
    throw_invalid_out("ns", out.n_elem, ns::get_n_basis());
#endif
  if(x < bspline.boundary_knots(0)) {
    if (ders==0){
      out  = tl1;
      out *= x - bspline.boundary_knots(0);
      out += tl0;

    } else if (ders == 1)
      out = tl1;
    else
      out.zeros();

    return;

  } else if (x > bspline.boundary_knots(1)) {
    if (ders==0){
      out  = tr1;
      out *= x - bspline.boundary_knots(1);
      out += tr0;

    } else if (ders==1)
      out = tr1;
    else
      out.zeros();

    return;
  }

  out = trans(bspline(x, ders));
}

vec ns::trans(const vec &x) const {
  vec out = q_matrix * (intercept ? x : x(span(1, x.n_elem - 1)));
  return out(span(2, out.size() - 1));
}

iSpline::iSpline(const vec &boundary_knots, const vec &interior_knots,
                 const bool intercept, const int order):
  intercept(intercept), order(order),
  bspline(boundary_knots, interior_knots, false, order + 1) { }


void iSpline::operator()(vec &out, double const x, const int der) const {
#ifdef DO_CHECKS
  if(out.n_elem != iSpline::get_n_basis())
    throw_invalid_out("iSpline", out.n_elem, iSpline::get_n_basis());
#endif
  if(x < 0){
    out.zeros();
    return;

  }
  else if(x <= 1){
    vec &b = wrk;
    bspline(b, x, der);
    int const js = (bspline.interior_knots.size()>0) ?
      std::lower_bound(
        bspline.knots.begin(),
        /* TODO: should this not be end and not -1? */
        bspline.knots.end() - 1L, x) -
        bspline.knots.begin() : order + 1;
    for(uword j = b.size(); j-- >0;)
      if ((int)j > js)
        b[j] = 0.0;
      else if(j!=b.size()-1)
        b[j] += b[j+1];
    if (der==0)
      for(uword j = b.size() - 1; j-- > 0;)
        if ((int)j < js - (order+1))
          b[j] = 1.0;

    if(intercept)
      out = b;
    else
      out = b.subvec(1, b.n_elem - 1);
    return;

  }
  else if(der > 0)
    out.zeros();
  else
    out.fill(1);
}


mSpline::mSpline(const vec &boundary_knots, const vec &interior_knots,
                 const bool intercept, const int order) :
  bspline(boundary_knots, interior_knots, true, order),
  intercept(intercept) { }

void mSpline::operator()(vec &out, double const x, const int der) const {
#ifdef DO_CHECKS
  if(out.n_elem != mSpline::get_n_basis())
    throw_invalid_out("mSpline", out.n_elem, mSpline::get_n_basis());
#endif
  bspline(wrk, x, der);
  for (uword j = 0; j < (uword)bspline.get_n_basis(); j++) {
    double denom = bspline.knots(j + bspline.order) - bspline.knots(j);
    wrk(j) *= denom > 0.0 ? bspline.order / denom : 0.0;
  }

  if(intercept)
    out = wrk;
  else
    out = wrk.subvec(1, wrk.size() - 1);
}

} // namespace splines


// [[Rcpp::export(rng = false)]]
Rcpp::XPtr<splines::ns>
  get_ns_ptr(const arma::vec &knots, arma::vec const &boundary_knots,
             bool const intercept){
    if(boundary_knots.n_elem < 2)
      throw std::invalid_argument("get_ns_ptr: invalid boundary_knots size");

    return Rcpp::XPtr<splines::ns>(
      new splines::ns(boundary_knots, knots, intercept), true);
  }

// [[Rcpp::export(rng = false)]]
arma::mat ns_cpp(arma::vec const &x, SEXP ns_ptr){
  Rcpp::XPtr<splines::ns> basis(ns_ptr);
  size_t const n_x = x.size(),
             n_col = basis->get_n_basis();

  arma::mat out(n_x, n_col);
  arma::vec wrk(n_col);

  for(unsigned i = 0; i < n_x; ++i){
    basis->operator()(wrk, x[i]);
    out.row(i) = wrk.t();
  }

  return out;
}
