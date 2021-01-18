#include <ggdmc.hpp>

using namespace Rcpp;

Prior::Prior (List & pprior)
{
  using namespace arma;

  std::vector<std::string> pnames = pprior.attr("names");
  m_npar = pnames.size();

  vec p0(m_npar), p1(m_npar), l(m_npar), u(m_npar);
  uvec d(m_npar), lg(m_npar);

  for (size_t i = 0; i < m_npar; i++) {
    List a_list = pprior[pnames[i]];
    unsigned int a_dist = a_list.attr("dist");

    d[i]  = a_dist;
    p0[i] = a_list[0];
    p1[i] = a_list[1];
    l[i]  = a_list[2];
    u[i]  = a_list[3];
    lg[i] = a_list[4];
  }
  m_d  = d;
  m_p0 = p0;
  m_p1 = p1;
  m_l  = l;
  m_u  = u;
  m_lg = lg;
}

Prior::~Prior()
{
  // Rcout << "Prior destructor\n";
}

void Prior::dprior(double * pvector, double * out)
{
  double x, l, u;

  for (size_t i = 0; i < m_npar; i++)
  {
    // NA go here; NA will be converted to 0 (unsigned int type)
    if ( ISNAN(m_p1[i]) || ISNAN(m_d[i]) ) {
      out[i] = m_lg[i] ? R_NegInf : 0;
    } else if ( m_d[i] == TNORM ) {

      l = ISNAN(m_l[i]) ? R_NegInf : m_l[i];
      u = ISNAN(m_u[i]) ? R_PosInf : m_u[i];

      tnorm * obj = new tnorm(m_p0[i], m_p1[i], l, u, m_lg[i]);
      out[i] = obj->d(pvector[i]);
      delete obj;

    } else if ( m_d[i] == BETA_LU ) {

      l = ISNAN(m_l[i]) ? 0 : m_l[i]; // In case the user enters NAs.
      u = ISNAN(m_u[i]) ? 1 : m_u[i];

      x = (pvector[i] - l) / (u -  l);

      // Note m_l differs from m_lg !!!
      out[i] = !m_lg[i] ? R::dbeta(x, m_p0[i], m_p1[i], false) / (u - l) :
                          R::dbeta(x, m_p0[i], m_p1[i], true) - std::log(u - l);

    } else if ( m_d[i] == GAMMA_L ) {

      l = ISNAN(m_l[i]) ? 0 : m_l[i];
      x = ( !R_FINITE(l) ) ? pvector[i] : pvector[i] - l;
      out[i] = R::dgamma(x, m_p0[i], m_p1[i], m_lg[i]);

    } else if ( m_d[i] == LNORM_L ) {

      l = ISNAN(m_l[i]) ? 0 : m_l[i];
      x = ( !R_FINITE(l) ) ? pvector[i] : pvector[i] - l;
      out[i] = R::dlnorm(x, m_p0[i], m_p1[i], m_lg[i]);

    } else if ( m_d[i] == 5 ) {

      out[i] = R::dunif(pvector[i], m_p0[i], m_p1[i], m_lg[i]);

    } else if ( m_d[i] == CONSTANT ) {

      out[i] = m_lg[i] ? R_NegInf : 0;

    } else if ( m_d[i] == TNORM_TAU ) {
      l = ISNAN(m_l[i]) ? R_NegInf : m_l[i];
      u = ISNAN(m_u[i]) ? R_PosInf : m_u[i];

      tnorm * obj = new tnorm(m_p0[i], m_p1[i], l, u, m_lg[i]);
      out[i] = obj->d2(pvector[i]);
      delete obj;

    } else {
      Rcpp::Rcout << "Distribution type undefined" << "\n";
      out[i] = m_lg[i] ? R_NegInf : 0;
    }
  }

}

arma::vec Prior::dprior(arma::vec pvector)
{
  double * pvec = new double[m_npar];
  double * tmp  = new double[m_npar];

  for (size_t i = 0; i < m_npar; i++) pvec[i] = pvector[i];

  dprior(pvec, tmp);

  arma::vec out(m_npar);
  for (size_t i = 0; i < m_npar; i++)
  {
    // Rcout << "tmp[i] " << tmp[i] << " - " << R_FINITE(tmp[i]) << "\n";

    if ( !R_FINITE(tmp[i]) )
    {
      out[i] = m_lg[i] ? -23.02585 : 1e-10; // critical to hierarchical
    }
    else
    {
      out[i] = tmp[i];
    }
  }

  delete [] pvec;
  delete [] tmp;

  return(out);

}

arma::vec Prior::rprior()
// Used by ininitlise & R's rprior
{
  // replace DMC modified r-function; used in initialise.cpp internally
  double l, u;
  arma::vec out(m_npar); out.fill(NA_REAL);

  // [p1 p2]: [mean sd]; [shape1 shape2]; [shape scale]; [meanlog sdlog]
  for (size_t i = 0; i < m_npar;  i++) {
    if ( ISNAN(m_d[i]) ) {
      out[i] = NA_REAL;

    } else if ( m_d[i] == 1 ) {         // tnorm
      l = ISNAN(m_l[i]) ? R_NegInf : m_l[i];
      u = ISNAN(m_u[i]) ? R_PosInf : m_u[i];

      tnorm * obj = new tnorm(m_p0[i], m_p1[i], l, u);
      out[i] = obj->r();
      delete obj;

    } else if ( m_d[i] == 2 ) {  // beta_ul
      l = ISNAN(m_l[i]) ? 0 : m_l[i];
      u = ISNAN(m_u[i]) ? 1 : m_u[i];
      out[i] = l + R::rbeta(m_p0[i], m_p1[i]) * (u - l);


    } else if ( m_d[i] == 3 ) {  // gamma_l
      l = ISNAN(m_l[i]) ? 0 : m_l[i];
      out[i] = R::rgamma(m_p0[i], m_p1[i]) + l;

    } else if ( m_d[i] == 4 ) {  // lnorm_l
      l = ISNAN(m_l[i]) ? 0 : m_l[i];
      out[i] = R::rlnorm(m_p0[i], m_p1[i]) + l;

    } else if ( m_d[i] == 5 ) {
      out[i] = R::runif(m_p0[i], m_p1[i]);
    } else if ( m_d[i] == 6 ){  // constant
      out[i] = m_p0[i];

    } else if ( m_d[i] == 7 ) { // tnorm2

      Rcout << "Distribution type not supported\n";

      // l = ISNAN(lower[i]) ? R_NegInf : lower[i];
      // u = ISNAN(upper[i]) ? R_PosInf : upper[i];
      // out[i] = rtn_scalar2(p1[i], p2[i], l, u);

    } else {
      Rcout << "Distribution type not supported\n";
      out[i] = NA_REAL;
    }
  }

  return out;
}

double Prior::sumlogprior(arma::vec pvector)
{
  arma::vec out = dprior(pvector);
  // den.replace(arma::datum::inf, 1e-10);
  // out.replace(R_PosInf, 1e-10);

  return arma::accu(out);
}

void Prior::print(std::string str) const
{
  Rcpp::Rcout << str << ":\n";
  Rcpp::Rcout << "[Location, scale, lower, upper]:\n";

  for (size_t i=0; i<m_npar; i++)
  {
    Rcpp::Rcout << "[" << m_p0[i] << ", " <<
      m_p1[i] << ", " << m_l[i] << ", " << m_u[i] << "]" << std::endl;
  }
}

// [[Rcpp::export]]
NumericMatrix rprior_mat(List prior, unsigned int n) {

  if (n < 1) stop("n must be greater or equal to 1");

  Prior * obj = new Prior(prior);
  CharacterVector pnames = prior.attr("names");
  unsigned int npar = pnames.size();

  NumericMatrix out(n, npar);
  for (size_t i=0; i<n; i++)
  {
    arma::vec tmp = obj->rprior();

    for (size_t j=0; j<npar; j++)
    {
      out(i,j) = tmp[j];
    }
  }

  Rcpp::colnames(out) = pnames;
  return out;
}


// [[Rcpp::export]]
double test_sumlogprior(arma::vec pvec, List prior)
{
  Prior      * p0 = new Prior (prior);
  double out = p0->sumlogprior(pvec);;
  delete p0;
  return out;
}

// [[Rcpp::export]]
arma::vec test_dprior(arma::vec pvec, List prior)
{
  Prior      * p0 = new Prior (prior);
  arma::vec out = p0->dprior(pvec);
  delete p0;
  return out;
}
