//  Copyright (c) 2020 Chaim Even-Zohar
//
//  This file is part of the R package "independence".
//
//  "independence" is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  "independence" is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with "independence".  If not, see <https://www.gnu.org/licenses/>.


#include <Rcpp.h>
#include <vector>
#include <cmath>

// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::plugins(cpp11)]]

#define RCPP_CHECK_INTERRUPT_EVERY 16384

//--- data types --------------------------------------------------------------

// 64-bit integers are always available
typedef signed long long int64;
typedef unsigned long long uint64;

// 128-bit integers are sometimes available
#ifdef __SIZEOF_INT128__
__extension__ typedef signed __int128 int128;
__extension__ typedef unsigned __int128 uint128;
#endif

// max data size, to avoid overflows
#define MAX_TAUSTAR_64 (102569)
#define MAX_TAUSTAR_128 (6721987087ULL)
#define MAX_HOEFFDING_64 (14081)
#define MAX_HOEFFDING_128 (100413509)

// used for entries of a permutations, positions, length, etc.
// try: unsigned long -> unsigned long long if data exceeds 2^32.
// change exported functions accordingly
typedef unsigned long entry;

//' Maximum data size for the T* statistic
//'
//' Data at most this size should not cause an overflow in the computation
//' of Tau*. This may depend on the availability of 64-bit or 128-bit words
//' to the compiler in use.
//'
//' @return
//' 6721987087 in 128-bit mode,
//' 102569 in 64-bit mode.
//'
//' @seealso
//' \code{\link{tau.star.test}}
//'
//' @export max_taustar
//'
//' @examples
//' max_taustar()
// [[Rcpp::export(name = "max_taustar")]]
double max_taustar()
#ifdef __SIZEOF_INT128__
{ return MAX_TAUSTAR_128; }
#else
{ return MAX_TAUSTAR_64; }
#endif

//' Maximum data size for Hoeffding's statistics
//'
//' Data at most this size should not cause an overflow in the computation
//' of Hoeffding's D statistic or its refinement. This may depend on the
//' availability of 64-bit or 128-bit words to the compiler in use.
//'
//' @return
//' 100413509 in 128-bit mode,
//' 14081 in 64-bit mode.
//'
//' @seealso
//' \code{\link{hoeffding.D.test}}
//' \code{\link{hoeffding.refined.test}}
//'
//' @export max_hoeffding
//'
//' @examples
//' max_hoeffding()
// [[Rcpp::export(name = "max_hoeffding")]]
double max_hoeffding()
#ifdef __SIZEOF_INT128__
{ return MAX_HOEFFDING_128; }
#else
{ return MAX_HOEFFDING_64; }
#endif

//--- sumtree -----------------------------------------------------------------

// Array with O(log n) insertion and prefix sum
template <typename element>
class sumtree
{
  entry leng;
  std::vector<element> tree;

public:

  sumtree(entry length) : leng(length), tree(length+1,0) { }
  ~sumtree() { }

  inline void add(entry i, const element what)
  {
    tree[i] += what;
    while (i) tree[i-=i&-i] += what;
  }

  inline element sum_suffix(entry i)
  {
    element sum = 0;
    do {sum += tree[i];} while (i && ((i+=i&-i) < leng));
    return sum;
  }

  inline element sum_prefix(entry i)
  {
    element sum = tree[0];
    do {sum -= tree[i];} while (i && ((i+=i&-i) < leng));
    return sum;
  }
};

//--- taustar -----------------------------------------------------------------

// Subroutine for taustar computation
template <typename element>
element taustar_subroutine(const std::vector<entry>& perm)
{
  entry i, n = perm.size();
  element u, d, ud, du, udu, sum = 0;
  sumtree<element> S(n);
  sumtree<element> Su(n);
  sumtree<element> Sd(n);
  sumtree<element> Sud(n);

  for (i=0; i<n; ++i)
  {
    u = S.sum_prefix(perm[i]);
    d = S.sum_suffix(perm[i]+1);
    du = Sd.sum_prefix(perm[i]);
    ud = Su.sum_suffix(perm[i]+1);
    udu = Sud.sum_prefix(perm[i]);
    S.add(perm[i], 1);
    Su.add(perm[i], u);
    Sd.add(perm[i], d);
    Sud.add(perm[i], ud);
    sum += 2*udu-d*du-u*ud+u*d*(i-1);
    if (!(i%RCPP_CHECK_INTERRUPT_EVERY)) Rcpp::checkUserInterrupt();
  }

  return sum;
}

// Count 2*concordant-discordant quadrauples
template <typename element>
element taustar_count(const std::vector<entry>& perm)
{
  entry i, n = perm.size();
  std::vector<entry> mirror(n);
  element sum = 0;

  sum += taustar_subroutine<element>(perm);
  for (i=0; i<n; ++i) mirror[i] = perm[n-1-i];
  sum += taustar_subroutine<element>(mirror);
  for (i=0; i<n; ++i) mirror[perm[i]] = i;
  sum += taustar_subroutine<element>(mirror);
  for (i=0; i<n; ++i) mirror[n-1-perm[i]] = i;
  sum += taustar_subroutine<element>(mirror);

  element total = element(n)*(n-1)/2*(n-2)/3*(n-3)/4;
  element tstar = total*2 - (sum/4)*3;

  return tstar;
}

//' Compute the Tau* statistic
//'
//' This is an internal CPP function, used by the R function
//' \code{\link{tau.star.test}}.
//'
//' Given (X1,Y1),...,(Xn,Yn), the Tau*_n statistic only depends on the
//' permutation P that satisfies rank Yi = P[rank Xi].
//' This function computes Tau*_n given P in O(n log n) time.
//'
//' @param perm
//'     An integer vector containing exactly 0,1,...,n-1 in any order.
//'
//'     The validity of the input is not checked by this function.
//'
//' @return
//'     The Tau* statistic of \code{perm}.
//'
//'     The normalization is such that \emph{-1/3 <= Tau* <= 2/3}.
//'
//'     The return value -1.0 indicates an error.
//'
//' @aliases calc.taustar
//'
//' @export
//'
//' @examples
//'
//' .calc.taustar(0:3)
//' ## [1] 0.6666667
//'
//' .calc.taustar(c(0,2,1,3))
//' ## [1] -0.3333333
//'
//' set.seed(397)
//' .calc.taustar(order(runif(1000))-1)
//' ## [1] 0.004392385
//'
// [[Rcpp::export(name = ".calc.taustar")]]
double calc_taustar(const std::vector<unsigned long>& perm)
{
  entry n = perm.size();

  if ((n >= 4) && (n <= MAX_TAUSTAR_64))
  {
    int64 count = taustar_count<uint64>(perm);
    return count * 8.0 / n / (n-1) / (n-2) / (n-3);
  }

#ifdef __SIZEOF_INT128__
  if ((n >= 4) && (n <= MAX_TAUSTAR_128))
  {
    int128 count = taustar_count<uint128>(perm);
    return count * 8.0 / n / (n-1) / (n-2) / (n-3);
  }
#endif

  return -1.0;
}

//--- hoeffding ---------------------------------------------------------------

// Count 2*concordant-discordant quintuples
template <typename element>
element hoeffding_count(const std::vector<entry>& perm)
{
  entry n = perm.size();
  element a, b, c, sum = 0;
  sumtree<element> S(n);

  for (entry i = 0; i < n; ++i)
  {
    S.add(perm[i],1);
    a = i;
    b = perm[i];
    c = S.sum_prefix(perm[i]);
    sum += a*(a-1)/2*b*(b-1) - (a-1)*(b-1)*c*(n-2) + c*(c-1)/2*(n-2)*(n-3);
    if (!(i%RCPP_CHECK_INTERRUPT_EVERY)) Rcpp::checkUserInterrupt();
  }

  return sum;
}

//' Compute Hoeffding's D statistic
//'
//' This is an internal CPP function, used by the R function
//' \code{\link{hoeffding.D.test}}.
//'
//' Given (X1,Y1),...,(Xn,Yn), Hoeffding's Dn only depends on the
//' permutation P that satisfies rank Yi = P[rank Xi].
//' This function computes Dn given P in O(n log n) time.
//'
//' @param perm
//'     An integer vector containing exactly 0,1,...,n-1 in any order.
//'
//'     The validity of the input is not checked by this function.
//'
//' @return
//'     Hoeffding's D statistic of \code{perm}.
//'
//'     The normalization is such that \emph{-1/60 <= D <= 1/30}.
//'
//'     The return value -1.0 indicates an error.
//'
//' @export
//'
//' @examples
//'
//' .calc.hoeffding(0:4)
//' ## [1] 0.03333333
//'
//' .calc.hoeffding(c(0,3,2,1,4))
//' ## [1] -0.01666667
//'
//' set.seed(397)
//' .calc.hoeffding(order(runif(1000))-1) * 36
//' ## [1] 0.004349087
// [[Rcpp::export(name = ".calc.hoeffding")]]
double calc_hoeffding(const std::vector<unsigned long>& perm)
{
  entry n = perm.size();

  if ((n >= 5) && (n <= MAX_HOEFFDING_64))
  {
    int64 count = hoeffding_count<uint64>(perm);
    return count * 2.0 / n / (n-1) / (n-2) / (n-3) / (n-4);
  }

#ifdef __SIZEOF_INT128__
  if ((n >= 5) && (n <= MAX_HOEFFDING_128))
  {
    int128 count = hoeffding_count<uint128>(perm);
    return count * 2.0 / n / (n-1) / (n-2) / (n-3) / (n-4);
  }
#endif

  return -1.0;
}

//' Compute the refined Hoeffding statistic
//'
//' This is an internal CPP function, used by the R function
//' \code{\link{hoeffding.refined.test}}.
//'
//' Given (X1,Y1),...,(Xn,Yn), the refined Hoeffding statistic Rn only depends
//' on the permutation P that satisfies rank Yi = P[rank Xi].
//' This function computes Rn given P in O(n log n) time.
//'
//' @param perm
//'     An integer vector containing exactly 0,1,...,n-1 in any order.
//'
//'     The validity of the input is not checked by this function.
//'
//' @return
//'     The refined Hoeffding statistic of \code{perm}.
//'
//'     The normalization is such that \emph{-1/180 <= R <= 1/90}.
//'
//'     The return value -1.0 indicates an error.
//'
//' @export
//'
//' @examples
//'
//' .calc.refined(0:4)
//' ## [1] 0.01111111
//'
//' .calc.refined(c(0,3,2,1,4))
//' ## [1] -0.005555556
//'
//' set.seed(397)
//' .calc.refined(order(runif(1000))-1) * 36
//' ## [1] 0.004414034
// [[Rcpp::export(name = ".calc.refined")]]
double calc_refined(const std::vector<unsigned long>& perm)
{
    double D = calc_hoeffding(perm);
    double T = calc_taustar(perm);
    if ((D == -1.0) || (T == -1.0)) return -1.0;
    return -D/2.0 + T/24.0;
}

//--- end ---------------------------------------------------------------------
