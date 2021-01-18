/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2016  Serge Iovleff, Université Lille 1, Inria

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as
    published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with this program; if not, write to the
    Free Software Foundation, Inc.,
    59 Temple Place,
    Suite 330,
    Boston, MA 02111-1307
    USA

    Contact : S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
*/

/* Project:  stkpp::STatistiK
 * Purpose:  Main pseudo-random uniform, Gaussian and exponential generators.
 * Author:   Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_RandBase.h
 *  @brief Declaration of the RandBase class.
 *  The RandBase class furnish the
 *  main pseudo random number generators (uniform, exponential, gauss).
 **/

#ifndef STK_RANDBASE_H
#define STK_RANDBASE_H

#include "STKernel/include/STK_Misc.h"

#ifdef IS_RTKPP_LIB /* is rtkpp lib */

#include <R_ext/Random.h>
#include <Rmath.h>

namespace STK
{
/** @ingroup Laws
 *  @brief class for the Base random generator.
 **/
class RandBase
{
  public:
    /** Default constructor */
    inline RandBase(){}
    /** destructor */
    inline ~RandBase(){}
    /** pseudo-random uniform generator.
     *  \f[
     *   f(x) = 1, \qquad 0< x <1
     *  \f]
     *  @return a uniform number in (0,1) using the R random generator
    **/
    inline Real randUnif() { GetRNGstate();
                             Real s = Rf_runif(0, 1);
                             PutRNGstate();
                             return s;
                           }
    /** @return same as randUnif().*/
    inline Real operator()() { return randUnif(); }
    /** Pseudo-random gaussian generator of the gaussian probability law:
     * \f[     f(x) = \frac{1}{\sqrt{2\pi}}
     *         \exp\left\{-\frac{x^2}{2}\right\}.
     * \f]
     * @param mu,sigma mean and standard deviation of the gaussian distribution
     * @return a real number from a normal (Gaussian) distribution.
    **/
    inline Real randGauss( Real const& mu = 0, Real const& sigma = 1)
                         { GetRNGstate();
                           Real s = Rf_rnorm(mu, sigma);
                           PutRNGstate();
                           return s;
                         }
    /** real number in (0,n) */
    inline Real rand( Real const& n ) { GetRNGstate();
                                        Real s = Rf_runif(0, n);
                                        PutRNGstate();
                                              return s;}
    /** DEPRECATED. real number in (0,n) */
    inline Real randDblExc( Real const& n ) { return rand(n);}
};
} // namespace STK

#else

// MersenneTwister header.
#include "MersenneTwister.h"

namespace STK
{

/** @ingroup Laws
 *  @brief class for the Base random generator.
 *
 * This class inherit from MTRand which should not be used directly.
 * Using RandBase, one get a Type safe generator for use in STK++
 * applications.
 *
 * This class furnish :
 * - A pseudo normalized uniform random variate generator
 * - A pseudo normalized Gaussian random variate generator
 * - A pseudo normalized exponential random variate generator
 *
 * @author George Marsaglia andWai Wan Tsang,
 *        "The Ziggurat Method for Generating Random Variables",
 *         Journal of Statistical Software,  Volume 5, 2000, Issue 8.
 *
 * @author Jurgen A Doornik,
 *        "An improved Ziggurat Method to generate Normal Random Sample"
 *        http://www.doornik.com/research.html (2005).
 *
 * For the exponential Law we remove the old method and use directly
 * the inverse pdf method.
 **/
class RandBase: protected MTRand
{
  public:
    /** Default constructor.
     * auto-initialize with /dev/urandom or time() and clock().
     * @param glimit maximal value of the boxes in the ziggourat method
     * @param gvol volume of each box in the ziggourat method
     * @param gsize number of boxes
     **/
      RandBase( Real const& glimit = 3.442619855899
              , Real const& gvol   = 9.91256303526217e-3
              , int gsize  = 128
              );
    /** Initialize with a simple int seed.
     * @param oneSeed seed of the generator
     * @param glimit maximal value of the boxes in the ziggourat method
     * @param gvol volume of each box in the ziggourat method
     * @param gsize number of boxes
     **/
    RandBase( int oneSeed
            , Real const& glimit = 3.442619855899
            , Real const& gvol   = 9.91256303526217e-3
            , int gsize  = 128
            );

    /** Initialize with a seed Array.
     * @param bigSeed seed of the generator
     * @param glimit maximal value of the boxes in the ziggourat method
     * @param gvol volume of each box in the ziggourat method
     * @param gsize number of boxes
     **/
    template< class TContainer1D>
    RandBase( TContainer1D const& bigSeed
            , Real const& glimit = 3.442619855899
            , Real const& gvol   = 9.91256303526217e-3
            , int gsize  = 128
            )
            : gsize_(gsize), glimit_(glimit), gvol_(gvol)
    {
      // dimension
      int first = bigSeed.begin(), size = bigSeed.size();
      uint32* arraySeed = new uint32[size];
      // cast int in uint32
      for (int i=first; i<bigSeed.end(); i++)
      { arraySeed[i-first] = uint32(bigSeed[i]);}
      // Re-seeding functions with same behavior as initializers
      seed(arraySeed, size);
      delete[] arraySeed;
      // initialize zigourat method for gaussian random generator
      gaussInit();
    }
   /** destructor. */
    ~RandBase();
   /** Pseudo-random int uniform generator.
    *  \f[
    *  f(x) = 1/(n+1), \quad 0\leq x \leq n
    *  \f]
    *  Return a [0,n] uniform integer number for n < 2^32 using the
    *  Mersenne Twister method. This is a wrapper of the MTRand class.
    *  @sa STK::MTRand
    **/
    inline int randDiscreteUnif() { return int(Real(randInt()));}
    /** pseudo-random uniform generator.
     *  This is a wrapper of the MTRand class.
     *  \f[
     *   f(x) = 1, \qquad 0< x <1
     *  \f]
     *  @return a uniform number in (0,1) using the Mersenne Twister method.
     *  @sa STK::MTRand
    **/
    inline Real randUnif() { return Real(randDblExc());}

    /** @return same as randUnif().*/
    inline Real operator()() { return (Real)MTRand::operator()(); }
    /** real number in [0,1] */
    inline Real rand() { return (Real)MTRand::rand(); }
    /** real number in [0,n] */
    inline Real rand( Real const& n ) { return (Real)MTRand::rand((double)n); }
    /** real number in [0,1) */
    inline Real randExc() { return (Real)MTRand::randExc(); }
    /** real number in [0,n) */
    inline Real randExc( Real const& n ) { return (Real)MTRand::randExc((double)n); }
    /** real number in (0,1) */
    inline Real randDblExc() { return (Real)MTRand::randDblExc(); }
    /** real number in (0,n) */
    inline Real randDblExc( Real const& n ) { return (Real)MTRand::randDblExc((double)n); }

    /** Pseudo-random gaussian generator of the gaussian probability law:
     * \f[     f(x) = \frac{1}{\sqrt{2\pi}}
     *         \exp\left\{-\frac{x^2}{2}\right\}.
     * \f]
     * @param mu mean of the gaussian distribution
     * @param sigma standard deviation of the gaussian distribution
     * @return a real number from a normal (Gaussian) distribution.
    **/
    Real randGauss(Real const& mu = 0, Real const& sigma = 1);
    /** Pseudo-random exponential generator.
     * \f[
     *    f(x) = \exp\left\{ -x \right\}.
     * \f]
     *  @return a real number from an exponential normalized
     *  distribution using the inverse pdf method.
    **/
    Real randExp();

  private:
    /*  Gauss parameters. */
    /** Number of box for the gaussian ziggourat method.
     * [0] is bottom box and [size_-1] is top box.
     **/
    const int gsize_;
    /** limit of the bottom box.    */
    Real const   glimit_;
    /** volume of each box and of the remaining tail. */
    Real const   gvol_;
    /** kn holds coordinates, such that each rectangle has same area.
     *  wn holds kn[i+1]/kn[i].
     *  fn holds exp(-0.5 * kn[i] * kn[i]).
     **/
    Real *kn, *wn, *fn;
    /** Initialization of the Zigourrat method. */
    void gaussInit();
};


/* Initialize with a simple int seed. */
inline RandBase::RandBase( int oneSeed, Real const& glimit, Real const& gvol, int gsize)
                         : MTRand( (uint32)oneSeed ), gsize_(gsize), glimit_(glimit), gvol_(gvol)
{ gaussInit();}

/* auto-initialize with /dev/urandom or time() and clock() */
inline RandBase::RandBase( Real const& glimit, Real const& gvol, int gsize)
                         : MTRand(), gsize_(gsize), glimit_(glimit), gvol_(gvol)
{ gaussInit();}

/* destructor. */
inline RandBase::~RandBase()
{
  delete [] kn;
  delete [] wn;
  delete [] fn;
}

/* Pseudo-random gaussian generator.
 * \f[     f(x) = \frac{1}{\sqrt{2\pi}}
 *         \exp\left\{-\frac{x^2}{2}\right\}
 * \f]
 *  Return a real number from a normal (Gaussian) normalized
 *  distribution.
**/
inline Real RandBase::randGauss(Real const& mu, Real const& sigma)
{
  while(1)
  {
    // uniforms number
    Real u = 2.0 * randUnif() - 1.0;
    // random box
    uint32 i = randInt() & 0x7F; //
    // squeeze step : try the rectangular boxes
    if (std::abs(u) < wn[i]) return mu + sigma * u * kn[i];
    // bottom box: we have to sample from the tail
    if (i == 0)
    {
      Real x;
      // loop
      do { x = randExp() / glimit_;}
      while ( 2.0 * randExp() < x * x);
      // result sign(u) * (x+dr)
      return mu + sigma * sign(u, x+glimit_);
    }
    else // other box
    { // x : result, f1 : intermediary result
      Real x = u * kn[i], f1 = fn[i+1];
      // reject step : is this a sample from the wedges ?
      if ( f1 + randUnif()*( fn[i]*exp(x*x/2.)-f1) < 1.0)
      { return mu + sigma * x;}
    }
  }
  // avoid warning at compilation
  return 0.;
}

/* Pseudo-random exponential generator.
 * \f[
 *      f(x) = \exp\left\{ -x \right\}
 * \f]
 *  Return a real number from an exponential normalized
 *  distribution.
**/
inline Real RandBase::randExp()
{ return -Real(std::log(randUnif()));}

/* Initialization (Zigourrat method).  */
inline void RandBase::gaussInit()
{
  kn = new Real[gsize_+1];
  wn = new Real[gsize_];
  fn = new Real[gsize_+1];

  Real f = exp(-0.5 * glimit_ * glimit_);
  kn[0]  = gvol_ / f; // [0] is bottom block: gvol/f(gbound)

  kn[gsize_] = 0;
  fn[gsize_] = 1;

  kn[1]  = glimit_;
  fn[1]  = f;

  for (int i = 2; i < gsize_; ++i)
  {
    kn[i] = sqrt(-2 * log(gvol_ / kn[i-1] + f));
    fn[i] = (f = exp(-0.5 * kn[i] * kn[i]));
  }
  for (int i = 0; i < gsize_; ++i)
    wn[i] = kn[i + 1] / kn[i];
}

} // namespace STK

#endif

#endif //STK_RANDBASE_H

