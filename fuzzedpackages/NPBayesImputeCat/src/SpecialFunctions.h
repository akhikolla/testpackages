#pragma once
#include <vector>
#include <climits>
#include <cstdio>
#include <ctime>
#include <cmath>
using namespace std;

class MTRand {
public:
  typedef unsigned long uint32;
  enum { N = 624 };
  enum { SAVE = N + 1 };

protected:
  enum { M = 397 };  // period parameter
  uint32 state[N];   // internal state
  uint32 *pNext;     // next value to get from state
  int left;          // number of values left before reload needed

  //Methods
public:
  MTRand( const uint32& oneSeed );
  MTRand( uint32 *const bigSeed, uint32 const seedLength = N );
  MTRand();

  // Access to 32-bit random numbers
  double rand();                          // real number in [0,1]
  double rand( const double& n );         // real number in [0,n]
  double randExc();                       // real number in [0,1)
  double randExc( const double& n );      // real number in [0,n)
  double randDblExc();                    // real number in (0,1)
  double randDblExc( const double& n );   // real number in (0,n)
  uint32 randInt();                       // integer in [0,2^32-1]
  uint32 randInt( const uint32& n );      // integer in [0,n] for n < 2^32
  double operator()() { return rand(); }  // same as rand()

  // Access to nonuniform random number distributions
  double randNorm( const double& mean = 0.0, const double& variance = 0.0 );

  // Re-seeding functions with same behavior as initializers
  void seed( const uint32 oneSeed );
  void seed( uint32 *const bigSeed, const uint32 seedLength = N );
  void seed();

protected:
  void initialize( const uint32 oneSeed );
  void reload();
  uint32 hiBit( const uint32& u ) const { return u & 0x80000000UL; }
  uint32 loBit( const uint32& u ) const { return u & 0x00000001UL; }
  uint32 loBits( const uint32& u ) const { return u & 0x7fffffffUL; }
  uint32 mixBits( const uint32& u, const uint32& v ) const
  { return hiBit(u) | loBits(v); }
  uint32 twist( const uint32& m, const uint32& s0, const uint32& s1 ) const
  { return m ^ (mixBits(s0,s1)>>1) ^ (-loBit(s1) & 0x9908b0dfUL); }
  static uint32 hash( time_t t, clock_t c );
};


inline MTRand::MTRand( const uint32& oneSeed )
{ seed(oneSeed); }

inline MTRand::MTRand( uint32 *const bigSeed, const uint32 seedLength )
{ seed(bigSeed,seedLength); }

inline MTRand::MTRand()
{ seed(); }

inline double MTRand::rand()
{ return double(randInt()) * (1.0/4294967295.0); }

inline double MTRand::rand( const double& n )
{ return rand() * n; }

inline double MTRand::randExc()
{ return double(randInt()) * (1.0/4294967296.0); }

inline double MTRand::randExc( const double& n )
{ return randExc() * n; }

inline double MTRand::randDblExc()
{ return ( double(randInt()) + 0.5 ) * (1.0/4294967296.0); }

inline double MTRand::randDblExc( const double& n )
{ return randDblExc() * n; }

inline double MTRand::randNorm( const double& mean, const double& variance )
{
  // Return a real number from a normal (Gaussian) distribution with given
  // mean and variance by Box-Muller method
  double r = sqrt( -2.0 * log( 1.0-randDblExc()) ) * variance;
  double phi = 2.0 * 3.14159265358979323846264338328 * randExc();
  return mean + r * cos(phi);
}

inline MTRand::uint32 MTRand::randInt()
{
  // Pull a 32-bit integer from the generator state
  // Every other access function simply transforms the numbers extracted here

  if( left == 0 ) reload();
  --left;

  uint32 s1;
  s1 = *pNext++;
  s1 ^= (s1 >> 11);
  s1 ^= (s1 <<  7) & 0x9d2c5680UL;
  s1 ^= (s1 << 15) & 0xefc60000UL;
  return ( s1 ^ (s1 >> 18) );
}

inline MTRand::uint32 MTRand::randInt( const uint32& n ) {
  uint32 used = n;
  used |= used >> 1;
  used |= used >> 2;
  used |= used >> 4;
  used |= used >> 8;
  used |= used >> 16;
  uint32 i;
  do
    i = randInt() & used;
  while( i > n );
  return i;
}

inline void MTRand::seed( const uint32 oneSeed )
{
  initialize(oneSeed);
  reload();
}


inline void MTRand::seed( uint32 *const bigSeed, const uint32 seedLength ) {
  initialize(19650218UL);
  int i = 1;
  uint32 j = 0;
  int k = ( N > seedLength ? N : seedLength );
  for( ; k; --k )
  {
    state[i] =
      state[i] ^ ( (sta