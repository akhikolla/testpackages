#ifndef RANDGEN_H
#define RANDGEN_H

// [AMR] Note that a C++ RNG is used in more recent ROPTLIB versions

#if !defined(R_BUILD)
#error This version of randgen.h only works with R builds
#endif

/* generates a random number on [0,1]-real-interval */
double genrand_real1();

/* generates a random number on [0,1)-real-interval */
double genrand_real2();

/* generate a random number following the standard normal distribution*/
double genrand_gaussian();

#endif
