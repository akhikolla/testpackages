/*
 * cblas_drotg.c
 *
 * The program is a C interface to drotg.
 *
 * Written by Keita Teranishi.  2/11/1998
 *
 */

#ifndef CBLAS_DROTG_H_
#define CBLAS_DROTG_H_

inline void cblas_drotg(double *a, double *b, double *c, double *s) {
  F77_NAME(drotg)(a, b, c, s);
}

#endif
