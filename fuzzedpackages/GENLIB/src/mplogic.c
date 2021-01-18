/*
    mplogic.c

    by Michael J. Fromberger <sting@linguist.dartmouth.edu>
    Copyright (C) 1998 Michael J. Fromberger, All Rights Reserved

    Bitwise logical operations on MPI values

    $Id: mplogic.c,v 2.0 2003/09/23 20:58:20 sebasLec Exp $
 */

#include "mplogic.h"
#include <stdlib.h>

/* Some things from the guts of the MPI library we make use of... */
extern mp_err   s_mp_lshd(mp_int *mp, mp_size p);
extern void     s_mp_rshd(mp_int *mp, mp_size p);

#define  s_mp_clamp(mp)\
   { while(USED(mp) > 1 && DIGIT((mp), USED(mp) - 1) == 0) USED(mp) -= 1; }

/* {{{ Lookup table for population count */

static unsigned char bitc[] = {
   0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4, 
   1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 
   1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 
   2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 
   1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 
   2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 
   2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 
   3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 
   1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 
   2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 
   2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 
   3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 
   2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 
   3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 
   3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 
   4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8
};

/* }}} */

/*------------------------------------------------------------------------*/
/*
  mpl_not(a, b)    - compute b = ~a
  mpl_and(a, b, c) - compute c = a & b
  mpl_or(a, b, c)  - compute c = a | b
  mpl_xor(a, b, c) - compute c = a ^ b
 */

/* {{{ mpl_not(a, b) */

mp_err mpl_not(mp_int *a, mp_int *b)
{
  mp_err   res;
  int      ix;

  ARGCHK(a != NULL && b != NULL, MP_BADARG);

  if((res = mp_copy(a, b)) != MP_OKAY)
    return res;

  /* This relies on the fact that the digit type is unsigned */
  for(ix = 0; ix < USED(b); ix++) 
    DIGIT(b, ix) = ~DIGIT(b, ix);

  s_mp_clamp(b);

  return MP_OKAY;

} /* end mpl_not() */

/* }}} */

/* {{{ mpl_and(a, b, c) */

mp_err mpl_and(mp_int *a, mp_int *b, mp_int *c)
{
  mp_int  *which, *other;
  mp_err   res;
  int      ix;

  ARGCHK(a != NULL && b != NULL && c != NULL, MP_BADARG);

  if(USED(a) <= USED(b)) {
    which = a;
    other = b;
  } else {
    which = b;
    other = a;
  }

  if((res = mp_copy(which, c)) != MP_OKAY)
    return res;

  for(ix = 0; ix < USED(which); ix++)
    DIGIT(c, ix) &= DIGIT(other, ix);

  s_mp_clamp(c);

  return MP_OKAY;

} /* end mpl_and() */

/* }}} */

/* {{{ mpl_or(a, b, c) */

mp_err mpl_or(mp_int *a, mp_int *b, mp_int *c)
{
  mp_int  *which, *other;
  mp_err   res;
  int      ix;

  ARGCHK(a != NULL && b != NULL && c != NULL, MP_BADARG);

  if(USED(a) >= USED(b)) {
    which = a;
    other = b;
  } else {
    which = b;
    other = a;
  }

  if((res = mp_copy(which, c)) != MP_OKAY)
    return res;

  for(ix = 0; ix < USED(which); ix++)
    DIGIT(c, ix) |= DIGIT(other, ix);

  return MP_OKAY;

} /* end mpl_or() */

/* }}} */

/* {{{ mpl_xor(a, b, c) */

mp_err mpl_xor(mp_int *a, mp_int *b, mp_int *c)
{
  mp_int  *which, *other;
  mp_err   res;
  int      ix;

  ARGCHK(a != NULL && b != NULL && c != NULL, MP_BADARG);

  if(USED(a) >= USED(b)) {
    which = a;
    other = b;
  } else {
    which = b;
    other = a;
  }

  if((res = mp_copy(which, c)) != MP_OKAY)
    return res;

  for(ix = 0; ix < USED(which); ix++)
    DIGIT(c, ix) ^= DIGIT(other, ix);

  s_mp_clamp(c);

  return MP_OKAY;

} /* end mpl_xor() */

/* }}} */

/*------------------------------------------------------------------------*/
/*
  mpl_rsh(a, b, d)     - b = a >> d
  mpl_lsh(a, b, d)     - b = a << d
 */

/* {{{ mpl_rsh(a, b, d) */

mp_err mpl_rsh(mp_int *a, mp_int *b, mp_digit d)
{
  mp_err   res;
  mp_digit dshift, bshift;

  ARGCHK(a != NULL && b != NULL, MP_BADARG);

  dshift = d / DIGIT_BIT;  /* How many whole digits to shift by */
  bshift = d % DIGIT_BIT;  /* How many bits to shift by         */

  if((res = mp_copy(a, b)) != MP_OKAY)
    return res;

  /* Shift over as many whole digits as necessary */
  if(dshift) 
    s_mp_rshd(b, dshift);

  /* Now handle any remaining bit shifting */
  if(bshift)
  {
    mp_digit  prev = 0, next, mask = (1 << bshift) - 1;
    int       ix;

    /*
      'mask' is a digit with the lower bshift bits set, the rest
      clear.  It is used to mask off the bottom bshift bits of each
      digit, which are then shifted on to the top of the next lower
      digit.
     */
    for(ix = USED(b) - 1; ix >= 0; ix--) {
      /* Take off the lower bits and shift them up... */
      next = (DIGIT(b, ix) & mask) << (DIGIT_BIT - bshift);

      /* Shift down the current digit, and mask in the bits saved
	 from the previous digit 
       */
      DIGIT(b, ix) = (DIGIT(b, ix) >> bshift) | prev;
      prev = next;
    }
  }

  s_mp_clamp(b);   

  return MP_OKAY;

} /* end mpl_rsh() */

/* }}} */

/* {{{ mpl_lsh(a, b, d) */

mp_err mpl_lsh(mp_int *a, mp_int *b, mp_digit d)
{
  mp_err   res;
  mp_digit dshift, bshift;

  ARGCHK(a != NULL && b != NULL, MP_BADARG);

  dshift = d / DIGIT_BIT;
  bshift = d % DIGIT_BIT;

  if((res = mp_copy(a, b)) != MP_OKAY)
    return res;

  if(dshift)
    if((res = s_mp_lshd(b, dshift)) != MP_OKAY)
      return res;

  if(bshift){ 
    int       ix;
    mp_digit  prev = 0, next, mask = (1 << bshift) - 1;

    for(ix = 0; ix < USED(b); ix--) {
      next = (DIGIT(b, ix) >> (DIGIT_BIT - bshift)) & mask;

      DIGIT(b, ix) = (DIGIT(b, ix) << bshift) | prev;
      prev = next;
    }
  }

  s_mp_clamp(b);

  return MP_OKAY;

} /* end mpl_lsh() */

/* }}} */

/*------------------------------------------------------------------------*/
/*
  mpl_num_set(a, num)

  Count the number of set bits in the binary representation of a.
  Returns MP_OKAY and sets 'num' to be the number of such bits, if
  possible.  If num is NULL, the result is thrown away, but it is
  not considered an error.

  mpl_num_clear() does basically the same thing for clear bits.
 */

/* {{{ mpl_num_set(a, num) */

mp_err mpl_num_set(mp_int *a, int *num)
{
  int            ix, db, nset = 0;
  mp_digit       cur;
  unsigned char  reg;

  ARGCHK(a != NULL, MP_BADARG);

  for(ix = 0; ix < USED(a); ix++) {
    cur = DIGIT(a, ix);
    
    for(db = 0; db < sizeof(mp_digit); db++) {
      reg = (cur >> (CHAR_BIT * db)) & UCHAR_MAX;

      nset += bitc[reg];
    }
  }

  if(num)
    *num = nset;

  return MP_OKAY;

} /* end mpl_num_set() */

/* }}} */

/* {{{ mpl_num_clear(a, num) */

mp_err mpl_num_clear(mp_int *a, int *num)
{
  int            ix, db, nset = 0;
  mp_digit       cur;
  unsigned char  reg;

  ARGCHK(a != NULL, MP_BADARG);

  for(ix = 0; ix < USED(a); ix++) {
    cur = DIGIT(a, ix);
    
    for(db = 0; db < sizeof(mp_digit); db++) {
      reg = (cur >> (CHAR_BIT * db)) & UCHAR_MAX;

      nset += bitc[UCHAR_MAX - reg];
    }
  }

  if(num)
    *num = nset;

  return MP_OKAY;


} /* end mpl_num_clear() */

/* }}} */

/*------------------------------------------------------------------------*/
/*
  mpl_parity(a)

  Determines the bitwise parity of the value given.  Returns MP_EVEN
  if an even number of digits are set, MP_ODD if an odd number are
  set.
 */

/* {{{ mpl_parity(a) */

mp_err mpl_parity(mp_int *a)
{
  int      ix, par = 0;
  mp_digit cur;

  ARGCHK(a != NULL, MP_BADARG);

  for(ix = 0; ix < USED(a); ix++) {
    int   shft = (sizeof(mp_digit) * CHAR_BIT) / 2;

    cur = DIGIT(a, ix);

    /* Compute parity for current digit */
    while(shft != 0) {
      cur ^= (cur >> shft);
      shft >>= 1;
    }
    cur &= 1;

    /* XOR with running parity so far   */
    par ^= cur;
  }

  if(par)
    return MP_ODD;
  else
    return MP_EVEN;

} /* end mpl_parity() */

/* }}} */

/*------------------------------------------------------------------------*/
/* HERE THERE BE DRAGONS                                                  */
