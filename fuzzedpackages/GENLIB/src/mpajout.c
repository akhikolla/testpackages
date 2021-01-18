/*
    mpajout.c

    by Sébastien 
    
    Arbitrary precision integer arithmetic library
	
	Set or Clear a specific bit...

    $Id: mpajout.c,v 2.0 2003/09/23 20:58:20 sebasLec Exp $
 */

#include "mplogic.h"
#include <stdlib.h>

/* Petite fonctions de l'intérieur de MPI que nous allons utilisé*/
extern mp_err   s_mp_pad(mp_int *mp, mp_size min);

#define  s_mp_clamp(mp)\
   { while(USED(mp) > 1 && DIGIT((mp), USED(mp) - 1) == 0) USED(mp) -= 1; }

mp_err mpl_bit_set(mp_int *a, int bit)
{
  unsigned int ddigit,bbit;
  mp_err   res;
  //int bob;

  ARGCHK(a != NULL, MP_BADARG);

  //bob=DIGIT_BIT;
  ddigit = bit / DIGIT_BIT;
  bbit	 = bit % DIGIT_BIT;

  if((res = s_mp_pad(a,ddigit+1)) != MP_OKAY)
    return res;

  DIGIT(a, ddigit) = DIGIT(a, ddigit) | (1 << bbit);  
  return MP_OKAY;

} 

mp_err mpl_bit_clear(mp_int *a, int bit)
{
  unsigned int ddigit,bbit;

  ARGCHK(a != NULL, MP_BADARG);

  ddigit = bit / DIGIT_BIT;
  bbit	 = bit % DIGIT_BIT;

  if( ddigit >= USED(a) ) //il faut clear un bit déjà clear...
	  return MP_OKAY;

  DIGIT(a, ddigit) = DIGIT(a, ddigit) &  ~(1 << bbit);
  
  //Enleve les espace inutiles si requis
  s_mp_clamp(a)

  return MP_OKAY;
} 

