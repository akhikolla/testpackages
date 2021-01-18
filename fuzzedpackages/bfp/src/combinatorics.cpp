#include "combinatorics.h"
#include <iostream>
#include <cstdlib>
#include <types.h>

using std::cout;

void ksub_next ( int n, int k, IntVector& a, bool *more, int &m, int &m2)

//****************************************************************************80
//
//  Purpose:
//
//    KSUB_NEXT generates the subsets of size K from a set of size N.
//
//  Modified:
//
//    29 May 2003
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms,
//    Academic Press, 1978, second edition,
//    ISBN 0-12-519260-6.
//
//  Parameters:
//
//    Input, int N, the size of the set from which subsets are drawn.
//
//    Input, int K, the desired size of the subsets.  K must
//    be between 0 and N.
//
//    Output, int A[K].  A[I] is the I-th element of the
//    subset.  Thus A[I] will be an integer between 1 and N.
//    Note that the routine will return the values in A
//    in sorted order: 1 <= A[0] < A[1] < ... < A[K-1] <= N
//
//    Input/output, bool *MORE.  Set MORE = FALSE before first call
//    for a new sequence of subsets.  It then is set and remains
//    TRUE as long as the subset computed on this call is not the
//    final one.  When the final subset is computed, MORE is set to
//    FALSE as a signal that the computation is done.
//
{
  int j;

  if ( k < 0 || n < k )
  {
      Rf_error("\nKSUB_NEXT - Fatal error!\nN = %d\nK = %d\nbut 0 <= K <= N is required!\n",
               n,
               k);
  }

  if ( !( *more ) )
  {
    m2 = 0;
    m = k;
  }
  else
  {
    if ( m2 < n-m )
    {
      m = 0;
    }
    m = m + 1;
    m2 = a[k-m];
  }

  for ( j = 1; j <= m; j++ )
  {
    a[k+j-m-1] = m2 + j;
  }

  *more = ( a[0] != (n-k+1) );

  return;
}

// ****************************************************************************************//

void comp_next ( int n, int k, IntVector& a, bool *more, int &h, int &t)

//****************************************************************************80
//
//  Purpose:
//
//    COMP_NEXT computes the compositions of the integer N into K parts.
//
//  Discussion:
//
//    The routine computes one composition on each call until there are no more.
//    For instance, one composition of 6 into 3 parts is
//    3+2+1, another would be 6+0+0.
//
//  Example:
//
//    The 28 compositions of 6 into three parts are:
//
//      6 0 0,  5 1 0,  5 0 1,  4 2 0,  4 1 1,  4 0 2,
//      3 3 0,  3 2 1,  3 1 2,  3 0 3,  2 4 0,  2 3 1,
//      2 2 2,  2 1 3,  2 0 4,  1 5 0,  1 4 1,  1 3 2,
//      1 2 3,  1 1 4,  1 0 5,  0 6 0,  0 5 1,  0 4 2,
//      0 3 3,  0 2 4,  0 1 5,  0 0 6.
//
//  Modified:
//
//    28 May 2003
//
//  Author:
//
//    Albert Nijenhuis, Herbert Wilf,
//
//    C++ translation by John Burkardt.
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms,
//    Academic Press, 1978, second edition,
//    ISBN 0-12-519260-6.
//
//  Parameters:
//
//    Input, int N, the integer whose compositions are desired.
//
//    Input, int K, the number of parts in the composition.
//
//    Input/output, int A[K], the parts of the composition.
//
//    Input/output, bool *MORE.
//    Set MORE = FALSE on first call.  It will be reset to TRUE on return
//    with a new composition.  Each new call returns another composition until
//    MORE is set to FALSE when the last composition has been computed
//    and returned.
//
{

  int i;

  if ( ! ( *more ) )
  {
    t = n;
    h = 0;
    a[0] = n;
    for ( i = 1; i < k; i++ )
    {
       a[i] = 0;
    }
  }
  else
  {
    if ( 1 < t )
    {
      h = 0;
    }

    h = h + 1;
    t = a[h-1];
    a[h-1] = 0;
    a[0] = t - 1;
    a[h] = a[h] + 1;

  }

  *more = ( a[k-1] != n );

  return;
}
