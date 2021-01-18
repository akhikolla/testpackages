#include "growfunctions.h"
using namespace Rcpp;
using namespace arma;
using namespace std;

SEXP lsqcluster(umat& S, ucolvec& ordscore, mat& phat, field<ucolvec>& bigS)
{
        BEGIN_RCPP
        //field<ucolvec> bigS;
        bigS.set_size(1,1);
        int n = S.n_rows;  int np = S.n_cols;
        //mat phat(np,np); 
        phat.zeros();
        //ucolvec ordscore(n); 
        ordscore.zeros();
        mat delta; delta.eye(np,np);
        mat diffsq(np,np); diffsq.zeros();
        urowvec s(np); colvec score(n);
	   unsigned long minscoreloc;
        int i, j, k, l;
        // create phat
        for(i = 0; i < n; i++)
        {
	       delta.eye();
            s = S.row(i);
            for(j = 0; j < np; j++)
            {
                for(k = j + 1; k < np; k++)
                {
                    if( s(j) == s(k) )
                    {
                        delta(j,k) = delta(k,j) = 1;
                    } /* end if */

                } /* end loop k */
            } /* end loop j */

            phat += delta;
        } /* end loop i*/
        phat *= (1/double(n));

	// compute least squares score
        for(i = 0; i < n; i++)
        {
	       delta.eye();
            s = S.row(i);
            for(j = 0; j < np; j++)
            {
                for(k = j + 1; k < np; k++)
                {
                    if( s(j) == s(k) )
                    {
                        delta(j,k) = delta(k,j) = 1;
                    } /* end if */

                } /* end loop k */
            } /* end loop j */
            diffsq = pow((delta - phat),2);
            score(i) = sum(sum(diffsq));
        } /* end loop i*/
        // find index of minimum score and return
        ordscore = sort_index(score);
        minscoreloc = ordscore(0);

        // extract cluster with minimum LS score (from truncated sampler of WB)
        urowvec s_min          = S.row(minscoreloc);
        urowvec cluster	      = unique( s_min ); 
        int M		           = cluster.n_elem;
        bigS.set_size(M,1);
        for( l = 0; l < M; l++ )
        {
          // entries in s.min not sequential
	     bigS(l,0)	= find( s_min == cluster(l) ) + 1; // The +1 to achieve R labels
        }

        END_RCPP
} /* end function lsqcluster */

