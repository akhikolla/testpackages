c------------------------------------------------------------
c     n             I  the # of subjects
c     delta1	    I  the censoring indicator for y_ij
c                        =1, if uncensored
c                        =0, if censored 
c     xmat          I  $x_{ij}$, first gaptimes
c     ymat          I  $y_{ij}$, second gaptimes
c                      gap = the time between two successive evnts
c     g1mat         I  Matrix with Kaplan-Meier estimators based on first pair of bivariate gap times
c     mstar         I  mstar_i = max(m_i -1, 1) where m_i is the number of episodes j for subject i
c     l1            I  limit of data for g1mat function such that max support of x_i < tau_c (support for g1)
c     amat          I  matrix of covariates
c------------------------------------------------------------
c  Author: Sandra Castro-Pearson
c  Based on R code from: Chi Hyun Lee
c  Division of Biostatistics, School of Public Health
c  University of Minnesota 
c  Latest update: July, 2018
c------------------------------------------------------------
      subroutine ymproee(n, nparams, di, xmati, ymati,
     +gmati, L, expA, subsum, kcount)

      implicit none
      integer n, nparams, kcount, j, k
      double precision di(kcount),xmati(kcount),ymati(kcount)
      double precision gmati(kcount), expA(n,2), subsum(n)
      double precision L, ts(2), minmax(2), trm, tmpsum
      
      minmax(2) = L

      do 12 j=1,n
	tmpsum = 0
	do 11 k=1,kcount
	  ts(1) = xmati(k) + ymati(k)
          ts(2) = xmati(k)*expA(j, 1)+ymati(k)*expA(j, 2)
          minmax(1) = maxval(ts)
          trm=di(k)*(log(minval(minmax))-log(L))/gmati(k)
          tmpsum = tmpsum + trm
 11     continue
      subsum(j) = tmpsum / kcount
 12   continue 

      return 
      end



