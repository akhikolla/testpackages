c------------------------------------------------------------
c     n             I  the # of subjects
c     dx	    I  the censoring indicator for x_ij
c                        =1, if uncensored
c                        =0, if censored 
c     dy	    I  the censoring indicator for y_ij
c                        =1, if uncensored
c                        =0, if censored 
c     xmat          I  $x_{ij}$, first gaptimes
c     ymat          I  $y_{ij}$, second gaptimes
c                      gap = the time between two successive evnts
c     gmatx         I  Matrix with Kaplan-Meier estimators based on first bivariate gap times
c     gmaty         I  Matrix with Kaplan-Meier estimators based on second bivariate gap times
c     mstar         I  mstar_i = max(m_i -1, 1) where m_i is the number of episodes j for subject i
c     l1            I  limit of data for g1mat function such that max support of x_i < tau_c (support for g1)
c     l2            I  limit of data for g2mat function such that max support of y_i < tau_c (support for g2)
c     amat          I  matrix of covariates
c     expAx
c     expAy
c------------------------------------------------------------
c  Author: Sandra Castro-Pearson
c  Based on R code from: Chi Hyun Lee
c  Division of Biostatistics, School of Public Health
c  University of Minnesota 
c  Latest update: July, 2018
c------------------------------------------------------------
      subroutine mprovar(n, nparams, xmat, ymat,
     +gmatx, gmaty, l1, l2, expAx, expAy, subsumx, subsumy,
     +dx, dy, mstar, mc)

      implicit none
      integer n, nparams, mc, kcount, j, k
      double precision gmatx(n, mc), gmaty(n, mc), l1, l2
      double precision dx(n, mc), dy(n, mc)
      double precision xmat(n, mc), ymat(n, mc)
      double precision mstar(n), expAx(n), expAy(n)
      double precision subsumx(n), subsumy(n)
      double precision tsx(2), xmax(2), xtrm, xtmpsum
      double precision tsy(2), ymax(2), ytrm, ytmpsum

      xmax(2) = l1
      ymax(2) = l2

      do 12 j=1,n
        kcount = int(mstar(j))
        xtmpsum = 0
        ytmpsum = 0

	do 11 k=1,kcount
	  tsx(1) = xmat(j, k)
          tsx(2) = xmat(j, k)/expAx(j)
	  tsy(1) = xmat(j, k) + ymat(j, k)
          tsy(2) = xmat(j, k)/expAx(j) + ymat(j, k)/expAx(j)
          xmax(1) = maxval(tsx)
	  ymax(1) = maxval(tsy)
          xtrm=dx(j, k)*(log(minval(xmax))-log(l1))/gmatx(j, k)
          ytrm=dy(j, k)*(log(minval(ymax))-log(l2))/gmaty(j, k)
          xtmpsum = xtmpsum + xtrm
          ytmpsum = ytmpsum + ytrm

 11     continue
      subsumx(j) = xtmpsum / kcount
      subsumy(j) = ytmpsum / kcount
 12   continue 

      return 
      end



