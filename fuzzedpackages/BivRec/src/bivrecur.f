c------------------------------------------------------------
c     n             I  the # of subjects
c     gtime(n,mc)   I  ragged array  - make N x max(count) for completeness
c                      jth combined survival time for ith subject (V+W)
c     markvar1      I  $V_{ij}$, first gaptimes
c     markvar2      I  $W_{ij}$, second gaptimes
c                      gap = the time between two successive evnts
c     ctime(n)      I  weight function
c     count(n)      I  the # of gap times for per subject
c     mc            I  max count
c     cen(n,mc)     I  the censoring indicator for gtime(k,n)
c                        cen(.,.)=1, if uncensored
c                        cen(.,.)=0, if censored
c     ucen(n)       I  count(n)-1
c     nd            I  number of unique support points
c     udt(nd)       I  *unique and ordered* support points
c     tmpindex      I  the maximal k such that $t_k^*$<=u1+u2
c     r(nd)         O  risk set at time t_k^*
c     d(nd,1)       O  risk mass at t_k^*, constraint by mark var
c     d(nd,2)       O  risk mass at t_k^*
c     sest(nd)      O  the estimated survival function
c     Fest(nd)      O  the estimated prob at (x, u1, u2)
c     prob          O  the estimated prob at (u1, u2)
c     std           O  the aymptotic standard error estimate
c------------------------------------------------------------
c  Author: Chiung-Yu Huang
c  Division of Biostatistics, School of Public Health
c  University of Minnesota
c  Latest update: Feb 12, 2004
c------------------------------------------------------------
      subroutine bivrecur(n,gtime,ctime,count,mc,m,
     +cen,ucen,nd,udt,tot,gap,event,
     +r,d,sest,var,markvar1,markvar2,mark1,mark2,u1,u2,Fest,
     +tmpindex,prob,std)

      integer n,mc,nd,cumni,curj,tot,m(n), tmpindex, flagk
      double precision gtime(n,mc),cen(n,mc),mark1(n,mc),mark2(n,mc)
      double precision gap(tot),event(tot),markvar1(tot),markvar2(tot)
      double precision ctime(n)
      double precision count(n),ucen(n)
      double precision udt(nd),r(nd),d(nd,2)
      double precision sest(nd),Fest(nd)
      double precision fai1, fai2
      double precision u1,u2,prob,std, var
      double precision w, psi

      cumni=0
      do 11 i=1,n
         do 10 j=1,m(i)
            curj = cumni + j
            gtime(i,j) = gap(curj)
            cen(i,j) = event(curj)
            mark1(i,j) = markvar1(curj)
            mark2(i,j) = markvar2(curj)
 10      continue
         cumni=cumni+m(i)
 11   continue

      do 30 k=1,nd
         r(k)=0.
         d(k,1)=0.
         d(k,2)=0.
         do 25 i=1,n
            if (count(i).gt.1.) then
               do 20 j=1,idint(ucen(i))
                  if(gtime(i,j).ge.udt(k)) then
                     r(k)=r(k)+ctime(i)/(ucen(i)*n)
                  endif
                  if(gtime(i,j).eq.udt(k)) then
                     d(k,2)=d(k,2)+ctime(i)/(ucen(i)*n)
                  endif

                  if((gtime(i,j).eq.udt(k)).and.(mark1(i,j).le.u1).
     +               and.(mark2(i,j).le.u2)) then
                     d(k,1)=d(k,1)+ctime(i)/( ucen(i)*n )
                  endif
 20            continue
            else
               if(gtime(i,1).ge.udt(k)) r(k)=r(k)+ctime(i)/n
               if((gtime(i,1).eq.udt(k)).and.(cen(i,1).gt.0))
     +              d(k,2)=d(k,2)+ctime(i)/n
     		if((gtime(i,1).eq.udt(k)).and.(cen(i,1).gt.0).and.
     +		   (mark1(i,1).le.u1).and.(mark2(i,1).le.u2))
     +              d(k,1)=d(k,1)+ctime(i)/n
            endif
 25      continue
 30   continue

c Calculate joint cdf estimate at time T (i.e. time nd)

      sest(1)=1.-d(1,2)/r(1)
      do 50 i=2,nd
         sest(i)=sest(i-1)*(1.-d(i,2)/r(i))
 50   continue


      Fest(1)=d(1,1)/r(1)

      do 60 i=2,nd
	 Fest(i)=Fest(i-1)+sest(i-1)*d(i,1)/r(i)
 60   continue


c     report the estimated bivaraite distribution function

      prob=Fest(tmpindex)

c------------------------------------------------------------------
c     calculate the standard deviation


      var=0.

      do 170 i=1,n
         psi=0.
         flagk=0


         do 190 k=1,nd

            flagk = flagk + 1
            fai1 = 0.
            fai2 = 0.
            w = 0.

            if ( udt(k).le.(u1+u2) ) then

               if ( count(i).gt.1 ) then
                  do 200 j1=1, idint( ucen(i) )

                     if ( gtime(i,j1) .eq. udt(k) ) then
                        fai2 = fai2 + ctime(i)/ucen(i)
                     endif

                     if ( (gtime(i,j1) .eq. udt(k)).and.
     +                    ( mark1(i,j1).le.u1 ).and.(mark2(i,j1).le.u2))
     +                    then
                        fai1 = fai1 + ctime(i)/ucen(i)
                     endif
                     if ( gtime(i,j1).ge. udt(k)) then
                        w = w + ctime(i)/ucen(i)
                     endif

 200              continue
               else
                  if ( (gtime(i,1).eq.udt(k)).and. (cen(i,1).eq.1)) then
                     fai2 = fai2 + ctime(i)
                  endif

                  if( (gtime(i,1).eq.udt(k)).and.(cen(i,1).eq.1).and.
     +                (mark1(i,1).le.u1) .and. (mark2(i,1).le.u2) ) then
                     fai1 = fai1 + ctime(i)
                  endif

                  if( gtime(i,1).ge.udt(k) ) then
                     w = w + ctime(i)
                  endif
               endif


               if ( flagk .eq. 1 ) then
                  psi = psi+ fai2 * ( Fest(1)-prob )/(sest(1)*r(1) )
     +                  - w*( Fest(1)-prob)*d(1,2)/(sest(1)*r(1)*r(1))
     +                  + fai1/r(1)
     +                  - w * d(1,1)/( r(1)*r(1) )
               else
                  psi =psi+fai2*(Fest(k)-prob)* sest(k-1)/(sest(k)*r(k))
     +                 -w*(Fest(k)-prob)*sest(k-1)*d(k,2)/(sest(k)
     +                         *r(k)*r(k))
     +                 +fai1 * sest(k-1) / r(k)
     +                 -w * sest(k-1) * d(k,1)/( r(k)*r(k) )
               endif

            endif

 190     continue

         var = var + psi * psi/(n*n)

 170  continue

         std = dsqrt(var)

      return
      end


