c------------------------------------------------------------
c     n             I  the # of subjects
c     gtime(n,mc)   I  ragged array  - make N x max(count) for completeness
c                      jth gap time for ith subject
c                      gap = the time between two successive evnts
c     ctime(n)      I  the final gap time of each subject (always censored)
c     count(n)      I  the # of gap times for per subject
c     mc            I  max count
c     cen(n,mc)     I  the censoring indicator for gtime(k,n)
c                        cen(.,.)=1, if uncensored
c                        cen(.,.)=0, if censored 
c     ucen(n)       I  the # of uncensored gap times per subject
c     nd            I  number of unique support points
c     udt(nd)       I  *unique and ordered* support points
c     r(nd)         O  risk set at time T
c     d(nd)         O  # events at time T
c     sest(nd)      O  the estimated survival function
c     std(nd)       O  the aymptotic standard error estimate
c------------------------------------------------------------

      subroutine onesamp(n,gtime,ctime,count,mc,m,
     +cen,ucen,nd,udt,tot,gap,event,
     +r,d,sest,std)
      integer n,mc,nd,cumni,curj,tot,m(n)
      double precision gtime(n,mc),cen(n,mc)
      double precision gap(tot),event(tot)
      double precision ctime(n)
      double precision count(n),ucen(n)
      double precision udt(nd),r(nd),d(nd)
      double precision sest(nd),std(nd)
      double precision hai,fai,phi,w

c     cumni is a counter : m(i-1) where m(0) = 0
c     curj is a counter within a counter

      cumni=0 
      do 11 i=1,n
         do 10 j=1,m(i)
            curj = cumni + j
            gtime(i,j) = gap(curj)
            cen(i,j) = event(curj)   
 10      continue
         cumni=cumni+m(i)
 11   continue 

      do 30 k=1,nd
         r(k)=0.
         d(k)=0.
         do 25 i=1,n
            if (count(i).gt.1.) then
               do 20 j=1,idint(ucen(i))
                  if(gtime(i,j).ge.udt(k)) then
                     r(k)=r(k)+ctime(i)/ucen(i)
                  endif
                  if(gtime(i,j).eq.udt(k)) then
                     d(k)=d(k)+ctime(i)/ucen(i)
                  endif
 20            continue
            else
               if(gtime(i,1).ge.udt(k)) r(k)=r(k)+ctime(i)
               if((gtime(i,1).eq.udt(k)).and.(cen(i,j).gt.0))
     +              d(k)=d(k)+ctime(i)
            endif
 25      continue
 30   continue

c Calculate survivor estimate at time T (i.e. time nd)
      sest(1)=1.-d(1)/r(1)
      do 50 i=2,nd
         sest(i)=sest(i-1)*(1.-d(i)/r(i))
 50   continue

c Calculate standard error
      do 180 k=1,nd
         phi=0.
         do 170 i=1,n
            fai=0.
            w=0.
            do 140 k2=1,k
               hai=0.
               if (count(i).gt.1.) then
                  do 130 j=1,idint(ucen(i))
                     if(gtime(i,j).ge.udt(k2)) then
                        hai=hai+ctime(i)/ucen(i)
                     endif
 130              continue
               elseif (gtime(i,1).ge. udt(k2)) then
                  hai=hai+ctime(i)
               endif
               w=w+hai*d(k2)/(r(k2)*r(k2))
 140        continue
            do 160 j2=1,idint(ucen(i))
               if (gtime(i,j2).lt. udt(k)) then
                  do 150 k3=1,nd
                     if (gtime(i,j2).eq. udt(k3)) then
                        fai=fai+ctime(i)/(ucen(i)*r(k3))
                     endif
 150              continue
               endif
 160        continue
            phi=phi+(w-fai)*(w-fai)
 170     continue
         std(k)=dsqrt(phi)*sest(k)
 180  continue

      return
      end


