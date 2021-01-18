c----------------------------------------------------------------------
c     Stanislav Nagy
c     nagy at karlin.mff.cuni.cz
c     27/06/2016
c
c     Fortran implementation of functional depths and kernel smoothing
c     - Nagy, Ferraty: Depth for Noisy Random Functions
c     - Nagy, Gijbels, Hlubinka: Depth-Based Recognition of Shape 
c           Outlying Functions
c
c----------------------------------------------------------------------

c----------------------------------------------------------------------
c     kernel smoothing of a set of functions
c----------------------------------------------------------------------

      SUBROUTINE KERNSM(T,X,G,M,N,H,K,R)
c     kernel smoothing of a single function
c     T vector of observation points (M)
c     X vector of observed values (M)
c     G grid to evaluate (N)
c     H bandwidth
c     K kernel code to use, see below
c     R resulting estimate (N)

      integer M,N,K
      integer I,J
      double precision T(M),X(M),G(N),R(N),H
      double precision DEN,P
      
      DO 10 I=1,N
            R(I)=0.0
            DEN = 0.0
            DO 5 J=1,M
                  CALL KERN((G(I)-T(J))/H,P,K)
                  R(I) = R(I) + X(J)*P
                  DEN = DEN + P
5           CONTINUE
            IF (DEN.GT.0) THEN
                  R(I) = R(I)/DEN
               ELSE
                  R(I) = 10**6
            ENDIF      
10    CONTINUE
      RETURN
      END
      
c----------------------------------------------------------------------

      SUBROUTINE KERN(T,R,K)
c     kernel functions
      
      double precision T,R
      INTEGER K
      
      IF(K.EQ.1) THEN
c     uniform kernel
            IF (ABS(T).LE.1.0) THEN
                  R = .5
                  ELSE
                  R = 0.0
            ENDIF
      ENDIF
      IF(K.EQ.2) THEN
c     triangular kernel
            IF (ABS(T).LE.1.0) THEN
                  R = 1-ABS(T)
                  ELSE
                  R = 0.0
            ENDIF
      ENDIF
      IF(K.EQ.3) THEN
c     Epanechnikov kernel
            IF (ABS(T).LE.1.0) THEN
                  R = 3.0/4.0*(1.0-T**2.0)
                  ELSE
                  R = 0.0
            ENDIF
      ENDIF
      IF(K.EQ.4) THEN
c     biweight kernel
            IF (ABS(T).LE.1.0) THEN
                  R = 15.0/16.0*(1.0-T**2.0)**2
                  ELSE
                  R = 0.0
            ENDIF
      ENDIF
      IF(K.EQ.5) THEN
c     triweight kernel
            IF (ABS(T).LE.1.0) THEN
                  R = 35.0/32.0*(1.0-T**2.0)**3
                  ELSE
                  R = 0.0
            ENDIF
      ENDIF
      IF(K.EQ.6) THEN
c     Gaussian kernel
            R = (2*4.D0*DATAN(1.D0))**(-0.5)*EXP(-0.5*T**2)
      ENDIF
      RETURN
      END
  
c----------------------------------------------------------------------
      
      SUBROUTINE CVKERNSM(T,X,G,M,N,H,NH,KER,TRE,XRE,TNRE,XNRE,NR,NT,R)
c     as in KERNSM, but with automated BW seletction
c     from the vector H(NH)
c     KER kernel code
c     TRE,XRE vectors from T,X for omission for the CV (NR*NT)
c     TNRE,XNRE supplementary vectors to TRE,NRE ((M-NR)*NT)
c     NR no of random elements for each H choice
c     NT no of random trials for each H choice
c     for each H we choose NR random points from X
c     compute the kernel smoother without these points
c     and evaluate the error. repeat NT times, then average
c     the errors to get the CV estimate

      integer M,N,NH,NR,NT,KER
      integer I,J,K,OM(NR),L
      double precision T(M),X(M),G(N),R(N),H(NH),TRE(NR*NT),XRE(NR*NT)
      double precision TNRE((M-NR)*NT),XNRE((M-NR)*NT)
      double precision TOM(NR),XOM(NR),TNOM(M-NR),XNOM(M-NR)
      double precision CV(NH),SR(NR),MCV
      
      DO 30 I=1,NH
            CV(I) = 0.0
            DO 20 J=1,NT
                  DO 5 K=1,NR
                        TOM(K) = TRE((J-1)*NR+K)
                        XOM(K) = XRE((J-1)*NR+K)
5                 CONTINUE
                  DO 6 K=1,(M-NR)
                        TNOM(K) = TNRE((J-1)*(M-NR)+K)
                        XNOM(K) = XNRE((J-1)*(M-NR)+K)
6                 CONTINUE
            CALL KERNSM(TNOM,XNOM,TOM,(M-NR),NR,H(I),KER,SR)
                  DO 10 K=1,NR
                  CV(I) = CV(I) + (XOM(K)-SR(K))**2
10                CONTINUE
20          CONTINUE            
30    CONTINUE
      K = 0
      MCV = CV(1)+1
      DO 40 I=1,NH
            IF (CV(I).LT.MCV) THEN
                  K = I
                  MCV = CV(I)
            ENDIF
40    CONTINUE
      CALL KERNSM(T,X,G,M,N,H(K),KER,R)
      RETURN
      END
      
c----------------------------------------------------------------------
c     fucntional data depth computation      
c----------------------------------------------------------------------

      SUBROUTINE funD1(A,B,M,N,D,funSDEP,funHDEP,fIsdep,fIhdep,IAsdep,
     +IAhdep)
c     univariate integrated depth
       
      INTEGER N,M,D
      INTEGER IAsdep(M),IAhdep(M)
      double precision A(M*D),B(N*D),funSDEP(M),funHDEP(M),
     +fIsdep(M),fIhdep(M) 
      double precision BH(N)
      double precision hSDEP,hHDEP
      INTEGER I,J,K

c     initializing
      DO 5 I=1,M
            funSDEP(I) = 0.0
            funHDEP(I) = 0.0
            fISDEP(I) = 2.0
            fIHDEP(I) = 2.0 
            IAsdep(I) = 0
            IAhdep(I) = 0
5     CONTINUE
c     essential loop, computing 1D depths at each point
      DO 30 I=1,D
            DO 10 K=1,N
                  BH(K) = B((I-1)*N+K)
10          CONTINUE
            DO 20 J=1,M
                  hSDEP = 0.0
                  hHDEP = 0.0
                  CALL fD1(A((I-1)*M+J),N,BH,hSDEP,hHDEP)
c     integrated depth evaluation
                  funSDEP(J) = funSDEP(J) + hSDEP
                  funHDEP(J) = funHDEP(J) + hHDEP
c     counting the area of the smallest depth for each function
                  IF (hSDEP.EQ.fISDEP(J)) THEN
                        IAsdep(J) = IAsdep(J)+1
                  ELSEIF (hSDEP.LT.fISDEP(J)) THEN 
                        IAsdep(J) = 1
                  ENDIF
                  IF (hHDEP.EQ.fIHDEP(J)) THEN
                        IAhdep(J) = IAhdep(J)+1
                  ELSEIF (hHDEP.LT.fIHDEP(J)) THEN 
                        IAhdep(J) = 1
                  ENDIF 
c     infimal depth evaluation                 
                  fISDEP(J) = min(fISDEP(J),hSDEP)
                  fIHDEP(J) = min(fIHDEP(J),hHDEP)
20          CONTINUE
30    CONTINUE
c     dividing the resulting depths by number of points d
      DO 40 I=1,M
            funSDEP(I) = funSDEP(I)/(D+0.0)
            funHDEP(I) = funHDEP(I)/(D+0.0)
40    CONTINUE
      RETURN
      END 
      
c----------------------------------------------------------------------

      SUBROUTINE fD1(U,N,X,SDEP,HDEP)
c     computes 1D simplicial and halfspace depth of U wrt X of size N
c     used in the univariate integrated depth

      double precision U,X(N),SDEP,HDEP
      INTEGER N,RB,RA,I
      
      RB = 0
      RA = 0
      DO 10 I=1,N
            IF (U .LE. X(I)) THEN 
                  RB = RB + 1
                  ENDIF
            IF (U .GE. X(I)) THEN 
                  RA = RA + 1
                  ENDIF
10    CONTINUE     
      HDEP = min(RA+0.0,RB+0.0)/(N+0.0)
      SDEP = (RA+0.0)*(RB+0.0)/(K(N,2)+0.0)
      RETURN
      END
      
c----------------------------------------------------------------------

      INTEGER FUNCTION K(M,J)
c     combination number (m choose j)

      integer m,j
      IF (M.LT.J) THEN
          K=0
      ELSE
          IF (J.EQ.1) K=M
          IF (J.EQ.2) K=(M*(M-1))/2
          IF (J.EQ.3) K=(M*(M-1)*(M-2))/6
      ENDIF
      RETURN
      END
      
c----------------------------------------------------------------------

      SUBROUTINE metrl2(A,B,M,N,D,METR)
c     computes a fast approximation of the L2 metric between A and B
c     supporting functions on regular grids only
c     used in the h-mode depth 
        
      INTEGER N,M,D
      double precision A(M*D),B(N*D),METR(M*N)
      INTEGER I,J,K

      DO 15 I=1,M
            DO 10 J=1,N
            METR((J-1)*M+I) = 0.0
                  DO 5 K=1,D
                  METR((J-1)*M+I) = METR((J-1)*M+I) + (A((K-1)*M+I)-
     +B((K-1)*N+J))**2
5                 CONTINUE
            METR((J-1)*M+I) = sqrt(METR((J-1)*M+I) - 
     +((A((0)*M+I)-B((0)*N+J))**2+(A((D-1)*M+I)-B((D-1)*N+J))**2)
     +/(2.0))
10          CONTINUE
15    CONTINUE
      RETURN
      END 
      
c----------------------------------------------------------------------

      SUBROUTINE metrC(A,B,M,N,D,METR)
c     computes a fast approximation of the C metric between A and B
c     supporting functions on regular grids only
c     used in the h-mode depth 
        
      INTEGER N,M,D
      double precision A(M*D),B(N*D),METR(M*N)
      INTEGER I,J,K

      DO 15 I=1,M
            DO 10 J=1,N
            METR((J-1)*M+I) = 0.0
                  DO 5 K=1,D
                  METR((J-1)*M+I) = MAX(METR((J-1)*M+I),A((K-1)*M+I)-
     +B((K-1)*N+J))
            METR((J-1)*M+I) = MAX(METR((J-1)*M+I),-A((K-1)*M+I)+
     +B((K-1)*N+J))
5                 CONTINUE
10          CONTINUE
15    CONTINUE
      RETURN
      END 

c----------------------------------------------------------------------

      double precision FUNCTION MAX(A,B)
c     maximum of two numbers A and B

      double precision A,B
      
      IF (A.LE.B) THEN
            MAX = B
      ELSE
            MAX = A
      ENDIF
      RETURN
      END 

      double precision FUNCTION AdjLPindicator(EVAL,J,B,V)
c     adjusted band depth core function, smoothing (exp(-u)), L2 metric
c     b is a vector of length eval
c     v is a matrix j*eval

      INTEGER EVAL, J
      double precision B(EVAL), V(J*EVAL)
      double precision MINI, MAXI, DIST, POWER
      INTEGER I, II

      DIST = 0.0
      POWER = 1.0
c     power in exp(-power*dist) for weighting
      DO 10 I=1,EVAL 
            MINI = V((I-1)*J+1)
            MAXI = V((I-1)*J+1)
            DO 5 II=1,J
            IF (MINI.GT.V((I-1)*J+II)) THEN
                  MINI = V((I-1)*J+II)
           ENDIF
            IF (MAXI.LT.V((I-1)*J+II)) THEN
                  MAXI = V((I-1)*J+II)
            ENDIF
5           CONTINUE
            IF  ((B(I).GE.MINI).AND.(B(I).LE.MAXI)) THEN
                  DIST = DIST + 0.0
            ELSE
                  IF(B(I).GT.MAXI) THEN
                  DIST = DIST + (B(I)-MAXI)**2
                  ENDIF
                  IF(B(I).LT.MINI) THEN
                  DIST = DIST + (MINI-B(I))**2
                  ENDIF
            ENDIF
10    CONTINUE
      DIST = DIST/(EVAL+0.0)
      DIST = EXP(-POWER*DIST)
      AdjLPindicator = DIST
      RETURN
      END

      SUBROUTINE AdjLP(EVAL,J,M,KOMB,COM,B,V,DJ)
c     adjusted band depth function, smoothing (exp(-u)), L2 metric
c     eval     no of time points for each functions
c     J        order of the depth
c     m        sample size
c     komb     m choose J number
c     com      vector of all the combinations of J elements from m
c     b        functional values of the x function, length eval
c     v        matrix of sample funct. values, dimension [m x eval]

      INTEGER EVAL,J,M,KOMB,COM(KOMB*J)
      double precision B(EVAL),V(M*EVAL),DJ
      double precision VPOM(J*EVAL)
      INTEGER C,I,II
      double precision ADJLPINDICATOR

      DJ = 0.0
      DO 10 C=1,KOMB
            DO 5 I=1,J
                  DO 3 II=1,EVAL
                        VPOM((II-1)*J+I) = V((II-1)*M+COM((C-1)*J+I))
3                 CONTINUE
5           CONTINUE 
            DJ = DJ + AdjLPindicator(EVAL,J,B,VPOM)     
10    CONTINUE
      DJ = DJ/(KOMB+0.0)
      RETURN
      END

c----------------------------------------------------------------------
      
      double precision FUNCTION AdjCindicator(EVAL,J,B,V)
c     adjusted band depth core function, smoothing (exp(-u)), C metric
c     b is a vector of length eval
c     v is a matrix j*eval

      INTEGER EVAL, J
      double precision B(EVAL), V(J*EVAL)
      double precision MINI, MAXI, DIST, POWER, MAX
      INTEGER I, II

      DIST = 0.0
      POWER = 1.0
c     power in exp(-power*dist) for weighting
      DO 10 I=1,EVAL 
            MINI = V((I-1)*J+1)
            MAXI = V((I-1)*J+1)
            DO 5 II=1,J
            IF (MINI.GT.V((I-1)*J+II)) THEN
                  MINI = V((I-1)*J+II)
           ENDIF
            IF (MAXI.LT.V((I-1)*J+II)) THEN
                  MAXI = V((I-1)*J+II)
            ENDIF
5           CONTINUE
            IF  ((B(I).GE.MINI).AND.(B(I).LE.MAXI)) THEN
                  DIST = DIST + 0.0
            ELSE
                  IF(B(I).GT.MAXI) THEN
                  DIST = MAX(DIST,(B(I)-MAXI))
                  ENDIF
                  IF(B(I).LT.MINI) THEN
                  DIST = MAX(DIST,(MINI-B(I)))
                  ENDIF
            ENDIF
10    CONTINUE
      DIST = EXP(-POWER*DIST)
      AdjCindicator = DIST
      RETURN
      END

      SUBROUTINE AdjC(EVAL,J,M,KOMB,COM,B,V,DJ)
c     adjusted band depth function, smoothing (exp(-u)), C metric
c     eval     no of time points for each functions
c     J        order of the depth
c     m        sample size
c     komb     m choose J number
c     com      vector of all the combinations of J elements from m
c     b        functional values of the x function, length eval
c     v        matrix of sample funct. values, dimension [m x eval]

      INTEGER EVAL,J,M,KOMB,COM(KOMB*J)
      double precision B(EVAL),V(M*EVAL),DJ
      double precision VPOM(J*EVAL)
      INTEGER C,I,II
      double precision ADJCINDICATOR

      DJ = 0.0
      DO 10 C=1,KOMB
            DO 5 I=1,J
                  DO 3 II=1,EVAL
                        VPOM((II-1)*J+I) = V((II-1)*M+COM((C-1)*J+I))
3                 CONTINUE
5           CONTINUE 
            DJ = DJ + AdjCindicator(EVAL,J,B,VPOM)     
10    CONTINUE
      DJ = DJ/(KOMB+0.0)
      RETURN
      END
      
c-----------------------------------
c     DiffDepth.f
c-----------------------------------

c-------------------------------------
c     first Difference Depth for 1D functions
c     for 1d functions, halfspace and simplicial depth
c-------------------------------------

      SUBROUTINE DiffD(A,B,M,N,D,REP,RN,funSDEP,funHDEP,funSDEPm,
     +funHDEPm,PSDEP,PHDEP,IAsdep,IAhdep)

c     using fdepth computes 2dimensional diff depth of functions in A
c     with respect to the functions in B
c     M     size of functions in A
c     N     size of functions in B
c     D     dimensionality of functions
c     REP   is the number of simulations to approximate the real value
c           of depth, if set to 0 full comutation is made
c     RN    random numbers, arrray of randoms from 1:D of size 2*REP
c     funSDEPm    DiffDepth when taking infimum of projected depths
c     funSDEP     DiffDepth when taking integral of projected depths
c     P.DEP pointwise depth, depth at each point of domain x domain
c     P.DEP is computed only if M=1, for simplicity
        
      INTEGER N,M,D,REP
      INTEGER IAsdep(M),IAhdep(M)
      double precision A(M*D),B(N*D),funSDEP(M),funHDEP(M),PSDEP(D*D),
     +PHDEP(D*D)
      double precision funSDEPm(M),funHDEPm(M),B1H(N),B2H(N)
      double precision hSDEP,hHDEP
      double precision hALPHA(N)
      INTEGER hF(N),I,J,K,RN(2*REP)

      DO 5 I=1,M
            funSDEP(I) = 0.0
            funHDEP(I) = 0.0
            funSDEPm(I) = 2.0
            funHDEPm(I) = 2.0
            IAsdep(I) = 0
            IAhdep(I) = 0
5     CONTINUE
c     initialization for easier recognition of unfilled elements
      IF (M.EQ.1) THEN 
            DO 8 I=1,(D*D)
                  PSDEP(I) = -1.0
                  PHDEP(I) = -1.0
8                 CONTINUE
            ENDIF
c     essential loop, computing 2D depths at each point
      IF (REP.EQ.0) THEN
c     if we compute the complete, non-approximated depth
      DO 30 I=1,D
            DO 25 J=(I+1),D
c     not taking into account diagonal elements, there is 0 depth
                  DO 10 K=1,N
                  B1H(K) = B((I-1)*N+K)
                  B2H(K) = B((J-1)*N+K)
10                CONTINUE
                  DO 20 L=1,M
                        hSDEP = 0.0
                        hHDEP = 0.0
                        hALPHA(1) = N+0.0
                        hF(1) = N
      CALL fD2(A((I-1)*M+L),A((J-1)*M+L),N,B1H,B2H,hALPHA,hF,
     +hSDEP,hHDEP)
                        funSDEP(L) = funSDEP(L) + hSDEP
                        funHDEP(L) = funHDEP(L) + hHDEP
c     counting the area of the smallest depth for each function
                        IF (hSDEP.EQ.funSDEPm(L)) THEN
                              IAsdep(L) = IAsdep(L)+1
                        ELSEIF (hSDEP.LT.funSDEPm(L)) THEN 
                              IAsdep(L) = 1
                        ENDIF
                        IF (hHDEP.EQ.funHDEPm(L)) THEN
                              IAhdep(L) = IAhdep(L)+1
                        ELSEIF (hHDEP.LT.funHDEPm(L)) THEN 
                              IAhdep(L) = 1
                        ENDIF 
                        funSDEPm(L) = min(funSDEPm(L),hSDEP)
                        funHDEPm(L) = min(funHDEPm(L),hHDEP)
                        IF (M.EQ.1) THEN
                              PSDEP((I-1)*D+J) = hSDEP
                              PHDEP((I-1)*D+J) = hHDEP
                              ENDIF
20                CONTINUE
25          CONTINUE
30    CONTINUE
c       dividing the resulting depths by number of points d
c     for infimal area, the result is *2+D because the diagonal 
c     has zero depth by default, and only the lower triangle 
c     of the matrix is actually computed
      DO 40 I=1,M
            funSDEP(I) = 2.0*funSDEP(I)/(D*(D-1.0))
c     diagonals are always zero, do not count into the sum 
c ((D+0.0)*(D+1.0)/2+0.0-D)
            funHDEP(I) = 2.0*funHDEP(I)/(D*(D-1.0))
c ((D+0.0)*(D+1.0)/2+0.0-D)
            IAhdep(I) = IAhdep(I)*2+D
            IAsdep(I) = IAsdep(I)*2+D
40    CONTINUE
      ELSE 
c     else of if(REP.EQ.0), that is if we approximate
      DO 70 I=1,REP
c     going through random elements
                  DO 50 K=1,N
                  B1H(K) = B((RN(2*I-1)-1)*N+K)
                  B2H(K) = B((RN(2*I)-1)*N+K)
50                CONTINUE
                  DO 60 L=1,M
                        hSDEP = 0.0
                        hHDEP = 0.0
                        hALPHA(1) = N+0.0
                        hF(1) = N
      CALL fD2(A((RN(2*I-1)-1)*M+L),A((RN(2*I)-1)*M+L),N,
     +B1H,B2H,hALPHA,hF,hSDEP,hHDEP)
                        funSDEP(L) = funSDEP(L) + hSDEP
                        funHDEP(L) = funHDEP(L) + hHDEP
c     counting the area of the smallest depth for each function
                        IF (hSDEP.EQ.funSDEPm(L)) THEN
                              IAsdep(L) = IAsdep(L)+1
                        ELSEIF (hSDEP.LT.funSDEPm(L)) THEN 
                              IAsdep(L) = 1
                        ENDIF
                        IF (hHDEP.EQ.funHDEPm(L)) THEN
                              IAhdep(L) = IAhdep(L)+1
                        ELSEIF (hHDEP.LT.funHDEPm(L)) THEN 
                              IAhdep(L) = 1
                        ENDIF 
                        funSDEPm(L) = min(funSDEPm(L),hSDEP)
                        funHDEPm(L) = min(funHDEPm(L),hHDEP)
                        IF (M.EQ.1) THEN
                              PSDEP((RN(2*I-1)-1)*D+RN(2*I)) = hSDEP
                              PHDEP((RN(2*I-1)-1)*D+RN(2*I)) = hHDEP
                              ENDIF
60                CONTINUE
70    CONTINUE
      DO 80 I=1,M
            funSDEP(I) = funSDEP(I)/(REP+0.0) 
            funHDEP(I) = funHDEP(I)/(REP+0.0)
80    CONTINUE      
      ENDIF
      RETURN
      END 
      
c----------------------------------------------------------------------

      SUBROUTINE fD2(U,V,N,X,Y,ALPHA,F,SDEP,HDEP)

c        Rousseuw, P.J., and Ruts, I. (1996), AS 307 : Bivariate location
c        depth, Applied Statistics (JRRSS-C), vol.45, 516-526

      double precision U,V,X(n),Y(n),ALPHA(n)
      double precision P,P2,EPS,D,XU,YU,ANGLE,ALPHK,BETAK,SDEP,HDEP
      INTEGER F(N),GI
      integer n,nums,numh,nt,i,nn,nu,ja,jb,nn2,nbad,nf,j,ki,k
      NUMS=0
      NUMH=0
      SDEP=0.0
      HDEP=0.0
      IF (N.LT.1) RETURN
      P=ACOS(-1.0)
      P2=P*2.0
      EPS=0.00000001
      NT=0
C
C  Construct the array ALPHA.
C
      DO 10 I=1,N
          D=SQRT((X(I)-U)*(X(I)-U)+(Y(I)-V)*(Y(I)-V))
          IF (D.LE.EPS) THEN
              NT=NT+1
          ELSE
              XU=(X(I)-U)/D
              YU=(Y(I)-V)/D
              IF (ABS(XU).GT.ABS(YU)) THEN
                  IF (X(I).GE.U) THEN
                      ALPHA(I-NT)=ASIN(YU)
                      IF(ALPHA(I-NT).LT.0.0) THEN
                          ALPHA(I-NT)=P2+ALPHA(I-NT)
                      ENDIF
                  ELSE
                      ALPHA(I-NT)=P-ASIN(YU)
                  ENDIF
              ELSE
                  IF (Y(I).GE.V) THEN
                      ALPHA(I-NT)=ACOS(XU)
                  ELSE
                      ALPHA(I-NT)=P2-ACOS(XU)
                  ENDIF
              ENDIF
              IF (ALPHA(I-NT).GE.(P2-EPS)) ALPHA(I-NT)=0.0
          ENDIF
  10  CONTINUE
      NN=N-NT
      IF (NN.LE.1) GOTO 60
C
C  Sort the array ALPHA.
C
      CALL SORT(ALPHA,NN)
C
C  Check whether theta=(U,V) lies outside the data cloud.
C
      ANGLE=ALPHA(1)-ALPHA(NN)+P2
      DO 20 I=2,NN
          ANGLE=MAX(ANGLE,(ALPHA(I)-ALPHA(I-1)))
  20  CONTINUE
      IF (ANGLE.GT.(P+EPS)) GOTO 60
C
C  Make smallest alpha equal to zero,
C  and compute NU = number of alpha < pi.
C
      ANGLE=ALPHA(1)
      NU=0
      DO 30 I=1,NN
          ALPHA(I)=ALPHA(I)-ANGLE
          IF (ALPHA(I).LT.(P-EPS)) NU=NU+1
  30  CONTINUE
      IF (NU.GE.NN) GOTO 60
C
C  Mergesort the alpha with their antipodal angles beta,
C  and at the same time update I, F(I), and NBAD.
C
      JA=1
      JB=1
      ALPHK=ALPHA(1)
      BETAK=ALPHA(NU+1)-P
      NN2=NN*2
      NBAD=0
      I=NU
      NF=NN
      DO 40 J=1,NN2
          IF ((ALPHK+EPS).LT.BETAK) THEN
              NF=NF+1
              IF (JA.LT.NN) THEN
                  JA=JA+1
                  ALPHK=ALPHA(JA)
              ELSE
                  ALPHK=P2+1.0
              ENDIF
          ELSE
              I=I+1
              IF (I.EQ.(NN+1)) THEN
                  I=1
                  NF=NF-NN
              ENDIF
              F(I)=NF
              NBAD=NBAD+K((NF-I),2)
              IF (JB.LT.NN) THEN
                  JB=JB+1
                  IF ((JB+NU).LE.NN) THEN
                      BETAK=ALPHA(JB+NU)-P
                  ELSE
                      BETAK=ALPHA(JB+NU-NN)+P
                  ENDIF
              ELSE
                  BETAK=P2+1.0
              ENDIF
          ENDIF
  40  CONTINUE
      NUMS=K(NN,3)-NBAD
C
C  Computation of NUMH for halfspace depth.
C
      GI=0
      JA=1
      ANGLE=ALPHA(1)
      NUMH=MIN(F(1),(NN-F(1)))
      DO 50 I=2,NN
          IF(ALPHA(I).LE.(ANGLE+EPS)) THEN
              JA=JA+1
          ELSE
              GI=GI+JA
              JA=1
              ANGLE=ALPHA(I)
          ENDIF
          KI=F(I)-GI
          NUMH=MIN(NUMH,MIN(KI,(NN-KI)))
   50 CONTINUE
C
C  Adjust for the number NT of data points equal to theta:
C
  60  NUMS=NUMS+K(NT,1)*K(NN,2)+K(NT,2)*K(NN,1)+K(NT,3)
      IF (N.GE.3) SDEP=(NUMS+0.0)/(K(N,3)+0.0)
      NUMH=NUMH+NT
      HDEP=(NUMH+0.0)/(N+0.0)
      RETURN
      END
c----------------------------------------------------------------------


c      INTEGER FUNCTION K(M,J)
c      integer m,j
c      IF (M.LT.J) THEN
c          K=0
c      ELSE
c          IF (J.EQ.1) K=M
c          IF (J.EQ.2) K=(M*(M-1))/2
c          IF (J.EQ.3) K=(M*(M-1)*(M-2))/6
c      ENDIF
c      RETURN
c      END

c----------------------------------------------------------------------

      SUBROUTINE SORT(B,N)
C  Sorts an array B (of length N<=1000) in O(NlogN) time.

      double precision B(N),x(n)
      integer q(n),i,n
      
      call indexx(n,b,q)
      
      do 10 i=1,n
        x(i)=b(i)
10    continue
 
      do 20 i=1,n
        b(i)=x(q(i))
20    continue
      end
c----------------------------------------------------------------------

      SUBROUTINE INDEXX(N,ARRIN,INDX)

      integer INDX(N)
      integer n,j,l,ir,indxt,i
      double precision arrin(n),q

      DO 11 J=1,N
        INDX(J)=J
11    CONTINUE
      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          INDXT=INDX(L)
          Q=ARRIN(INDXT)
        ELSE
          INDXT=INDX(IR)
          Q=ARRIN(INDXT)
          INDX(IR)=INDX(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            INDX(1)=INDXT
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(ARRIN(INDX(J)).LT.ARRIN(INDX(J+1)))J=J+1
          ENDIF
          IF(Q.LT.ARRIN(INDX(J)))THEN
            INDX(I)=INDX(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        INDX(I)=INDXT
      GO TO 10
      END
      
cccccccccccccccccc
c
c     newdepth
c
cccccccccccccccccc

c----------------------------------------------------------------------
      SUBROUTINE funMD(A,B,M,N,D,Q,funDEP)
      
      INTEGER N,M,D,I,J      
      double precision A(M*D),B(N*D),METR1(N*N),METR2(N*M),funDEP(M)
      double precision Q,H,PI
c     Q is quantile, default 0.15, Cuevas 0.2     
      
      CALL metrl2(B,B,N,N,D,METR1)
      CALL metrl2(A,B,M,N,D,METR2)
      CALL SORT(METR1,N*N)
      H = METR1(INT(Q*(N*N+0.0)))
      PI=4.D0*DATAN(1.D0)
      DO 10 I=1,M*N
            METR2(I) = EXP(-(METR2(I)/H)**2/2.0)/SQRT(2.0*PI)
10    CONTINUE
      DO 20 I=1,M
            funDEP(I)=0.0
            DO 15 J=1,N
                  funDEP(I) = funDEP(I)+METR2((J-1)*M+I)
15          CONTINUE
20    CONTINUE      
      RETURN
      END
      

c----------------------------------------------------------------------
      SUBROUTINE funRPD1(A,B,M,N,D,NPROJ,NPROJ2,V,funSDEP,
     +funHDEP,funRDEP)
      
c     computes fast random projection depth for 1D functions using 
c     NPROJ simple iid gaussian processes projection and then 1D
c     simplicial or halfspace depth 
c     RDEP random halfspace depth (Cuesta-Albertos, Nieto-Reyes 2008)
c     NPROJ2 nr of projections for RDEP

      INTEGER M,N,D,NPROJ,NPROJ2
      double precision A(M*D),B(N*D),funSDEP(M),funHDEP(M),funRDEP(M)
      double precision V(D*NPROJ),AP,BP(N),VNORM
      INTEGER I,J,K
      double precision hSDEP,hHDEP

      DO 5 I=1,M
            funSDEP(I) = 0.0
            funHDEP(I) = 0.0
            funRDEP(I) = 2.0
5     CONTINUE      
c     main loop
      DO 40 K=1,NPROJ
c     generating a single gaussian processes
      VNORM = 0.0
      DO 10 I=1,D
             Vnorm = VNORM + V((K-1)*D+I)**2
10    CONTINUE
      VNORM = sqrt(VNORM - (V((K-1)*D+1)**2+V((K-1)*D+D)**2)/(2.0))
c     loop projecting j-th function of the random sample B on v
      DO 18 J=1,N
            BP(J) = 0.0
            DO 15 I=1,D
                  BP(J) = BP(J)+ B((I-1)*N+J)*V((K-1)*D+I)/VNORM
15          CONTINUE
            BP(J) = BP(J)/(D+0.0)
18    CONTINUE
c     loop with projecting j-th function of the random sample A on v
      DO 30 J=1,M
            AP = 0.0
c     compute projection of j-th function A
            DO 20 I=1,D
                  AP = AP + A((I-1)*M+J)*V((K-1)*D+I)/VNORM
20          CONTINUE
            AP = AP/(D+0.0)
            hSDEP = 0.0
            hHDEP = 0.0
            CALL fD1(AP,N,BP,hSDEP,hHDEP)
            funSDEP(J) = funSDEP(J) + hSDEP
            funHDEP(J) = funHDEP(J) + hHDEP
            IF (K .LE. NPROJ2) THEN 
                  funRDEP(J) = min(funRDEP(J),hHDEP)
                  ENDIF
30    CONTINUE
40    CONTINUE
c     averaging projection depth over all projections
      DO 50 I=1,M
            funSDEP(I) = funSDEP(I)/(NPROJ+0.0)
            funHDEP(I) = funHDEP(I)/(NPROJ+0.0)
50    CONTINUE
      RETURN
      END
      
c-------------------------------------
c Fraiman Muniz type integrated depths
c for 1 and 2d functions, halfspace and simplicial depth
c-------------------------------------

      SUBROUTINE funD2(A1,A2,B1,B2,M,N,D,funSDEP,funHDEP,
     +fIsdep,fIhdep,IAsdep,IAhdep)

c     using fdepth computes 2dimensional depth of functions in A1,A2
c     with respect to the functions in B1,B2
c     M size of functions in A
c     N size of functions in B
c     D dimensionality of functions
        
      INTEGER N,M,D
      INTEGER IAsdep(M),IAhdep(M)
      double precision A1(M*D),A2(M*D),B1(N*D),B2(N*D),funSDEP(M),
     +funHDEP(M),fIsdep(M),fIhdep(M) 
      double precision B1H(N),B2H(N)
      double precision hSDEP,hHDEP
      double precision hALPHA(N)
      INTEGER hF(N),I,J,K

      DO 5 I=1,M
            funSDEP(I) = 0.0
            funHDEP(I) = 0.0
            fISDEP(I) = 2.0
            fIHDEP(I) = 2.0 
            IAsdep(I) = 0
            IAhdep(I) = 0           
5     CONTINUE
c       essential loop, computing 2D depths at each point
      DO 30 I=1,D
            DO 10 K=1,N
                  B1H(K) = B1((I-1)*N+K)
                  B2H(K) = B2((I-1)*N+K)
10          CONTINUE
            DO 20 J=1,M
                  hSDEP = 0.0
                  hHDEP = 0.0
                  hALPHA(1) = N+0.0
                  hF(1) = N
      CALL fD2(A1((I-1)*M+J),A2((I-1)*M+J),N,B1H,B2H,hALPHA,hF,
     +hSDEP,hHDEP)
                  funSDEP(J) = funSDEP(J) + hSDEP
                  funHDEP(J) = funHDEP(J) + hHDEP
c     counting the area of the smallest depth for each function
                  IF (hSDEP.EQ.fISDEP(J)) THEN
                        IAsdep(J) = IAsdep(J)+1
                  ELSEIF (hSDEP.LT.fISDEP(J)) THEN 
                        IAsdep(J) = 1
                  ENDIF
                  IF (hHDEP.EQ.fIHDEP(J)) THEN
                        IAhdep(J) = IAhdep(J)+1
                  ELSEIF (hHDEP.LT.fIHDEP(J)) THEN 
                        IAhdep(J) = 1
                  ENDIF 
c     infimal depth evaluation
                  fISDEP(J) = min(fISDEP(J),hSDEP)
                  fIHDEP(J) = min(fIHDEP(J),hHDEP)                  
20          CONTINUE
30    CONTINUE
c       dividing the resulting depths by number of points d
      DO 40 I=1,M
            funSDEP(I) = funSDEP(I)/(D+0.0)
            funHDEP(I) = funHDEP(I)/(D+0.0)
40    CONTINUE
      RETURN
      END 
      
c----------------------------------------------------------------------
c     2D RANDOM PROJECTION DEPTH
c     PROJECTIONS ON WHITE NOISE
c     DEPENDS ON THE RANDOM NUMBER GENERATOR VERSION
c----------------------------------------------------------------------
      
      SUBROUTINE funRPD2(A1,A2,B1,B2,M,N,D,NPROJ,V,Q,funSDEP,
     +funHDEP,funMDEP,funSDDEP,funHDDEP)
      
c     computes fast random projection depth for 2D functions using 
c     NPROJ simple iid gaussian processes projection and then 2D
c     simplicial or halfspace depth 
c     Q is quantile used in MDEP
c     funDDEP is double random projection (2*D-dim)->(2-dim)->(1-dim)
c     and then appying Halfspace depth (equiv to quantile as Cuevas 2007)

      INTEGER M,N,D,NPROJ
      double precision A1(M*D),A2(M*D),B1(N*D),B2(N*D),funSDEP(M),
     +funHDEP(M),funMDEP(M),Q,funSDDEP(M),funHDDEP(M)
      double precision V(D*NPROJ+2*NPROJ),A1P,A2P,B1P(N),B2P(N),
     +VNORM,VNORMF,AP(2*M),BP(2*N),VF(2),AFP,BFP(N)
c     last 2*NPROJ of V are for the second random projection 2d->1d
c     VF(2) is the final projection
      INTEGER I,J,K
      double precision hSDEP,hHDEP,hMDEP(M),hHDDEP,hSDDEP
      double precision hALPHA(N)
      INTEGER hF(N)

      DO 5 I=1,M
            funSDEP(I) = 0.0
            funHDEP(I) = 0.0
            funMDEP(I) = 0.0
            funSDDEP(I) = 0.0
            funHDDEP(I) = 0.0
5     CONTINUE      
c     main loop
      DO 40 K=1,NPROJ
c     generating a single gaussian processes
      VNORM = 0.0
      VNORMF = sqrt(V(NPROJ*D+(K-1)*2+1)**2+V(NPROJ*D+(K-1)*2+2)**2)
      VF(1) = V(NPROJ*D+(K-1)*2+1)/VNORMF
      VF(2) = V(NPROJ*D+(K-1)*2+2)/VNORMF
c     final projection is now normed
      DO 10 I=1,D
             VNORM = VNORM + V((K-1)*D+I)**2
10    CONTINUE
      VNORM = sqrt(VNORM - (V((K-1)*D+1)**2+V((K-1)*D+D)**2)/(2.0))
c     loop projecting j-th function of the random sample B on v
      DO 18 J=1,N
            B1P(J) = 0.0
            B2P(J) = 0.0
            DO 15 I=1,D
                  B1P(J) = B1P(J)+ B1((I-1)*N+J)*V((K-1)*D+I)/VNORM
                  B2P(J) = B2P(J)+ B2((I-1)*N+J)*V((K-1)*D+I)/VNORM
15          CONTINUE
            BFP(J) = B1P(J)*VF(1)+B2P(J)*VF(2)
18    CONTINUE
      DO 19 I=1,N
            BP(I)=B1P(I)
            BP(N+I)=B2P(I)
19    CONTINUE
c     loop with projecting j-th function of the random sample A on v
      DO 30 J=1,M
            A1P = 0.0
            A2P = 0.0
c     compute projection of j-th function A
            DO 20 I=1,D
                  A1P = A1P + A1((I-1)*M+J)*V((K-1)*D+I)/VNORM
                  A2P = A2P + A2((I-1)*M+J)*V((K-1)*D+I)/VNORM
20          CONTINUE
            AP(J) = A1P
            AP(J+M) = A2P
            hSDEP = 0.0
            hHDEP = 0.0
            hMDEP(J) = 0.0
            hALPHA(1) = N+0.0
            hF(1) = N
            CALL fD2(A1P,A2P,N,B1P,B2P,hALPHA,hF,
     +hSDEP,hHDEP)
            funSDEP(J) = funSDEP(J) + hSDEP
            funHDEP(J) = funHDEP(J) + hHDEP
c     and now the second projection of A
            AFP = A1P*VF(1)+A2P*VF(2)
            hHDDEP = 0.0
            hSDDEP = 0.0
            CALL fD1(AFP,N,BFP,hSDDEP,hHDDEP)
            funSDDEP(J) = funSDDEP(J) + hSDDEP
            funHDDEP(J) = funHDDEP(J) + hHDDEP
30    CONTINUE
      CALL funMD(AP,BP,M,N,2,Q,hMDEP)
      DO 35 I=1,M
            funMDEP(I) = funMDEP(I) + hMDEP(I)
35    CONTINUE
40    CONTINUE
c     averaging projection depth over all projections
      DO 50 I=1,M
            funSDEP(I) = funSDEP(I)/(NPROJ+0.0)
            funHDEP(I) = funHDEP(I)/(NPROJ+0.0)
            funMDEP(I) = funMDEP(I)/(NPROJ+0.0)
            funSDDEP(I) = funSDDEP(I)/(NPROJ+0.0)
            funHDDEP(I) = funHDDEP(I)/(NPROJ+0.0)
50    CONTINUE
      RETURN
      END
      
c
c
c     halfregiondepth.f
c
c

      SUBROUTINE HRD(A,B,M,N,D,FD)
c     computes fast half-region depth of Lopez-Pintado and Romo 2011
      INTEGER M,N,D
      double precision A(M*D),B(N*D),FD(M)
      INTEGER I,J,K, U, L, UI, LI

      DO 30 I=1,M
            FD(I) = 0.0
            U = 0
            L = 0
            DO 20 J=1,N
                  UI = 0
                  LI = 0
                  K = 0
                  DO 10 WHILE ((K .LT. D) 
     +                  .AND. (UI .EQ. 0 .OR. LI .EQ. 0)) 
                        K = K + 1
                        IF(A((K-1)*M+I) .GT. B((K-1)*N+J)) UI = UI + 1
                        IF(A((K-1)*M+I) .LT. B((K-1)*N+J)) LI = LI + 1
10                CONTINUE
                  IF (UI .EQ. 0) U = U+1
                  IF (LI .EQ. 0) L = L+1
20          CONTINUE
            FD(I) = (MIN(U,L)+0.0)/(N+0.0)             
30    CONTINUE
      RETURN
      END 

c----------------------------------------------------------------------      
      SUBROUTINE BD(A,B,M,N,D,FD)
c     computes fast 2nd order band depth of Lopez-Pintado and Romo 2009
      INTEGER M,N,D
      double precision A(M*D),B(N*D),FD(M),L,U
      INTEGER I,J,JJ,K,W,WI

      DO 30 I=1,M
            FD(I) = 0.0
            W = 0
            DO 20 J=1,(N-1)
                  DO 15 JJ=(J+1),N
                        WI = 0
                        K = 0
                  DO 10 WHILE ( WI .GE. 0) 
c     -1 means TRUE, the function is in the band of J and JJ
c     -2 means FALSE, the function crosses
                        WI = WI + 1
                        K = K+1
                        L = MIN(B((K-1)*N+J),B((K-1)*N+JJ))
                        U = MAX(B((K-1)*N+J),B((K-1)*N+JJ))
                        IF (A((K-1)*M+I) .LT. L .OR. 
     +  A((K-1)*M+I) .GT. U) WI = -2
                        IF (WI .NE. -2 .AND. WI .EQ. D) WI = -1
10                CONTINUE
                  IF (WI .EQ. -1) W = W+1
15                CONTINUE
20          CONTINUE
            FD(I) = (W+0.0)/((N*(N-1))/2 + 0.0)          
30    CONTINUE
      RETURN
      END 
      
c-------------------------------------
c bivariate halfspace and simplicial depth
c-------------------------------------

      SUBROUTINE DPTH2(A1,A2,B1,B2,M,N,SDEP,HDEP)

c     computes 2dimensional depth of A1,A2
c     with respect to B1,B2
c     M size of A
c     N size of B
        
      INTEGER N,M
      double precision A1(M),A2(M),B1(N),B2(N),SDEP(M),HDEP(M)
      double precision hSDEP,hHDEP
      double precision hALPHA(N)
      INTEGER hF(N),I,J,K

      DO 5 I=1,M
            SDEP(I) = 0.0
            HDEP(I) = 0.0
5     CONTINUE
      DO 10 J=1,M
            hSDEP = 0.0
            hHDEP = 0.0
            hALPHA(1) = N+0.0
            hF(1) = N
      CALL fD2(A1(J),A2(J),N,B1,B2,hALPHA,hF,
     +hSDEP,hHDEP)
            SDEP(J) = hSDEP
            HDEP(J) = hHDEP
10    CONTINUE
      RETURN
      END 
      
c-------------------------------------
c univariate halfspace and simplicial depth
c-------------------------------------

      SUBROUTINE DPTH1(A1,B1,M,N,SDEP,HDEP)

c     computes 1dimensional depth of A1
c     with respect to B1
c     M size of A
c     N size of B
        
      INTEGER N,M
      double precision A1(M),B1(N),SDEP(M),HDEP(M)
      double precision hSDEP,hHDEP
      double precision hALPHA(N)
      INTEGER hF(N),I,J,K

      DO 5 I=1,M
            SDEP(I) = 0.0
            HDEP(I) = 0.0
5     CONTINUE
      DO 10 J=1,M
            hSDEP = 0.0
            hHDEP = 0.0
            hALPHA(1) = N+0.0
            hF(1) = N
      CALL fD1(A1(J),N,B1,hSDEP,hHDEP)
            SDEP(J) = hSDEP
            HDEP(J) = hHDEP
10    CONTINUE
      RETURN
      END 
