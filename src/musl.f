      SUBROUTINE MUSL ( FLIN, FDIF, N, IHOM, A, B, MA, MB, BCV, AMP,
     &  ER, NRTI, TI, NTI, X, U, NU, Q, D, KPART, PHIREC, W, LW, IW,
     &  LIW, IERROR )

c*********************************************************************72
c
cc MUSL solves a linear two point boundary value problem.
c
c  Modified:
c
c    10 March 2013
c
c  Author:
c
c    Robert Mattheij, G.W.M. Staarink.
c
c  Reference:
c
c    Uri Ascher, Robert Mattheij, Robert Russell,
c    Numerical Solution of Boundary Value Problems for 
c    Ordinary Differential Equations,
c    Prentice Hall, 1988,
c    ISBN: 0-13-627266-5,
c    LC: QA379.A83.
c
c  Parameters:
c
c    Input, external FLIN, the name of a subroutine which evaluates the
c    homogeneous part of the differential equation, A(X)*Y(X), placing
c    the result in F.  FLIN has the form
c      subroutine flin ( x, y, f )
c      double precision x, y(n), f(n)
c
c    Input, external FDIF, the name of a subroutine which evaluates the
c    inhomogeneous part of the differential equation, A(X)*Y(X)+Q(X), placing
c    the result in F.  FDIF has the form
c      subroutine fdif ( x, y, f )
c      double precision x, y(n), f(n)
c
c    Input, integer N, the order of the system.
c
c    Input, integer IHOM, indicates whether the system is homogeneous or not:
c    0, homogeneous;
c    1, nonhomogeneous.
c
c    Input, double precision A, B, the endpoints of the interval.
c
c    Input, double precision MA(N,N), MB(N,N), boundary condition arrays.
c
c    Input, double precision BCV(N), the boundary values.  The boundary 
c    conditions are defined by MA*Y(A)+MB*Y(B)=BCV.
c
c    Input, double precision AMP, the allowed incremental factor of the
c    homogeneous solutions between successive output points.  If AMP <= 1,
c    the defaults are:
c    NRTI = 0: max(ER(1),ER(2)/ER(3));
c    NRTI nonzero: infinity.
c
c    Input/output, double precision ER(5).  On input, ER(1) and ER(2) contain
c    the relative and absolute error tolerances, respectively, and ER(3)
c    contains the machine precision constant.  On output, ER(4) contains
c    an estimate of the condition of the boundary value problem, and ER(5)
c    contains an estimated error amplification factor.
c
c    Input/output, integer NRTI, determines the output points.
c    0, the subroutine determines the output points using AMP;
c    1, the user supplies NTI output points in TI.
c    otherwise, the program computes NRTI+1 output points, equally spaced
c    between A and B, inclusive.
c    On output NRTI contains the number of output points.
c
c    Input/output, double precision TI(NTI).  If NRTI = 1, then on input,
c    TI contains the points in [A,B] where the solution is desired, in 
c    increasing order.  Otherwise, the values of TI are determined by 
c    the program.
c
c    Input, integer NTI, the dimension of TI, and one of the dimensions of
c    Y, U, Q, D, and PHI.  NTI must be greater than the total number of
c    output points plus 3.
c
c    Output, double precision Y(N,NTI), contains the solution of the 
c    boundary value problem at each output point TI().
c
c    Output, double precision U(NU,NTI), contains the relevant elements of
c    the upper triangular matrix Uk.
c
c    Input, integer NU, used as one of the dimensions of U and PHI.
c    NU must be at least N*(N+1)/2.
c
c    Output, double precision Q(N,N,NTI), contains the orthogonal matrix
c    Qk associated with each shooting point.
c
c    Output, double precision D(N,NTI).  If IHOM = 1, then on exit, for
c    D(I,K) contains the I-th inhomogeneous term at the K-th
c    shooting points. 
c
c    Output, integer KPART, contains the order or size of the left upper
c    blocks of the matrices Uk which correspond to the globally increasing modes.
c
c    Output, double precision PHI(NU,NTI), an upper triangular fundamental
c    solution of the multiple shooting recursion.
c
c    Workspace, double precision W(LW).
c
c    Input, integer LW, the dimension of W, which must be at least
c    8*N+2*N*N.
c
c    Workspace, integer IW(LIW).
c
c    Input, integer LIW, the dimension of IW, which must be at least 3*N.
c
c    Output, integer IERROR, error indicator.  A return value of 0 means
c    no error was detected.
c
      implicit double precision (A-H,O-Z)

      integer liw
      integer lw
      integer n
      integer nti
      integer nu

      double precision bcv(n)
      double precision D(N,NTI)
      double precision er(5)
      external fdif
      external flin
      integer i
      integer i2
      integer i3
      integer ierror
      integer ihom
      integer iw(liw)
      integer j
      integer j2
      integer j3
      integer j4
      double precision MA(N,N)
      double precision MB(N,N)
      double precision PHIREC(NU,NTI)
      double precision Q(N,N,NTI)
      double precision TI(NTI)
      double precision U(NU,NTI)
      double precision W(LW)
      double precision X(N,NTI)

      I2 = 1 + N
      I3 = I2 + N
      J2 = 1 + N
      J3 = J2 + 7*N
      J4 = J3 + N*N
      IW(I2) = 0
      IERROR = 0
      I = N * (N+1) / 2

      IF ((N.LT.1).OR.(IHOM.LT.0).OR.(NRTI.LT.0).OR.(NTI.LT.5).OR.
     1 (NU.LT.I).OR.(A.EQ.B)) THEN
        IERROR = 100
        WRITE(*,100) IERROR
        CALL ERRHAN(IERROR,ER,0)
        RETURN
      end if

      IF ((ER(1).LT.0.D0).OR.(ER(2).LT.0.D0).OR.(ER(3).LT.0.D0)) THEN
        IERROR = 101
        WRITE(*,100) IERROR
        CALL ERRHAN(IERROR,ER,0)
        RETURN
      end if

      IF ((LW.LT.8*N+2*N*N).OR.(LIW.LT.3*N)) THEN
        IERROR = 103
        WRITE(*,100) IERROR
        CALL ERRHAN(IERROR,ER,0)
        RETURN
      end if
c
C  SET Q(I) = I.
c
      DO I = 1, N
        DO J = 1, N
          Q(I,J,1) = 0.0D+00
        end do
        Q(I,I,1) = 1.0D+00
      end do

      CALL DDUR(FLIN,FDIF,N,IHOM,A,B,NRTI,AMP,TI,NTI,ER,Q,U,NU,D,D,
     1          KPART,W(1),W(J2),W(J3),W(J4),IW(1),IW(I2),IW(I3),IERROR)

      IF (IERROR.eq. 200 .or. ierror .eq. 213 ) THEN
        WRITE(*,110) IERROR
        CALL ERRHAN(IERROR,ER,0)
      else IF ( IERROR.NE.0 ) then
        WRITE(*,100) IERROR
        CALL ERRHAN(IERROR,ER,0)
        RETURN
      end if

      IW(I2) = 0
      J5 = J2+N

      CALL DGTUR(N,IHOM,NRTI,U,NU,NTI,Q,D,ER,IW(1),KPART,IW(I2),
     1           IW(I3),W(1),W(J2),W(J5),W(J3),W(J4),IERROR)

      CALL DKPCH(N,U,NU,NTI,NRTI,KPART,IW(1),ER,IERROR)

      IF (IERROR.NE.0) THEN
        WRITE(*,110) IERROR
        CALL ERRHAN(IERROR,ER,0)
      end if

      CALL DFUNRC(N,KPART,NRTI,U,NU,NTI,IW(1),PHIREC,W(1),W(J2),IERROR)

      IF (IERROR.NE.0) then
        WRITE(*,100) IERROR
        CALL ERRHAN(IERROR,ER,0)
        RETURN
      end if

      IF ( IHOM .NE. 0 ) THEN

        CALL DPSR(N,KPART,NRTI,U,NU,NTI,D,IW(1),X,W(1),W(J2),IERROR)

        IF (IERROR.NE.0) then
          WRITE(*,100) IERROR
          CALL ERRHAN(IERROR,ER,0)
          RETURN
        end if

      end if

      CALL DSBVP(N,IHOM,MA,MB,NRTI,Q,NTI,BCV,PHIREC,NU,X,X,ER(3),
     1           ER(4),IW(1),IW(I2),W(J3),W(J4),W(1),W(J2),IERROR)

      IF (IERROR.NE.0) then
        WRITE(*,100) IERROR
        CALL ERRHAN(IERROR,ER,0)
        RETURN
      end if

      RETURN

  100 FORMAT(' TERMINAL ERROR IN MUSL : IERROR =',I4)
  110 FORMAT(' WARNING ERROR IN MUSL : IERROR =',I4)
      END
      SUBROUTINE DCSHPO(FDIF,X0T,N,A,B,ER,TI,NTI,NRTI,ALINC,X,S,Q,U,NU,
     1                  W,WTI,WF,WF0,KKK,IP,WGR,LWGR,IERROR)

c*********************************************************************72
c
cc DCSHPO determines the shooting points.
C
      implicit double precision (A-H,O-Z)

      DIMENSION ER(5),TI(NTI),ALINC(4),X(N,NTI),S(N,NTI),Q(N,N,NTI),
     1          U(NU,NTI),W(N,7),WTI(NTI),WF(N,N),WF0(N,N),WGR(LWGR)
      INTEGER KKK(N),IP(N)
      EXTERNAL FDIF,X0T

      IERROR = 0
      EPS = DSQRT(ER(3))
      TI(1) = A
      WTI(1) = A
      NOTI = 2
      IF (NRTI.EQ.0) TI(2) = B
      IF (NRTI.EQ.1) THEN
        IF (A.LE.B) THEN
          DO 1100 I = 2 , NTI
            IF (TI(I-1).GE.TI(I)) GOTO 5200
            IF (TI(I).GE.B) GOTO 1200
 1100     CONTINUE
 1200     NOTI=I
        ELSE
          DO 1300 I = 2 , NTI
            IF (TI(I-1).LE.TI(I)) GOTO 5200
            IF (TI(I).LE.B) GOTO 1400
 1300     CONTINUE
 1400     NOTI=I
        end if
        IF (TI(NOTI).NE.B) GOTO 5300
      end if
      IF (NRTI.GT.1) THEN
        NOTI = NRTI + 1
        ST = (B-A) / NRTI
        TI(NRTI+1) = B
        DO I = 2 , NRTI
          TI(I) = TI(I-1) + ST
        end do
      end if
      IF (LWGR.LT.NOTI) GOTO 5550
C     DETERMINATION OF THE SHOOTING POINTS. GIVEN OUTPUTPOINTS WILL
C     BE SHOOTING POINTS. IF THE INCREMENT BETWEEN TWO OUTPUTPOINTS
C     IS GREATER THAN ALINC(1), A NEW SHOOTING POINT IS INSERTED.
C     DURING THE DETERMINATION THE SHOOTING POINTS WILL BE STORED IN
C     THE ARRAY WTI.
      DO I = 1 , N
        KKK(I) = I*(I+1) / 2
        DO J = 1 , N
          IF (I.EQ.J) THEN
            Q(I,I,1) = 1.D0
          ELSE
            Q(I,J,1) = 0.D0
          end if
          Q(I,J,NTI) = Q(I,J,1)
        end do
      end do
      WTI(1) = TI(1)
      JTI = 2
      NHI = 5
      NRTI = 2
      MWGR = 0
      DO 2600 K = 2 , NTI
        KM1 = K - 1
        T1 = WTI(KM1)
        T2 = TI(JTI)
 1650   CALL X0T(T1,X(1,KM1))
        IF (T1.EQ.B) GOTO 2700
        X1 = 0.D0
        DO I = 1 , N
          ST = DABS(X(I,KM1))
          S(I,K) = X(I,KM1)
          IF (ST.GT.X1) X1 = ST
          DO J = 1 , N
            WF(J,I) = X(J,KM1) + EPS * Q(J,I,KM1)
          end do
        end do

        IF (X1.EQ.0.D0) X1 = 1.D0
        IWGR = 0
        IER = 1
        IF (K.GT.2) IER=-1
 1900   CALL DRKFMS(FDIF,N,S(1,K),T1,T2,ER(1),ER(2),ER(3),WF,N,HI,NHI,
     1              W,IER)
        IF (IER.GT.3) GOTO 5400
C     CHECK WHETHER T2 SHOULD BE A SHOOTING POINT
        X2 = 0.D0
        X3 = 0.D0
        DO I = 1 , N
          ST = DABS(S(I,K))
          IF (ST.GT.X2) X2 = ST
          DO J = 1 , N
            ST = DABS(WF(J,I)-S(J,K))
            IF (ST.GT.X3) X3 = ST
          end do
        end do

        IF ((X2/X1.LT.ALINC(1)).AND.(X3.LT.ALINC(3))) THEN
          IWGR = IWGR + 1
          MWGR = MWGR + 1
          IF (MWGR.GT.LWGR) GOTO 5600
          T1 = T2
          WGR(MWGR) = T2
          IER = -1
          IF (T2.EQ.TI(JTI)) GOTO 2040
          T2 = TI(JTI)
          GOTO 1900
        end if
        IF ((X2/X1.LE.ALINC(2)).AND.(X3.LE.ALINC(4))) THEN
          MWGR = MWGR + 1
          IF (MWGR.GT.LWGR) GOTO 5600
          IWGR = IWGR + 1
          WGR(MWGR) = T2
          IF (T2.EQ.TI(JTI)) GOTO 2040
          WTI(K) = T2
          GOTO 2050
        end if
C     INCREMENT GREATER THAN ALINC(2)
        IF (IWGR.EQ.0) THEN
          T1 = WTI(KM1)
          T2 = (T2 + T1) / 2.D0
          GOTO 1650
        ELSE
          WTI(K) = WGR(MWGR)
          GOTO 2050
        end if
C     A SHOOTING POINT IS DETERMINED.
C     COMPUTATION OF JK, QK AND UK
 2040   JTI = JTI + 1
        WTI(K) = T2
 2050   DO I = 1 , N
          DO J = 1 , N
            WF(J,I) = (WF(J,I) - S(J,K)) / EPS
            IF (K.EQ.2) WF0(J,I) = WF(J,I)
          end do
        end do

        CALL DQUDEC(WF,N,N,N,N,W(1,1),W(1,2),W(1,3))
        IF (K.EQ.2) THEN
          ISORT = 0
          DO 2400 ITEL = 1 , N
            CALL DSORTD(W(1,1),N,N,KP,IP,ISORT,IER)
            IF (IER.EQ.0) GOTO 2500
            DO I = 1 , N
              IK = IP(I)
              DO J = 1 , N
                WF(J,I) = WF0(J,IK)
                Q(J,I,1) = Q(J,IK,NTI)
              end do
            end do

            DO I = 1 , N
              DO J = 1 , N
                WF0(I,J) = WF(I,J)
                Q(I,J,NTI) = Q(I,J,1)
              end do
            end do
            CALL DQUDEC(WF,N,N,N,N,W(1,1),W(1,2),W(1,3))
 2400     CONTINUE
        end if
 2500   CALL DFQUS(WF,N,N,N,W(1,1),W(1,2),K,KKK,U,NU,NTI,Q,WF0,W(1,3))
        NRTI = K
 2600 CONTINUE
 2700 DO 2800 I = 1 , NRTI
        TI(I) = WTI(I)
 2800 CONTINUE
      IF (TI(NRTI).NE.B) GOTO 5500
      RETURN
 5200 IERROR = 120
      RETURN
 5300 IERROR = 121
      RETURN
 5400 IERROR = 210 + IER
      RETURN
 5500 IERROR = 122
      RETURN
 5550 IERROR = 123
      IP(1) = NOTI
      RETURN
 5600 IERROR = 219
      IP(1) = int ( (B - A) / (T1 - A) * (LWGR+1) + 1 )
      RETURN
      END
      SUBROUTINE DCSAOJ(FDIF,N,ER,TI,NTI,NRTI,X,JA,NG,KKK,S,Q,U,NU,W,
     1                  WF,WF0,IP,WGR,LWGR,IERROR)

c*********************************************************************72
c
cc DCSAOJ
C
      implicit double precision (A-H,O-Z)

      DIMENSION ER(5),TI(NTI),X(N,NTI),S(N,NTI),Q(N,N,NTI),U(NU,NTI),
     1          W(N,7),WF(N,N),WF0(N,N),WGR(LWGR)
      INTEGER KKK(N),IP(N)
      DIMENSION HS(5)
      EXTERNAL FDIF
C
C     INPUT ERROR 122
      IF (NTI-NRTI.LT.1) GOTO 5000
      IERROR = 0
      NHI = 0
      EPS = DSQRT(ER(3))
      IF (JA.EQ.1) THEN
        DO I = 1 , N
          DO J = 1 , N
            Q(I,J,NTI) = Q(I,J,1)
          end do
        end do
      end if
      MWGR = 1
      IF (NG.EQ.1) MWGR = 0
      DO 1900 K = 2 , NRTI
        KM1 = K - 1
        T1 = TI(KM1)
        T2 = TI(K)
        DO 1200 I = 1 , N
          S(I,K) = X(I,KM1)
 1200   CONTINUE
        IF (JA.EQ.1) THEN
          DO I = 1 , N
            DO J = 1 , N
              WF(J,I) = X(J,KM1) + EPS * Q(J,I,KM1)
            end do
          end do
          NF = N
        ELSE
          NF = 0
        end if
        IF (NG.EQ.1) THEN
         NHI = 5
         IER = -1
         IF (K.EQ.2) IER = 1
 1320    CALL DRKFMS(FDIF,N,S(1,K),T1,T2,ER(1),ER(2),ER(3),WF,NF,HI,NHI,
     1               W,IER)
         IF (IER.GT.3) GOTO 5100
         MWGR = MWGR + 1
         IF (MWGR.GT.LWGR) GOTO 5200
         WGR(MWGR) = T2
         IF (T2.NE.TI(K)) THEN
           T1 = T2
           T2 = TI(K)
           IER = -1
           GOTO 1320
         end if
        ELSE
          T0 = T1
 1340     ST = (WGR(MWGR) - T0) / 5.D0
          DO 1360 I = 1 , 5
            HS(I) = ST
 1360     CONTINUE
          CALL DRKFGS(FDIF,N,S(1,K),T1,HS,5,W)
          IF (JA.EQ.1) THEN
            DO 1380 I = 1 , N
              T1 = T0
              CALL DRKFGS(FDIF,N,WF(1,I),T1,HS,5,W)
 1380       CONTINUE
          end if
            T0 = WGR(MWGR)
            MWGR = MWGR + 1
            IF (T0.NE.T2) GOTO 1340
        end if
C     IF JA=1 COMPUTATION OF JK, QK AND UK
        IF (JA.EQ.1) THEN
          DO I = 1 , N
            DO J = 1 , N
              WF(J,I) = (WF(J,I) - S(J,K)) / EPS
              IF (K.EQ.2) WF0(J,I) = WF(J,I)
            end do
          end do
          CALL DQUDEC(WF,N,N,N,N,W(1,1),W(1,2),W(1,3))
          IF (K.EQ.2) THEN
            ISORT = 0
            DO 1700 ITEL = 1 , N
              CALL DSORTD(W(1,1),N,N,KP,IP,ISORT,IER)
              IF (IER.EQ.0) GOTO 1800

              DO I = 1 , N
                IK = IP(I)
                DO J = 1 , N
                  WF(J,I) = WF0(J,IK)
                  Q(J,I,1) = Q(J,IK,NTI)
                end do
              end do

              DO I = 1 , N
                DO J = 1 , N
                  WF0(I,J) = WF(I,J)
                  Q(I,J,NTI) = Q(I,J,1)
                end do
              end do

              CALL DQUDEC(WF,N,N,N,N,W(1,1),W(1,2),W(1,3))
 1700       CONTINUE
          end if
 1800     CALL DFQUS(WF,N,N,N,W(1,1),W(1,2),K,KKK,U,NU,NTI,Q,WF0,W(1,3))
        end if
 1900 CONTINUE
      RETURN
 5000 IERROR = 122
      RETURN
 5100 IERROR = 210 + IER
      RETURN
 5200 IERROR = 219
      IP(1) = int ( 
     &  (TI(NRTI) - TI(1)) / (T1 - TI(1)) * (LWGR+1) + 1 )
      RETURN
      END
      SUBROUTINE DJINGX(G,N,NRTI,X,S,NTI,Q,U,NU,ER,KKK,JAC,D,KPART,
     1                  PHIREC,W,WF,WF0,WBMA,WBMB,IP1,IP2,IERROR)

c*********************************************************************72
c
cc DJINGX
C
      implicit double precision (A-H,O-Z)

      DIMENSION X(N,NTI),S(N,NTI),Q(N,N,NTI),U(NU,NTI),ER(5),D(N,NTI),
     1          PHIREC(NU,NTI),W(N,3),WF(N,N),WF0(N,N),WBMA(N,N),
     2          WBMB(N,N)
      INTEGER KKK(N),IP1(N),IP2(N)
      EXTERNAL G
C
      IERROR = 0
C     COMPUTATION OF THE INHOMOGENEOUS TERMS OF THE RECURSION.
      DO 1400 I = 2 , NRTI
        DO 1100 J = 1 , N
           W(J,1) = S(J,I) - X(J,I)
 1100   CONTINUE
        DO 1300 J = 1 , N
          SOM = 0.D0
          DO 1200 K = 1 , N
            SOM = SOM + Q(K,J,I) * W(K,1)
 1200     CONTINUE
          D(J,I) = SOM
 1300   CONTINUE
 1400 CONTINUE
      IF (JAC.NE.0) THEN
      IP1(1) = 0
      CALL DGTUR(N,1,NRTI,U,NU,NTI,Q,D,ER,KKK,KPART,IP1,IP2,W(1,1),
     1           W(1,2),W(1,3),WF,WF0,IERROR)
      CALL DKPCH(N,U,NU,NTI,NRTI,KPART,KKK,ER,IERROR)
      CALL DFUNRC(N,KPART,NRTI,U,NU,NTI,KKK,PHIREC,W(1,1),W(1,2),IERROR)
      IF (IERROR.NE.0) RETURN
      end if
      CALL DPSR(N,KPART,NRTI,U,NU,NTI,D,KKK,S,W(1,1),W(1,2),IERROR)
      IF (IERROR.NE.0) RETURN
C     COMPUTATION OF THE BOUNDARY CONDITIONS.
      IF (JAC.EQ.0) THEN
        CALL G(N,X(1,1),X(1,NRTI),W(1,1),WF,WF0)
      ELSE
        CALL G(N,X(1,1),X(1,NRTI),W(1,1),WBMA,WBMB)
      end if

      DO I = 1 , N
        W(I,1) = -W(I,1)
      end do

      CALL DSBVP(N,1,WBMA,WBMB,NRTI,Q,NTI,W(1,1),PHIREC,NU,S,S,
     1           ER(3),ER(4),KKK,IP1,WF,WF0,W(1,2),W(1,3),IERROR)
      RETURN
      END
      SUBROUTINE DCHINC(X,S,N,NTI,U,NU,NRTI,ALINCR,KKK,IPG)

c*********************************************************************72
c
cc DCHINC
C
      implicit double precision (A-H,O-Z)

      DIMENSION X(N,NTI),S(N,NTI),U(NU,NTI)
      INTEGER KKK(N),IPG(NTI)
C
      IPG(1) = 0
      DO 1500 K = 2 , NRTI
        KM1 = K - 1
        E1 = 0.D0
        E2 = 0.D0
        DO 1000 J = 1 , N
          E3 = DABS(X(J,KM1))
          IF (E3.GT.E1) E1 = E3
          E3 = DABS(S(J,K))
          IF (E3.GT.E2) E2 = E3
 1000   CONTINUE
        E3 = E2 / E1
        IF (E3.GT.ALINCR) THEN
          IPG(1) = 1
          IPG(K) = 1
        ELSE
          IPG(K) = 0
        end if
 1500 CONTINUE
      DO 2500 K = 2 , NRTI
        E1 = 0.D0
        DO 2000 J = 1 , N
          E2 = U(KKK(J),K)
          IF (E2.GT.E1) E1 = E2
 2000   CONTINUE
        IF (E1.GT.ALINCR) THEN
          IPG(1) = 1
          IPG(K) = 1
        end if
 2500 CONTINUE
      RETURN
      END
      SUBROUTINE DNEWPO(FDIF,N,TOL,EPSMA,TI,NTI,NRTI,ALINC,X,S,KKK,Q,U,
     1                  NU,IPG,WG,LWG,W,WF,WF0,IERROR)

c*********************************************************************72
c
cc DNEWPO
C
      implicit double precision (A-H,O-Z)

      DIMENSION TI(NTI),ALINC(4),X(N,NTI),S(N,NTI),Q(N,N,NTI),
     1          U(NU,NTI),WG(LWG),W(N,7),WF(N,N),WF0(N,N)
      INTEGER KKK(N),IPG(NTI)
      DIMENSION HS(10)
      EXTERNAL FDIF
C
      IERROR = 0
      IF (IPG(1).EQ.0) RETURN
C     INPUT ERROR 122
      IF (NTI-NRTI.LT.1) GOTO 5000
      LHI = 0
      NF = 0
      JAC = 0
      IF (IPG(1).EQ.-1) JAC = 1
      B = TI(NRTI)
      EPS = DSQRT(EPSMA)
      MWGR = 1
      K = 2
C     COMPUTING TOTAL NUMBER OF CURRENT GRID POINTS AND EXPECTED NEEDED
C     NUMBER OF GRIDPOINTS
      NEG = 0
      IK = 0
      IJ = 1
      DO 1300 I = 2 , NRTI
        IM1 = I - 1
        DO 1100 J = IJ , LWG
          IF (WG(J).EQ.TI(IM1)) IK = J
          IF (WG(J).EQ.TI(I)) GOTO 1200
 1100   CONTINUE
 1200   IJ = J
        IF (IPG(I).NE.0) NEG = NEG + IJ - IK
 1300 CONTINUE
      NTOTGR = IJ
      NEG = NEG + NTOTGR
 1500 IF ((K.GT.NTI).OR.(NRTI.GE.NTI)) GOTO 5000
      NSHPO = 0
      KM1 = K - 1
      KPL1 = K + 1
      T1 = TI(KM1)
      IF (T1.GE.B) GOTO 4900
      T2 = TI(K)
C     IF IPG(K)=0 NO SHOOTING POINT HAS TO BE INSERTED IN THE INTERVAL
C     (TI(K-1),TI(K)).
C     IF JAC=0 NO UPDATE FOR THE JACOBIAN IS NEEDED.
      IF (IPG(K).EQ.0) THEN
        DO 1600 I = 1 , N
          S(I,K) = X(I,KM1)
          IF (JAC.NE.0) THEN
            DO 1550 J = 1 , N
              WF(J,I) = X(J,KM1) + EPS * Q(J,I,KM1)
 1550       CONTINUE
          end if
 1600   CONTINUE
        T0 = T1
 1700   ST = (WG(MWGR) - T0) / 5.D0
        DO 1800 I = 1 , 5
          HS(I) = ST
 1800   CONTINUE
        CALL DRKFGS(FDIF,N,S(1,K),T1,HS,5,W)
        IF (JAC.NE.0) THEN
          DO 1900 I = 1 , N
            T1 = T0
            CALL DRKFGS(FDIF,N,WF(1,I),T1,HS,5,W)
 1900     CONTINUE
        end if
        T0 = WG(MWGR)
        T1 = T0
        MWGR = MWGR + 1
        IF (T1.NE.T2) GOTO 1700
        IF (JAC.NE.0) THEN
          DO I = 1 , N
            DO J = 1 , N
              WF(J,I) = (WF(J,I) - S(J,K)) / EPS
            end do
          end do
          CALL DQUDEC(WF,N,N,N,N,W(1,1),W(1,2),W(1,3))
          CALL DFQUS(WF,N,N,N,W(1,1),W(1,2),K,KKK,U,NU,NTI,Q,WF0,W(1,3))
        end if
        K = K + 1
        GOTO 1500
      end if
C     A NEW SHOOTING POINT IS NEEDED AND FROM NOW ON JACOBIAN UPDATES
C     ARE NECESSARY.
      JAC = 1
      NSHPO = 0
      KWGR = MWGR
      KWG = MWGR
      X1 = 0.D0
      DO I = 1 , N
        X2 = DABS(X(I,KM1))
        IF (X2.GT.X1) X1 = X2
        S(I,K) = X(I,KM1)
        S(I,KPL1) = X(I,KM1)
        Q(I,1,NTI) = 1.D0
        DO J = 1 , N
          WF(J,I) = X(J,KM1) + EPS * Q(J,I,KM1)
          WF0(J,I) = WF(J,I)
        end do
      end do

      IWGR = 0
      IWG = 0
 2100 CALL DRKFGG(FDIF,N,T1,WG(MWGR),S(1,K),WF,HS,5,W)
C     CHECK INCREMENT
      X2 = 0.D0
      X3 = 0.D0
      DO 2300 I = 1 , N
        ST = DABS(S(I,K))
        IF (ST.GT.X2) X2 = ST
        X4 = 0.D0
        DO 2200 J = 1 , N
          ST = DABS((WF(J,I)-S(I,K))/Q(I,1,NTI))
          IF (ST.GT.X3) X3 = ST
          IF (ST.GT.X4) X4 = ST
 2200   CONTINUE
        Q(I,2,NTI) = X4
 2300 CONTINUE
      IF ((X2/X1.LT.ALINC(1)).AND.(X3.LT.ALINC(3))) THEN
        IWGR = IWGR + 1
        T1 = WG(MWGR)
        MWGR = MWGR + 1
        IF (T1.EQ.T2) GOTO 3500
        GOTO 2100
      end if
      NSHPO = 1
      IF ((IWGR.EQ.0).AND.(WG(MWGR).EQ.T2)) GOTO 2800
      IF (KWGR.EQ.1) THEN
        T1 = TI(KM1)
      ELSE
        T1 = WG(KWGR-1)
      end if
      DO 2400 I = KWGR , KWGR+IWGR
        CALL DRKFGG(FDIF,N,T1,WG(I),S(1,KPL1),WF0,HS,10,W)
        T1 = WG(I)
 2400 CONTINUE
      X4 = 0.D0
      DO 2600 I = 1 , N
        ST = DABS(S(I,KPL1))
        IF (ST.GT.X4) X4 = ST
        X5 = 0.D0
        DO 2500 J = 1 , N
          ST = DABS((WF0(J,I)-S(J,KPL1))/Q(I,1,NTI))
          IF (ST.GT.X5) X5 = ST
 2500   CONTINUE
        Q(I,2,NTI) = Q(I,2,NTI) / X5
 2600 CONTINUE
      X5 = X2 / X4
      DO 2700 I = 1 , N
        IF (Q(I,2,NTI).GT.X5) X5 = Q(I,2,NTI)
 2700 CONTINUE
      IF ((X5.GE.0.9D0).AND.(X5.LE.1.12D0)) GOTO 3200
C     NEW GRID NEEDED.
 2800 JWGR = IWGR + 1
      IF (NTOTGR+JWGR.GT.LWG) GOTO 5200
      DO 2900 I = NTOTGR , MWGR , -1
        WG(I+JWGR) = WG(I)
 2900 CONTINUE
      J = KWGR - 1
      DO 3000 I = JWGR , 1 , -1
        WG(J+2*I) = WG(J+I)
 3000 CONTINUE
      J = KWGR - 2
      IF (KWGR.EQ.1) THEN
        WG(KWGR) = (WG(KWGR+1) + TI(1)) / 2.D0
        IK = 2
      ELSE
        IK = 1
      end if
      DO 3100 I = IK , JWGR
        I1 = J + 2*I
        WG(I1) = (WG(I1+1) + WG(I1-1)) / 2.D0
 3100 CONTINUE
      MWGR = MWGR + JWGR
 3200 IF (WG(MWGR).NE.T2) THEN
        MWGR = MWGR + 1
        KWGR = MWGR
        IWGR = 0
        X1 = 0.D0
        DO 3400 I = 1 , N
          ST = DABS(S(I,K))
          IF (ST.GT.X1) X1 = ST
          S(I,KPL1) = S(I,K)
          X4 = 0.D0
          DO 3300 J = 1 , N
            ST = DABS(WF(J,I))
            IF (ST.GT.X4) X4 = ST
            WF0(J,I) = WF(J,I)
 3300     CONTINUE
          Q(I,1,NTI) = X4
 3400   CONTINUE
        GOTO 2100
      ELSE
        MWGR = MWGR + 1
      end if
 3500 IF (NSHPO.EQ.0) GOTO 4100
C     DETERMINATION OF THE NEW SHOOTING POINT.
      T1 = TI(KM1)
      T2 = TI(K)
      TP = (T2 + T1) / 2.D0
      E1 = DABS(T2-T1)
      DO 3600 I = KWG , LWG
        ST = DABS(TP-WG(I))
        IF (ST.LT.E1) THEN
          E1 = ST
          JP = I
        ELSE
          GOTO 3700
        end if
 3600 CONTINUE
 3700 TP = WG(JP)
C     INSERT NEW SHOOTING POINT
      NRTI = NRTI + 1
      DO I = NRTI , KPL1 , -1
        IM1 = I-1
        TI(I) = TI(IM1)
        IPG(I) = IPG(IM1)
        DO J = 1 , N
          X(J,I) = X(J,IM1)
        end do
      end do

      TI(K) = TP
C     COMPUTATION OF S(.,K) AND WF
      DO I = 1 , N
        S(I,K) = X(I,KM1)
        DO J = 1 , N
          WF(J,I) = X(J,KM1) + EPS * Q(J,I,KM1)
        end do
      end do

      DO 4000 I = KWG , JP
        CALL DRKFGG(FDIF,N,T1,WG(I),S(1,K),WF,HS,5,W)
        T1 = WG(I)
 4000 CONTINUE
 4100 DO I = 1 , N
        DO J = 1 , N
          WF(J,I) = (WF(J,I) - S(J,K)) / EPS
        end do
      end do
 
      CALL DQUDEC(WF,N,N,N,N,W(1,1),W(1,2),W(1,3))
      CALL DFQUS(WF,N,N,N,W(1,1),W(1,2),K,KKK,U,NU,NTI,Q,WF0,W(1,3))
      IF (NSHPO.EQ.0) THEN
        K = K + 1
        GOTO 1500
      end if
C     COMPUTATION NEW STARTING VALUE FOR X(.,K).
      X2 = U(1,K)
      DO 4200 I = 2 , N
        IK = KKK(I)
        IF (U(IK,K).LT.X2) X2 = U(IK,K)
 4200 CONTINUE
      ST = DLOG(X2) * (TI(K) - TI(KPL1)) / (TI(K) - TI(KM1))
      ST = DEXP(ST)
      AE = DMIN1(TOL,TOL/ST)
      RE=1.1D-12 + 3.D0 * EPSMA
      T1 = TI(KPL1)
      T2 = TI(K)
      IER = 1
      DO 4300 I = 1 , N
        X(I,K) = X(I,KPL1)
 4300 CONTINUE
      CALL DRKFMS(FDIF,N,X(1,K),T1,T2,RE,AE,EPSMA,WF,NF,HI,LHI,W,IER)
      IF (IER.GT.3) GOTO 5100

      DO I = 1 , N
        S(I,KPL1) = X(I,K)
        DO J = 1 , N
          WF(J,I) = X(J,K) + EPS * Q(J,I,K)
        end do
      end do

      T1 = TP
      DO 4500 I = JP+1 , MWGR-1
        CALL DRKFGG(FDIF,N,T1,WG(I),S(1,KPL1),WF,HS,5,W)
        T1 = WG(I)
 4500 CONTINUE

      DO I = 1 , N
        DO J = 1 , N
          WF(J,I) = (WF(J,I) - S(J,KPL1)) / EPS
        end do
      end do

      CALL DQUDEC(WF,N,N,N,N,W(1,1),W(1,2),W(1,3))
      CALL DFQUS(WF,N,N,N,W(1,1),W(1,2),KPL1,KKK,U,NU,NTI,Q,WF0,W(1,3))
      K = K + 2
      GOTO 1500
 4900 IF (TI(NRTI).NE.B) GOTO 5000
      RETURN
 5000 IERROR = 122
      RETURN
 5100 IERROR = IER + 210
      RETURN
 5200 IERROR = 219
      IPG(1) = NEG
      RETURN
      END
      SUBROUTINE DBCMAV(N,MA,MB,BCV,IHOM,NRTI,Q,NTI,PHIREC,NU,Z,KKK,
     1                  R,BB,WS)

c*********************************************************************72
c
cc DBCMAV
C
      implicit double precision (A-H,O-Z)

      double precision MA(N,N),MB(N,N)
      DIMENSION BCV(N),Q(N,N,NTI),PHIREC(NU,NTI),Z(N,NTI),
     1          R(N,N),BB(N),WS(N,N)
      INTEGER KKK(N)
      LOGICAL HOM0,HOM1
C
      HOM0 = IHOM.EQ.0
      HOM1 = IHOM.NE.0
      DO 1000 I = 1 , N
        KKK(I) = I*(I+1)/2
 1000 CONTINUE
      DO 1010 I = 1 , N
        BB(I) = BCV(I)
 1010 CONTINUE
C     COMPUTATION OF MA.Q(1) STORED IN WS
      DO 1120 K = 1 , N
        CALL DMATVC(MA,N,N,N,N,Q(1,K,1),WS(1,K))
 1120 CONTINUE
      IF (HOM1) THEN
C     COMPUTATION OF MA.Q(1).Z(1) AND BCV-MA.Q(1).Z(1) STORED IN BB.
        DO 1210 I = 1 , N
          SOM = 0.D0
          DO 1220 J = 1 , N
            SOM = SOM + WS(I,J) * Z(J,1)
 1220     CONTINUE
          BB(I) = BB(I) - SOM
 1210   CONTINUE
      end if
C     COMPUTATION OF MA.Q(1).PHIREC(1) STORED IN R
      DO 1420 K = 1 , N
        NK = KKK(K) - K
        DO 1410 I = 1 , N
          SOM = 0.D+0
          DO 1400 J = 1 , K
            NJ = NK + J
            SOM = SOM + WS(I,J) * PHIREC(NJ,1)
 1400     CONTINUE
          R(I,K) = SOM
 1410   CONTINUE
 1420 CONTINUE
C     COMPUTATION OF MB.Q(N) STORED IN WS
      DO 1520 K = 1 , N
        CALL DMATVC(MB,N,N,N,N,Q(1,K,NRTI),WS(1,K))
 1520 CONTINUE
      IF (HOM1) THEN
C     COMPUTATION OF MB.Q(N).Z(N) AND BCV-MB.Q(N).Z(N) STORED IN BB
        DO 1610 I = 1 , N
          SOM = 0.D0
          DO 1600 J = 1 , N
            SOM = SOM + WS(I,J) * Z(J,NRTI)
 1600     CONTINUE
          BB(I) = BB(I) - SOM
 1610   CONTINUE
      end if
C     COMPUTATION OF MB.Q(N).PHIREC(N) AND
C     MA.Q(1).PHIREC(1) + MB.Q(N).PHIREC(N) STORED IN R
      DO 1820 K = 1 , N
        NK = KKK(K) - K
        DO 1810 I = 1 , N
          SOM = 0.D+0
          DO 1800 J = 1 , K
            NJ = NK + J
            SOM = SOM + WS(I,J) * PHIREC(NJ,NRTI)
 1800     CONTINUE
        R(I,K) = R(I,K) + SOM
 1810   CONTINUE
 1820 CONTINUE
      RETURN
      END
      SUBROUTINE DFQUS(PHI,N,IRM,IRN,DIAG,BET,K,KKK,U,NU,NTI,Q,
     1                 WPHI0,WSUM)

c*********************************************************************72
c
cc DFQUS
C
      implicit double precision (A-H,O-Z)

      DIMENSION PHI(IRM,IRN),DIAG(IRN),BET(IRN),U(NU,NTI),
     1          Q(IRM,IRM,NTI),WPHI0(IRM,IRM),WSUM(IRN)
      INTEGER KKK(N)
C
      CALL DQEVAL(N,N,PHI,IRM,IRN,BET,WPHI0,WSUM)

      DO I = 1 , N
        DO J = 1 , N
          Q(I,J,K) = WPHI0(I,J)
        end do
      end do

      U(1,K) = DIAG(1)
      IF (N.LT.2) GOTO 2000
      DO 1200 I = 2 , N
      NI = KKK(I)
      IMIN1 = I - 1
      NJ = NI - I
      U(NI,K) = DIAG(I)
      DO 1100 J = 1 , IMIN1
      NK = NJ + J
      U(NK,K) = PHI(J,I)
 1100 CONTINUE
 1200 CONTINUE
 2000 DO 2300 I = 1 , N
        NI = KKK(I)
        IF (U(NI,K).GE.0.D0) GOTO 2300
        DO 2100 J = 1 , N
          Q(J,I,K) = -Q(J,I,K)
 2100   CONTINUE
        DO 2200 J = I , N
          NJ = KKK(J) - J + I
          U(NJ,K) = -U(NJ,K)
 2200   CONTINUE
 2300 CONTINUE
      RETURN
      END
      SUBROUTINE DFUNRC(N,KPART,NRTI,U,NU,NTI,KKK,PHIREC,W1,W2,
     1                  IERROR)

c*********************************************************************72
c
cc DFUNRC
C
      implicit double precision (A-H,O-Z)

      DIMENSION PHIREC(NU,NTI),U(NU,NTI),W1(N),W2(N)
      INTEGER KKK(N)

      IERROR = 0
      IF (KPART.LT.0) KPART=0
      IF (KPART.GT.N) KPART=N
      KPPL1 = KPART + 1
C     COMPUTATION OF THE PHI1'S FOR I <= KPART
      IF (KPART.GT.0) THEN
      DO 1400 K = 1 , KPART
      KK = KKK(K)
      NK = KK - K
      NK1 = NK + 1
      DO 1150 I = 1 , K
      PHIREC(NK + I,NRTI) = 0.0D+0
 1150 CONTINUE
      PHIREC(KK,NRTI) = 1.0D+0
      DO 1350 I = 2 , NRTI
      IN = NRTI + 1 - I
      IN1 = IN + 1
      DO 1250 IJ = 1 , K
      W2(IJ) = PHIREC(NK + IJ,IN1)
 1250 CONTINUE
      CALL DSOLUP(K,U(1,IN1),W2,W1,IERROR)
      IF (IERROR.GT.0) GOTO 3000
      DO 1300 IJ = 1 , K
      PHIREC(NK + IJ,IN) = W1(IJ)
 1300 CONTINUE
 1350 CONTINUE
 1400 CONTINUE
      end if
C     THE REQUIRED PART OF PHIREC1 IS NOW COMPUTED.
C     THE REQUIRED PART OF PHIREC2 WILL BE COMPUTED NOW.
      IF (KPPL1.LE.N) THEN
C     COMPUTATION OF THE PHI2'S FOR I > KPART
      DO 1700 K = KPPL1 , N
      KK = KKK(K)
      NK = KK - K
      NK1 = NK + KPPL1
      DO 1500 I = NK1 , KK
      PHIREC(I,1) = 0.0D+0
 1500 CONTINUE
      PHIREC(KK,1) = 1.0D+0
      DO 1650 I = 2 , NRTI
      IN1 = I - 1
      DO 1600 IJ = KPPL1 , K
      SOM = 0.0D+0
      DO 1550 IN = IJ , K
      MK = KKK(IN) - IN
      SOM = SOM + PHIREC(NK + IN , IN1) * U(MK + IJ , I)
 1550 CONTINUE
      PHIREC(NK + IJ , I) = SOM
 1600 CONTINUE
 1650 CONTINUE
 1700 CONTINUE
C     COMPUTATION OF THE PHI1'S FOR I > KPART, IF KPART > 0
      IF (KPART.GT.0) THEN
      DO 2050 K = KPPL1 , N
      KK = KKK(K)
      NK = KK - K
      DO 1750 I = 1, KPART
      PHIREC(NK + I, NRTI) = 0.0D+0
 1750 CONTINUE
      DO 2000 I = 2 , NRTI
      IN = NRTI + 1 - I
      IN1 = IN + 1
      DO 1850 IJ = 1 , KPART
      SOM = 0.0D+0
      DO 1800 NK1 = KPPL1 , K
      MK = KKK(NK1 - 1)
      SOM = SOM + U(MK + IJ,IN1) * PHIREC(NK + NK1 ,IN)
 1800 CONTINUE
      W2(IJ) = - SOM + PHIREC(NK + IJ , IN1)
 1850 CONTINUE
      CALL DSOLUP(KPART,U(1,IN1),W2,W1,IERROR)
      IF (IERROR.GT.0) GOTO 3000
      DO 1950 IJ = 1 , KPART
      PHIREC(NK + IJ , IN) = W1(IJ)
 1950 CONTINUE
 2000 CONTINUE
 2050 CONTINUE
      end if
      end if
      RETURN
 3000 IERROR = 250
      RETURN
      END
      SUBROUTINE DGTUR(N,IHOM,NRTI,U,KU,NTI,Q,D,ER,KKK,KPART,IP1,IP2,
     1                 WDIA,WBET,WSUM,WPHI,WPHI0,IERROR)

c*********************************************************************72
c
cc DGTUR
C
      implicit double precision (A-H,O-Z)

      DIMENSION U(KU,NTI),Q(N,N,NTI),D(N,NTI),ER(5),WDIA(N),
     1          WBET(N),WSUM(N),WPHI(N,N),WPHI0(N,N)
      INTEGER KKK(N),IP1(N),IP2(N)
      LOGICAL NOSORT
C
      IERROR = 0
      ISORT = IP1(1)
      NKP = N
      IF (ISORT.NE.0) NKP = IP1(2)
      IF (NKP.LE.0) THEN
        KPART = 0
        RETURN
      end if
      NTI1 = NTI - 1
      NOSORT = .FALSE.
C.....CHECK WHETHER TRANSFORMATION IS ALLOWED.
      TOL = DMAX1(ER(1),ER(2))
      TOL1 = TOL / ER(3)
      DO 1100 I = 1 , N
        IP1(I) = I
 1100 CONTINUE
      DO 1300 I = 1 , NKP
        IK = KKK(I)
        WDIA(I) = 0.0D+0
        DO 1200 J = 2 , NRTI
          SOM = DABS(U(IK,J))
          IF (SOM.GE.TOL1) NOSORT = .TRUE.
          WDIA(I) = WDIA(I) + DLOG10(SOM)
 1200 CONTINUE
 1300 CONTINUE
      KPART = 0
      DO 1400 I = 1 , NKP
        IF (WDIA(I).GT.0) KPART = KPART + 1
 1400 CONTINUE
C     PARTITIONING IS DETERMINED
      IF ((KPART.EQ.0).OR.(KPART.EQ.NKP)) RETURN
      IF (NOSORT) RETURN
      IFLAG = 0
      K = 0
      N1 = NKP - 1
C     DETERMINATION OF THE PERMUTATION WHICH GLOBALLY ORDERS THE MODES
      DO 1700 I = 1 , N1
        SORT = WDIA(I)
        IF (SORT.LE.0.D0) THEN
          IPL1 = I + 1
          DO 1500 J = IPL1 , NKP
            SOM = WDIA(J)
            IF ((SOM.GT.0.D0).AND.(SOM-SORT.GT.1.D0)) THEN
              IFLAG = 1
              WDIA(J) = SORT
              WDIA(I) = SOM
              JM = IP1(I)
              IP1(I) = IP1(J)
              IP1(J) = JM
              GOTO 1600
            end if
 1500     CONTINUE
 1600     CONTINUE
        end if
 1700 CONTINUE
      IF (IFLAG.EQ.0) RETURN
C.....THE PERMUTATION P IS NOW DETERMINED
C.....COMPUTATION OF Q1P,U2P,INV(P)D1,INV(P)W1
      IERROR = 1
      DO I = 1 , N
        DO J = 1 , N
          Q(I,J,NTI) = Q(I,J,1)
          WPHI(I,J) = 0.0D+0
        end do
      end do

      DO 2500 I = 1 , N
        IK = IP1(I)
        IP1(I) = IP2(IK)
        IK1 = KKK(IK) - IK
        DO 2300 J = 1 , N
           Q(J,I,1) = Q(J,IK,NTI)
 2300   CONTINUE
        DO 2400 J = 1 , IK
          WPHI(J,I) = U(IK1+J,2)
 2400   CONTINUE
 2500 CONTINUE
      DO 2600 I = 1 , N
        IP2(I) = IP1(I)
 2600 CONTINUE
C.....TRANSFORMATION OF THE U'S,Q'S,D'S AND WI'S
      DO 3400 I = 2 , NRTI
        CALL DQUDEC(WPHI,N,N,N,N,WDIA,WBET,WSUM)
        CALL DFQUS(WPHI,N,N,N,WDIA,WBET,NTI,KKK,U,KU,NTI,Q,
     1             WPHI0,WSUM)
        DO 2700 J = 1 , KU
          U(J,I) = U(J,NTI)
 2700   CONTINUE
        IF (IHOM.EQ.0) GOTO 2950
        CALL DTAMVC(Q(1,1,NTI),N,N,N,N,D(1,I),WBET)
        DO 2900 J = 1 , N
          D(J,I) = WBET(J)
 2900   CONTINUE
 2950   DO 3000 J = 1 , N
          CALL DMATVC(Q(1,1,I),N,N,N,N,Q(1,J,NTI),WPHI0(1,J))
 3000   CONTINUE

        DO J = 1 , N
          DO L = 1 , N
            Q(J,L,I) = WPHI0(J,L)
          end do
        end do

        I1 = I + 1
        IF (I1.GT.NRTI) GOTO 3400
C.....COMPUTATION OF UI+1.OI
        DO J = 1 , N
          DO L = 1 , N
            SOM = 0.0D+0
            DO K = L , N
              IK = KKK(K) - K
              SOM = SOM + U(IK+L,I1) * Q(K,J,NTI)
            end do
            WPHI(L,J) = SOM
          end do
        end do

 3400 CONTINUE

      RETURN
      END
      SUBROUTINE DINPRO(A,NA,B,NB,N,C)

c*********************************************************************72
c
cc DINPRO
C
C      DINPRO IS THE double precision VERSION OF SUBROUTINE
C      VVIPP FROM THE ACCU LIBRARY
C
       double precision A(*),B(*),C
       INTEGER NN,JL,JJ

       C = 0.0D+0
       IF (N.LE.0) GOTO 2000
       NN = 1 - NB
       JL = NA * N
       DO 1000 JJ = 1 , JL , NA
        NN = NN + NB
        C = C + A(JJ) * B(NN)
 1000  CONTINUE
 2000  CONTINUE
      RETURN
      END
      SUBROUTINE DKPCH(N,U,NU,NTI,NRTI,KPART,KKK,ER,IER)

c*********************************************************************72
c
cc DKPCH
C
      implicit double precision (A-H,O-Z)

      DIMENSION U(NU,NTI),ER(5)
      INTEGER KKK(N)

      IER = 0
      PSI = DMAX1(ER(1),ER(2)) / ER(3)
      ER(5) = 0.0D+0
      BLOWUQ = 0.0D+0
      KPPL1 = KPART + 1
      IF (KPART.GE.N) GOTO 1300
      DO 1200 K = KPPL1 , N
        PROD = 1.0D+0
        KK = KKK(K)
        DO 1100 I = 2 , NRTI
          PROD = 1.0D+0 + PROD * DABS(U(KK,I))
          IF (BLOWUQ.LT.PROD) BLOWUQ = PROD
 1100   CONTINUE
 1200 CONTINUE
 1300 IF (KPART.LT.1) GOTO 2000
      DO 1500 K = 1 , KPART
        PROD = 1.0D+0
        KK = KKK(K)
        DO 1400 I = 2 , NRTI
          IN = NRTI + 2 - I
          PROD = 1.0D+0 + PROD / DABS(U(KK,IN))
          IF (ER(5).LT.PROD) ER(5) = PROD
 1400   CONTINUE
 1500 CONTINUE
      IF (KPART.GE.N) GOTO 2000
      BLOW1 = 0.0D+0
      BLOW2 = 0.0D+0
      DO 1700 K = 1 , KPART
        KK = KKK(K)
        AMX = 1.0D+0 / DABS(U(KK,NRTI))
        AMX2 = AMX
        IF (NRTI.LT.3) GOTO 1650
        DO 1600 I = 3 , NRTI
        IN = NRTI + 2 - I
        PROD = 1.0D+0 / DABS(U(KK,IN))
        AMX1 = PROD * AMX2
        AMX2 = DMAX1(AMX1,PROD)
        AMX = DMAX1(AMX,AMX2)
 1600   CONTINUE
 1650   IF (AMX.GT.BLOW1) BLOW1 = AMX
 1700 CONTINUE
      DO 1900 K = KPPL1 , N
        KK = KKK(K)
        AMX = DABS(U(KK,2))
        AMX2 = AMX
        DO 1800 I = 3 , NRTI
        PROD = DABS(U(KK,I))
        AMX1 = PROD * AMX2
        AMX2 = DMAX1(PROD,AMX1)
        AMX = DMAX1(AMX,AMX2)
 1800   CONTINUE
        IF (AMX.GT.BLOW2) BLOW2 = AMX
 1900 CONTINUE
      BLOW = BLOW1 * BLOW2
      ER(5) = ER(5) + BLOW
 2000 ER(5) = DMAX1(ER(5),BLOWUQ)
      IF (ER(5).LT.PSI) RETURN
      IER = 240
      RETURN
      END
      SUBROUTINE DMATVC(A,M,N,IRM,IRN,Y,X)

c*********************************************************************72
c
cc DMATVC
C
      implicit double precision (A-H,O-Z)

      DIMENSION A(IRM,IRN),Y(IRN),X(IRM)

      IF (M.LE.0.OR.N.LE.0) GOTO 2000
      DO 1100 I=1,M
      S = 0.0D+0
      DO 1000 J=1,N
      S = S + A(I,J) * Y(J)
 1000 CONTINUE
      X(I) = S
 1100 CONTINUE
 2000 CONTINUE
      RETURN
      END
      SUBROUTINE DPSR(N,KPART,NRTI,U,NU,NTI,D,KKK,V,W1,W2,IERROR)

c*********************************************************************72
c
cc DPSR
C
      implicit double precision (A-H,O-Z)

      DIMENSION U(NU,NTI),D(N,NTI),V(N,NTI),W1(N),W2(N)
      INTEGER KKK(N)

      IERROR = 0
      KPPL1 = KPART + 1
C     COMPUTATION OF VI2
      IF (KPART.GE.N) GOTO 1500
      DO 1300 I = KPPL1 , N
      V(I,1) = 0.0D+0
 1300 CONTINUE
      DO 1450 I = 2 , NRTI
      IN = I - 1
      DO 1400 K = KPPL1 , N
      SOM = 0.0D+0
      DO 1350 L = K , N
      MK = KKK(L) - L
      SOM = SOM + U(MK+K,I) * V(L,IN)
      MK1 = MK + K
 1350 CONTINUE
      V(K,I) = SOM + D(K,I)
 1400 CONTINUE
 1450 CONTINUE
C     COMPUTATION OF VI1
 1500 IF (KPART.LT.1) RETURN
      DO 1600 I = 1 , KPART
      V(I,NRTI) = 0.0D+0
 1600 CONTINUE
      DO 1850 I = 2 , NRTI
      IN = NRTI + 1 - I
      IN1 = IN + 1
      DO 1700 K = 1 , KPART
      SOM = 0.0D+0
      IF (KPPL1.GT.N) GOTO 1660
      DO 1650 L = KPPL1 , N
      MK = KKK(L) - L
      SOM = SOM + U(MK+K,IN1) * V(L,IN)
 1650 CONTINUE
 1660 W2(K) = - SOM - D(K,IN1) + V(K,IN1)
 1700 CONTINUE
      CALL DSOLUP(KPART,U(1,IN1),W2,W1,IERROR)
      IF (IERROR.NE.0) GOTO 3000
      DO 1800 K = 1 , KPART
      V(K,IN) = W1(K)
 1800 CONTINUE
 1850 CONTINUE
      RETURN
 3000 IERROR = 250
      RETURN
      END
      SUBROUTINE DQEVAL(M,N,QR,IRM,IRN,BET,QREV,WD)

c*********************************************************************72
c
cc DQEVAL
C
      implicit double precision (A-H,O-Z)

      DIMENSION QR(IRM,IRN),QREV(IRM,IRM),BET(IRN),WD(IRM)

      DO 1100 I=N,M
      DO 1050 J=I,M
        QREV(I,J) = -QR(I,N) * QR(J,N) * BET(N)
        QREV(J,I) = QREV(I,J)
 1050 CONTINUE
      QREV(I,I) = QREV(I,I) + 1.0D+0
 1100 CONTINUE
      IF (N.LT.2) RETURN
      DO 1300 IK = 2 , N
      I = N + 2 - IK
C     THE NECESARY INNERPRODUCTS OF THE UK'S ARE STORED IN THE
C     ARRAY WD(1:IRM)
C
      IMIN1 = I - 1
      WD(IMIN1) = QR(IMIN1,IMIN1)
      IM = M - I + 1
      DO 1150 J = I , M
        CALL DINPRO(QR(I,IMIN1),1,QREV(I,J),1,IM,WD(J))
C
C     COMPUTATION OF (I - BET(I-1).U(I-1)') * PRODUCT
C     ( ( -BET(K).U(K).U(K)') K= I TO N)
C
 1150 CONTINUE

      DO K = I , M
        DO L = I , M
          QREV(K,L) = QREV(K,L) - BET(IMIN1) * QR(K,IMIN1) * WD(L)
        end do
      end do

      DO 1250 J = I , M
        QREV(J,IMIN1) = - BET(IMIN1) * QR(J,IMIN1) * WD(IMIN1)
        QREV(IMIN1,J) = - BET(IMIN1) * QR(IMIN1,IMIN1) * WD(J)
 1250 CONTINUE
      QREV(IMIN1,IMIN1) = 1.0D+0 - BET(IMIN1)*WD(IMIN1)*WD(IMIN1)
 1300 CONTINUE
      RETURN
      END
      SUBROUTINE DQUDEC(A,M,N,IRM,IRN,DIAG,BET,WY)

c*********************************************************************72
c
cc DQUDEC
C
      implicit double precision (A-H,O-Z)

      DIMENSION A(IRM,IRN),DIAG(IRN),BET(IRN),WY(IRN)

      IF (M.GT.IRM) GOTO 3000
      IF (N.GT.IRN) GOTO 3010

      DO 2000 K = 1 , N
        K1 = K + 1
        IK = M - K + 1
        CALL DINPRO(A(K,K),1,A(K,K),1,IK,SIGMA)
        IF (SIGMA.EQ.0.D0) GOTO 1400
        QRKK = A(K,K)
        ALPHAK = DSQRT(SIGMA)
        IF (QRKK.GE.0.0D+0) ALPHAK = - ALPHAK
        DIAG(K) = ALPHAK
        BETA = 1.0D+0 / (SIGMA - QRKK * ALPHAK)
        BET(K) = BETA
        A(K,K) = QRKK - ALPHAK
        IF (K1.GT.N) GOTO 2000
        DO 1200 J = K1 , N
          CALL DINPRO(A(K,K),1,A(K,J),1,IK,WY(J))
          WY(J) = WY(J) * BETA
 1200   CONTINUE

        DO J = K1 , N
          DO I = K  , M
           A(I,J) = A(I,J) - A(I,K) * WY(J)
          end do
        end do
        GOTO 2000
 1400   BET(K) = 0.D0
        DIAG(K) = 0.D0
 2000 CONTINUE
      RETURN
 3000 WRITE(*,100)
      RETURN
 3010 WRITE(*,110)
      RETURN
C
  100 FORMAT(' ERROR: M > IRM ; NO DECOMPOSITION DONE')
  110 FORMAT(' ERROR: N > IRN ; NO DECOMPOSITION DONE')
      END
      SUBROUTINE DSBVP(N,IHOM,M0,MN,NRTI,Q,NTI,BCV,PHIREC,NU,V,
     1                 WI,EPS,COND,KKK,IP,WQ,WS,W1,W2,IER)

c*********************************************************************72
c
cc DSBVP
C
      implicit double precision (A-H,O-Z)

      double precision M0(N,N),MN(N,N)
      DIMENSION Q(N,N,NTI),BCV(N),PHIREC(NU,NTI),V(N,NTI),
     1          WI(N,NTI),WQ(N,N),WS(N,N),W1(N),W2(N)
      INTEGER KKK(N),IP(N)

      CALL DBCMAV(N,M0,MN,BCV,IHOM,NRTI,Q,NTI,PHIREC,NU,WI,KKK,WQ,
     1            W1,WS)

      CALL DCROUT(WQ,N,N,W1,EPS,IP,W2,IER)
 
      IF (IER.NE.0) THEN
       IER = 260
       RETURN
      end if
c
C     COMPUTATION OF THE CONDITION NUMBER
c
      COND = 0.0D+0
      DO I = 1 , N
        DO J = 1 , N
          W2(J) = 0.0D+0
        end do
        W2(I) = 1.0D+0
        CALL DSOLDE(WQ,N,N,IP,W2)
        SOMIN = 0.0D+0
        DO J = 1 , N
          SOMIN = SOMIN + DABS(W2(J))
        end do
        IF (SOMIN.GT.COND) COND = SOMIN
      end do
c
C     COMPUTATION OF EI'=PHIREC(I)Y' AND EI=Q(I)EI'
c
      DO K=1,NRTI
        DO I=1,N
          SOMIN=0.0D+0
          DO J=I,N
            NJ=KKK(J) - J
            NI= NJ + I
            SOMIN = SOMIN + PHIREC(NI,K) * W1(J)
          end do
          IF (IHOM.EQ.0) THEN
            WI(I,K) = SOMIN
          ELSE
            WI(I,K) = WI(I,K) + SOMIN
          end if
        end do
      end do
c
C     COMPUTATION OF Q(I)XI'
c
      DO K=1,NRTI
        CALL DMATVC(Q(1,1,K),N,N,N,N,WI(1,K),W2)
        DO J=1,N
          V(J,K) = W2(J)
        end do
      end do

      RETURN
      END
      SUBROUTINE DSOLUP(M,U,B,X,IER)

c*********************************************************************72
c
cc DSOLUP
C
      implicit double precision (A-H,O-Z)

      DIMENSION U(*),X(*),B(*)

      IER = 0
      IF (M.LE.0) RETURN
      M1 = M * (M+1) / 2
      IF (U(M1).EQ.0.0D+0) GOTO 1200
      X(M) = B(M) / U(M1)
      IF (M.EQ.1) RETURN
      MMIN1 = M - 1
      DO 1100 I = 1 , MMIN1
      MI = M - I
      SOMIN = 0.0D+0
      MJ = (MI + 1) * MI / 2
      MIJ = MJ
      DO 1000 J = 1 , I
      MIJ = MIJ + MI - 1 + J
      SOMIN = SOMIN + U(MIJ) * X(MI+J)
 1000 CONTINUE
      IF (U(MJ).EQ.0.0D+0) GOTO 1200
      X(MI) = (B(MI) - SOMIN) / U(MJ)
 1100 CONTINUE
      RETURN
 1200 IER = 1
      RETURN
      END
      SUBROUTINE DUPUP(N,U1,U2,U3)

c*********************************************************************72
c
cc DUPUP
C
      implicit double precision (A-H,O-Z)

      DIMENSION U1(*),U2(*),U3(*)
C
      DO I = 1 , N
      IK = I * (I+1) / 2 - I
        DO J = 1 , I
          SOM = 0.0D+0
          DO L = J , I
            IL = L * (L+1) / 2 - L
            SOM = SOM + U1(IL+J) * U2(IK+L)
          end do
          U3(IK+J) = SOM
        end do
      end do

      RETURN
      END
      SUBROUTINE DAMTES(AMP,WDIA,N,IFLAG)

c*********************************************************************72
c
cc DAMTES
C
      implicit double precision (A-H,O-Z)

      DIMENSION WDIA(*)

      AMAXU = 0.0D+0
      DO 1000 I = 1 , N
      ABSDIA = DABS(WDIA(I))
      IF (ABSDIA.GT.AMAXU) AMAXU = ABSDIA
 1000 CONTINUE
      IF (AMAXU.GT.AMP) GOTO 1020
      IFLAG = 0
      RETURN
 1020 IFLAG = 1
      RETURN
      END
      SUBROUTINE DTAMVC(A,M,N,IM,IN,Y,X)

c*********************************************************************72
c
cc DTAMVC
C
      implicit double precision (A-H,O-Z)

      DIMENSION A(IM,IN),Y(IM),X(IN)

      IF (M.LE.0 .OR. N.LE.0) RETURN
      DO 1100 I = 1 , N
        SOM = 0.D0
        DO 1000 J = 1 , M
          SOM = SOM + A(J,I) * Y(J)
 1000   CONTINUE
        X(I) = SOM
 1100 CONTINUE

      RETURN
      END
      SUBROUTINE DCDI(D,N,NTI,Q,K,WSUM)

c*********************************************************************72
c
cc DCDI
C
      implicit double precision (A-H,O-Z)

      DIMENSION D(N,NTI),Q(N,N,NTI),WSUM(N)

      CALL DTAMVC(Q(1,1,K),N,N,N,N,D(1,K),WSUM)

      DO 1000 I = 1 , N
        D(I,K) = WSUM(I)
 1000 CONTINUE

      RETURN
      END
      SUBROUTINE DCNRHS ( DIAG, N, NRHS )

c*********************************************************************72
c
cc DCNRHS
C
c  Modified:
c
c    09 January 2013
c
      implicit none

      integer n

      double precision a
      double precision b
      double precision DIAG(N)
      integer i
      integer nrhs

      A = DABS(DIAG(1))
      NRHS = 1

      DO I = 2, N

        B = DABS(DIAG(I))

        IF ( B .LT. 1.D0 ) then
          B = 1.D0 / B
        end if

        IF ( A .lt. B ) then
          A = B 
          NRHS = I
        end if

      end do

      RETURN
      END
      SUBROUTINE DCPHIS(FLIN,FDIF,N,A,B,ER,K,IHOM,NRHS,WY,WI,Q,NTI,W,
     &                  WPHI,WPHI0,IER)

c*********************************************************************72
c
cc DCPHIS
C
      implicit double precision (A-H,O-Z)

      DIMENSION ER(5),WY(N),WI(N,NTI),Q(N,N,NTI),
     1 W(N,7),WPHI(N,N),WPHI0(N,N)
      EXTERNAL FLIN,FDIF
C     LOCAL ARRAYS
      DIMENSION HI(5),HI2(5)
C
C     COMPUTATION OF THE FUNDAMENTAL SOLUTION
      KMIN1 = K - 1
      NTI1 = NTI - 1
      IF ((NRHS.NE.0).OR.(IHOM.NE.0)) GOTO 1600
C     NRHS=0 AND IHOM=0 , DETERMONATION OF NRHS ON THE FIRST MINOR
C     SHOOTING INTERVAL
      DO 1500 L = 1 , N
        T1 = A
        T2 = B
        DO 1100 I = 1 , N
         WY(I) = Q(I,L,1)
 1100   CONTINUE
        IER = 1
        CALL DRKFSM(FLIN,N,WY,T1,T2,ER(1),ER(2),ER(3),HI2,NHI2,W,IER)
        IF (IER.GT.3) RETURN
        DEL2 = (T2 - A) / NHI2
        IF (L.EQ.1) GOTO 1200
        IF (DABS(DEL2).GT.DABS(DEL1)) GOTO 1500
 1200   DEL1 = DEL2
        NRHS = L
        NHI = NHI2
        DO 1300 I = 1 , NHI
          HI(I) = HI2(I)
 1300   CONTINUE
        T3 = T2
        DO 1400 I = 1 , N
          WPHI(I,NRHS) = WY(I)
          WPHI0(I,NRHS)= WY(I)
 1400   CONTINUE
 1500 CONTINUE
      T1 = A
      T2 = B
      B = T3
      GOTO 2400
 1600 IF (IHOM.EQ.0) GOTO 1900
C     INHOMOGENEOUS SYSTEM, COMPUTATION OF THE PARTICULAR SOLUTION
      DO 1700 I = 1 , N
        WY(I) = 0.D0
 1700 CONTINUE
      IER = 1
      T1 = A
      CALL DRKFSM(FDIF,N,WY,T1,B,ER(1),ER(2),ER(3),HI,NHI,W,IER)
      IF (IER.GT.3) RETURN
      DO 1800 I = 1 , N
        WI(I,K) = WY(I)
 1800 CONTINUE
      NRHS = 0
      GOTO 2400
 1900 DO 2000 I = 1 , N
        WY(I) = Q(I,NRHS,KMIN1)
 2000 CONTINUE
      IER = 1
      T1 = A
      CALL DRKFSM(FLIN,N,WY,T1,B,ER(1),ER(2),ER(3),HI,NHI,W,IER)
      IF (IER.GT.3) RETURN
      IF (KMIN1.GT.1) GOTO 2200
      DO 2100 I = 1 , N
        WPHI(I,NRHS) = WY(I)
        WPHI0(I,NRHS)= WY(I)
 2100 CONTINUE
      GOTO 2400
 2200 DO 2300 I = 1 , N
        WPHI(I,NRHS) = WY(I)
 2300 CONTINUE
      GOTO 2400
 2400 CONTINUE
C     COMPUTATION OF THE BASIC HOMOGENEOUS SOLUTIONS, WHICH HAVE NOT
C     BEEN COMPUTED YET
      DO 2900 I = 1 , N
        IF (I.EQ.NRHS) GOTO 2900
        DO 2500 J = 1 , N
          WY(J) = Q(J,I,KMIN1)
 2500   CONTINUE
        T1 = A
        CALL DRKFGS(FLIN,N,WY,T1,HI,NHI,W)
        IF (KMIN1.GT.1) GOTO 2700
        DO 2600 J = 1 , N
          WPHI(J,I) = WY(J)
          WPHI0(J,I) = WY(J)
 2600   CONTINUE
        GOTO 2900
 2700   DO 2800 J = 1 , N
          WPHI(J,I) = WY(J)
 2800   CONTINUE
 2900 CONTINUE
      CALL DQUDEC(WPHI,N,N,N,N,W(1,1),W(1,2),WY)
      RETURN
      END
      SUBROUTINE DCROUT(A,N,M,B,EPS,P,V,IER)

c*********************************************************************72
c
cc DCROUT
C
C
C     DCROUT IS AN double precision FORTRAN VERSION OF
C     T.DEKKER'S SOLDEC
c
      double precision A(M,*),B(*),V(*),EPS
      INTEGER P(*)

      IER = 0
      CALL DLUDEC(A,N,M,EPS,P,V,IER)
      IF (IER.NE.0) RETURN
      CALL DSOLDE(A,N,M,P,B)
      RETURN
      END
      SUBROUTINE DINTCH(A,NA,B,NB,N)

c*********************************************************************72
c
cc DINTCH
C
C     DINTCH IS A double precision VERSION OF
C     J.M. VAN KATS SUBROUTINE INTCHA
c
      double precision A(*),B(*),C

      IF (N.LE.0) RETURN
      NN = 1 - NB
      JL = NA * N
      DO 1000 JJ = 1 , JL , NA
        NN = NN + NB
        C = B(NN)
        B(NN) = A(JJ)
        A(JJ) = C
 1000 CONTINUE
      RETURN
      END
      SUBROUTINE DLUDEC(A,N,M,EPS,P,V,IER)

c*********************************************************************72
c
cc DLUDEC
C
C     DLUDEC IS AN FORTRAN VERSION OF T. DEKKER'S DEC
c
      implicit double precision (A-H,O-Z)

      DIMENSION A(M,*),V(*)
      INTEGER P(*),PK

      IER = 0
      R = -1.0D+0
      Z = 1.0D+0
      DO 1000 I = 1 , N
        CALL DINPRO(A(I,1),M,A(I,1),M,N,S)
        S = DSQRT(S)
        IF (S.GT.R) R = S
        V(I) = 1.0D+0 / S
 1000 CONTINUE
      E = EPS * R
      DO 5000 K = 1 , N
        R = -1.0D+0
        K1 = K - 1
        DO 2000 I = K , N
          CALL DINPRO(A(I,1),M,A(1,K),1,K1,S)
          A(I,K) = A(I,K) - S
          S = DABS(A(I,K)) * V(I)
          IF (S.LE.R) GOTO 2000
          R = S
          PK = I
 2000   CONTINUE
        P(K) = PK
        V(PK) = V(K)
        Z = A(PK,K)
        IF (DABS(Z).LT.E) GOTO 7000
        IF (PK.EQ.K) GOTO 3000
        CALL DINTCH(A(K,1),M,A(PK,1),M,N)
 3000   KK = K + 1
        IF (K.EQ.N) RETURN
        DO 4000 I = KK , N
          CALL DINPRO(A(K,1),M,A(1,I),1,K1,S)
          A(K,I) = (A(K,I) - S) / Z
 4000   CONTINUE
 5000 CONTINUE
 7000 IER = 1
      RETURN
      END
      SUBROUTINE DSOLDE(A,N,M,P,B)

c*********************************************************************72
c
cc DSOLDE
C
C     DSOLDE IS double precision VERSION OF H.ZWEERUS-VINK'S
C     AND J.M. VAN KAT'S SUBROUTINE SOLDEC
c
      implicit double precision (A-H,O-Z)

      DIMENSION A(M,*),B(*)
      INTEGER P(*),PK

      DO 1000 K = 1 , N
        R = B(K)
        PK = P(K)
        CALL DINPRO(A(K,1),M,B,1,K-1,C)
        B(K) = (B(PK) - C) / A(K,K)
        IF (PK.NE.K) B(PK) = R
 1000 CONTINUE
      IF (N.EQ.1) RETURN
      K = N
 2000 K = K - 1
      CALL DINPRO(A(K,K+1),M,B(K+1),1,N-K,C)
      B(K) = B(K) - C
      IF (K.GT.1) GOTO 2000
      RETURN
      END
      SUBROUTINE DDUR(FLIN,FDIF,N,IHOM,A,B,NRTI,AMP,TI,NTI,ER,Q,U,KU,WI,
     1                D,KPART,WY,W,WF,WF0,KKK,IP1,IP2,IERROR)

c*********************************************************************72
c
cc DDUR
C
      implicit double precision (A-H,O-Z)

      DIMENSION TI(NTI),ER(5),Q(N,N,NTI),U(KU,NTI),WI(N,NTI),
     1 D(N,NTI),WY(N),W(N,7),WF(N,N),WF0(N,N)

      INTEGER KKK(N),IP1(N),IP2(N)

      LOGICAL HOM0,HOM1,SORT,INC

      EXTERNAL FLIN,FDIF
C
C.....SETTING STARTING PARAMETERS
c
      IERROR = 0
      JFLAG = 0
      ISORT = IP1(1)
      NKP = N
      IF (ISORT.EQ.1) NKP = IP1(2)
      JSORT = ISORT + 2
      HOM0 = IHOM.EQ.0
      HOM1 = IHOM.NE.0
      SORT = .TRUE.
      ITH = IHOM
      TI(1) = A
      NTI2 = NTI - 3
      NTI1 = NTI - 1
      IF (NRTI.GT.NTI2) GOTO 7200
      TOL = DMAX1(ER(1),ER(2))
      AMP3 = TOL / ER(3)
c
C.....DETERMINE OPTION FOR THE OUTPUT POINTS
      IF (NRTI.EQ.0) THEN
C.....NRTI = 0
        IF (AMP.LE.1.D0) AMP = TOL / (10.D0*ER(3))
        NOTI = 2
        TI(2) = B
      end if
      IF (NRTI.EQ.1) THEN
C.....NRTI = 1
        IF (A.LT.B) THEN
          DO 1100 I = 2 , NTI
            IF (TI(I-1).GE.TI(I)) GOTO 7000
            IF (TI(I).GE.B) GOTO 1300
 1100     CONTINUE
        ELSE
          DO 1200 I = 2 , NTI
            IF (TI(I-1).LE.TI(I)) GOTO 7000
            IF (TI(I).LE.B) GOTO  1300
 1200     CONTINUE
        end if
 1300   IF (TI(I).NE.B) GOTO 7100
        NOTI = I
      end if

      IF (NRTI.GT.1) THEN
C.....NRTI > 1
        SOM = (B-A) / NRTI
        NOTI = NRTI + 1
        TI(NOTI) = B
        DO 1400 I = 2 , NRTI
          TI(I) = TI(I-1) + SOM
 1400   CONTINUE
      end if
      IF (AMP.GT.1.D0) THEN
C.....AMP > 1
        INC = .TRUE.
        AMP1 = AMP
        AMP2 = 2.D0 * AMP
      ELSE
        INC = .FALSE.
      end if
      IF (NOTI.GT.NTI2) GOTO 7200
      DO I = 1 , N
        KKK(I) = I * (I+1) / 2
        IP2(I) = I
        DO J = 1 , N
          Q(J,I,NTI) = Q(J,I,1)
        end do
      end do
      IU = 0
c
C.....Q1 = I , QNTI = Q1 FOR MAKING IT POSSIBLE TO PERMUTATE Q1
C.....LATERON
C.....COMPUTATION OF Q2 AND U2
      T1 = TI(1)
      T2 = TI(2)
      NRHS = 1
      IF (HOM0) NRHS = 0

      CALL DCPHIS(FLIN,FDIF,N,T1,T2,ER,2,IHOM,NRHS,WY,WI,Q,NTI,W,
     1            WF,WF0,IER)
      IF (IER.GT.3) GOTO 7300
      IF (IER.EQ.3) IERROR = 213
      CALL DAMTES(AMP3,W(1,1),N,IFLAG)
      IF (IFLAG.NE.0) SORT = .FALSE.
c
C.....ORDERING OF U2
c
      IF (SORT) THEN
        CALL DSORTD(W(1,1),N,NKP,KPART,IP1,ISORT,IFLAG)
        IF (IFLAG.NE.0) THEN
          JFLAG = 1
          DO I = 1 , N
            IK = IP1(I)
            IP1(I) = IP2(IK)
            DO J = 1 , N
              WF(J,I) = WF0(J,IK)
            end do
          end do
          DO 1700 I = 1 , N
            IP2(I) = IP1(I)
 1700     CONTINUE
          CALL DQUDEC(WF,N,N,N,N,W(1,1),W(1,2),WY)
        end if
      end if
      IF (HOM0) CALL DCNRHS(W(1,1),N,NRHS)
      CALL DFQUS(WF,N,N,N,W(1,1),W(1,2),2,KKK,U,KU,NTI,Q,
     1           WF0,W(1,3))
      IF (HOM1) CALL DCDI(D,N,NTI,Q,2,W(1,3))
      IFLAG = 0
      IF (INC) THEN
        CALL DAMTES(AMP2,W(1,1),N,IFLAG)
        IF (IFLAG.EQ.1) IERROR = 200
      end if
      T1 = T2

C.....COMPUTATION OF V2
 1800 IF (T1.EQ.TI(2)) GOTO 3500
      IFLAG = 0
      IF (INC) CALL DAMTES(AMP1,W(1,1),N,IFLAG)
      T2 = TI(2)
      CALL DCPHIS(FLIN,FDIF,N,T1,T2,ER,3,IHOM,NRHS,WY,WI,Q,NTI,W,WF,
     1            WF0,IER)
      IF (IER.GT.3) GOTO 7300
      IF (HOM0) CALL DCNRHS(W(1,1),N,NRHS)
      CALL DFQUS(WF,N,N,N,W(1,1),W(1,2),3,KKK,U,KU,NTI,Q,WF0,W(1,3))
 1850 IF (IFLAG.EQ.1) THEN
        NOTI = NOTI + 1
        IF (NOTI.GT.NTI2) GOTO 7200
        DO 1900 I = NOTI , 3 , -1
          TI(I) = TI(I-1)
 1900   CONTINUE
        TI(2) = T1
        IU = 1
        T1 = T2
        GOTO 3500
      end if
      CALL DUPUP(N,U(1,3),U(1,2),U(1,4))
      DO 2000 I = 1 , N
        IK = KKK(I)
        W(I,4) = U(IK,4)
 2000 CONTINUE
      IF (INC) THEN
        CALL DAMTES(AMP2,W(1,4),N,IFLAG)
        IF (IFLAG.EQ.1) GOTO 1850
      end if

      IF (HOM1) THEN
C.....COMPUTATION OF U3.D2
      DO 2200 I = 1 , N
        SOM = 0.0D+0
        DO 2100 J = I , N
          IK = KKK(J) - J
          SOM = SOM + U(IK+I,3) * D(J,2)
 2100   CONTINUE
        D(I,2) = SOM
 2200 CONTINUE
 
C.....COMPUTATION OF D'I+1, STORED IN D3
      CALL DCDI(D,N,NTI,Q,3,W(1,3))
      DO 2300 I = 1 , N
        D(I,2) = D(I,2) + D(I,3)
 2300 CONTINUE
      end if

c
C.....CHECK ORDERING OF U4, IF NOT ORDERED U4.P = Q4.U4' IS
C.....COMPUTED, STORED IN Q4 AND U4
      IFLAG = 0
      IF (SORT) THEN
        CALL DAMTES(AMP3,W(1,4),N,IK)
        IF (IK.NE.0) THEN
          SORT = .FALSE.
          GOTO  3100
        end if
        DO I = 1 , N
          IK = KKK(I) - I
          DO J = 1 , N
            IF (J.LE.I) THEN
              WF0(J,I) = U(IK+J,4)
            ELSE
              WF0(J,I) = 0.D0
            end if
          end do
        end do
        CALL DSORTD(W(1,4),N,NKP,KPART,IP1,JSORT,IFLAG)
        IF (IFLAG.NE.0) THEN
          JFLAG = 1
          DO I = 1 , N
            IK = IP1(I)
            IP1(I) = IP2(IK)
            DO J = 1 , N
              WF(J,I) = WF0(J,IK)
            end do
          end do
          DO I = 1 , N
            IP2(I) = IP1(I)
          end do
          CALL DQUDEC(WF,N,N,N,N,W(1,1),W(1,2),WY)
C.....U4 WAS NOT ORDERED, A NEW U4' IS COMPUTED, STORED IN U4
          CALL DFQUS(WF,N,N,N,W(1,1),W(1,2),4,KKK,U,KU,NTI,Q,WF0,W(1,3))
C.....COMPUTATION OF NEW Q3' = Q3.Q4, STORED IN Q2
          DO 2700 J = 1 , N
            CALL DMATVC(Q(1,1,3),N,N,N,N,Q(1,J,4),Q(1,J,2))
 2700     CONTINUE
          IF (HOM0) CALL DCNRHS(W(1,1),N,NRHS)
          IF (HOM1) THEN
C.....COMPUTATION OF NEW D2' = INV(Q4).D2, STORED IN D2
            CALL DTAMVC(Q(1,1,4),N,N,N,N,D(1,2),W(1,3))
            DO 3000 I = 1 , N
              D(I,2) = W(I,3)
 3000       CONTINUE
          end if
          GOTO 3250
        end if
      end if

c
C.....U4 WAS ALREADY ORDERED SO Q2 = Q3
 3100 DO J = 1 , N
        DO I = 1 , N
          Q(J,I,2) = Q(J,I,3)
        end do
      end do
C.....U2 := U4
 3250 DO 3300 I = 1 , KU
       U(I,2) = U(I,4)
 3300 CONTINUE
      T1 = T2
      GOTO 1800
 3500 DO I = 1 , N
        IK = IP2(I)
        DO J = 1 , N
          Q(J,I,1) = Q(J,IK,NTI)
        end do
      end do

      IF (TI(2).EQ.B) GOTO 6000
      DO 5000 K = 3 , NTI2
        KMIN1 = K - 1
        KPL1 = K + 1
        KPL2 = K + 2
        IF (TI(KMIN1).EQ.B) GOTO 6000
        IF (IU.EQ.1) THEN
          CALL DAMTES(AMP2,W(1,1),N,IFLAG)
          IF (IFLAG.EQ.1) IERROR = 200
          IU = 0
          T2 = TI(K)
        ELSE
          T1 = TI(KMIN1)
          T2 = TI(K)
          CALL DCPHIS(FLIN,FDIF,N,T1,T2,ER,K,IHOM,NRHS,WY,WI,Q,NTI,W,
     1                WF,WF0,IER)
          IF (IER.GT.3) GOTO 7300
          IF (HOM0) CALL DCNRHS(W(1,1),N,NRHS)
          CALL DFQUS(WF,N,N,N,W(1,1),W(1,2),K,KKK,U,KU,NTI,Q,WF0,W(1,3))
          T1 = T2
          T2 = TI(K)
        end if
        IF (HOM1) CALL DCDI(D,N,NTI,Q,K,W(1,3))
 4000   IF (T1.EQ.TI(K)) GOTO 5000
        IFLAG = 0
        IF (INC) CALL DAMTES(AMP,W(1,1),N,IFLAG)
        CALL DCPHIS(FLIN,FDIF,N,T1,T2,ER,KPL1,IHOM,NRHS,WY,WI,Q,NTI,W,
     1              WF,WF0,IER)
        IF (IER.GT.3) GOTO 7300
        IF (HOM0) CALL DCNRHS(W(1,1),N,NRHS)
        CALL DFQUS(WF,N,N,N,W(1,1),W(1,2),KPL1,KKK,U,KU,NTI,Q,
     1             WF0,W(1,3))
 4100   IF (IFLAG.EQ.1) THEN
          NOTI = NOTI + 1
          IF (NOTI.GT.NTI2) GOTO 7200
          DO 4200 I = NOTI , KPL1 , -1
            TI(I) = TI(I-1)
 4200     CONTINUE
          TI(K) = T1
          IU = 1
          T1 = T2
          GOTO 5000
        end if
        CALL DUPUP(N,U(1,KPL1),U(1,K),U(1,KPL2))
        DO 4300 I = 1 , N
          IK = KKK(I)
          W(I,4) = U(IK,KPL2)
 4300   CONTINUE
        IFLAG = 0
        IF (INC) THEN
          CALL DAMTES(AMP2,W(1,4),N,IFLAG)
          IF (IFLAG.EQ.1) GOTO 4100
        end if
        IF (HOM1) THEN
C.....COMPUTATION OF UK+1.DK
          DO 4500 I = 1 , N
            SOM = 0.0D+0
            DO 4400 J = I , N
              IK = KKK(J) - J
              SOM = SOM + U(IK+I,KPL1) * D(J,K)
 4400       CONTINUE
           D(I,K) = SOM
 4500     CONTINUE
          CALL DCDI(D,N,NTI,Q,KPL1,W(1,3))
          DO I = 1 , N
            D(I,K) = D(I,K) + D(I,KPL1)
          end do
        end if
        DO I = 1 , N
          DO J = 1 , N
            Q(I,J,K) =Q(I,J,KPL1)
          end do
        end do
        DO I = 1 , KU
          U(I,K) = U(I,KPL2)
        end do
        T1 = T2
        T2 = TI(K)
        GOTO 4000
 5000 CONTINUE

 6000 continue

      NRTI = NOTI
      IF (TI(NRTI).NE.B) GOTO 7200
      RETURN
 7000 IERROR = 120
      RETURN
 7100 IERROR = 121
      RETURN
 7200 IERROR = 122
      RETURN
 7300 IERROR = IER + 210
      RETURN
      END
      SUBROUTINE DSORTD ( DIAG, N, NKP, KPART, IP, ISORT, IFLAG )

c*********************************************************************72
c
cc DSORTD carries out sorting.
C
c  Modified:
c
c    09 January 2013
c
      implicit none

      integer n

      double precision a
      double precision b
      double precision diag(n)
      integer i
      integer iflag
      integer ip(n)
      integer isort
      integer j
      integer jm
      integer k
      integer kpart
      integer n1
      integer nkp
      integer nn
      double precision sort

      iflag = 0
      DO I = 1 , N
        IP(I) = I
      end do

      NN = N
c
C  ISORT = 1 ; CHECK WHETHER PERMUTATION IS NECESSARY.
c
      IF ( ISORT .EQ. 1 ) THEN
        J = 0
        DO I = 1 , NKP
          IF (DABS(DIAG(I)).GT.1.D0) J = J + 1
        end do
        IF (J.EQ.NKP) THEN
          KPART = J
          RETURN
        end if
        NN = NKP
      end if
c
c  Monotonic ordering.
c
      IF (ISORT.LT.2) THEN
        K = 0
        N1 = NN - 1
        DO I = 1 , N1
          SORT = 0.D0
          DO J = I , NN
            A = DABS(DIAG(J))
            IF (A.GT.SORT) THEN
              SORT = A
              JM = J
            end if
          end do
          IF (SORT.GT.1.D0) K = K+1
          IF (JM.NE.I) THEN
            IFLAG = IFLAG + 1
            A = DIAG(JM)
            DIAG(JM) = DIAG(I)
            DIAG(I) = A
            J = IP(JM)
            IP(JM) = IP(I)
            IP(I) = J
          end if
        end do
        KPART = K
        IF (DABS(DIAG(NN)).GT.1.D0) then
          KPART = KPART + 1
        end if
        RETURN
      end if
c
C  ISORT = 2 OR 3; NO MONOTONIC ORDENING
c
      IF (ISORT.EQ.3) then
        NN = NKP
      end if

      N1 = NN - 1

      DO I = 1 , N1
        A = DABS(DIAG(I))
        IF (A.LT.1.D0) THEN
          DO J = I+1 , NN
            B = DABS(DIAG(J))
            IF ((B.GT.1.D0).AND.(B/A.GE.10.D0)) THEN
              IFLAG = IFLAG + 1
              SORT = DIAG(J)
              DIAG(J) = DIAG(I)
              DIAG(I) = SORT
              JM = IP(J)
              IP(J) = IP(I)
              IP(I) = JM
              GOTO 2100
            end if
          end do
 2100     CONTINUE
        end if
      end do

      KPART = 0
      DO I = 1 , NN
        IF (DABS(DIAG(I)).GT.1.D0) then
          KPART = KPART + 1
        end if
      end do

      RETURN
      END
      SUBROUTINE DRKFGS(F,NEQN,Y,T,HI,NHI,W)

c*********************************************************************72
c
cc DRKFGS
C
      implicit double precision (A-H,O-Z)

      DIMENSION Y(NEQN),HI(5),W(NEQN,7)
      EXTERNAL F

      CALL F(T,Y,W(1,1))

      DO I = 1 , NHI
        CALL DRKF1S(F,NEQN,Y,T,HI(I),W)
        DO J = 1 , NEQN
          Y(J) = W(J,7)
        end do
        T = T + HI(I)
        IF (I.LT.NHI) CALL F(T,Y,W(1,1))
      end do

      RETURN
      END
      SUBROUTINE DRKFSM(F,NEQN,Y,T,TOUT,RE,AE,EPS,HI,NHI,W,IFLAG)

c*********************************************************************72
c
cc DRKFSM implements a Fehlberg 4th/5th order Runge-Kutta method.
C
C     FEHLBERG FOURTH-FIFTH ORDER RUNGE-KUTTA METHOD
C
      implicit double precision (A-H,O-Z)

      DIMENSION Y(NEQN),HI(5),W(NEQN,7)
      LOGICAL HFAILD,OUTPUT
      EXTERNAL F
C
C     REMIN IS THE MINIMUM ACCEPTABLE VALUE OF RELERR.  ATTEMPTS
C     TO OBTAIN HIGHER ACCURACY WITH THIS SUBROUTINE ARE USUALLY
C     VERY EXPENSIVE AND OFTEN UNSUCCESSFUL.
C
      DATA REMIN/1.D-12/

      NHI = 0
      KFLAG = 0
C
C
C     CHECK INPUT PARAMETERS
C
C
      IF (NEQN .LT. 1) GO TO 10
      IF ((RE .LT. 0.0D0)  .OR.  (AE .LT. 0.0D0)) GO TO 10
      MFLAG=IABS(IFLAG)
      IF (IFLAG.NE.1) GOTO 10

      U26 = 26.0D0*EPS
      GO TO 50
C
C     INVALID INPUT
   10 IFLAG=8
      RETURN
C
   50 RER=2.0D0*EPS+REMIN
      IF (RE .GE. RER) GO TO 55
C
C     RELATIVE ERROR TOLERANCE TOO SMALL
      RE=RER
      KFLAG=3

   55 DT=TOUT-T

      A=T
      CALL F(A,Y,W(1,1))
      IF (T .NE. TOUT) GO TO 65
      IFLAG=2
      IF (KFLAG.EQ.3) IFLAG = 3
      RETURN

   65 H=DABS(DT)
      TOLN=0.
      DO 70 K=1,NEQN
        S = RE * DABS(Y(K))
        TOL = S + DMAX1(RE,AE)
        IF (TOL .LE. 0.D0) GO TO 70
        TOLN=TOL
        YPK=DABS(W(K,1))
        IF (YPK*H**5 .GT. TOL) H=(TOL/YPK)**0.2D0
   70 CONTINUE
      IF (TOLN .LE. 0.0D0) H=0.0D0
      H=DMAX1(H,U26*DMAX1(DABS(T),DABS(DT)))
C
C     SET STEPSIZE FOR INTEGRATION IN THE DIRECTION FROM T TO TOUT
C
      H=DSIGN(H,DT)
C
C     INITIALIZE OUTPUT POINT INDICATOR
C
      OUTPUT= .FALSE.
C
C     TO AVOID PREMATURE UNDERFLOW IN THE ERROR TOLERANCE FUNCTION,
C     SCALE THE ERROR TOLERANCES
C
      SCALE=2.0D0/RE
      AAE=SCALE*AE
C
C     STEP BY STEP INTEGRATION
C
  100 HFAILD= .FALSE.
C
C     SET SMALLEST ALLOWABLE STEPSIZE
C
      HMIN=U26*DABS(T)
C
C     ADJUST STEPSIZE IF NECESSARY TO HIT THE OUTPUT POINT.
C     LOOK AHEAD TWO STEPS TO AVOID DRASTIC CHANGES IN THE STEPSIZE AND
C     THUS LESSEN THE IMPACT OF OUTPUT POINTS ON THE CODE.
C
      DT=TOUT-T
      IF (DABS(DT) .GE. 2.0D0*DABS(H)) GO TO 200
      IF (DABS(DT) .GT. DABS(H)) GO TO 150
C
C     THE NEXT SUCCESSFUL STEP WILL COMPLETE THE INTEGRATION TO THE
C     OUTPUT POINT
C
      OUTPUT= .TRUE.
      H=DT
      GO TO 200
C
  150 H=0.5D0*DT
C
C     ADVANCE AN APPROXIMATE SOLUTION OVER ONE STEP OF LENGTH H
C
  200 CALL DRKF1S(F,NEQN,Y,T,H,W)
C
C     COMPUTE AND TEST ALLOWABLE TOLERANCES VERSUS LOCAL ERROR ESTIMATES
C     AND REMOVE SCALING OF TOLERANCES. NOTE THAT RELATIVE ERROR IS
C     MEASURED WITH RESPECT TO THE AVERAGE OF THE MAGNITUDES OF THE
C     SOLUTION AT THE BEGINNING AND END OF THE STEP.
C
      EEOET=0.0D0
      DO K=1,NEQN
        ET=DABS(Y(K))+DABS(W(K,7))+AAE
        IF (ET .GT. 0.0D0) GO TO 240
C
C       INAPPROPRIATE ERROR TOLERANCE
        IFLAG=5
        RETURN

  240   EE=DABS((-2090.0D0*W(K,1)+(21970.0D0*W(K,4)-15048.0D0*W(K,5)))+
     1                        (22528.0D0*W(K,3)-27360.0D0*W(K,6)))
        EEOET=DMAX1(EEOET,EE/ET)

      end do

      ESTTOL=DABS(H)*EEOET*SCALE/752400.0D0

      IF (ESTTOL .LE. 1.0D0) GO TO 260
C
C     UNSUCCESSFUL STEP
C                       REDUCE THE STEPSIZE , TRY AGAIN
C                       THE DECREASE IS LIMITED TO A FACTOR OF 1/10
C
      HFAILD= .TRUE.
      OUTPUT= .FALSE.
      S=0.1D0
      IF (ESTTOL .LT. 59049.0D0) S=0.9D0/ESTTOL**0.2D0
      H=S*H
      IF (DABS(H) .GT. HMIN) GO TO 200
C
C     REQUESTED ERROR UNATTAINABLE AT SMALLEST ALLOWABLE STEPSIZE
      IFLAG=6
      KFLAG=6
      RETURN
C
C
C     SUCCESSFUL STEP
C                        STORE SOLUTION AT T+H
C                        AND EVALUATE DERIVATIVES THERE
C
  260 T=T+H
      NHI = NHI + 1
      HI(NHI) = H
      DO K=1,NEQN
        Y(K)=W(K,7)
      end do
      IF (NHI.EQ.5) GOTO 320
      A=T
      CALL F(A,Y,W(1,1))
C
C
C                       CHOOSE NEXT STEPSIZE
C                       THE INCREASE IS LIMITED TO A FACTOR OF 5
C                       IF STEP FAILURE HAS JUST OCCURRED, NEXT
C                          STEPSIZE IS NOT ALLOWED TO INCREASE
C
      S=5.0D0
      IF (ESTTOL .GT. 1.889568D-4) S=0.9D0/ESTTOL**0.2D0
      IF (HFAILD) S=DMIN1(S,1.0D0)
      H=DSIGN(DMAX1(S*DABS(H),HMIN),H)
C
C     END OF CORE INTEGRATOR
C
C
C     SHOULD WE TAKE ANOTHER STEP
C
      IF (OUTPUT) GO TO 300
      IF (IFLAG .GT. 0) GO TO 100
C
C
C     INTEGRATION SUCCESSFULLY COMPLETED
C
C     INTERVAL MODE
  300 T=TOUT
      IFLAG=2
      IF (KFLAG.EQ.3) IFLAG = 3
      RETURN

  320 TOUT = T
      IFLAG = 2
      IF (KFLAG.EQ.3) IFLAG = 3
      RETURN
      END
      SUBROUTINE DRKF1S(F,NEQN,Y,T,H,W)

c*********************************************************************72
c
cc DRKF1S
C
C     FEHLBERG FOURTH-FIFTH ORDER RUNGE-KUTTA METHOD
C
      implicit double precision (A-H,O-Z)

      DIMENSION Y(NEQN),W(NEQN,7)
      EXTERNAL F

      CH=H/4.0D0
      DO K=1,NEQN
        W(K,6)=Y(K)+CH*W(K,1)
      end do
      CALL F(T+CH,W(1,6),W(1,2))

      CH=3.0D0*H/32.0D0
      DO K=1,NEQN
        W(K,6)=Y(K)+CH*(W(K,1)+3.0D0*W(K,2))
      end do
      CALL F(T+3.0D0*H/8.0D0,W(1,6),W(1,3))

      CH=H/2197.0D0
      DO K=1,NEQN
        W(K,6) = Y(K) + CH * (1932.0D0 * W(K,1)
     1                + (7296.0D0 * W(K,3) - 7200.0D0 * W(K,2)))
      end do
      CALL F(T+12.0D0*H/13.0D0,W(1,6),W(1,4))

      CH=H/4104.0D0
      DO K=1,NEQN
        W(K,6)=Y(K)+CH*((8341.0D0*W(K,1)-845.0D0*W(K,4))+
     1                    (29440.0D0*W(K,3)-32832.0D0*W(K,2)))
      end do
      CALL F(T+H,W(1,6),W(1,5))

      CH=H/20520.0D0
      DO K=1,NEQN
        W(K,2)=Y(K)+CH*((-6080.0D0*W(K,1)+(9295.0D0*W(K,4)-
     1         5643.0D0*W(K,5)))+(41040.0D0*W(K,2)-28352.0D0*W(K,3)))
      end do
      CALL F(T+H/2.0D0,W(1,2),W(1,6))
C
C     COMPUTE APPROXIMATE SOLUTION AT T+H
C
      CH=H/7618050.0D0
      DO K=1,NEQN
        W(K,7)=Y(K)+CH*((902880.0D0*W(K,1)+(3855735.0D0*W(K,4)-
     1        1371249.0D0*W(K,5)))+(3953664.0D0*W(K,3)+
     2        277020.0D0*W(K,6)))
      end do

      RETURN
      END
      SUBROUTINE DRKFMS(F,NEQN,Y,T,TOUT,RE,AE,EPS,PHI,NPHI,
     1                  HI,NHI,W,IFLAG)

c*********************************************************************72
c
cc DRKFMS
C
C     FEHLBERG FOURTH-FIFTH ORDER RUNGE-KUTTA METHOD
C
      implicit double precision (A-H,O-Z)

      DIMENSION Y(NEQN),PHI(NEQN,NEQN),W(NEQN,7)
      LOGICAL HFAILD,OUTPUT
      EXTERNAL F
C
C     REMIN IS THE MINIMUM ACCEPTABLE VALUE OF RELERR.  ATTEMPTS
C     TO OBTAIN HIGHER ACCURACY WITH THIS SUBROUTINE ARE USUALLY
C     VERY EXPENSIVE AND OFTEN UNSUCCESSFUL.
C
      DATA REMIN/1.D-12/

      LHI = 0
      IST = 0
      IUST = 0
      KFLAG = 0
C
C
C     CHECK INPUT PARAMETERS
C
C
      IF (NEQN .LT. 1) GO TO 10
      IF ((RE .LT. 0.0D0)  .OR.  (AE .LT. 0.0D0)) GO TO 10
      MFLAG=IABS(IFLAG)
      IF (MFLAG.NE.1) GOTO 10

      U26 = 26.0D0*EPS
      GO TO 50
C
C     INVALID INPUT
   10 IFLAG=8
      RETURN

   50 RER=2.0D0*EPS+REMIN
      IF (RE .GE. RER) GO TO 55
C
C     RELATIVE ERROR TOLERANCE TOO SMALL
      RE=RER
      KFLAG=3

   55 DT=TOUT-T

      A=T
      CALL F(A,Y,W(1,1))
      IF (T .NE. TOUT) GO TO 65
      IFLAG=2
      IF (KFLAG.EQ.3) IFLAG = 3
      RETURN

   65 IF (IFLAG.EQ.1) THEN
      H=DABS(DT)
      TOLN=0.
      DO 70 K=1,NEQN
        S = RE * DABS(Y(K))
        TOL = S + DMAX1(RE,AE)
        IF (TOL .LE. 0.D0) GO TO 70
        TOLN=TOL
        YPK=DABS(W(K,1))
        IF (YPK*H**5 .GT. TOL) H=(TOL/YPK)**0.2D0
   70 CONTINUE
      IF (TOLN .LE. 0.0D0) H=0.0D0
      H=DMAX1(H,U26*DMAX1(DABS(T),DABS(DT)))
      HI = H
      ELSE
      H = DMIN1(DABS(HI),DABS(DT))
      end if
C
C     SET STEPSIZE FOR INTEGRATION IN THE DIRECTION FROM T TO TOUT
C
      H=DSIGN(H,DT)
C
C     INITIALIZE OUTPUT POINT INDICATOR
C
      OUTPUT= .FALSE.
C
C     TO AVOID PREMATURE UNDERFLOW IN THE ERROR TOLERANCE FUNCTION,
C     SCALE THE ERROR TOLERANCES
C
      SCALE=2.0D0/RE
      AAE=SCALE*AE
C
C     STEP BY STEP INTEGRATION
C
  100 HFAILD= .FALSE.
C
C     SET SMALLEST ALLOWABLE STEPSIZE
C
      HMIN=U26*DABS(T)
C
C     ADJUST STEPSIZE IF NECESSARY TO HIT THE OUTPUT POINT.
C     LOOK AHEAD TWO STEPS TO AVOID DRASTIC CHANGES IN THE STEPSIZE AND
C     THUS LESSEN THE IMPACT OF OUTPUT POINTS ON THE CODE.
C
      DT=TOUT-T
      IF (DABS(DT) .GE. 2.0D0*DABS(H)) GO TO 200
      IF (DABS(DT) .GT. DABS(H)) GO TO 150
C
C     THE NEXT SUCCESSFUL STEP WILL COMPLETE THE INTEGRATION TO THE
C     OUTPUT POINT
C
      OUTPUT= .TRUE.
      H=DT
      GO TO 200

  150 H=0.5D0*DT
C
C     ADVANCE AN APPROXIMATE SOLUTION OVER ONE STEP OF LENGTH H
C
  200 CALL DRKF1S(F,NEQN,Y,T,H,W)
      IST = IST + 1
C
C     COMPUTE AND TEST ALLOWABLE TOLERANCES VERSUS LOCAL ERROR ESTIMATES
C     AND REMOVE SCALING OF TOLERANCES. NOTE THAT RELATIVE ERROR IS
C     MEASURED WITH RESPECT TO THE AVERAGE OF THE MAGNITUDES OF THE
C     SOLUTION AT THE BEGINNING AND END OF THE STEP.
C
      EEOET=0.0D0
      DO K = 1, NEQN

        ET=DABS(Y(K))+DABS(W(K,7))+AAE
        IF (ET .GT. 0.0D0) GO TO 240
C
C       INAPPROPRIATE ERROR TOLERANCE
        IFLAG=5
        RETURN

  240   EE=DABS((-2090.0D0*W(K,1)+(21970.0D0*W(K,4)-15048.0D0*W(K,5)))+
     1                        (22528.0D0*W(K,3)-27360.0D0*W(K,6)))
        EEOET=DMAX1(EEOET,EE/ET)

      end do

      ESTTOL=DABS(H)*EEOET*SCALE/752400.0D0

      IF (ESTTOL .LE. 1.0D0) GO TO 260
C
C
C     UNSUCCESSFUL STEP
C                       REDUCE THE STEPSIZE , TRY AGAIN
C                       THE DECREASE IS LIMITED TO A FACTOR OF 1/10
C
      IUST = IUST + 1
      HFAILD= .TRUE.
      OUTPUT= .FALSE.
      S=0.1D0
      IF (ESTTOL .LT. 59049.0D0) S=0.9D0/ESTTOL**0.2D0
      H=S*H
      IF (DABS(H) .GT. HMIN) GO TO 200
C
C     REQUESTED ERROR UNATTAINABLE AT SMALLEST ALLOWABLE STEPSIZE
      IFLAG=6
      KFLAG=6
      RETURN
C
C
C     SUCCESSFUL STEP
C                        STORE SOLUTION AT T+H
C                        AND EVALUATE DERIVATIVES THERE
C
  260 LHI = LHI + 1
      DO 270 K = 1 , NEQN
        Y(K) = W(K,7)
  270 CONTINUE
      IF (NPHI.GT.0) THEN
      DO 290 K = 1 , NPHI
        CALL F(T,PHI(1,K),W(1,1))
        CALL DRKF1S(F,NEQN,PHI(1,K),T,H,W)
        DO 280 L = 1 , NEQN
          PHI(L,K) = W(L,7)
  280   CONTINUE
  290 CONTINUE
      end if
      T = T + H
      IF (LHI.EQ.NHI) GOTO 320
      A=T
      CALL F(A,Y,W(1,1))
C
C
C                       CHOOSE NEXT STEPSIZE
C                       THE INCREASE IS LIMITED TO A FACTOR OF 5
C                       IF STEP FAILURE HAS JUST OCCURRED, NEXT
C                          STEPSIZE IS NOT ALLOWED TO INCREASE
C
      S=5.0D0
      IF (ESTTOL .GT. 1.889568D-4) S=0.9D0/ESTTOL**0.2D0
      IF (HFAILD) S=DMIN1(S,1.0D0)
      H=DSIGN(DMAX1(S*DABS(H),HMIN),H)
      HI = H
C
C     END OF CORE INTEGRATOR
C
C     SHOULD WE TAKE ANOTHER STEP
C
      IF (OUTPUT) GO TO 300
      GO TO 100
C
C     INTEGRATION SUCCESSFULLY COMPLETED
C
C     INTERVAL MODE
  300 T=TOUT
      IFLAG=2
      IF (KFLAG.EQ.3) IFLAG = 3
      RETURN
C
  320 TOUT = T
      IFLAG = 2
      IF (KFLAG.EQ.3) IFLAG = 3

      RETURN
      END
      SUBROUTINE DRKFGG(FDIF,N,T1,T2,X,PHI,HS,NHS,W)

c*********************************************************************72
c
cc DRKFGG
C
      implicit double precision (A-H,O-Z)

      DIMENSION X(N),PHI(N,N),HS(NHS),W(N,7)
      EXTERNAL FDIF

      T0 = T1
      ST = (T2-T1) / NHS
      DO 1000 I = 1 , NHS
        HS(I) = ST
 1000 CONTINUE
      CALL DRKFGS(FDIF,N,X,T1,HS,NHS,W)
      DO 1100 I = 1 , N
        T1 = T0
        CALL DRKFGS(FDIF,N,PHI(1,I),T1,HS,NHS,W)
 1100 CONTINUE
      RETURN
C     END OF DRKFGG
      END
      SUBROUTINE ERRHAN(IERROR,ER,NEG)

c*********************************************************************72
c
cc ERRHAN handles errors.
C
      implicit double precision (A-H,O-Z)

      DIMENSION ER(5)
C
      IF (IERROR.EQ.100) WRITE(*,100)
      IF (IERROR.EQ.101) WRITE(*,101)
      IF (IERROR.EQ.103) WRITE(*,103)
      IF (IERROR.EQ.105) WRITE(*,105)
      IF (IERROR.EQ.106) WRITE(*,106)
      IF (IERROR.EQ.120) WRITE(*,120)
      IF (IERROR.EQ.121) WRITE(*,121)
      IF (IERROR.EQ.122) WRITE(*,122)
      IF (IERROR.EQ.123) WRITE(*,123) NEG
      IF (IERROR.EQ.200) WRITE(*,200)
      IF (IERROR.EQ.213) WRITE(*,213) ER(1)
      IF (IERROR.EQ.215) WRITE(*,215)
      IF (IERROR.EQ.218) WRITE(*,218)
      IF (IERROR.EQ.219) WRITE(*,219) NEG
      IF (IERROR.EQ.230) WRITE(*,230)
      IF (IERROR.EQ.231) WRITE(*,231)
      IF (IERROR.EQ.240) WRITE(*,240) ER(5)
      IF (IERROR.EQ.250) WRITE(*,250)
      IF (IERROR.EQ.260) WRITE(*,260)
      RETURN
  100 FORMAT(' INPUT ERROR: N<1 OR IHOM<0 OR NRTI<0 OR NTI<5 OR A=B OR',
     1 ' NU<N*(N+1)/2')
  101 FORMAT(' INPUT ERROR: ER(1), ER(2) OR ER(3) IS NEGATIVE')
  103 FORMAT(' INPUT ERROR: LW<8*N+2*N*N OR LIW<3*N')
  105 FORMAT(' INPUT ERROR: N<1 OR NRTI<0 OR NTI<3 OR A=B OR',
     1 ' NU<N*(N+1)/2')
  106 FORMAT(' INPUT ERROR: LW<7*N+3*N*NTI+4*N*N OR LIW<3*N+NTI')
  120 FORMAT(' INPUT ERROR: NRTI=1, BUT SUPPLIED OUTPUT POINTS IN TI ',
     1 /,' ARE NOT IN STRICT MONOTONE ORDER')
  121 FORMAT(' INPUT ERROR: NRTI=1, BUT A OR B ARE NOT INCLUDED IN THE',
     1 /,' GIVEN OUTPUTPOINTS IN TI')
  122 FORMAT(' INPUT ERROR: NTI TOO SMALL !')
  123 FORMAT(' INPUT ERROR: LWG IS LESS THAN THE NUMBER OF OUTPUT  ',
     1 'POINTS.',/,' LWG MUST BE GREATER THAN:',I10)
  200 FORMAT(' INCREMENT ON MINOR SHOOTING INTERVAL IS',
     1 ' GREATER THAN AMP')
  213 FORMAT(' ER(1) TOO SMALL, CHANGED INTO:',D16.9)
  215 FORMAT(' PURE RELATIVE ERROR TEST IMPOSSIBLE; MUST USE NON-ZERO',
     1 /,' ABSOLUTE TOLERANCE')
  218 FORMAT(' N, ER(1) OR ER(2) IS NEGATIVE')
  219 FORMAT(' ARRAY WG IS TOO SMALL. ESTIMATED VALUE FOR LWG:',I10)
  230 FORMAT(' NECESSARY DAMPING FACTOR IS TOO SMALL; YOU MAY TRY MORE',
     1 /,' OUTPUT POINTS OR START WITH A BETTER APPROXIMATION OF',/,
     2 ' THE SOLUTION.')
  231 FORMAT(' NUMBER OF ITERATIONS > ITLIM')
  240 FORMAT(' BAD DICHOTOMY; AMPLIFICATION FACTOR = ',D12.5)
  250 FORMAT(' ONE OF THE UPPER TRIANGULAR MATRIX OF THE MULTIPLE',
     1 /,' SHOOTING RECURSION IS SINGULAR')
  260 FORMAT(' PROBLEM TOO ILL-CONDITIONED WITH RESPECT TO THE BC')
      END
      subroutine timestamp ( )

c*********************************************************************72
c
cc TIMESTAMP prints out the current YMDHMS date as a timestamp.
c
c  Licensing:
c
c    This code is distributed under the MIT license.
c
c  Modified:
c
c    12 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    None
c
      implicit none

      character * ( 8 ) ampm
      integer d
      character * ( 8 ) date
      integer h
      integer m
      integer mm
      character * ( 9 ) month(12)
      integer n
      integer s
      character * ( 10 ) time
      integer y

      save month

      data month /
     &  'January  ', 'February ', 'March    ', 'April    ', 
     &  'May      ', 'June     ', 'July     ', 'August   ', 
     &  'September', 'October  ', 'November ', 'December ' /

      call date_and_time ( date, time )

      read ( date, '(i4,i2,i2)' ) y, m, d
      read ( time, '(i2,i2,i2,1x,i3)' ) h, n, s, mm

      if ( h .lt. 12 ) then
        ampm = 'AM'
      else if ( h .eq. 12 ) then
        if ( n .eq. 0 .and. s .eq. 0 ) then
          ampm = 'Noon'
        else
          ampm = 'PM'
        end if
      else
        h = h - 12
        if ( h .lt. 12 ) then
          ampm = 'PM'
        else if ( h .eq. 12 ) then
          if ( n .eq. 0 .and. s .eq. 0 ) then
            ampm = 'Midnight'
          else
            ampm = 'AM'
          end if
        end if
      end if

      write ( *, 
     &  '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) 
     &  d, month(m), y, h, ':', n, ':', s, '.', mm, ampm

      return
      end

