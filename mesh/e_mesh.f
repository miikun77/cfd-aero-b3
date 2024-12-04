C======= <<< E_MESH : MESH GENERATION BY ELLIPTIC EQUATIONS >>> =========
C=======================================================================
      PARAMETER ( JDIM=501, KDIM=200 )
      DIMENSION X(JDIM,KDIM), Y(JDIM,KDIM)
      COMMON/BASE/MTIN,MTOUT

C                                           << INPUT & SURFACE MESH >>
      CALL EINPUT
C                                           << READ INITIAL MESH. >>
      CALL MESH_IO(1, JDIM,KDIM, JMAX,KMAX, X,Y, MTIN)
C                                           << SOLVE POISSON EQS. >>
      CALL ECORE(JDIM,KDIM, JMAX,KMAX, X,Y)
C                                           << WRITE MESH         >>
      CALL MESH_IO(2, JDIM,KDIM, JMAX,KMAX, X,Y, MTOUT)

      STOP
      END
C=======================================================================
      SUBROUTINE EINPUT
C=======================================================================
      COMMON/CSPACE/MAXIT,DSB(2),A1,A2,B1,B2,OMGP(2),OMGQ(2),OMGA
      COMMON/BASE/MTIN,MTOUT
      COMMON/CJFREE/JSFREE,JEFREE
CIN   NAMELIST/INPUT/MTIN,MTOUT,MAXIT,DS1,DS2,A1,A2,B1,B2,
CIN  &               OMGP,OMGQ,OMGA,JSFREE,JEFREE
C................................ << READ INPUT DATA >>
      MTIN   = 11
      MTOUT  = 12
      MAXIT  = 500
      DS1    = 0.01
      DS2    = 2.0
      A1     = 0.7
      A2     = 0.7
      B1     = 0.7
      B2     = 0.7
      JSFREE = 1
      JEFREE = 1
      DO 10 I=1,2
      OMGP(I) = 0.06
  10  OMGQ(I) = 0.06
      OMGA    = 1.0
CIN   READ(5,INPUT)
          DSB(1)  = DS1
          DSB(2)  = DS2
      WRITE(6,*)' INPUT  FILE NO. : MT1   =',MT1
      WRITE(6,*)' OUTPUT FILE NO. : MT2   =',MT2
      WRITE(6,*)' SPACING AT WALL : DS1   =',DS1
      WRITE(6,*)' SPACING AT OUTER: DS2   =',DS2
      RETURN
      END
C=======================================================================
      SUBROUTINE ECORE(JDIM,KDIM, JMAX,KMAX, X,Y)
C=======================================================================
      PARAMETER (LD=301)
      DIMENSION X(JDIM,KDIM), Y(JDIM,KDIM)
     & ,XXB(LD,2), YXB(LD,2), XEB(LD,2), YEB(LD,2), RJB(LD,2)
     & ,FB1(LD,2), FB2(LD,2), FB3(LD,2),  PS(LD,2),  QS(LD,2)
     & ,XX(LD),   YX(LD), XEE(LD), YEE(LD), XXE(LD), YXE(LD) 
     & ,XYJS(LD), XR(LD),  YR(LD),   A(LD),   B(LD),   C(LD), WORK(LD)
      COMMON/CSPACE/MAXIT,DSB(2),A1,A2,B1,B2,OMGP(2),OMGQ(2),OMGA
      COMMON/CJFREE/JSFREE,JEFREE
C...............................................
      DO 50     I = 1, 2
                K = 1
      IF(I.EQ.2)K = KMAX
        DO 10    J  = 2, JMAX-1
        XXB(J,I)    = 0.5 *(X(J+1,K)-X(J-1,K))
  10    YXB(J,I)    = 0.5 *(Y(J+1,K)-Y(J-1,K))
        XXB(1,I)    = X(   2,K) -X(   1,K)
        YXB(1,I)    = Y(   2,K) -Y(   1,K)
        XXB(JMAX,I) = X(JMAX,K) -X(JMAX-1,K)
        YXB(JMAX,I) = Y(JMAX,K) -Y(JMAX-1,K)
          DO 20    J = 1,JMAX
          XEB(J,I) = -DSB(I)*YXB(J,I)/SQRT(XXB(J,I)**2+YXB(J,I)**2)
          YEB(J,I) =  DSB(I)*XXB(J,I)/SQRT(XXB(J,I)**2+YXB(J,I)**2)
  20      RJB(J,I) = 1.0/(XXB(J,I)*YEB(J,I)-XEB(J,I)*YXB(J,I))
            DO 30  J = 2,JMAX-1
            XXXB     = X(J+1,K) -2.*X(J,K) +X(J-1,K)
            YXXB     = Y(J+1,K) -2.*Y(J,K) +Y(J-1,K)
            XXEB     = 0.5*(XEB(J+1,I)-XEB(J-1,I))
            YXEB     = 0.5*(YEB(J+1,I)-YEB(J-1,I))
            ALPB     = XEB(J,I)**2+YEB(J,I)**2
            GAMB     = XXB(J,I)**2+YXB(J,I)**2
            BETB     = XXB(J,I)*XEB(J,I)+YXB(J,I)*YEB(J,I)
            FB1(J,I) = -GAMB*RJB(J,I)**2
            FB2(J,I) = ( -ALPB*XXXB+2.0*BETB*XXEB )*RJB(J,I)**2
  30        FB3(J,I) = ( -ALPB*YXXB+2.0*BETB*YXEB )*RJB(J,I)**2
              DO 40   J = 1,JMAX
              PS(J,I) = 0.
  40          QS(J,I) = 0.
  50  CONTINUE
C........................................... MAIN ITERATION LOOP ......
      IF(MAXIT.NE.1 .AND. MOD(MAXIT,2).EQ.1) MAXIT=MAXIT+1
      DO 200 ITER=1,MAXIT
       IF(MOD(ITER,2).EQ.1) RESI = 0.
       DO 90 I=1,2
       IF(I.EQ.1)THEN
        DO 60  J = 2,JMAX-1
        XEE(J)   = .5*(-7.*X(J,1)+8.*X(J,2)-X(J,3)) -3.*XEB(J,1)
  60    YEE(J)   = .5*(-7.*Y(J,1)+8.*Y(J,2)-Y(J,3)) -3.*YEB(J,1)
       ELSE
        DO 70 J= 2,JMAX-1
        XEE(J)=.5*(-7.*X(J,KMAX)+8.*X(J,KMAX-1)-X(J,KMAX-2))+3.*XEB(J,2)
  70    YEE(J)=.5*(-7.*Y(J,KMAX)+8.*Y(J,KMAX-1)-Y(J,KMAX-2))+3.*YEB(J,2)
       ENDIF
        DO 80   J = 2,JMAX-1
          RB1     = FB2(J,I)+FB1(J,I)*XEE(J)
          RB2     = FB3(J,I)+FB1(J,I)*YEE(J)
          DP      = ( YEB(J,I)*RB1-XEB(J,I)*RB2)*RJB(J,I) -PS(J,I)
          DQ      = (-YXB(J,I)*RB1+XXB(J,I)*RB2)*RJB(J,I) -QS(J,I)
          DPP     = OMGP(I)*ABS(DP)/(ABS(PS(J,I))+0.001)
          DQQ     = OMGQ(I)*ABS(DQ)/(ABS(QS(J,I))+0.001)
          PS(J,I) = PS(J,I) + DP *OMGP(I)/MAX(1.,MIN(DPP,10.))
  80      QS(J,I) = QS(J,I) + DQ *OMGQ(I)/MAX(1.,MIN(DQQ,10.))
  90  CONTINUE
      DO 150 KK=2,KMAX-1
        IF(MOD(ITER,2).EQ.0)THEN
          K    = KMAX +1 -KK
        ELSE
          K    = KK
        ENDIF
          KP   = K+1
          KM   = K-1
        DO 100 J = 2,JMAX-1
          JP     = J+1
          JM     = J-1
          XX(J)  = 0.5*(X(JP,K)-X(JM,K))
          YX(J)  = 0.5*(Y(JP,K)-Y(JM,K))
          XXE(J) = .25*(X(JP,KP)-X(JP,KM)-X(JM,KP)+X(JM,KM))
 100      YXE(J) = .25*(Y(JP,KP)-Y(JP,KM)-Y(JM,KP)+Y(JM,KM))
        XX(   1)  = X(   2,K)-X(   1,  K)
        YX(   1)  = Y(   2,K)-Y(   1,  K)
        XX(JMAX)  = X(JMAX,K)-X(JMAX-1,K)
        YX(JMAX)  = Y(JMAX,K)-Y(JMAX-1,K)
        XXE(1)    = .5*(X(2,KP)-X(2,KM)-X(1,KP)+X(1,KM))
        YXE(1)    = .5*(Y(2,KP)-Y(2,KM)-Y(1,KP)+Y(1,KM))
        XXE(JMAX) = .5*(X(JMAX,KP)-X(JMAX,KM)-X(JMAX-1,KP)+X(JMAX-1,KM))
        YXE(JMAX) = .5*(Y(JMAX,KP)-Y(JMAX,KM)-Y(JMAX-1,KP)+Y(JMAX-1,KM))
          DO 110  J = 1,JMAX
            XE      = 0.5*(X(J,KP)-X(J,KM))
            YE      = 0.5*(Y(J,KP)-Y(J,KM))
            XYJS(J) = (XX(J)*YE-XE*YX(J))**2
            ALP     = XE**2 + YE**2
            BET     = XX(J)*XE + YX(J)*YE
            GAM     = XX(J)**2 + YX(J)**2
            A(J)    =  ALP
            B(J)    = -2.*ALP-2.*GAM
            C(J)    =  ALP
            XR(J)   = 2.*BET*XXE(J) -GAM*(X(J,KP)+X(J,KM))
 110        YR(J)   = 2.*BET*YXE(J) -GAM*(Y(J,KP)+Y(J,KM))
          DO 120 J=1,JMAX
            P      = PS(J,1)*EXP(-A1*(K-1)) + PS(J,2)*EXP(-A2*(KMAX-K))
            Q      = QS(J,1)*EXP(-B1*(K-1)) + QS(J,2)*EXP(-B2*(KMAX-K))
            B(J)   = B(J) - XYJS(J)*ABS(P)
            IF(P.GT.0.) THEN
              C(J)   = C(J)+XYJS(J)*P
            ELSE
              A(J)   = A(J)-XYJS(J)*P
            ENDIF
              B(J)   = B(J)-XYJS(J)*ABS(Q)
              IF(Q.GT.0.) THEN
                XR(J)  = XR(J) - XYJS(J)*Q*X(J,KP)
                YR(J)  = YR(J) - XYJS(J)*Q*Y(J,KP)
              ELSE
                XR(J)  = XR(J) + XYJS(J)*Q*X(J,KM)
                YR(J)  = YR(J) + XYJS(J)*Q*Y(J,KM)
              ENDIF
 120      CONTINUE
          XR(   1) = XR(   1) - A(   1)*X(   1,K)
          XR(JMAX) = XR(JMAX) - C(JMAX)*X(JMAX,K)
          YR(   1) = YR(   1) - A(   1)*Y(   1,K)
          YR(JMAX) = YR(JMAX) - C(JMAX)*Y(JMAX,K)
        CALL TRIDSL(1,JMAX,A,B,C,XR,WORK)
        CALL TRIDSL(1,JMAX,A,B,C,YR,WORK)
          DO 130 J=2,JMAX-1
            DX     = XR(J)-X(J,K)
            DY     = YR(J)-Y(J,K)
            X(J,K) = X(J,K) + OMGA*DX
            Y(J,K) = Y(J,K) + OMGA*DY
 130        RESI   = RESI + DX**2 + DY**2
 150  CONTINUE

      IF(JSFREE.EQ.1) THEN
C       ........... REDESTRIBUTE ON J=1
        DXS  = X( 2,KMAX)-X( 2,1)
        DYS  = Y( 2,KMAX)-Y( 2,1)
        DO 160 K=1,KMAX
        IF(ABS(DXS).GT.1.E-5)THEN
          X(   1,K) = X(   1,1)
     &              + (X(1,KMAX)-X(1,1))*(X(2,K)-X(2,1))/DXS
        ENDIF
        IF(ABS(DYS).GT.1.E-5)THEN
          Y(   1,K) = Y(   1,1)
     &               + (Y(1,KMAX)-Y(1,1))*(Y(2,K)-Y(2,1))/DYS
         ENDIF
 160    CONTINUE
      ENDIF

      IF(JEFREE.EQ.1) THEN
C       ........... REDESTRIBUTE ON J=JMAX
        JM   = JMAX-1
        DXE  = X(JM,KMAX)-X(JM,1)
        DYE  = Y(JM,KMAX)-Y(JM,1)
        DO 170 K=1,KMAX
        IF(ABS(DXE).GT.1.E-5)THEN
           X(JMAX,K) = X(JMAX,1)
     &             + (X(JMAX,KMAX)-X(JMAX,1))*(X(JM,K)-X(JM,1))/DXE
        ENDIF
        IF(ABS(DYE).GT.1.E-5)THEN
          Y(JMAX,K) = Y(JMAX,1)
     &             + (Y(JMAX,KMAX)-Y(JMAX,1))*(Y(JM,K)-Y(JM,1))/DYE
        ENDIF
 170    CONTINUE
      ENDIF

        IF(MOD(ITER,2).EQ.0)THEN
          RESI   = SQRT(RESI/FLOAT(JMAX*KMAX))
          IF(RESI .GT. 1.E5) GO TO 900
          WRITE(6,*)'ITERATION NO. & RESIDUAL =',ITER,RESI
          IF( RESI .LT. 0.1*MIN(DSB(1),DSB(2)) ) GO TO 220
        ENDIF
 200  CONTINUE
      WRITE(6,*)'? NOT CONVERGED'
 220  RETURN
 900  WRITE(6,*)'?? FAILED : CHANGE OMGA ETC.'
      STOP
      END
C================== TRI-DIAGONAL SYSTEM OF LINEAR EQUATIONS ============
      SUBROUTINE TRIDSL(L,M,A,B,C,D,W)
C=======================================================================
      DIMENSION A(1),B(1),C(1),D(1),W(1)
      W(L)    = C(L)/B(L)
      D(L)    = D(L)/B(L)
      DO 10 I = L+1,M
      TEMP    = 1./(B(I)-A(I)*W(I-1))
      W(I)    = C(I)*TEMP
  10  D(I)    = (D(I)-A(I)*D(I-1))*TEMP
      DO 20 J = L+1,M
      I       = M+L-J
  20  D(I)    = D(I)-W(I)*D(I+1)
      RETURN
      END
