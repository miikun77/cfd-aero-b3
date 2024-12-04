C==============<< DISTRIBUTE MESH POINTS ON A 2D CURVE >>===============
      SUBROUTINE GRCV2D(NIN,XIN,YIN, N,K,DL1,DL2, XOUT,YOUT)
C=======================================================================
      DIMENSION XIN(NIN),YIN(NIN),XOUT(N),YOUT(N)
     &         ,SN(301),ARCS(301),CSPX(301),CSPY(301)
C................................. XIN, YIN ---> ARCS (ARC LENGTH)
      ARCS(1) = 0.
      DO 10  J=2,NIN
  10  ARCS(J) =ARCS(J-1)+SQRT((XIN(J)-XIN(J-1))**2+(YIN(J)-YIN(J-1))**2)
      IF(NIN.GT.2) THEN
        CALL SPLINE(NIN,ARCS,XIN,CSPX)
        CALL SPLINE(NIN,ARCS,YIN,CSPY)
      ENDIF
C..................................<<CLUSTERING>>
      SL1   = DL1/ARCS(NIN)
      SL2   = DL2/ARCS(NIN)
      IF(K.NE.0) THEN
         SL1   = DL1/FLOAT(N-1)
         SL2   = DL2/FLOAT(N-1)
      ENDIF
      IF(SL2.LE.0.) THEN
        CALL CLST1(N,SL1,SN)
      ELSE
        CALL CLST2(N,SL1,SL2,SN)
      ENDIF
C.......................................... ARCS ---> XOUT, YOUT
      IF (NIN.EQ.2) THEN
        DO 20 L=1,N
         XOUT(L)  = (1.-SN(L))*XIN(1)+SN(L)*XIN(2)
  20     YOUT(L)  = (1.-SN(L))*YIN(1)+SN(L)*YIN(2)
      ELSE
        DO 30 L=1,N
         SND   = SN(L)*ARCS(NIN)
         CALL SPEVAL(NIN,ARCS,XIN,CSPX,SND,XOUT(L))
  30     CALL SPEVAL(NIN,ARCS,YIN,CSPY,SND,YOUT(L))
      ENDIF
      RETURN
      END
C==============<< DISTRIBUTE MESH POINTS ON A 3D CURVE >>==============
      SUBROUTINE GRCV3D(NIN,XIN,YIN,ZIN, N,K,DL1,DL2, XOUT,YOUT,ZOUT)
C=======================================================================
      DIMENSION XIN(NIN),YIN(NIN),ZIN(NIN),XOUT(N),YOUT(N),ZOUT(N)
     &         ,SN(301),ARCS(301),CSPX(301),CSPY(301),CSPZ(301)
C................................. XIN, YIN, ZIN ---> ARCS (ARC LENGTH)
      ARCS(1) = 0.
      DO 10  J=2,NIN
  10  ARCS(J) =ARCS(J-1)+SQRT((XIN(J)-XIN(J-1))**2+(YIN(J)-YIN(J-1))**2 
     &                       +(ZIN(J)-ZIN(J-1))**2)
      IF(NIN.GT.2) THEN
        CALL SPLINE(NIN,ARCS,XIN,CSPX)
        CALL SPLINE(NIN,ARCS,YIN,CSPY)
        CALL SPLINE(NIN,ARCS,ZIN,CSPZ)
      ENDIF
C..................................<<CLUSTERING>>
      SL1   = DL1/ARCS(NIN)
      SL2   = DL2/ARCS(NIN)
      IF(K.NE.0) THEN
         SL1   = DL1/FLOAT(N-1)
         SL2   = DL2/FLOAT(N-1)
      ENDIF
      IF(SL2.LE.0.) THEN
        CALL CLST1(N,SL1,SN)
      ELSE
        CALL CLST2(N,SL1,SL2,SN)
      ENDIF
C....................................... ARCS ---> XOUT, YOUT, ZOUT
      IF (NIN.EQ.2) THEN
        DO 20 L=1,N
         XOUT(L)  = (1.-SN(L))*XIN(1)+SN(L)*XIN(2)
         YOUT(L)  = (1.-SN(L))*YIN(1)+SN(L)*YIN(2)
  20     ZOUT(L)  = (1.-SN(L))*ZIN(1)+SN(L)*ZIN(2)
      ELSE
        DO 30 L=1,N
         SND   = SN(L)*ARCS(NIN)
         CALL SPEVAL(NIN,ARCS,XIN,CSPX,SND,XOUT(L))
         CALL SPEVAL(NIN,ARCS,YIN,CSPY,SND,YOUT(L))
  30     CALL SPEVAL(NIN,ARCS,ZIN,CSPZ,SND,ZOUT(L))
      ENDIF
      RETURN
      END
C=======================================================================
      SUBROUTINE SPLINE(NIN,X,Y,FDP)
C=======================================================================
      DIMENSION X(NIN),Y(NIN),FDP(NIN),A(301),B(301),C(301),R(301)
      N1     = NIN - 1
      C(1)   = X(2) - X(1)
      DO 10 I=2,N1
        C(I)  = X(I+1) - X( I )
        A(I)  = C(I-1)
        B(I)  = 2.*(A(I) + C(I))
  10    R(I)  = 6.*((Y(I+1) - Y(I))/C(I) - (Y(I) - Y(I-1))/C(I-1))
      B(2)   = B( 2) + C( 1)
      B(N1)  = B(N1) + C(N1)
      DO 20 I = 3,N1
        T       = A(I)/B(I-1)
        B(I)    = B(I) - T * C(I-1)
  20    R(I)    = R(I) - T * R(I-1)
      FDP(N1)   = R(N1)/B(N1)
      DO 30    I = 2,NIN-2
        NMI      = NIN - I
  30    FDP(NMI) = (R(NMI) - C(NMI)*FDP(NMI+1))/B(NMI)
      FDP(  1)   = FDP( 2)
      FDP(NIN)   = FDP(N1)
      RETURN
      END
C=======================================================================
      SUBROUTINE SPEVAL(NIN,X,Y,FDP,XX,F)
C=======================================================================
      DIMENSION X(NIN),Y(NIN),FDP(NIN)
      DO 10 I=1,NIN-1
  10  IF (XX.LE.X(I+1)) GO TO 20
  20  DXM   = XX - X(I)
      DXP   = X(I+1) - XX
      DEL   = X(I+1) - X(I)
      F     = FDP( I )*DXP*(DXP*DXP/DEL - DEL)/6. + Y( I )*DXP/DEL
     &       +FDP(I+1)*DXM*(DXM*DXM/DEL - DEL)/6. + Y(I+1)*DXM/DEL
      RETURN
      END
C====<< DESTRIBUTE POINTS ON UNIT-LENGTH LINE BY ONE-SIDED SPACING >>===
!  長さ１の線上でN個の点を、左端の幅( d )だけ指定して指数的に分布さえるサブルーチン

      SUBROUTINE CLST1(N,DX1,X)
      
C X(1)=0., X(2)=DX1, X(3)=X(2)+DX1*E, , , X(N)=X(N-1)+DX1*E**(N-2)=1.0

!  入力 n（点の数）とd（左端の幅）から、1.0 = d*(1.-e**(n-1))/(1.-e) を満たすeを求める。
!  そのeを用いて以下の点分布を出力する。
!    x(1)=0.,  x(2)=d,  x(3)=x(2)+d*e, , ,  x(N)=x(N-1)+d*e**(N-2)=1.0

!	Input : N = number of points
!			d = x(2)-x(1): Spacing on the LHS of a unit-length line.
!   Output : coordinates:		x(i), i=1,N
!            streaching factor: e
C=======================================================================
      DIMENSION X(N)
        N1   = N-1
        DX   = 1./FLOAT(N1)
          IF (ABS(DX-DX1).LT.1.E-5) THEN
           E  = 1.
           GO TO 50
          ENDIF
            E1  = 1.2
  10        E   = E1
            IF (DX1.GT.DX) E=1./E
      ITER  = 0
  30  ITER  = ITER+1
        EP   = E
        E    = EP-(DX1*EP**N1-EP+1.-DX1)/(DX1*FLOAT(N1)*E**(N1-1)-1.)
        IF ( ABS(E-EP) .LT. 1.E-5 ) GO TO 40
        IF ( ITER .LT. 100 ) GO TO 30
        WRITE(6,*)'? ITERATION FAILED IN SUB.CLST1'
  40      IF (ABS(E-1.).LT.1.E-5) THEN
           E1  = E1**2
           GO TO 10
          ENDIF
  50        X(1)   = 0.
            DX     = DX1
            DO 60 I=2,N
            X(I)   = X(I-1)+DX
  60        DX     = DX*E
      RETURN
      END
C=<< DISTRIBUTE N-POINTS ON A UNIT LENGTH LINE WITH BOTH END SPACINGS>>=
C INPUTS  : N       = NUMBER OF MESH POINTS,
C           SP1,SP2 = END SPACINGS (0.<SP1,SP2<0.3)
C OUTPUTS : XOUT(N) = COORDINATE

      SUBROUTINE CLST2(N,SP1,SP2,XOUT)
C=======================================================================
      DIMENSION XIN(1000),XOUT(N)
      IF(SP1.LE.0. .OR. SP2.LE.0.) THEN
        WRITE(6,*)'? WRONG VALUES IN SP1, SP2  :   SP1, SP2=',SP1,SP2
        WRITE(6,*)'? SP1, SP2 MUST BE 0.<SP1, SP2< 0.3'
        STOP
      ENDIF
      SMAL   = 0.01*MIN(SP1,SP2)
      N1     = N-1
      DO 10 I= 1,N
  10  XIN(I) = FLOAT(I-1)/FLOAT(N1)
      S11   = 1./(SP1*FLOAT(N1))
      S21   = 1./(SP2*FLOAT(N1))
      CALL STRET (XIN,XOUT,S11,S21,N)
      DX1   = XOUT(2)
      DX2   = XOUT(N)-XOUT(N1)
      IF( ABS(SP1-DX1).LE.SMAL .AND. ABS(SP2-DX2).LE.SMAL ) GO TO 30
C <<< IMPROVE ACCURACY OF END-SPACING CONTROL. >>>
      ITR    = 0
  20  ITR    = ITR+1
      DS1    = S11*0.1
      DS2    = S21*0.1
      S12    = S11+DS1
      S22    = S21+DS2
      CALL STRET (XIN,XOUT,S12,S21,N)
      DX1DS1 = (XOUT(2)-DX1)/DS1
      DX2DS1 = (XOUT(N)-XOUT(N1)-DX2)/DS1
      CALL STRET (XIN,XOUT,S11,S22,N)
      DX1DS2 = (XOUT(2)-DX1)/DS2
      DX2DS2 = (XOUT(N)-XOUT(N1)-DX2)/DS2
      DF1    = SP1-DX1
      DF2    = SP2-DX2
      D      = DX1DS1*DX2DS2-DX1DS2*DX2DS1
      DS1    = (DF1*DX2DS2-DF2*DX1DS2)/D
      DS2    = (DF2*DX1DS1-DF1*DX2DS1)/D
      S12    = S11+DS1
      S22    = S21+DS2
      CALL STRET (XIN,XOUT,S12,S22,N)
      DX1    = XOUT(2)
      DX2    = XOUT(N)-XOUT(N1)
      IF(ABS(SP1-DX1).LE.SMAL .AND. ABS(SP2-DX2).LE.SMAL) GO TO 30
      S11    = S12
      S21    = S22
      IF(ITR.LE.20) GO TO 20
  30  XOUT(1) = 0.
      XOUT(N) = 1.
      RETURN
      END
C== << TWO-SIDED STRETCHING FUNCTION OF VINOKUR : JCP 50, 2 (1983)>> ===
      SUBROUTINE STRET(X,Y,S1,S2,N)
C ======================================================================
      DIMENSION X(N),Y(N)
      B      = SQRT(S1*S2)
      A      = B/S2
      IF(B .LT. 0.999) THEN
        DZ     = FASIN(B)
        DO 10 J= 1,N
        TANX   = TAN(DZ*X(J))
  10    Y(J)   = TANX/(A*SIN(DZ)+(1.-A*COS(DZ))*TANX)
      ELSEIF(B .GT. 1.001) THEN
        DZ     = FASINH(B)
        DO 20 J= 1,N
        TANHX  = TANH(DZ*X(J))
  20    Y(J)   = TANHX/(A*SINH(DZ)+(1.-A*COSH(DZ))*TANHX)
      ELSE
        DO 30 J=1,N
        U      = X(J)*(1.+2.*(B-1.)*(X(J)-.5)*(1.-X(J)))
  30    Y(J)   = U/(A+(1.-A)*U)
      ENDIF
      RETURN
      END
C=================== << INVERSION OF Y=SIN(X)/X >> =====================
      FUNCTION FASIN(Y)
C ======================================================================
      DATA A3,A4,A5,A6/-2.6449341,6.794732,-13.205501,11.726095/
      DATA B2,B3,B4,B5/.057321429,.048774238,-.053337753,.075845134/
      DATA PI/3.1415927/
      IF(Y .LE. 0.26938972) THEN
        FASIN = PI*((((((A6*Y+A5)*Y+A4)*Y+A3)*Y+1.)*Y-1.)*Y+1.)
      ELSE
        YB    = 1.-Y
        FASIN = SQRT(6.*YB)*(((((B5*YB+B4)*YB+B3)*YB+B2)*YB+.15)*YB+1.)
      ENDIF
      RETURN
      END
C=================== << INVERSION OF Y=SINH(X)/X >> ====================
      FUNCTION FASINH(Y)
C ======================================================================
      DATA A2,A3,A4,A5/.057321429,-.024907295,.0077424461,-.0010794123/
      DATA B0,B1,B2,B3,B4/-.02041793,.24902722,1.9496443,-2.6294547
     &                   ,8.56795911/
      IF(Y .LE. 2.7829681) THEN
        YB     = Y-1.
        FASINH = SQRT(6.*YB)*(((((A5*YB+A4)*YB+A3)*YB+A2)*YB-.15)*YB+1.)
      ELSE
        V      = LOG(Y)
        W      = 1./Y-0.028527431
        FASINH = V+LOG(2.*V)*(1.+1./V)+(((B4*W+B3)*W+B2)*W+B1)*W+B0
      ENDIF
      RETURN
      END
