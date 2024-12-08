C======= <<< TMESH : MESH GENERATION BY TRANSFINITE INTERPOLATION >>====
C=======================================================================
      PARAMETER ( JDIM=501, KDIM=200 )
      DIMENSION X(JDIM,KDIM),Y(JDIM,KDIM)
     &         ,XTR(JDIM,KDIM),YTR(JDIM,KDIM)
      COMMON/BASE/JMAX,KMAX,JWAKE,ROUT,XDW,DS1,DS2,MTOUT

C                                           << INPUT & SURFACE MESH >>
      CALL TINPUT(X,Y)
C                                           << VOLUME MESH >>
      CALL TCORE(JDIM,KDIM, JMAX,KMAX, X,Y, XTR,YTR)
C                                           << OUTPUT >>
      CALL MESH_IO(2, JDIM,KDIM, JMAX,KMAX, X,Y, MTOUT)

      STOP
      END
C=======================================================================
      SUBROUTINE TINPUT(X,Y)
C=======================================================================
      PARAMETER ( JDIM=501, KDIM=200, LD=501 )
      DIMENSION X(JDIM,KDIM),Y(JDIM,KDIM)
      DIMENSION X1(LD),Y1(LD),X2(LD),Y2(LD)
      COMMON/BASE/JMAX,KMAX,JWAKE,ROUT,XDW,DS1,DS2,MTOUT
CIN   NAMELIST/INPUTG/JMAX,KMAX,JWAKE,ROUT,XDW,DS1,DS2,MTOUT
C............................................................. INPUT
CDEFALULT
      JMAX  = 201
      KMAX  = 50
      JWAKE = 32
      DS1   = 0.03
      DS2   = 0.8
      ROUT  = 5.
      XDW   = 5.
      IAF   = 1
C
      MTOUT = 11
CIN   READ (5,INPUTG)
          IF(MOD(JMAX,2).EQ.0) THEN
            JMAX   = JMAX + 1
            WRITE(6,*)'JMAX MUST BE ODD NO.,   MODIFIED TO ',JMAX
          ENDIF
          JLEAD  = (JMAX+1)/2
          JTAIL1 = JWAKE+1
          JTAIL2 = JMAX+1-JTAIL1
          NWALL  = JMAX-(JTAIL1-1)*2
      WRITE(6,*)'<<< INPUT DATA FOR AIRFOIL-MESH >>>...........'
      WRITE(6,*)'  NUMBER OF MESH POINTS : JMAX =',JMAX,' , KMAX =',KMAX
      IF(JWAKE.GE.1) THEN
        WRITE(6,*)' GENERATE C-MESH'
        WRITE(6,*)' NUMBER OF POINTS ON THE WAKE-CUT: JWAKE =',JWAKE
      ELSE
        WRITE(6,*)' GENERATE O-MESH'
      ENDIF
      WRITE(6,*)' LEADING EDGE MESH NO.  : JLEAD        =',JLEAD
      WRITE(6,*)' TRAILING EDGE MESH NO. : JTAIL1,JTAIL2=',JTAIL1,JTAIL2
      WRITE(6,*)' NO. OF MESH ON THE WALL: NWALL        =',NWALL
      WRITE(6,*)' OUTER RADIUS           : ROUT         =',ROUT
      IF(JWAKE.GE.1)WRITE(6,*)' DOWNSTREAM LOCATION: XDW=',XDW
      WRITE(6,*)' SPACING AT WALL & OUTER BOUNDARY: DS1, DS2 =',DS1,DS2
      WRITE(6,*)' OUTPUT FILE NO.      : MTOUT =',MTOUT
C
C................................ << MESH ON INNER BOUNDARY : AIRFOIL >>

      IF (IAF.EQ.1)THEN
CDS     ** arbitrary airfoil reader
        call airfoil(JDIM,KDIM,X,Y)
      ELSE
CDS     ** NACA 4-digit airfoil generation
        call naca4(JDIM,KDIM,X,Y,IAF,0,0,12)
      ENDIF

C................................ << MESH ON INNER BOUNDARY : WAKE >>
      IF(JWAKE.GE.1) THEN
        DSX    = SQRT( (X(JTAIL1,1)-X(JTAIL1+1,1))**2 +
     &                 (Y(JTAIL1,1)-Y(JTAIL1+1,1))**2  )
        X1(2)  = XDW
        Y1(2)  = 0.
        X1(1)  = 1.
        Y1(1)  = 0.
        CALL GRCV2D(2,X1,Y1,JTAIL1,0, DSX, 0., X2,Y2)
        DO 50 J =1,JTAIL1-1
          JJ      = JTAIL1+1-J
          X(J       ,1) = X2(JJ)
          X(JMAX+1-J,1) = X(J,1)
          Y(J       ,1) = 0.
  50      Y(JMAX+1-J,1) = Y(J,1)
      ENDIF
C......................................... << MESH ON OUTER-BOUNDARY >>
      DATA PI/3.1415927/, NWK2/50/
      IF(JWAKE.GE.1) THEN
        X1(1)  =  XDW
        Y1(1)  = -ROUT
        X1(2)  =  1.
        Y1(2)  = -ROUT
        DO 60 J = 3,NWK2-2
          TH      = 1.5*PI-PI/FLOAT(NWK2-2-1)*FLOAT(J-2)
          X1(J)  = 1.+ROUT*COS(TH)
  60      Y1(J)  =    ROUT*SIN(TH)
            X1(NWK2-1) =  1.
            Y1(NWK2-1) =  ROUT
            X1(NWK2)   =  XDW
            Y1(NWK2)   =  ROUT
      ELSE
        DO 62 J = 1,NWK2
          TH     = -2.*PI/FLOAT(NWK2-1)*FLOAT(J-1)
          X1(J)  = ROUT*COS(TH) + 0.5
  62      Y1(J)  = ROUT*SIN(TH)
      ENDIF
      CALL GRCV2D(NWK2,X1,Y1,JMAX,1, 1., 1., X2,Y2)
      DO 70   J = 1,JMAX
        X(J,KMAX) = X2(J)
  70    Y(J,KMAX) = Y2(J)
C.................................................. 
      DO 80 J=1,JMAX,JMAX-1
      X1(1)  = X(J,1)
      Y1(1)  = Y(J,1)
      X1(2)  = X(J,KMAX)
      Y1(2)  = Y(J,KMAX)
      CALL GRCV2D(2,X1,Y1,KMAX,0, DS1, DS2, X2,Y2)
      DO 80 K = 2,KMAX-1
      X(J,K)  = X2(K)
  80  Y(J,K)  = Y2(K)
      RETURN
      END
C=======================================================================
      FUNCTION YNACA(T,XX,IAF)
C=======================================================================
C T     (INPUT)  : THICKNESS OF NACA00?? (=.12 FOR NACA0012)
C XX    (INPUT)  : X COORDINATE   0.<=XX<=1.
C IAF    (INPUT) : IF IAF=12 THEN NACA12 for NASA TURBULENCE CASE
C YNACA (OUTPUT) : OUTPUT OF Y COORDINATE OF NACA00?? AT X=XX
      IF (IAF.EQ.12) THEN
CDS nakahashi 
        X     = XX*1.00893
CDS opened TE (finite thickness at TE)
        a4    = 0.1015 
      ELSE
        X     = XX
CDS closed TE
        a4    = 0.1036
      ENDIF

CDS
CDS      YNACA=T/.2*(.2969*SQRT(X)-.126*X-.3516*X**2+.2843*X**3-.1015*X**4)
CDS
      YNACA=T/.2*(.2969*SQRT(X)-.126*X-.3516*X**2+.2843*X**3-a4*X**4)
CDS
CDS   scaled down if NACA0012 case
      if (IAF.EQ.12)YNACA=YNACA/1.00893
      
CDS      write(*,*) X, XX, YNACA
      RETURN
      END

C=======================================================================
      SUBROUTINE CNACA(CM,CP,X,YT,XN,YN)
C=======================================================================
C CM    (INPUT)  : MAX CAMBER OF NACA?000 (=.02 FOR NACA2412)
C CP    (INPUT)  : MAX CAMBER POSITION OF NACA0?00 (=.4 FOR NACA2412)
C X     (INPUT)  : X COORDINATE   0.<=XX<=1.
C YT    (INPUT)  : THICKNESS (Y)
C XN (OUTPUT) : OUTPUT OF X COORDINATE OF NACA???? 
C YN (OUTPUT) : OUTPUT OF Y COORDINATE OF NACA????
CDS      X     = XX
      IF (X.LT.CP) THEN
        YC = CM/CP**2 * (2*CP*X-X**2)
        DYC = 2*CM/CP**2 * (CP-X)
      ELSE
        YC = CM/(1.-CP)**2 * (1.-2*CP+2*CP*X-X**2)
        DYC = 2*CM/(1.-CP)**2 * (CP-X)
      ENDIF

      theta = atan(DYC)

CDS YT>0 for upper surface
CDS YT<0 for lower surface
      XN = X  - YT*sin(theta)
      YN = YC + YT*cos(theta)

CDS      CNACA = YT + C
CDS      write(*,*) X,YT,CNACA

      RETURN
      END

C==<< TRANSFINITE INTERPOLATION USING ARC-LENGTH BLENDING FUNCTION >>==
      SUBROUTINE TCORE(JD,KD,JMAX,KMAX,X,Y,X1,Y1)
C=======================================================================
      DIMENSION X(JD,KD),Y(JD,KD),X1(JD,KD),Y1(JD,KD), S(300)
      S(1)   = 0.
      DO 10 K=2,KMAX
  10  S(K)   =S(K-1)+SQRT((X(1,K)-X(1,K-1))**2+(Y(1,K)-Y(1,K-1))**2)
      DO 20 J=1,JMAX
      DO 20 K=1,KMAX
      BETA   = 1.-S(K)/S(KMAX)
      X1(J,K)= BETA*X(J,1) + (1.-BETA)*X(J,KMAX)
  20  Y1(J,K)= BETA*Y(J,1) + (1.-BETA)*Y(J,KMAX)
      S(1)   = 0.
      DO 30 J=2,JMAX
  30  S(J)   =S(J-1)+SQRT((X(J,1)-X(J-1,1))**2+(Y(J,1)-Y(J-1,1))**2)
      DO 40 J=1,JMAX
      DO 40 K=1,KMAX
      ALPHA  = 1.-S(J)/S(JMAX)
      X(J,K) = X1(J,K) + ALPHA*(X(1,K)-X1(1,K)) +
     &                   (1.-ALPHA)*(X(JMAX,K)-X1(JMAX,K))
  40  Y(J,K) = Y1(J,K) + ALPHA*(Y(1,K)-Y1(1,K)) +
     &                   (1.-ALPHA)*(Y(JMAX,K)-Y1(JMAX,K))
      RETURN
      END

C==<< NACA 4digit airfoil generator >>==
      SUBROUTINE NACA4(JDIM,KDIM,X,Y,IAF,NM,NP,NT)
C=======================================
C IF IAF = 12 THEN NACA0012 for Comparison (NASA Case)
C ELSE General NACA 4-digit airfoil
C  
C NACAMPTT (M=NM, P=NP, TT=NT)
C CMNACA: maximum camber, NM/100
C CPNACA: maximum camber location, NP/10
C TNACA: thickness, NT/100 
C
C EXAMPLE: NACA2412 
C CMNACA=0.02
C CPNACA=0.4
C TNACA =0.12
C
      PARAMETER ( LD=501 )
      DIMENSION X(JDIM,KDIM),Y(JDIM,KDIM)
      DIMENSION X1(LD),Y1(LD),X2(LD),Y2(LD)
      DIMENSION XL(LD),YL(LD),XU(LD),YU(LD)
      COMMON/BASE/JMAX,KMAX,JWAKE,ROUT,XDW,DS1,DS2,MTOUT

      JLEAD  = (JMAX+1)/2
      JTAIL1 = JWAKE+1
      JTAIL2 = JMAX+1-JTAIL1
      NWALL  = JMAX-(JTAIL1-1)*2

C       < NACA0012 >
C        CMNACA = 0.0
C        CPNACA = 0.0
C        TNACA  = 0.12

C       < NACA2412 >
C        CMNACA = 0.02
C        CPNACA = 0.4
C        TNACA  = 0.12

C       < NACA2315 >
C        CMNACA = 0.02
C        CPNACA = 0.3
C        TNACA  = 0.15

CDS
CDS     < NACA0012 for validation >
      IF (IAF.EQ.12)THEN
        WRITE(*,*)'NACA0012' 
        WRITE(*,*)' -- turbulence modeling resource case'
        NM = 0
        NP = 0
        NT = 12
C        CMNACA = 0.0
C        CPNACA = 0.0
C        TNACA  = 0.12
      ENDIF

      CMNACA = NM/100.
      CPNACA = NP/10.
      TNACA  = NT/100.

        JLEADW = 301
C                            AIRFOIL LOWER SURFACE
        DO I=1,JLEADW
          X1(I)  = 1.-FLOAT(I-1)/FLOAT(JLEADW-1)
          YT     = -YNACA(TNACA,X1(I),IAF)
          CALl CNACA(CMNACA,CPNACA,X1(I),YT,XL(I),YL(I))
C          write(12,*)XL(I),YL(I)
        ENDDO
        NWALLH = (NWALL+1)/2
        CALL GRCV2D(JLEADW,XL,YL,NWALLH,1, 0.5, 0.5, X2,Y2)
        DO 20 J=JTAIL1,JLEAD
        JJ     = J-JTAIL1+1
        X(J,1) = X2(JJ)
 20     Y(J,1) = Y2(JJ)
C                            AIRFOIL UPPER SURFACE
        DO I=1,JLEADW
          X2(I)  = X1(JLEADW-I+1)
          YT  = YNACA(TNACA,X2(I),IAF) 
          CALl CNACA(CMNACA,CPNACA,X2(I),YT,XU(I),YU(I))
C          write(12,*)XU(I),YU(I)
        ENDDO
        CALL GRCV2D(JLEADW,XU,YU,NWALLH,1, 0.5, 0.5, X1,Y1)
        DO 40 J=JLEAD+1,JTAIL2
        JJ     = J-JLEAD+1
        X(J,1) = X1(JJ)
 40     Y(J,1) = Y1(JJ)

      RETURN
      END

C==<< arbitrary airfoil reader >>==
      SUBROUTINE AIRFOIL(JDIM,KDIM,X,Y)
C=======================================
      PARAMETER ( LD=501 )
      DIMENSION X(JDIM,KDIM),Y(JDIM,KDIM)
      DIMENSION X1(LD),Y1(LD),X2(LD),Y2(LD)
      DIMENSION XN(LD),YN(LD)
      COMMON/BASE/JMAX,KMAX,JWAKE,ROUT,XDW,DS1,DS2,MTOUT

      JLEAD  = (JMAX+1)/2
      JTAIL1 = JWAKE+1
      JTAIL2 = JMAX+1-JTAIL1
      NWALL  = JMAX-(JTAIL1-1)*2
      NWALLH = (NWALL+1)/2

      mt=10

      open(mt,file='work/airfoil-snowed.dat')

CDS   JAIRF: ODD NUMBER (starting from TE)
CDS   LOWER TE TO LE
CDS   UPPER LE TO TE
CDS   1 LE point, but 2 TE points (can be same)
CDS
      read(mt,*)JAIRF
      write(*,*)'reading airfoil.dat, # points=',JAIRF

      JLEADW=(JAIRF+1)/2

      DO I=1,JLEADW
        read(mt,*)X1(I),Y1(I)
CDS        write(13,*)X1(I),Y1(I)
      enddo
      X2(1)=X1(JLEADW)
      Y2(1)=Y1(JLEADW)
      DO I=2,JLEADW
        read(mt,*)X2(I),Y2(I)
CDS        write(13,*)X2(I),Y2(I)
      enddo

      close(mt)

C                AIRFOIL LOWER SURFACE
      CALL GRCV2D(JLEADW,X1,Y1,NWALLH,1, 0.5, 0.5, XN,YN)
      DO J=JTAIL1,JLEAD
        JJ     = J-JTAIL1+1
        X(J,1) = XN(JJ)
        Y(J,1) = YN(JJ)
CDS       write(14,*)J,X(J,1),Y(J,1)
      ENDDO
C                AIRFOIL UPPER SURFACE
      CALL GRCV2D(JLEADW,X2,Y2,NWALLH,1, 0.5, 0.5, XN,YN)
      DO J=JLEAD+1,JTAIL2
        JJ     = J-JLEAD+1
        X(J,1) = XN(JJ)
        Y(J,1) = YN(JJ)
CDS       write(14,*)J,X(J,1),Y(J,1)
      ENDDO

      RETURN
      END
