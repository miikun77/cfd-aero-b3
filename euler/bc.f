C === Initial condition and setup of SHOCK-REFLECTION (i-problem=1)====
      SUBROUTINE setup_SL(JDIM,KDIM, Q)
C =====================================================================
      DIMENSION Q(JDIM,KDIM,4)
      COMMON/CNUMB/JMAX, KMAX
      COMMON/CFLOW/GAMMA, FSMACH, ALPHA

!    upper surface condition for a oblique shock generation.
!    Freestream Mach =2 is assumed.

      GAMI   = GAMMA-1.
      PIN    = 1.0/GAMMA
      RHOIN  = 1.0

      BETA   = 0.5688
      FSBETA = (FSMACH*SIN(BETA))**2
      P2     = PIN * (1. + 2.*GAMMA/(GAMMA+1.)*(FSBETA-1.) )
      PFAC   = P2/PIN
      THETA  = ATAN( 2./TAN(BETA)*( FSBETA-1.)/
     &                            (FSMACH**2*(GAMMA+COS(2.*BETA))+2.) )
      RHOUP  = RHOIN * (GAMMA+1.)*FSBETA/(GAMI*FSBETA +2.)
      FSM2   = SQRT( (1.+GAMI*.5*FSBETA)/(GAMMA*FSBETA-GAMI*.5) )
     &        /SIN(BETA-THETA)
      UUP    = FSM2*COS(THETA)
      VUP    =-FSM2*SIN(THETA)
      EUP    = 0.5*RHOUP*FSM2**2 + P2/GAMI

      DO 150 J=1,JMAX
        Q(J,KMAX,1) = RHOUP
        Q(J,KMAX,2) = UUP*RHOUP
        Q(J,KMAX,3) = VUP*RHOUP
        Q(J,KMAX,4) = EUP
 150  CONTINUE

      return
      end
C === BOUNDARY CONDITION FOR SHOCK-REFLECTION (i-problem=1) ===========
      SUBROUTINE BC_shock(JDIM,KDIM, Q,XYJ,XY, PP)
C =====================================================================
      DIMENSION Q(JDIM,KDIM,4) , XYJ(JDIM,KDIM), XY(JDIM,KDIM,4)
     &         ,PP(JDIM,KDIM)
      COMMON/CNUMB/JMAX, KMAX
      COMMON/CFLOW/GAMMA, FSMACH, ALPHA
C
      GAMI   = GAMMA-1.

!----- slip surface conditions for Euler computation --------------
!  Pressure and velocity are extraporated and density is determined by total enthalpy.
!  More sophisticated treatment is to use normal momentum (see T. Pulliam's papers)

      HSTFS  = 1./GAMI + 0.5*FSMACH**2	! =  THE STAGNATION FREESTREAM ENTHALPY

      DO 10 J=2,JMAX-1
C
        DRXY  = 1./(Q(J,2,1)*SQRT( XY(J,2,3)**2 + XY(J,2,4)**2 ))
        VN2   = ( XY(J,2,3)*Q(J,2,2) + XY(J,2,4)*Q(J,2,3) )*DRXY
        VT2   = ( XY(J,2,4)*Q(J,2,2) - XY(J,2,3)*Q(J,2,3) )*DRXY
        P2    = GAMI*( Q(J,2,4)-.5*(Q(J,2,2)**2+Q(J,2,3)**2)/Q(J,2,1) )
     &         *XYJ(J,2)
C
        VT    = VT2
        VN    = 0.
        PW    = P2
        DXY   = 1./SQRT( XY(J,1,3)**2 + XY(J,1,4)**2 )
        U1    = ( XY(J,1,4)*VT + XY(J,1,3)*VN )*DXY
        V1    = (-XY(J,1,3)*VT + XY(J,1,4)*VN )*DXY
C
        RJ       = 1./XYJ(J,1)
        Q(J,1,1) = GAMMA*PW/(GAMI* (HSTFS -0.5*(U1**2+V1**2) ) )
        Q(J,1,1) = Q(J,1,1)*RJ
        Q(J,1,2) = U1*Q(J,1,1)
        Q(J,1,3) = V1*Q(J,1,1)
        Q(J,1,4) = PW*RJ/GAMI + .5*(Q(J,1,2)**2+Q(J,1,3)**2)/Q(J,1,1)
C
 10   CONTINUE
C
C
C...<< BOUNDARY CONDITION ON EXIT PLANE >>...
C
      DO 20 K=1,KMAX
      TEMP        = XYJ(JMAX-1,K)/XYJ(JMAX,K)
      Q(JMAX,K,1) = Q(JMAX-1,K,1)*TEMP
      Q(JMAX,K,2) = Q(JMAX-1,K,2)*TEMP
      Q(JMAX,K,3) = Q(JMAX-1,K,3)*TEMP
      Q(JMAX,K,4) = Q(JMAX-1,K,4)*TEMP
 20   CONTINUE
C
C   FOR SUBSONIC FLOW, PRESSURE IS FIXED AT EXIT PLANE.
C
      IF(FSMACH.LE.1.) THEN
        DO 22 K=1,KMAX
         PEXIT       = 1./GAMMA
         Q(JMAX,K,4) = PEXIT/((GAMMA-1.)*XYJ(JMAX,K)) +
     &                .5*(Q(JMAX,K,2)**2+Q(JMAX,K,3)**2)/Q(JMAX,K,1)
 22     CONTINUE
      ENDIF
C
C
C...<< ON UPPER BOUNDARY >>...


C    --> FIXED TO INITIAL CONDITION
C
C
C...<< INFLOW BOUNDARY CONDITION >>...
C
C    --> FIXED TO INITIAL CONDITION
C
C
      RETURN
      END
C === BOUNDARY CONDITION FOR GAM-PROBLEM (i-problem=2) ================
      SUBROUTINE BC_GAM(JDIM,KDIM, Q,XYJ,XY,PP)
C =====================================================================
      DIMENSION Q(JDIM,KDIM,4) , XYJ(JDIM,KDIM), XY(JDIM,KDIM,4)
     &         ,PP(JDIM,KDIM)
      COMMON/CNUMB/JMAX, KMAX
      COMMON/CFLOW/GAMMA, FSMACH, ALPHA

      GAMI   = GAMMA-1.

!----- slip surface conditions for Euler computation --------------
!  Pressure and velocity are extraporated and density is determined by total enthalpy.
!  More sophisticated treatment is to use normal momentum (see T. Pulliam's papers)

      HSTFS  = 1./GAMI + 0.5*FSMACH**2	! =  THE STAGNATION FREESTREAM ENTHALPY

      DO 2 K=1,3
      DO 2 J=1,JMAX
      PP(J,K) = GAMI*( Q(J,K,4)-.5*(Q(J,K,2)**2+Q(J,K,3)**2)/Q(J,K,1) )
     &         *XYJ(J,K)
  2   CONTINUE

!   OBTAIN U AND V ON BODY VIA TANGENCY AND EXTRAPOLATION
!   p_wall is extraporated from inner pressure 
!                   by 0th (c_wall=0) or 1st order (c_wall=1).
      c_wall = 1.
      c1 = 1. +c_wall
      c2 = 0. -c_wall

      k   = 1
      kp1 = k+1
      kp2 = k+2

      DO 10 J=2,JMAX-1
        UEXT = 2.*( XY(J,2,1)*Q(J,2,2)/Q(J,2,1)
     &             +XY(J,2,2)*Q(J,2,3)/Q(J,2,1) )
     &          - ( XY(J,3,1)*Q(J,3,2)/Q(J,3,1)
     &             +XY(J,3,2) *Q(J,3,3)/Q(J,3,1) )
        V = 0.
        U1 = (  XY(J,1,4)*UEXT - XY(J,1,2)*V )/XYJ(J,1)
        V1 = ( -XY(J,1,3)*UEXT + XY(J,1,1)*V )/XYJ(J,1)

        pw   = pp(j,kp1)*c1 + pp(j,kp2)*c2
        Q(J,K,1) = GAMMA*PW/(GAMI* (HSTFS - 0.5*(U1**2 + V1**2) ) )
        Q(J,K,1) = Q(J,K,1)/XYJ(J,K)
        Q(J,K,2) = U1*Q(J,K,1)
        Q(J,K,3) = V1*Q(J,K,1)
        Q(J,K,4) = PW/(GAMI*XYJ(J,K))
     &                 +.5 * ( Q(J,K,2)**2+Q(J,K,3)**2) /Q(J,K,1)
   10 CONTINUE

C...<< ON UPPER BOUNDARY >>...   SYMMETRY CONDITION
! Symmetry condition should be given by using the dummy mesh,
! but here a zero-th extraporation for slip condition is used for simplicity.

      DO 28 K=KMAX-2,KMAX
      DO 28 J=1,JMAX
      PP(J,K) = GAMI*( Q(J,K,4)-.5*(Q(J,K,2)**2+Q(J,K,3)**2)/Q(J,K,1) )
     &         *XYJ(J,K)
  28  CONTINUE

      c_wall = 0.		! c_wall=1.0 for 1st-oder slip condition
!						! zero-th order is used here for robustness.
      c1 = 1. +c_wall
      c2 = 0. -c_wall
 
      K   = KMAX
      KM1 = KMAX-1
      KM2 = KMAX-2

      DO 30 J=2,JMAX-1
        UEXT = 2.*( XY(J,KM1,1)*Q(J,KM1,2)/Q(J,KM1,1)
     &             +XY(J,KM1,2)*Q(J,KM1,3)/Q(J,KM1,1) )
     &          - ( XY(J,KM2,1)*Q(J,KM2,2)/Q(J,KM2,1)
     &             +XY(J,KM2,2)*Q(J,KM2,3)/Q(J,KM2,1) )
        V = 0.
        U1 = (  XY(J,K,4)*UEXT - XY(J,K,2)*V )/XYJ(J,K)
        V1 = ( -XY(J,K,3)*UEXT + XY(J,K,1)*V )/XYJ(J,K)

        pw   = pp(j,km1)*c1 + pp(j,km2)*c2
        Q(J,K,1) = GAMMA*PW/(GAMI* (HSTFS - 0.5*(U1**2 + V1**2) ) )
        Q(J,K,1) = Q(J,K,1)/XYJ(J,K)
        Q(J,K,2) = U1*Q(J,K,1)
        Q(J,K,3) = V1*Q(J,K,1)
        Q(J,K,4) = PW/(GAMI*XYJ(J,K))
     &             +.5 * ( Q(J,K,2)**2+Q(J,K,3)**2) /Q(J,K,1)
 30   CONTINUE

C...<< BOUNDARY CONDITION ON EXIT PLANE >>...

      DO 40 K=1,KMAX
      TEMP        = XYJ(JMAX-1,K)/XYJ(JMAX,K)
      Q(JMAX,K,1) = Q(JMAX-1,K,1)*TEMP
      Q(JMAX,K,2) = Q(JMAX-1,K,2)*TEMP
      Q(JMAX,K,3) = Q(JMAX-1,K,3)*TEMP
      Q(JMAX,K,4) = Q(JMAX-1,K,4)*TEMP
 40   CONTINUE

C   FOR SUBSONIC FLOW, PRESSURE IS FIXED AT EXIT PLANE.

      IF(FSMACH.LE.1.) THEN
        DO 42 K=1,KMAX
         PEXIT       = 1./GAMMA
         Q(JMAX,K,4) = PEXIT/((GAMMA-1.)*XYJ(JMAX,K)) +
     &                .5*(Q(JMAX,K,2)**2+Q(JMAX,K,3)**2)/Q(JMAX,K,1)
 42     CONTINUE
      ENDIF

C...<< INFLOW BOUNDARY CONDITION >>...

      TINF   = 1.
      TZ1    = TINF*(1.+0.5*GAMI*FSMACH**2)

      DO 50 K=1,KMAX
      P2     = GAMI*( Q(2,K,4)-.5*(Q(2,K,2)**2+Q(2,K,3)**2)/Q(2,K,1) )
     &         *XYJ(2,K)
      P3     = GAMI*( Q(3,K,4)-.5*(Q(3,K,2)**2+Q(3,K,3)**2)/Q(3,K,1) )
     &         *XYJ(3,K)
CC    P1     = 2.*P2 - P3
      P1     = P2
      RHO1   = (GAMMA*P1)**(1./GAMMA)
      T1     = GAMMA*P1/RHO1
      U1     = SQRT( 2./GAMI*(TZ1-T1) )
      V1     = 0.
      E1     = P1/GAMI + 0.5*RHO1*U1**2
      Q(1,K,1) = RHO1/XYJ(1,K)
      Q(1,K,2) = RHO1*U1/XYJ(1,K)
      Q(1,K,3) = 0.
      Q(1,K,4) = E1/XYJ(1,K)
  50  CONTINUE

      RETURN
      END
C === BOUNDARY CONDITIONS FOR AIRFOIL TOPOLOGY (C-MESH TYPE only) =====
!                                                     (i-problem = 3) =
      SUBROUTINE BC_airfoil(JDIM,KDIM, Q,XYJ,XY,PP)
C =====================================================================
      DIMENSION Q(JDIM,KDIM,4) , XYJ(JDIM,KDIM), XY(JDIM,KDIM,4)
     &         ,PP(JDIM,KDIM)
      COMMON/CNUMB/JMAX, KMAX
      COMMON/CFLOW/GAMMA, FSMACH, ALPHA
      common/cairfoil/jtail

      GAMI   = GAMMA-1.

!  xi coordinates OF START AND END OF BODY
      jtail1 = jtail
      jtail2 = jmax-jtail1+1

!----- slip surface conditions for Euler computation --------------
!  Pressure and velocity are extraporated and density is determined by total enthalpy.
!  More sophisticated treatment is to use normal momentum (see T. Pulliam's papers)

      HSTFS  = 1./GAMI + 0.5*FSMACH**2	! =  THE STAGNATION FREESTREAM ENTHALPY

      DO 2 K=1,3
      DO 2 J=1,JMAX
      PP(J,K) = GAMI*( Q(J,K,4)-.5*(Q(J,K,2)**2+Q(J,K,3)**2)/Q(J,K,1) )
     &         *XYJ(J,K)
  2   CONTINUE

!   OBTAIN U AND V ON BODY VIA TANGENCY AND EXTRAPOLATION
!   p_wall is extraporated from inner pressure 
!                   by 0th (c_wall=0) or 1st order (c_wall=1).
      c_wall = 1.
      c1 = 1. +c_wall
      c2 = 0. -c_wall

      K = 1
      kp1 = k+1
      kp2 = k+2

      DO 10 J=jtail1, jtail2
        UEXT = 2.*( XY(J,2,1)*Q(J,2,2)/Q(J,2,1)
     &             +XY(J,2,2)*Q(J,2,3)/Q(J,2,1) )
     &          - ( XY(J,3,1)*Q(J,3,2)/Q(J,3,1)
     &             +XY(J,3,2)*Q(J,3,3)/Q(J,3,1) )

        V = 0.
        U1 = (  XY(J,1,4)*(UEXT) - XY(J,1,2)*V )/XYJ(J,1)
        V1 = ( -XY(J,1,3)*(UEXT) + XY(J,1,1)*V )/XYJ(J,1)

        pw   = pp(j,kp1)*c1 + pp(j,kp2)*c2
        Q(J,K,1) = GAMMA*PW/(GAMI* (HSTFS - 0.5*(U1**2 + V1**2) ) )
        Q(J,K,1) = Q(J,K,1)/XYJ(J,K)
        Q(J,K,2) = U1*Q(J,K,1)
        Q(J,K,3) = V1*Q(J,K,1)
        Q(J,K,4) = PW/(GAMI*XYJ(J,K))
     &                 +.5 * ( Q(J,K,2)**2+Q(J,K,3)**2) /Q(J,K,1)
   10 CONTINUE

!-------------  WAKE -----------------

         DO 20 M=1,4
         DO 20 J=1,jtail1
            I = JMAX +1 - J
            QSAV = 0.5*( Q(J,2,M)*XYJ(J,2) + Q(I,2,M)*XYJ(I,2) )
            Q(J,1,M) = QSAV/XYJ(J,1)
            Q(I,1,M) = QSAV/XYJ(I,1)
   20    CONTINUE


C... << FAR FIELD BC >>......................
!  Far field values are fixed for simplicity.
!  More sophisticated one is to use Riemann invariant (see T. Pulliam's paper)

C    --> FIXED TO INITIAL CONDITION


C...<< BOUNDARY CONDITION ON EXIT PLANE >>...

      DO 40 K=1,KMAX
      TEMP        = XYJ(JMAX-1,K)/XYJ(JMAX,K)
      Q(JMAX,K,1) = Q(JMAX-1,K,1)*TEMP
      Q(JMAX,K,2) = Q(JMAX-1,K,2)*TEMP
      Q(JMAX,K,3) = Q(JMAX-1,K,3)*TEMP
      Q(JMAX,K,4) = Q(JMAX-1,K,4)*TEMP
      
      TEMP        = XYJ(2,K)/XYJ(1,K)
      Q(1,K,1) = Q(2,K,1)*TEMP
      Q(1,K,2) = Q(2,K,2)*TEMP
      Q(1,K,3) = Q(2,K,3)*TEMP
      Q(1,K,4) = Q(2,K,4)*TEMP
 40   CONTINUE

C   FOR SUBSONIC FLOW, PRESSURE IS FIXED AT EXIT PLANE.

      IF(FSMACH.LE.1.) THEN
        DO 42 K=1,KMAX
         PEXIT       = 1./GAMMA
         Q(JMAX,K,4) = PEXIT/((GAMMA-1.)*XYJ(JMAX,K)) +
     &                .5*(Q(JMAX,K,2)**2+Q(JMAX,K,3)**2)/Q(JMAX,K,1)
         Q(1   ,K,4) = PEXIT/((GAMMA-1.)*XYJ(1,K)) +
     &                .5*(Q(1,K,2)**2+Q(1,K,3)**2)/Q(1,K,1)
 42     CONTINUE
      ENDIF

      RETURN
      END
