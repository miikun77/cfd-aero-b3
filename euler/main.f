C =====================================================================
C --- MACCORMACK METHOD FOR THE 2-D EULER EQUATIONS
C               GENERALIZED COORDINATES

C  SEMINAR VERTION (91/3/5)              ??? SOME BUGS ???
!  Modified (2010/5/24)
C =====================================================================
      PARAMETER (JDIM=501,KDIM=101)

      DIMENSION X(JDIM,KDIM), Y(JDIM,KDIM), XYJ(JDIM,KDIM)
     &         ,Q(JDIM,KDIM,4), QP(JDIM,KDIM,4) ,XY(JDIM,KDIM,4)
     &         ,E(JDIM,KDIM,4), F(JDIM,KDIM,4), PP(JDIM,KDIM)
      DIMENSION RES(10000)
      COMMON/CNUMB/JMAX, KMAX
      COMMON/CFLOW/GAMMA, FSMACH, ALPHA
      COMMON/CTIME/DT,CDT,n_steps,ITIME,ITIMEL,TIME
      COMMON/CIOA/mesh_file, q_file
      character(32) :: input_filename, mesh_file, q_file

C..............................................

      WRITE(*,*) '------  MACCORMACK TWO-STEP EXPLICIT METHOD ------'
      WRITE(*,*) '------       FOR 2-D, EULER EQUATIONS       ------'
      WRITE(*,*) '------      IN GENERALIZED COORDINATES      ------'
      WRITE(*,*) '--------------------------------------------------'
      WRITE(*,*) ' '

C.......<< Initialization and Setup >> .....................

      CALL setup(JDIM,KDIM,X,Y,XYJ,Q,QP,XY)


C.......<< TIME-STEP LOOP >>......................

  100 CONTINUE
      IF(ITIMEL.GE.n_steps) GO TO 120

      ITIME   = ITIME  + 1
      ITIMEL  = ITIMEL + 1

C...<< COMPUTE TIME STEP AT EVERY 10 STEPS >>...

      IF(MOD(ITIMEL,100).EQ.1) CALL DTCAL(JDIM,KDIM,XYJ,Q,XY)

      TIME    = TIME   + DT

C...<< TIME INTEGRATION >>...

      CALL STEP(JDIM,KDIM, XYJ,XY, Q,E,F,QP,PP,NRES,RES)

      GO TO 100

  120 CONTINUE        ! << END TIME-STEP LOOP >>....................

C..<< OUTPUT >>

      DO 200 I=1,4
      DO 200 K=1,KMAX
      DO 200 J=1,JMAX
      Q(J,K,I) = Q(J,K,I)*XYJ(J,K)
 200  CONTINUE

      CALL q_io(2,q_file,JDIM,KDIM,Q)
      CALL surf_io(JDIM,KDIM,X,Y,Q)
      CALL res_io(NRES,RES)

      STOP
      END
C =====================================================================
      SUBROUTINE STEP(JDIM,KDIM, XYJ,XY, Q,E,F,QP,PP,NRES,RES)
C =====================================================================
      DIMENSION XYJ(JDIM,KDIM), XY(JDIM,KDIM,4)
     &         ,Q(JDIM,KDIM,4), QP(JDIM,KDIM,4)
     &         ,E(JDIM,KDIM,4), F(JDIM,KDIM,4), PP(JDIM,KDIM)
      DIMENSION RES(10000)
      COMMON/CNUMB/JMAX, KMAX
      COMMON/CTIME/DT,CDT,n_steps,ITIME,ITIMEL,TIME
      common/cprobem/i_problem

C...<< PREDICTOR STEP >>...

      CALL FLUXD(JDIM,KDIM, XYJ,XY, Q,E,F,-1.)

      DO 10 IQ=1,4
      DO 10 K=2,KMAX-1
      DO 10 J=2,JMAX-1
      QP(J,K,IQ) = Q(J,K,IQ) - DT*( E(J  ,K,IQ)-E(J-1,K  ,IQ)
     &                             +F(J  ,K,IQ)-F(J  ,K-1,IQ) )
 10   CONTINUE

C...<< BOUNDARY CONDITION >>

      if(i_problem.eq.1) call bc_shock(JDIM,KDIM, QP,XYJ,XY, PP)
      if(i_problem.eq.2) call bc_gam(JDIM,KDIM, QP,XYJ,XY, PP)
      if(i_problem.eq.3) call bc_airfoil(JDIM,KDIM, QP,XYJ,XY, PP)

C...<< CORRECTOR STEP >>...

      CALL FLUXD(JDIM,KDIM, XYJ,XY, QP,E,F,1.)

      DO 50 IQ=1,4
      DO 50 K =2,KMAX-1
      DO 50 J =2,JMAX-1
      QP(J,K,IQ) = .5*( Q(J,K,IQ)+QP(J,K,IQ)
     &                 -DT*( E(J+1,K  ,IQ)-E(J  ,K,IQ)
     &                      +F(J  ,K+1,IQ)-F(J  ,K,IQ) ) )
 50   CONTINUE


C...<< BOUNDARY CONDITION >>

      if(i_problem.eq.1) call bc_shock(JDIM,KDIM, QP,XYJ,XY, PP)
      if(i_problem.eq.2) call bc_gam(JDIM,KDIM, QP,XYJ,XY, PP)
      if(i_problem.eq.3) call bc_airfoil(JDIM,KDIM, QP,XYJ,XY, PP)

C...<< COMPUTE L2 NORM OF S  , THE RESIDUAL >>.................

      if(mod(itime,10).eq.0) then
        RESI = 0.
        DO 300 K=2,KMAX-1 
        DO 300 J=2,JMAX-1 
 300    RESI    = RESI + (QP(J,K,1)-Q(J,K,1))**2
        RESI = RESI/((JMAX-2)*(KMAX-2))
        RESI = SQRT( RESI ) / DT
CDS RESIDUAL
        NRES = ITIME/10
        RES(NRES) = RESI
        WRITE(6,310) ITIME, RESI
 310    FORMAT(1H ,'ITIME = ',I5,':',3X,'RESIDUAL=',1PE9.3)
      endif

!.......................
      DO 400 IQ=1,4
      DO 400 K =1,KMAX
      DO 400 J =1,JMAX
 400  Q(J,K,IQ) = QP(J,K,IQ)

      RETURN
      END
C ======================================================================
      SUBROUTINE FLUXD(JDIM,KDIM, XYJ,XY, Q,E,F, SIGN)
C =====================================================================
      DIMENSION XYJ(JDIM,KDIM), XY(JDIM,KDIM,4)
     &         ,Q(JDIM,KDIM,4), E(JDIM,KDIM,4), F(JDIM,KDIM,4)
      COMMON/CNUMB/JMAX, KMAX
      COMMON/CFLOW/GAMMA, FSMACH, ALPHA
      COMMON/CDISS/EPS4

      DO 10 K=1,KMAX
      DO 10 J=1,JMAX

         RR  = 1./Q(J,K,1)
         U   = Q(J,K,2)*RR
         V   = Q(J,K,3)*RR
         QU  = XY(J,K,1)*U + XY(J,K,2)*V
         QV  = XY(J,K,3)*U + XY(J,K,4)*V

C       PJ IS PRESSURE / J
         PJ       = (GAMMA-1.)*(Q(J,K,4)-.5*Q(J,K,1)*(U*U+V*V))

         E(J,K,1) = Q(J,K,1)*QU
         E(J,K,2) = Q(J,K,2)*QU + XY(J,K,1)*PJ
         E(J,K,3) = Q(J,K,3)*QU + XY(J,K,2)*PJ
         E(J,K,4) = QU*( Q(J,K,4)+PJ )

         F(J,K,1) = Q(J,K,1)*QV
         F(J,K,2) = Q(J,K,2)*QV + XY(J,K,3)*PJ
         F(J,K,3) = Q(J,K,3)*QV + XY(J,K,4)*PJ
         F(J,K,4) = QV*( Q(J,K,4)+PJ )

      IF(EPS4.NE.0.)THEN

C  ...<< DISSIPATION >>

        J2  = J
        IF(J2.EQ.1   ) J2 = 2
        IF(J2.EQ.JMAX) J2 = JMAX-1
        JP  = J2+1
        JM  = J2-1
        DX1  = 0.125*EPS4*( Q(JP,K,1)-2.*Q(J2,K,1)+Q(JM,K,1) )
        DX2  = 0.125*EPS4*( Q(JP,K,2)-2.*Q(J2,K,2)+Q(JM,K,2) )
        DX3  = 0.125*EPS4*( Q(JP,K,3)-2.*Q(J2,K,3)+Q(JM,K,3) )
        DX4  = 0.125*EPS4*( Q(JP,K,4)-2.*Q(J2,K,4)+Q(JM,K,4) )
        E(J,K,1) = E(J,K,1)+DX1*SIGN
        E(J,K,2) = E(J,K,2)+DX2*SIGN
        E(J,K,3) = E(J,K,3)+DX3*SIGN
        E(J,K,4) = E(J,K,4)+DX4*SIGN
C
        K2  = K
        IF(K2.EQ.1   ) K2 = 2
        IF(K2.EQ.KMAX) K2 = KMAX-1
        KP  = K2+1
        KM  = K2-1
        DE1  = 0.125*EPS4*( Q(J,KP,1)-2.*Q(J,K2,1)+Q(J,KM,1) )
        DE2  = 0.125*EPS4*( Q(J,KP,2)-2.*Q(J,K2,2)+Q(J,KM,2) )
        DE3  = 0.125*EPS4*( Q(J,KP,3)-2.*Q(J,K2,3)+Q(J,KM,3) )
        DE4  = 0.125*EPS4*( Q(J,KP,4)-2.*Q(J,K2,4)+Q(J,KM,4) )
        F(J,K,1) = F(J,K,1)+DE1*SIGN
        F(J,K,2) = F(J,K,2)+DE2*SIGN
        F(J,K,3) = F(J,K,3)+DE3*SIGN
        F(J,K,4) = F(J,K,4)+DE4*SIGN

      ENDIF

 10   CONTINUE

      RETURN
      END
