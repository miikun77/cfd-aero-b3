C =====================================================================
      SUBROUTINE setup(JDIM,KDIM,X,Y,XYJ,Q,QP,XY)
C =====================================================================
      DIMENSION X(JDIM,KDIM), Y(JDIM,KDIM), XYJ(JDIM,KDIM)
     &         ,Q(JDIM,KDIM,4), QP(JDIM,KDIM,4), XY(JDIM,KDIM,4)
      COMMON/CNUMB/JMAX, KMAX
      COMMON/CFLOW/GAMMA, FSMACH, ALPHA
      COMMON/CTIME/DT,CDT,n_steps,ITIME,ITIMEL,TIME
      LOGICAL    restart
      COMMON/CIOA/mesh_file, q_file
      COMMON/CDISS/EPS4
      common/cprobem/i_problem
      common/cairfoil/jtail
      character(32) :: mesh_file, q_file

      namelist/input/ i_problem, fsmach, alpha, mesh_file, q_file,
     &                cdt, n_steps, restart, eps4, jtail

! setting of default values
      i_problem = 3			! 1=shock-reflection, 2=gam problem, 3=airfoil
      FSMACH    = 0.675		! FREESTREAM MACH NUMBER
      alpha     = 0.		! angle of attack in degree(only for problem 3)
      mesh_file = 'gam.xyz'	! file name of mesh data 
      q_file    = 'gam.q'	! file name of Q data (computed data)
      CDT       = 0.9		! Courant number (cdt < 1.0)
      n_steps   = 5000		! Number of steps of computation
      restart   = .false.	! logical data to define new or restart computation
      EPS4      = fsmach*2.	! Coefficient of 4th order dissipation
      jtail     = 14		! lower TE's J #  of airfoil c-mesh (only for problem 3)

! other valuables
!     DT         : TIME INCREMENT determined using cdt
!     JMAX, KMAX : NUumbers of mes points in xi- and eta-directions
      GAMMA     = 1.4			! RATIO OF SPECIFIC HEAT


      WRITE(*,*) ' '
      WRITE(*,*) 'The 2D Euler solver using MacCormack method'

      mtin    = 11
      open(mtin,file="mesh.d",form='formatted')

        READ(mtin,INPUT)
        write(*,INPUT)

      close(mtin)

      COSANG = COS(3.141593*alpha/180.)
      SINANG = SIN(3.141593*alpha/180.)

C ................................................... GRID I/O

      CALL mesh_io(1,mesh_file,JDIM,KDIM,X,Y,xyj,xy)

C.........................................FLOW DATA I/O

      TIME   = 0.    ! 
      ITIME  = 0     ! NUMBER OF TIME STEPS FROM THE BEGINNING.
      ITIMEL = 0     ! NUMBER OF TIME STEPS IN THE CURRENT COMPUTATION.


C....<< READ IN INITIAL CONDITIONS IF ANY >>....

      GAMI   = GAMMA-1.
      PIN    = 1.0/GAMMA
      RHOIN  = 1.0
      UIN    = FSMACH*COSANG
      VIN    = FSMACH*SINANG
      EIN    = 0.5*RHOIN*FSMACH**2 + PIN/GAMI

      DO 140  K=1,KMAX
      DO 140  J=1,JMAX
      Q(J,K,1) = RHOIN
      Q(J,K,2) = UIN*RHOIN
      Q(J,K,3) = VIN*RHOIN
      Q(J,K,4) = EIN
 140  CONTINUE

      if(i_problem.eq.1) call setup_SL(jdim,kdim,q)

      if(restart) then
           WRITE(6,*)'This is a restart computation.'
           write(*,*)'Computed result will be overwritten on this file.'

           CALL q_IO(1,q_file,JDIM,KDIM,Q)
      endif

      DO 160 IQ=1,4
      DO 160  K=1,KMAX
      DO 160  J=1,JMAX
      Q(J,K,IQ)  = Q(J,K,IQ)/XYJ(J,K)
      QP(J,K,IQ) = Q(J,K,IQ)
 160  CONTINUE

      WRITE(6,*) '>> INITIALIZATION COMPLETED.'
      WRITE(6,*) ' '

      RETURN
      END
C =====================================================================
      SUBROUTINE DTCAL(JDIM,KDIM,XYJ,Q,XY)
C =====================================================================
      DIMENSION XYJ(JDIM,KDIM), Q(JDIM,KDIM,4), XY(JDIM,KDIM,4)
      COMMON/CNUMB/JMAX, KMAX
      COMMON/CFLOW/GAMMA, FSMACH, ALPHA
      COMMON/CTIME/DT,CDT,n_steps,ITIME,ITIMEL,TIME

      RCFLMX = 0.

      DO 10 K=1,KMAX
      DO 10 J=1,JMAX

         RR    = 1./Q(J,K,1)
         U     = Q(J,K,2)*RR
         V     = Q(J,K,3)*RR
         QU    = XY(J,K,1)*U + XY(J,K,2)*V
         QV    = XY(J,K,3)*U + XY(J,K,4)*V
         PJ    = (GAMMA-1.)*(Q(J,K,4)-Q(J,K,1)*(U*U+V*V)*.5)
         SONIC = SQRT(GAMMA*PJ*RR)

         RCFL  =  ABS(QU) +ABS(QV) +
     &            SONIC*SQRT( XY(J,K,1)**2+XY(J,K,2)**2
     &                       +XY(J,K,3)**2+XY(J,K,4)**2 )

         IF (RCFL.GT.RCFLMX) RCFLMX = RCFL

 10   CONTINUE

      DT    = CDT/RCFLMX

      WRITE(6,*) '>> TIME STEP WAS RENEWED. : DT =',DT

      RETURN
      END
C=======================================================================
      SUBROUTINE mesh_io(ID,mesh_file,JDIM,KDIM,X,Y,xyj,xy)
C=======================================================================
      DIMENSION X(JDIM,KDIM), Y(JDIM,KDIM),
     &          XYJ(JDIM,KDIM), XY(JDIM,KDIM,4)
      COMMON/CNUMB/JMAX, KMAX
      character(32) :: mesh_file

      MT     = 1
      open(mt,file=mesh_file,form='formatted')

      WRITE(6,*) ' '
      GO TO (10,20),ID
C.......................... READ GRID .................
 10   CONTINUE
        REWIND MT
        READ(MT,*) JMAX,KMAX
        WRITE(6,*) 'JMAX = ',JMAX,',   KMAX = ',KMAX
        IF (JMAX.GT.JDIM.OR.KMAX.GT.KDIM) THEN
          WRITE(6,*) '??? DIMENSION OVER IN meshio ???'
          WRITE(6,*) 'JDIM, KDIM = ',JDIM,KDIM
          STOP
        ENDIF

        READ(MT,*) ( ( X(J,K), J=1,JMAX ) ,K=1,KMAX )
     &           , ( ( Y(J,K), J=1,JMAX ) ,K=1,KMAX )
        REWIND MT
        WRITE(6,*) '>> GRID WAS READ FROM FILE-',MT
        WRITE(6,*) ' '

C....<< COMPUTE  XYJ: JACOBIANS OF METRIC TRANSFORMATIONS >>

      CALL JACOB(JDIM,KDIM,X,Y,XYJ)
      CALL XYMETS(JDIM,KDIM,X,Y,XYJ,XY)

      go to 90

 20   CONTINUE
C........................... WRITE GRID .................
        REWIND MT
          WRITE(MT,*) JMAX,KMAX
          WRITE(MT,100) ( ( X(J,K) ,J=1,JMAX ) ,K=1,KMAX )
     &              ,   ( ( Y(J,K) ,J=1,JMAX ) ,K=1,KMAX )
        REWIND MT

        WRITE(6,*) '>> GRID WAS WRITTEN ON FILE-',MT
        WRITE(6,*) ' '

 90   continue
      close(MT)
      
 100  format(10e14.5)

      RETURN
      end
C=======================================================================
      SUBROUTINE q_io(ID,q_file,JDIM,KDIM,Q)
C=======================================================================
      DIMENSION Q(JDIM,KDIM,4)
      COMMON/CNUMB/JMAX, KMAX
      COMMON/CTIME/DT,CDT,n_steps,ITIME,ITIMEL,TIME
      COMMON/CFLOW/GAMMA, FSMACH, ALPHA
      character(32) :: q_file

      mt    = 2
      open(mt,file=q_file,form='formatted')

      WRITE(6,*) ' '
      GO TO (30,40),ID

C........................... READ FLOW .................
 30   CONTINUE
        REWIND MT
          READ(MT,*) JMAXQ,KMAXQ

        IF (JMAXQ.NE.JMAX .OR. KMAXQ.NE.KMAX) THEN
          WRITE(6,*) 'MISMATCH IN INPUTTED GRID AND Q'
          WRITE(6,*) 'JMAX, KMAX FOR GRID ARE ',JMAX,KMAX
          WRITE(6,*) 'JMAX, KMAX FOR Q    ARE ',JMAXQ,KMAXQ
          STOP
        ENDIF

          read(mt,*) fsmach, alpha, re, time
          read(MT,*) ((( Q(J,K,I),j=1,jmax),k=1,kmax),I=1,4 )

        REWIND MT
C
        WRITE(6,*)'>> THIS IS A RESTART RUN'
        WRITE(6,*)'>> COMPUTED RESULTS WERE READ FROM FILE-',MT
        WRITE(6,*)'      ITIME, TIME = ',ITIME,TIME
        WRITE(6,*)' '

       go to 90
                         
C........................... WRITE Q .................
 40   CONTINUE
CDS        alpha = 0.
        REWIND MT

          WRITE(MT,*) JMAX,KMAX
          write(mt,100) fsmach, alpha, re, time
          WRITE(MT,100) ((( Q(J,K,I),j=1,jmax),k=1,kmax),I=1,4 )

        REWIND MT

        WRITE(6,*) '>> COMPUTED RESULTS WERE WRITTEN ON FILE-',MT
        WRITE(6,*)'      ITIME, TIME = ',ITIME,TIME
        WRITE(6,*) ' '

  90  continue
      close(MT)

 100  format(10e14.5)

      RETURN
      END
C====================<< METRIC CALCULATIONS >>==========================
      SUBROUTINE JACOB(JDIM,KDIM,X,Y,XYJ)
C=======================================================================
      COMMON/CNUMB/JMAX, KMAX
      DIMENSION X(JDIM,KDIM), Y(JDIM,KDIM), XYJ(JDIM,KDIM), XYT(4)
C
C.......................................................................
C  COMPUTE  XYJ:       JACOBIANS OF METRIC TRANSFORMATIONS
C                    = ( DXI/DX * DETA/DY ) - ( DXI/DY * DETA/DX )
C
C       X:       CARTESIAN X COORDINATES OF THE GRID.
C       Y:       CARTESIAN Y COORDINATES OF THE GRID.
C
C ....... << ASSUME GRID IS FIXED IN TIME >>
C
C         XIT = 0.    : DXI /DT
C         ETT = 0.    : DETA/DT
C.......................................................................
C
C
      DO 10 K=1,KMAX
      DO 10 J=1,JMAX
C
      CALL XYDIF(JDIM,KDIM,X,Y,XYT,J,K)
C
C......... << FORM JACOBIAN >>
C
      DINV = 1./ ( XYT(4) * XYT(1) - XYT(3) * XYT(2) )
      XYJ(J,K) = DINV
C
C......... << CHECK NEGATIVE JACOBIAN >>
C
      IF(DINV .LE. 0.0)THEN
        WRITE(6,*)'NEGATIVE JACOBIAN AT J,K = ',J,K,',; XYJ =',XYJ(J,K)
      ENDIF
C
   10 CONTINUE
C
      WRITE(6,*) '>> METRIC JACOBIAN WAS COMPUTED'
      WRITE(6,*) ' '
C
      RETURN
      END
C====================<< METRIC CALCULATIONS >>==========================
      SUBROUTINE XYMETS(JDIM,KDIM,X,Y,XYJ,XY)
C=======================================================================
      DIMENSION X(JDIM,KDIM), Y(JDIM,KDIM), XYJ(JDIM,KDIM), XYT(4)
      DIMENSION XY(JDIM,KDIM,4)

      COMMON/CNUMB/JMAX, KMAX
C
C... << CALCULATE METRICS >>
C                           D(XI)/DX = XY(1),     D(ETA)/DX = XY(3)
C                           D(XI)/DY = XY(2),     D(ETA)/DY = XY(4)
C
      DO 10 K=1,KMAX
      DO 10 J=1,JMAX

        CALL XYDIF(JDIM,KDIM,X,Y,XYT,J,K)

         DINV       =   XYJ(J,K)
         XY(J,K,1)  =   XYT(1)*DINV
         XY(J,K,2)  = - XYT(2)*DINV
         XY(J,K,3)  = - XYT(3)*DINV
         XY(J,K,4)  =   XYT(4)*DINV

  10  CONTINUE

      RETURN
      END
C====================<< METRIC CALCULATIONS >>==========================
      SUBROUTINE XYDIF(JDIM,KDIM,X,Y,XY,J,K)
C=======================================================================
      DIMENSION X(JDIM,KDIM), Y(JDIM,KDIM), XY(4)
      COMMON/CNUMB/JMAX, KMAX
C
C........ << COMPUTE XI DERIVATIVES OF X, Y >>
C
C  XY4 =D(X)/D(XI), XY3 = D(Y)/D(XI)
C
      IF    (J.EQ.   1)THEN
            XY(4) = .5*( -3.*X(J,K) +4.*X(J+1,K) - X(J+2,K) )
            XY(3) = .5*( -3.*Y(J,K) +4.*Y(J+1,K) - Y(J+2,K) )
      ELSEIF(J.EQ.JMAX)THEN
            XY(4) = ( 3.*X(J,K) -4.*X(J-1,K) + X(J-2,K) )*.5
            XY(3) = ( 3.*Y(J,K) -4.*Y(J-1,K) + Y(J-2,K) )*.5
      ELSE
            XY(4) = ( X(J+1,K) - X(J-1,K) )*.5
            XY(3) = ( Y(J+1,K) - Y(J-1,K) )*.5
      ENDIF
C
C........ << COMPUTE ETA DERIVATIVES OF X, Y >>
C
C  XY2 = D(X)/D(ETA),   XY1 = D(Y)/D(ETA)
C
      IF    (K.EQ.   1)THEN
            XY(2) = ( -3.*X(J,K) +4.*X(J,K+1) - X(J,K+2) )*.5
            XY(1) = ( -3.*Y(J,K) +4.*Y(J,K+1) - Y(J,K+2) )*.5
      ELSEIF(K.EQ.KMAX)THEN
            XY(2) = ( 3.*X(J,K) -4.*X(J,K-1) + X(J,K-2) )*.5
            XY(1) = ( 3.*Y(J,K) -4.*Y(J,K-1) + Y(J,K-2) )*.5
      ELSE
            XY(2) = ( X(J,K+1) - X(J,K-1) )*.5
            XY(1) = ( Y(J,K+1) - Y(J,K-1) )*.5
      ENDIF

      RETURN
      END

C=======================================================================
      SUBROUTINE surf_io(JDIM,KDIM,X,Y,Q)
C=======================================================================
      DIMENSION X(JDIM,KDIM), Y(JDIM,KDIM),
     &          XYJ(JDIM,KDIM), XY(JDIM,KDIM,4)
      DIMENSION Q(JDIM,KDIM,4)
      DIMENSION P(JDIM),CP(JDIM)
      COMMON/CNUMB/JMAX, KMAX
      COMMON/CFLOW/GAMMA, FSMACH, ALPHA
      COMMON/CTIME/DT,CDT,n_steps,ITIME,ITIMEL,TIME
      common/cairfoil/jtail

CDS      write(*,*)'gamma=',gamma
CDS      write(*,*)'fsmach=',fsmach
CDS      write(*,*)'alpha=',alpha

      COSANG = COS(3.141593*alpha/180.)
      SINANG = SIN(3.141593*alpha/180.)

CDS      write(*,*)'COSANG=',COSANG
CDS      write(*,*)'SINANG=',SINANG

      GAMI    = GAMMA-1.
      PINF    = 1.0/GAMMA
      RHOINF  = 1.0
CDS freestream speed of sound
      AINF    = 1.0
      UINF    = FSMACH*COSANG
      VINF    = FSMACH*SINANG
      VELINF  = FSMACH*AINF

      jtail1 = jtail
      jtail2 = jmax-jtail1+1

      do j=1,jmax
        p(j) =0.0
        cp(j)=0.0
      enddo

      mt    = 1
      open(mt,file='work/surf.dat',form='formatted')

CDS P = (gamma-1)*rho*[e0-0.5*VEL**2]
CDS pinf = 1/gamma
CDS CP = (p-pinf)/(0.5*rhoinf*VELinf**2)
CDS E0 = Q4/rho  (stagnation energy)
CDS E0inf = [1/(gamma-1)]*(pinf/rhoinf)+0.5*VELinf**2
CDS M = vel/c
CDS c = sqrt(gamma*p/rho)

      do j=jtail1,jtail2
        rho=q(j,1,1)
        u=q(j,1,2)/rho
        v=q(j,1,3)/rho
        vel=sqrt(u**2+v**2)
        e0=Q(j,1,4)/rho
        p(j)=GAMI*rho*(e0-0.5*vel**2)
        CP(j)=(p(j)-pinf)/(0.5*rhoinf*velinf**2)
        surfM = vel/sqrt(gamma*p(j)/rho)
        write(mt,100)x(j,1),y(j,1),rho,u,v,surfM,P(j),CP(j),-CP(j)
      enddo

      close(mt)

      cl =0.0
      clp=0.0
      cdp=0.0
      fx =0.0
      fy =0.0

      jj=0
CDs      do j=jtail1+1,(jmax+1)/2-1
      do j=jtail1,(jmax+1)/2-1
        jjl=(jmax+1)/2-jj
        jju=(jmax+1)/2+jj
        cpl=0.5*(cp(jjl)+cp(jjl-1))
        cpu=0.5*(cp(jju)+cp(jju+1))
CDS        clu=clu+(cpl-cpu)*(x(jju+1,1)-x(jju,1))
CDS        cll=cll+(cpl-cpu)*(x(jjl-1,1)-x(jjl,1))
CDS        write(*,*)j,jj,jjl,jju
        dxl= x(jjl-1,1)-x(jjl,1)
        dxu= x(jju+1,1)-x(jju,1)
        cl=cl+cpl*dxl
        cl=cl-cpu*dxu
CDS
CDS  ** cl and cd from surface pressure **
CDS
        dyl= y(jjl-1,1)-y(jjl,1)
        dyu= y(jju+1,1)-y(jju,1)
        dsl= sqrt(dxl**2+dyl**2)
        dsu= sqrt(dxu**2+dyu**2)
        thl= atan(dyl/dxl)
        thu= atan(dyu/dxu)
        pl =  0.5*(p(jjl)+p(jjl-1))
        pu = -0.5*(p(jju)+p(jju+1))
        fxl= -pl*dsl*sin(thl)
        fyl=  pl*dsl*cos(thl)
        fxu= -pu*dsu*sin(thu)
        fyu=  pu*dsu*cos(thu)
        fx =  fx + fxl + fxu
        fy =  fy + fyl + fyu
CDS     <back up cl and cd computation>
CDS        fliftl = -fxl*SINANG + fyl*COSANG
CDS        fdragl =  fxl*COSANG + fyl*SINANG
CDS        fliftu = -fxu*SINANG + fyu*COSANG
CDS        fdragu =  fxu*COSANG + fyu*SINANG
CDS        clp = clp + (fliftl+fliftu)*2/VELINF**2
CDS        cdp = cdp + (fdragl+fdragu)*2/VELINF**2
        jj=jj+1
      enddo

CDS      write(*,*)'fx=',fx
CDS      write(*,*)'fy=',fy
CDS      clift=(-fx*SINANG+fy*COSANG)*2/VELINF**2
CDS      cdrag=( fx*COSANG+fy*SINANG)*2/VELINF**2
CDS      write(*,*)'cl1=',clift
CDS      write(*,*)'cd1=',cdrag

      clp=(-fx*SINANG+fy*COSANG)*2/VELINF**2
      cdp=( fx*COSANG+fy*SINANG)*2/VELINF**2

      write(*,*)'Cl (integration of Cp_lower - Cp_upper) =',cl
      write(*,*)'Cl (surface integration of pressure) =',clp
      write(*,*)'Cd (surface integration of pressure) =',cdp
CDS      write(*,*)'Cl (computed from Cp) =',cll
CDS      write(*,*)'Cl (computed from Cp) =',clu

 100  format(10e14.5)

      RETURN
      END

C=======================================================================
      SUBROUTINE res_io(NRES,RES)
C=======================================================================
      DIMENSION RES(10000)

      mt    = 10
      open(mt,file='work/residual.dat',form='formatted')

      do i=1,nres
        write(mt,110)i*10,res(i)
      enddo

      close(mt)
 110  format(i7,e14.5)

      RETURN
      END
