C============== I/O ====================================================
      SUBROUTINE MESH_IO(ID, JDIM,KDIM, JMAX,KMAX, X,Y, MT)
C=======================================================================
      DIMENSION X(JDIM,KDIM),Y(JDIM,KDIM)

      OPEN (mt, FILE='work/mesh.xyz')	! file name of plot3d-format output

      REWIND MT
      IF(ID.EQ.1)THEN 
        READ(MT,*) JMAX,KMAX
        READ(MT,*) ((X(J,K),J=1,JMAX),K=1,KMAX),
     &             ((Y(J,K),J=1,JMAX),K=1,KMAX)
      ELSEIF(ID.EQ.2)THEN
        WRITE(MT,*) JMAX,KMAX
        WRITE(MT,100) ( ( X(J,K),J=1,JMAX ) , K=1,KMAX ),
     &                ( ( Y(J,K),J=1,JMAX ) , K=1,KMAX )
        WRITE(6,*) ' '
        WRITE(6,*) 'Generated mesh was written on (','mesh.xyz',')'
      ENDIF
      REWIND MT
      CLOSE (mt)

 100  format(10e14.5)

      RETURN
      END