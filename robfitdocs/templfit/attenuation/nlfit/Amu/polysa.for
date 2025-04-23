      PROGRAM POLY
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(2),PF(91),CONS(91),MPPF(64,1)
      open( unit=10, file='fapfa.tmp',access='sequential',
     2      form='unformatted',recordtype='fixed',SHARE='DENYNONE' )
      READ(10)NV,ND,MP,NDAT
      READ(10)(CONS(I),I=1,NV)
C      PRINT*,' I     CONS(I)'
C      DO I=1,NV
C        PRINT*,I,CONS(I)
C      ENDDO   
      REWIND 10
      open( unit=9, file='xdat.tmp',access='sequential',
     2      form='unformatted',recordtype='fixed',SHARE='DENYNONE' )
      DO J=1,NDAT
        READ(9)(X(I),I=1,ND)
        CALL POLY(X,PF,CONS,FA,MPPF,NV,ND,MP)
        WRITE(10)FA,(PF(I),I=1,NV)
      ENDDO     
      STOP 
      END
