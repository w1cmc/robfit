      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*256 RESOUT,RESOUTORD
      CHARACTER*80 CTEMP
      DIMENSION NORDER(:),Egy(:),RES(:),Z(:)
      NUMARGS=NARGS()
      IF(NUMARGS.GT.1)THEN
        CALL BGETARG(1,RESOUT)
      ENDIF
      ILCML=ILB(RESOUT)
      OPEN(1,FILE=RESOUT)
      J=0
      DO I=1,10000
        READ(1,'(A20)',END=20,ERR=20)CTEMP
        IF(CTEMP(:20).NE.'                    ')J=J+1
      ENDDO
20    JM=J      
      ALLOCATE(NORDER(JM),EGY(JM),RES(JM),Z(JM))
      REWIND(1)
      J=0
      DO I=1,10000
        READ(1,'(A)',END=20,ERR=20)CTEMP
        IF(CTEMP(:20).NE.'                    ')THEN
          J=J+1
          READ(CTEMP,*)Egy(J),RES(J),Z(J)
        ENDIF  
        IF(J.EQ.JM)GOTO 30
      ENDDO
30    CONTINUE 
      CLOSE(1) 
      PRINT*,'ENTER NPOW'
      READ(*,*)NPOW
      CALL SORT(NORDER,RES,JM)
      RESOUTORD=RESOUT(:ILCML-4)//'.ORD'
      OPEN(1,FILE=RESOUTORD)
      DO I=1,JM
        WRITE(1,'(1P,3E15.6)')RES(NORDER(I))**(1/(2D0*NPOW)),
     2   Egy(NORDER(I)),Z(NORDER(I))
      ENDDO
      CLOSE(1)
      END
          
        
        
        
        
        
      
      
      
      SUBROUTINE SORT(NORDER,RES,NTOTAL)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION NORDER(NTOTAL),RES(NTOTAL)
      NORDER(1)=1
      DO I=2,NTOTAL
        CALL TLOC(NORDER,RES,RES(I),I-1,J)
        DO K=I-1,J+1,-1
          NORDER(K+1)=NORDER(K)
        ENDDO
        NORDER(J+1)=I
      ENDDO
      RETURN
      END
      SUBROUTINE TLOC(NORDER,ENERGY,X,NMAX,J)
C CODE RETURNS J SUCH THAT ENERGY(NORDER(J))<=X<=ENERGY(NORDER(J+1)
C J=0 ==> X < ENERGY(NORDER(1)
C J=NMAX ==> X > ENERGY(NORDER(NMAX))
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION NORDER(NMAX),ENERGY(NMAX)
      JL=0
      JU=NMAX+1
10    IF(JU-JL.LE.1)GOTO 20
      JM=(JU+JL)/2
      IF(X.GT.ENERGY(NORDER(JM)))THEN
      JL=JM
      ELSE
      JU=JM
      ENDIF
      GOTO 10
20    J=JL
      RETURN
      END      
      
