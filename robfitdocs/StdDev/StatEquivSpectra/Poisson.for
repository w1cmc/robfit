      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/RAN/IX,IY,ITAB(128)
      REAL*4 FIN(:),FOUT(:)
      DATA N/100/
      ALLOCATE (FIN(N),FOUT(N))
      IX=131312345
      IY=757532145
      CALL RSEED(IX,IY)
      print*,'enter XB,X0,XE,HGT,W'
      READ(*,*)XB,X0,XE,HGT,W
      H=(XE-XB)/N
      OPEN(1,FILE='InSpect.out')
      DO I=1,N
        X=XB+(I-0.5D0)*H
        FIN(I)=HGT*H*(EXP(-((X-X0)/W)**2)/2)
        WRITE(1,'(1P,2E17.6)')X,FIN(I)
      ENDDO
      CLOSE(1)
      OPEN(2,FILE='OutSpect.out')
      CALL GENPOS(FIN,FOUT,N)  
      DO I=1,N
        X=XB+(I-0.5D0)*H
        WRITE(2,'(1P,2E17.6)')X,FOUT(I)
      ENDDO
      CLOSE(2)
      CALL FSYSTEM('..\..\GPLOT InSpect.out OutSpect.out')
      END
        
      
      

