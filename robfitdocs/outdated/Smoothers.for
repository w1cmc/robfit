      SUBROUTINE SM17(DAT,NDAT)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*4 DAT(NDAT),TDAT(-8:8)
      C0=0
      DO I=-8,8
        TDAT(I)=DAT(I+9)
        C0=C0+TDAT(I)
      ENDDO
      C0=C0/17
      DO I=1,8
        DAT(I)=C0
      ENDDO  
      DO I=9,NDAT-8
        DO J=-8,7
          TDAT(J)=TDAT(J+1)
        ENDDO
        TDAT(8)=DAT(I+8)
        C0=0 
        DO J=-8,8
          C0=C0+TDAT(J)
        ENDDO
        C0=C0/17
        DAT(I)=C0
      ENDDO
      DO I=1,8
        DAT(NDAT-8+I)=C0
      ENDDO  
      RETURN
      END

      SUBROUTINE SMHW(DAT,NDAT)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*4 DAT(NDAT),TDAT(-12:12)
      C0=0
      DO I=-12,12
        TDAT(I)=DAT(I+13)
c        C0=C0+TDAT(I)*HOFXw(1D0*I)
c        write(27,*)1d0*i,hofxw(1d0*i)
      ENDDO
      DO I=1,12
        DAT(I)=C0
      ENDDO  
      DO I=13,NDAT-12
        DO J=-12,11
          TDAT(J)=TDAT(J+1)
        ENDDO
        TDAT(12)=DAT(I+12)
        C0=0 
        DO J=-12,12
c          C0=C0+TDAT(J)*HOFXw(1D0*J)
        ENDDO
        DAT(I)=C0
      ENDDO
      DO I=1,12
        DAT(NDAT-12+I)=C0
      ENDDO  
      RETURN
      END      

      Function hofxw(x)
c \public_html\Fittery\robfit\CnmFit\Hoff.docx
C function extends from -12 to 12 and is 0.5 at +- 3
C the integral of h(x) form -infinity to infinity is W
C which is divided out here
      implicit real*8 (a-h,o-z)
      Data W/2D0/
      hofxw=(aigau(3+x/W)*aigau(3-x/W))/(6*W)
      return
      end

