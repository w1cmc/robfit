      FUNCTION AIGAU(XT)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION CON1(18),CON2(15)
      DATA SQRPI/1.772453850905516D0/
      DATA CON1/0.5D0,-0.9452049927406015D0,0.8901349154719547D0,
     2 -0.5265323590674972D0,0.2112033098149304D0,
     3 -0.5851110587267722D-01,0.1093104379246567D-01,
     4 -0.1264586529819479D-02,0.7007918572074021D-04,
     5 -3.018789152576802D0,4.186608620562856D0,
     6 -3.510610291883386D0,1.967990042365595D0,
     7 -0.7679453864061583D0,0.2096579949268288D0,
     8 -0.3887374800121166D-01,0.4482842710646635D-02,
     9 -0.2484242423263732D-03/  
      DATA CON2/-0.6931471805599452D0,1.995915713776837D0,
     2 -2.582310623341017D0,1.953688092571124D0,-0.9387916412290680D0,
     3 0.2892445861593335D0,-0.5400784330829852D-01,
     4 0.5078739222313374D-02,-0.1061251766837922D-03,      
     5 -1.251590673686981D0,0.7695652871670581D0,-0.2680046704303274D0,
     6 0.5332127596813570D-01,-0.5079508799796978D-02,
     7 0.1060850657795260D-03/
      IF(XT.LT.-27)THEN
C Tests of this limit are in Welcome.doc      
        AIGAU=0
        RETURN
      ELSEIF(XT.LT.-6.2D0)THEN
C Discussion in public_html/nlfit/testdata/genint/Expanding the integrand.htm      
        BR=ABRACK(XT)
        EXPP=BEXP(-XT*XT)
        AIGAU=-EXPP*BR/(2*XT*SQRPI)
      ELSEIF(XT.LT.0)THEN
C DISCUSION IN public_html\nlfit\ExpTimesPade.htm
        ANUM=CON1(1)+XT*(CON1(2)+XT*(CON1(3)+XT*(CON1(4)+XT*(CON1(5)
     2  +XT*(CON1(6)+XT*(CON1(7)+XT*(CON1(8)+XT*CON1(9))))))))
        ADEN=1+XT*(CON1(10)+XT*(CON1(11)+XT*(CON1(12)+XT*(CON1(13)
     2  +XT*(CON1(14)+XT*(CON1(15)+XT*(CON1(16)+XT*(CON1(17)+
     3  XT*CON1(18)))))))))
        BR=ANUM/ADEN
        AIGAU=BEXP(-XT*XT)*BR
      ELSEIF(XT.LT.7)THEN
C DISCUSSION IN ..\public_html\nlfit\MPNlfitAigau.htm
        XN=-XT
        ANUM=CON2(1)+XN*(CON2(2)+XN*(CON2(3)+XN*(CON2(4)+XN*(CON2(5)
     2  +XN*(CON2(6)+XN*(CON2(7)+XN*(CON2(8)+XN*CON2(9))))))))
        ADEN=1+XN*(CON2(10)+XN*(CON2(11)+XN*(CON2(12)+XN*(CON2(13)
     2  +XN*(CON2(14)+XN*CON2(15))))))
        EXPAD=ANUM/ADEN
        AIGAU=1-BEXP(EXPAD)
      ELSE
        AIGAU=1
      ENDIF
      RETURN
      END

      FUNCTION ABRACK(X)
      IMPLICIT REAL*8 (A-H,O-Z)
      ABRACK=1
      XP=1
      X2=X*X
      ANUM=1
      IS=-1
      M=1
      DEN=2*X2
      ALT=1
5     CONTINUE
      ALT=ALT*M/DEN
      IF(M.GT.1001)THEN
        STOP ' ABRACK NOT CONVERGING '
      ENDIF
      M=M+2
      ABRACK=ABRACK+IS*ALT
      IS=-IS
      IF(ALT.GT.1D-16)GOTO 5
      RETURN
      END

      Function hofxg(x)
c \public_html\Fittery\robfit\CnmFit\Hoff.docx
C function extends from -12 to 12 and is 0.5 at +- 3
C the integral of h(x) form -infinity to infinity is W
C which is divided out here
      implicit real*8 (a-h,o-z)
      Data W/2D0/
      hofxg=(aigau(3+x/W)*aigau(3-x/W))/(6*W)
      return
      end
      
      SUBROUTINE SMH(DAT,NDAT)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*4 DAT(NDAT),TDAT(-12:12)
      C0=0
      DO I=-12,12
        TDAT(I)=DAT(I+13)
        C0=C0+TDAT(I)*HOFXg(1D0*I)
        write(27,*)1d0*i,hofxg(1d0*i)
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
          C0=C0+TDAT(J)*HOFXg(1D0*J)
        ENDDO
        DAT(I)=C0
      ENDDO
      DO I=1,12
        DAT(NDAT-12+I)=C0
      ENDDO  
      RETURN
      END
