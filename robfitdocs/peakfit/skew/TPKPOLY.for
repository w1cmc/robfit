        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION ASKEW(6),AJUMPV(6)
        CHARACTER*1 CSERIES(6)
        DATA CSERIES/'1','2','3','4','5','6'/
        DATA ASKEW/.4D0,.6d0,.8d0,1d0,1.2d0,1.4d0/
        DATA AJUMPV/6*0D0/
        DATA DELTA/1D-4/
        GOTO 10
5       PRINT*,' Enter x,skew,ajump '
        READ(*,*)X,SKEW,AJUMP
        IF(SKEW.LT.0)GOTO 10
        CALL PKPOLM(X,SKEW,AJUMP,PX,PSKEW,PJUMP,FA1)
        PRINT*,' AFTER FIRST PKPOLYM  ',FA1
        XT=X-DELTA
        CALL PKPOLM(XT,SKEW,AJUMP,PXT,PSKEWT,PJUMPT,FAM)
        XT=X+DELTA
        CALL PKPOLM(XT,SKEW,AJUMP,PXT,PSKEWT,PJUMPT,FAP)
        dPKPdxn=(FAP-FAM)/(2*DELTA)
        PRINT*,'   PX   ',PX
        PRINT*,'dPKPdxn ',dPKPdxn
        SKEWT=SKEW-DELTA
        CALL PKPOLM(X,SKEWT,AJUMP,PXT,PSKEWT,PJUMPT,FAM)
        SKEWT=SKEW+DELTA
        CALL PKPOLM(X,SKEWT,AJUMP,PXT,PSKEWT,PJUMPT,FAP)
        PSKEWN=(FAP-FAM)/(2*DELTA)
        PRINT*,' PSKEW  ',PSKEW
        PRINT*,' PSKEWN ',PSKEWN
        AJUMPT=AJUMP-DELTA
        CALL PKPOLM(X,SKEW,AJUMPT,PXT,PSKEWT,PJUMPT,FAM)
        AJUMPT=AJUMP+DELTA
        CALL PKPOLM(X,SKEW,AJUMPT,PXT,PSKEWT,PJUMPT,FAP)
        PJUMPN=(FAP-FAM)/(2*DELTA)
        PRINT*,' PJUMP  ',PJUMP
        PRINT*,' PJUMPN ',PJUMPN
        GOTO 5
10      N=1000
        OPEN(1,FILE='pkpol2.grf')
        WRITE(1,*)'<Giraffe>'
        WRITE(1,*)'<ChartTitle>'
        WRITE(1,*)'skew 0.4, 0.6, 0.8, 1.0 1.2, 1.4' 
        WRITE(1,*)'</ChartTitle>'
        WRITE(1,*)'<LeftAxisCaption>'
        WRITE(1,*) 'log(100-2(X+40)+1000*Peak/(1+jump))'
        WRITE(1,*) '</LeftAxisCaption>'
        WRITE(1,*)'<BottomAxisCaption>X=E/W</BottomAxisCaption>'

        XL=-40
        XU=+10
        H=(XU-XL)/1000
        
        
        DO J=1,6
        WRITE(1,*)'<Series'//CSERIES(J)//'>'
        DO I=1,N+1
          X=XL+H*(I-1)
          BKG=0.01
          CALL PKPOLM(X,ASKEW(J),AJUMPV(J),PX,PSKEW,PJUMP,FA)
c          FA=1000*(FA-AJUMPV(J))+BKG
          FA=1000*FA/(1+AJUMPV(J))+BKG
          WRITE(1,'(2E15.6)')X,FA
        ENDDO
        WRITE(1,*)'</Series'//CSERIES(J)//'>'
        ENDDO

        H=(XU-XL)/100
        DO J=1,3
        WRITE(1,*)'<DataPoints'//CSERIES(J)//'>'
        DO I=1,101
          X=XL+H*(I-1)
          BKG=0.01
          FA=EXP(-X*X/2)
          IF(X.LT.-askew(2*j-1)**2)FA=exp(askew(2*j-1)**2*X+
     2     askew(2*j-1)**4/2)
          FA=1000*FA/(1+AJUMPV(J))+BKG
          WRITE(1,'(2E15.6)')X,FA
        ENDDO
        WRITE(1,*)'</DataPoints>'//CSERIES(J)//'>'
        ENDDO

        WRITE(1,*)'</Giraffe>'
        END

