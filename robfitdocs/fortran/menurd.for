      SUBROUTINE MENURD(FFCOM,FFVAR,IR,ILT)
      CHARACTER*256 FFVAR(*)
      CHARACTER*40 FFCOM(*)
      CHARACTER*140 CTEMP
      CHARACTER*256 NAPMLSP
      CHARACTER*1 CANS(140)
      DIMENSION IR(20)
      EQUIVALENCE (CTEMP(1:1),CANS(1))
      COMMON /NAMES/NAPMLSP      
      DATA CTEMP/' '/
      ICHANG=0
      ILR=0
1000  CONTINUE
      DO I=1,IR(1)
        IF(I.GT.2)THEN
          IWRT=LENTRIM(FFVAR(I))
          IWRT=MAX(1,IWRT)
          WRITE(*,'(I3,1X,A40,A)')I-2,FFCOM(I),FFVAR(I)(1:IWRT)
        ELSE
          WRITE(*,'(3X,A40)')FFCOM(I)
        ENDIF
      ENDDO
      DO I=1,2
        PRINT*,' '
      ENDDO
25    PRINT*,' ENTER MENU #; A TO ABORT; OR <CR> TO RUN'
      CTEMP=' '
      READ(*,'(A)',END=3000)CTEMP
      IF(CTEMP(:1).EQ.'A'.OR.CTEMP(:1).EQ.'a')STOP
      IF(CTEMP.EQ.' ')RETURN
      READ(CTEMP,*)LLINE
      LLINE=LLINE+1
2000  CONTINUE

        IOFF=0
        DO I=1,LLINE-1
          IOFF=IOFF+IR(I)
        ENDDO
        IOFFST=IOFF+2
      DO I=1,IR(LLINE)
          IOFF=IOFF+1
        IF(I.GT.2)THEN
          IWRT=MAX(1,LENTRIM(FFVAR(IOFF)))
          WRITE(*,'(I3,1X,A40,A)')I-2,FFCOM(IOFF),
     2      FFVAR(IOFF)(1:IWRT)
        ELSE
          WRITE(*,'(3X,A40)')FFCOM(IOFF)
        ENDIF
      ENDDO
      DO I=1,2
        PRINT*,' '
      ENDDO
      WRITE(*,5300)
5300  FORMAT(' ENTER LINE NUMBER TO CHANGE A LINE Y <--> N
     1'/'  <CR> TO RETURN TO THE MAIN MENU; LN # OR
     1'/'   LN #,VAR TO CHANGE A LINE; LN #,&D TO DELETE A LINE;
     1'/'   OR LN #,&A TO ADD A COMPLETE NEW LINE')
      DO I=1,2
        PRINT*,' '
      ENDDO
      CTEMP=' '
      READ(*,'(A)',ERR=1000,END=1000)CTEMP
      IF(CTEMP.EQ.' ')GOTO 1000
      READ(CTEMP,*)LSLINE
      ICOMMA=INDEX(CTEMP,',')
      IF(ICOMMA.EQ.0)ICOMMA=INDEX(CTEMP,' ')-1
      ICP1=ICOMMA+1
      IT=IOFFST+LSLINE
      IF(FFVAR(IT).EQ.'N '.AND.CTEMP(ICP1:ICP1).EQ.' ')THEN
        FFVAR(IT)='Y '
        ICHANG=1
        GOTO 2000
      ENDIF
      IF(FFVAR(IT).EQ.'Y '.AND.CTEMP(ICP1:ICP1).EQ.' ')THEN
        FFVAR(IT)='N '
        ICHANG=1
        GOTO 2000
      ENDIF
      IF(FFVAR(IT).EQ.'Log '.AND.CTEMP(ICP1:ICP1).EQ.' ')THEN
        FFVAR(IT)='Linear '
        ICHANG=1
        GOTO 2000
      ENDIF
      IF(FFVAR(IT).EQ.'Linear '.AND.CTEMP(ICP1:ICP1).EQ.' ')THEN
        FFVAR(IT)='Log '
        ICHANG=1
        GOTO 2000
      ENDIF
        
      IF(CTEMP(ICP1:ICP1+1).EQ.'&a'.OR.CTEMP(ICP1:ICP1+1).EQ.'&A')THEN
        IR(LLINE)=IR(LLINE)+1
        DO I=ILT,IT,-1
          FFCOM(I+1)=FFCOM(I)
          FFVAR(I+1)=FFVAR(I)
        ENDDO
        ILT=ILT+1
        PRINT*,' ENTER THE NEW COMMENT'
          READ(*,'(A)')CTEMP
        FFCOM(IT)=CTEMP
        PRINT*,' ENTER THE NEW VARIABLE'
          READ(*,'(A)')CTEMP
        FFVAR(IT)=CTEMP
          ICHANG=1
        GOTO 2000
      ENDIF
      IF(CTEMP(ICP1:ICP1+1).EQ.'&d'.OR.CTEMP(ICP1:ICP1+1).EQ.'&D')THEN
        IR(LLINE)=IR(LLINE)-1
        IT=IOFFST+LSLINE
        DO I=IT,ILT-1
          FFCOM(I)=FFCOM(I+1)
          FFVAR(I)=FFVAR(I+1)
        ENDDO
        ILT=ILT-1
        ICHANG=1
        GOTO 2000
      ENDIF
      IF(LSLINE.NE.0)THEN
        ICHANG=1
        IT=IOFFST+LSLINE
        IF(CTEMP(ICP1:).NE.' ')THEN
          FFVAR(IT)=CTEMP(ICP1:)
        ELSE
          WRITE(*,'('' ENTER THE NEW VARIABLE'')')
          READ(*,'(A)')FFVAR(IT)
        ENDIF
        GOTO 2000
      ENDIF
3000  RETURN      
      END

      SUBROUTINE CTON(CI,NMISS,ANUM)
      CHARACTER*64ci
      CHARACTER*1 KNUM(10),C
      DATA KNUM/'1','2','3','4','5','6','7','8','9','0'/
      POIFLG=0.
      SFLAG=0.
      NDP=0
      NUM=0
      NMISS=0
      DO 10 I=1,80
      NONUM=0
      C=CI(I:I+1)
      IF(C.EQ.' ')GOTO 100
        IF(C.EQ.'.')POIFLG=1.
        IF(C.EQ.'-')SFLAG=1.
        IF(C.EQ.'.')GOTO 10
        IF(C.EQ.'-')GOTO 10
        IF(POIFLG.EQ.1.)NDP=NDP+1
        DO 20 J=1,10
        IF(C.EQ.KNUM(J))THEN
          IF(J.NE.10)THEN
            NUM=NUM*10+J
          ELSE
            NUM=NUM*10
          ENDIF
          NONUM=1
        ENDIF
 20     CONTINUE
        IF(NONUM.NE.1)THEN
          NMISS=NMISS+1
        ENDIF
 10   CONTINUE
100   CONTINUE
      ANUM=NUM
      IF(SFLAG.EQ.1.)ANUM=-ANUM
      DO 200 J=1,NDP
      ANUM=ANUM/10
200   CONTINUE
      IF(NMISS.NE.0)THEN
        PRINT*,'THIS IS NOT AN INTEGER NMISS=',NMISS
      ENDIF
        CLOSE(17)
      RETURN
      END
      
      SUBROUTINE CHANGE(STR1,L1,STR2,L2,STR3,L3)
      CHARACTER*40 STR1,STR2,STR3,CTEMP
      IF(STR1(1:L1).NE.STR3(1:L1))RETURN
      CTEMP=STR3(L1+1:L3)
      STR3=STR2(1:L2)//CTEMP
      RETURN
      END
                        
      SUBROUTINE CCHAR(NTAPE,NC)
      CHARACTER*64 NTAPE
      NC=64
      DO 10 I=1,64
      IF(NTAPE(NC:NC).NE.' '.AND.NTAPE(NC:NC).NE.'\')GOTO 20
      NC=NC-1
10    CONTINUE
20    RETURN
      END
      SUBROUTINE CLBS(NTAPE,NC)
      CHARACTER*64 NTAPE
      NC=64
      DO 10 I=1,64
      IF(NTAPE(NC:NC).EQ.'\\')GOTO 20
      NC=NC-1
10    CONTINUE
20    RETURN
      END
      SUBROUTINE FLASTS(NTAPE,NC,SY)
      CHARACTER*(*) NTAPE
      CHARACTER*1 SY
      NC=64
      DO 10 I=1,64
      IF(NTAPE(NC:NC).EQ.SY)GOTO 20
      NC=NC-1
10    CONTINUE
20    RETURN
      END
