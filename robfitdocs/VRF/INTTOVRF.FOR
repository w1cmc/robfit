      SUBROUTINE INTTOVRF(NUMARGS
      IF(NUMARGS.EQ.3)THEN
        CALL BGETARG(1,NACMLSP)
        ILCML=LENTRIM(NACMLSP)
        NATEMP=NASTR('*.UNF',ILW,EX)
        CALL RDATA(NACMLSP) !makes namedat.unf and puts ibegc,iendc into 
C        intvals
        CALL RCONS(NACMLSP) !SETS DEFAULT IBEG'S AND IEND'S makes com.unf
        DUM(1)=CHNTOEGY(1D0*IDUM(14))
        DUM(2)=CHNTOEGY(1D0*IDUM(15))
        LDUM(4)=.FALSE. !lAGRANGE INERPOLATION
        LDUM(8)=.TRUE. !TURNS OFF MOST ROBMIN PRINTS.
        CALL PMLTOPGS(NACMLSP)
        IF(IDUM(18).EQ.2)CALL RPKS(NACMLSP,EX)
        CALL BGETARG(2,CPASTICHE)
        IF(CPASTICHE(1:8).EQ.'templfit')then
          CALL COMTOUNF()
          CLOSE(1)
c          COMMAND='start /MIN /WAIT '//
          COMMAND=
     2     'Templfit.exe '//NACMLSP(:ILCML)//' mlbat'
          ISTATUS = FSYSTEM(COMMAND)
          CALL NuclidesToPml(NACMLSP) !NANUCFIT UPDATED AND REPLACED
          OPEN(1,FILE='RfitFin.out')
          CALL REMFILES() 
          WRITE(1,*)CPASTICHE,' DONE'
          CLOSE(1)
          IF(.NOT.LDUM(23))THEN
            PRINT*,'ENTER ANYTHING TO STOP'
            READ(*,'(A)')ANS
          ENDIF  
          STOP
        ENDIF  
        GOTO 10 
      ENDIF  
      RETURN
      END
