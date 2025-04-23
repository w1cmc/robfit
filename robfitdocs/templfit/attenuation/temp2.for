        ELSEIF(NASP(ISP-3:ISP).EQ.'.+ml'.OR.NASP(ISP-3:ISP).EQ.
     2   '.+ML')THEN
          CALL CMLREAD('<Data >',NSP(:ISP),'temp.txt',IBNDAT)
          IF(IBNDAT.EQ.0)THEN
            PRINT*,' FILE ',NASP(:ISP),' CONTAINS NO DATA'
            STOP            
          ENDIF 
          IBEGC=1
          IENDC=IBNDAT
          ALLOCATE(DAT(IENDC))
          OPEN(1,FILE='temp.txt')
          READ(1,*)(DAT(I),I=1,IENDC)  
          CLOSE(1)
          OPEN(2,FILE=NADAT,FORM='UNFORMATTED')
          WRITE(2)IBNDAT
          WRITE(2)(DAT(I),I=1,IBNDAT)
          CLOSE(2)
          DEALLOCATE(DAT)
          IDUM(10)=IBEGC
          IDUM(11)=IENDC
          IDUM(12)=IBEGC
          IDUM(13)=IENDC
          RETURN

