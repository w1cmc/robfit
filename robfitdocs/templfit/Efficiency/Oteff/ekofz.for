        FUNCTION EKOFZ(Z)
        IMPLICIT REAL*8 (A-H,O-Z)
c The curve has a statistical error of less than 0.05 KeV.  The errors on
c Tl(80) and Bi(83) are + and - 0.12 KeV.        
        DIMENSION CONSN(4),CONSD(3)
        DATA CONSN/ 0.8004085118672191D-02,0.1633113527307094D-03,
     2   -0.4223372767326372D-05,0.2265957795031611D-07/
        DATA CONSD/0.2834374736285381E-02,-0.2402736923270774E-03,
     2   0.1455896647602572E-05/
       
        ANUM=Z*Z*(CONSN(1)+Z*(CONSN(2)+Z*(CONSN(3)+Z*CONSN(4))))
        ADEN=1+Z*(CONSD(1)+Z*(CONSD(2)+Z*CONSD(3)))
        EKOFZ=ANUM/ADEN
        RETURN
        END  
          
        FUNCTION EKOFZDER(Z,DEKDZ)
        IMPLICIT REAL*8 (A-H,O-Z)
c The curve has a statistical error of less than 0.05 KeV.  The errors on
c Tl(80) and Bi(83) are + and - 0.12 KeV.        
        DIMENSION CONSN(4),CONSD(3)
        DATA CONSN/ 0.8004085118672191D-02,0.1633113527307094D-03,
     2   -0.4223372767326372D-05,0.2265957795031611D-07/
        DATA CONSD/0.2834374736285381E-02,-0.2402736923270774E-03,
     2   0.1455896647602572E-05/
       
        ANUM=Z*Z*(CONSN(1)+Z*(CONSN(2)+Z*(CONSN(3)+Z*CONSN(4))))
        DNUM=Z*(2*CONSN(1)+Z*(3*CONSN(2)+Z*(4*CONSN(3)+Z*5*CONSN(4))))
        ADEN=1+Z*(CONSD(1)+Z*(CONSD(2)+Z*CONSD(3)))
        DDEN=CONSD(1)+Z*(2*CONSD(2)+Z*3*CONSD(3))
        EKOFZDER=ANUM/ADEN
        DEKDZ=(DNUM-EKOFZDER*DDEN)/ADEN
        RETURN
        END    

