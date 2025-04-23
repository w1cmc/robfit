      SUBROUTINE CxyToX(CIN,COUT,Y,N)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION CIN(N*N+3*N+1),COUT(2*N+2)
C *** N IS HIGHEST POWER OF LNX AND LNY 
C IN 0-->1, 1-->5, 2-->11 3--> 19 4--> 29 5--> 41 6-->55 
C OUT 0 -->2 1-->4 2-->6  3-->8   4--> 10 5--> 12 6-->14     
      AL=LOG(Y)
      IF(N.EQ.0)THEN
        COUT(1)=CIN(1)
        COUT(2)=1
      ELSEIF(N.EQ.1)THEN
C (C1+C2*ALX+C3*ALY)/(1+C4*ALX+C5*ALY)
C -->((C1+C3*ALY)+C2*ALX)/((1+C5*ALY)+C4*ALX)      
        COUT(1)=CIN(1)+AL*CIN(3)
        COUT(2)=CIN(2)
        COUT(3)=1+AL*CIN(5)
        COUT(4)=CIN(4)
      ELSEIF(N.EQ.2)THEN
        COUT(1)=CIN(1)+AL*(CIN(3)+AL*CIN(6))
        COUT(2)=CIN(2)+AL*CIN(5)
        COUT(3)=CIN(4)
        COUT(1)=1+AL*(CIN(8)+AL*CIN(11))
        COUT(2)=CIN(7)+AL*CIN(10)
        COUT(3)=CIN(9)
      ELSEIF(N.EQ.3)THEN
        COUT(1)=CIN(1)+AL*(CIN(3)+AL*(CIN(6)+AL*CIN(10)))
        COUT(2)=CIN(2)+AL*(CIN(5)+AL*CIN(9))
        COUT(3)=CIN(4)+AL*CIN(8)
        COUT(4)=CIN(7)
        COUT(5)=1+AL*(CIN(12)+AL*(CIN(15)+AL*CIN(19)))
        COUT(6)=CIN(11)+AL*(CIN(14)+AL*CIN(18))
        COUT(7)=CIN(13)+AL*CIN(17)
        COUT(8)=CIN(16)
      ELSEIF(N.EQ.4)THEN
        COUT(1)=CIN(1)+AL*(CIN(3)+AL*(CIN(6)+AL*(CIN(10)+AL*CIN(15))))
        COUT(2)=CIN(2)+AL*(CIN(5)+AL*(CIN(9)+AL*CIN(14)))
        COUT(3)=CIN(4)+AL*(CIN(8)+AL*CIN(13))
        COUT(4)=CIN(7)+AL*CIN(12)
        COUT(5)=CIN(11)
        COUT(6)=1+AL*(CIN(17)+AL*(CIN(20)+AL*(CIN(24)+AL*CIN(29))))
        COUT(7)=CIN(16)+AL*(CIN(19)+AL*(CIN(23)+AL*CIN(28)))
        COUT(8)=CIN(18)+AL*(CIN(22)+AL*CIN(27))
        COUT(9)=CIN(21)+AL*CIN(26)
        COUT(10)=CIN(25)
      ELSEIF(N.EQ.5)THEN
        COUT(1)=CIN(1)+AL*(CIN(3)+AL*(CIN(6)+AL*(CIN(10)+AL*
     2   (CIN(15)+AL*CIN(21)))))
        COUT(2)=CIN(2)+AL*(CIN(5)+AL*(CIN(9)+AL*(CIN(14)+AL*CIN(20))))
        COUT(3)=CIN(4)+AL*(CIN(8)+AL*(CIN(13)+AL*CIN(19)))
        COUT(4)=CIN(7)+AL*(CIN(12)+AL*CIN(18))
        COUT(5)=CIN(11)+AL*CIN(17)
        COUT(6)=CIN(16)
        COUT(7)=1+AL*(CIN(23)+AL*(CIN(26)+AL*(CIN(30)+AL*
     2   (CIN(35)+AL*CIN(41)))))
        COUT(8)=CIN(22)+AL*(CIN(25)+AL*(CIN(29)+AL*(CIN(34)+
     2   AL*CIN(40))))
        COUT(9)=CIN(24)+AL*(CIN(28)+AL*(CIN(33)+AL*CIN(39)))
        COUT(10)=CIN(27)+AL*(CIN(32)+AL*CIN(38))
        COUT(11)=CIN(31)+AL*CIN(37)
        COUT(12)=CIN(36)
      ELSEIF(N.EQ.6)THEN
        COUT(1)=CIN(1)+AL*(CIN(3)+AL*(CIN(6)+AL*(CIN(10)+AL*
     2   (CIN(15)+AL*(CIN(21)+AL*CIN(28))))))
        COUT(2)=CIN(2)+AL*(CIN(5)+AL*(CIN(9)+AL*(CIN(14)+AL*(CIN(20)+
     2   AL*CIN(27)))))
        COUT(3)=CIN(4)+AL*(CIN(8)+AL*(CIN(13)+AL*(CIN(19)+AL*CIN(26))))
        COUT(4)=CIN(7)+AL*(CIN(12)+AL*(CIN(18)+AL*CIN(25)))
        COUT(5)=CIN(11)+AL*(CIN(17)+AL*CIN(24))
        COUT(6)=CIN(16)+AL*CIN(23)
        COUT(7)=CIN(22)
        COUT(8)=1+AL*(CIN(30)+AL*(CIN(33)+AL*(CIN(37)+AL*
     2   (CIN(42)+AL*(CIN(48)+AL*CIN(55))))))
        COUT(9)=CIN(29)+AL*(CIN(32)+AL*(CIN(36)+AL*(CIN(41)+AL*(CIN(47)
     2   +AL*CIN(54)))))
        COUT(10)=CIN(31)+AL*(CIN(35)+AL*(CIN(40)+AL*(CIN(46)
     2   +AL*CIN(53))))
        COUT(11)=CIN(34)+AL*(CIN(39)+AL*(CIN(45)+AL*CIN(52)))
        COUT(12)=CIN(38)+AL*(CIN(44)+AL*CIN(51))
        COUT(13)=CIN(43)+AL*CIN(50)
        COUT(14)=CIN(49)
      ELSE
        PRINT*,'N=',N,' HAS NOT BEEN CODED IN'
        STOP 17
      ENDIF    
      RETURN
      END
