      VERR2=.FALSE.
      IF(FFVAR(49)(1:1).EQ.'Y'.OR.FFVAR(49)(1:1).EQ.'y')VERR2=.TRUE.
      WGT='D'
      IF(MAX(INDEX(FFVAR(50),'E'),INDEX(FFVAR(50),'e')).GT.0)WGT='E'
      WSM='N'
      IF(INDEX(FFVAR(51),'5').GT.0)WSM='5'
      IF(INDEX(FFVAR(51),'17').GT.0)WSM='7'
      IF(MAX(INDEX(FFVAR(51),'H'),INDEX(FFVAR(51),'h')).GT.0)WSM='H'
      VPTPLE=.FALSE.
      IF(FFVAR(52)(1:1).EQ.'Y'.OR.FFVAR(52)(1:1).EQ.'y')VPTPLE=.TRUE.
      DEALLOCATE(FFCOM,FFVAR)
      print*,'g3link LDUM(28) EgyEffPairs ',LDUM(28)
      PRINT*,'g3link LDUM(33) SIMPLIFIED EFF ',LDUM(33)
      IF(.NOT.LDUM(37))call dtodrho
      RETURN
