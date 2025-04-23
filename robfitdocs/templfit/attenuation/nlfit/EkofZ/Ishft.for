C *** The function ISHFT, LSHFTB, IRSHFP, IRSHFM should
C *** be used for compilers whose internal ISHFT does not 
C *** appear to work properly.  For some compilers
C *** -2147483684 --> Z80000000

      FUNCTION ISHFT(IT,N)
      IF(N.GT.0)THEN
        ISHFT=LSHFTB(IT,N)
      ELSE
        IF(IT.GT.0)THEN
          ISHFT=IRSHFP(IT,-N)
        ELSE
          ISHFT=IRSHFM(IT,-N)
        ENDIF
      ENDIF
      RETURN
      END


      FUNCTION LSHFTB(IT,N)
      LSHFTB=IT*2**N
      RETURN
      END

      FUNCTION IRSHFP(IT,N)
      IRSHFP=IT/2**N
      RETURN
      END


      FUNCTION IRSHFM(IT,N)
      DATA IAA/-2147483648/
      ITT=IT+IAA
      IRSHFM=ITT/2**N+2**(31-N)
      RETURN
      END
