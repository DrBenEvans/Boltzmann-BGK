      SUBROUTINE RFILLM(A,NA,MA,C) 
      REAL A(NA,MA),C 
      DO 1000 J=1,MA 
      DO 2000 I=1,NA 
      A(I,J)=C 
 2000 CONTINUE 
 1000 CONTINUE 
      RETURN 
      END 