      SUBROUTINE RFILLV(A,NA,C) 
	  INTEGER NA 
      REAL A(NA),C 
      DO 1000 I=1,NA 
      A(I)=C 
 1000 CONTINUE 
      RETURN 
      END 
